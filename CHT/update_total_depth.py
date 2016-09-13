#!/bin/env python
#
# Copyright 2013-14 Graham McVicker and Bryce van de Geijn
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
#
"""
This script updates the expected read depths of the input files 
for the Combined Haplotype Test. GC content biases and peakiness 
of the data (potentially caused by differences in ChIP efficiency 
or the expression of a few highly expressed genes) can vary greatly 
between individuals, sequencing libraries, and lanes. To correct 
for this, for each individual, we fit quartic splines across regions 
for GC content and combined (across all individuals) read depth.
The fitted splines provide an expected read depth for each
region for each individual. This expected read depth is used 
as an input to the CHT.
"""

import numpy as np
import sys
import gzip
import util
import argparse
import math
import tables
from scipy.optimize import *



def parse_options():
    parser=argparse.ArgumentParser(
        description="Adjusts the expected read depths for each "
        "region in the CHT input files. Expected read depths take "
        "into account the GC content and total read depth of each "
        "region (accross individuals). The adjusted read depths "
        "are computed by fitting quartic functions for these values "
        "for each individual.")

    parser.add_argument("infile_list", action='store', default=None)

    parser.add_argument("outfile_list", action='store', default=None)

    parser.add_argument("--seq", required=True,
                        help="Path to HDF5 file containing "
                        "genome sequence. Used to calculate GC content "
                        "of each region. Can be created "
                        "using fasta2h5 program.",
                        metavar="SEQ_H5_FILE")
    
    parser.add_argument("-i", "--fit_in_file", action='store',
                        dest='fit_in_file', default=None,
                        help="Read coefficients from specified file "
                        "rather than estimating them.")
    
    parser.add_argument("-o", "--fit_out_file", action='store',
                        dest='fit_out_file', default=None,
                        help="Estimate coefficients and write them "
                        "to specified file, but do not adjust read counts.")
    
    parser.add_argument("-m", "--min_counts", action='store', 
                        type=int, dest='min_counts', default=0,
                        help="only use rows with at least min_counts for fitting")
    
    parser.add_argument("--skip", action='store', type=int,
                        dest='skips', default=0,
                        help="specify a number of rows to skip between each row "
                        "used for estimating coefficients.")

    dflt_sample = 10000
    parser.add_argument("--sample", action='store', type=int,
                        default=dflt_sample, help="randomly sample this many rows "
                        "and use them for fitting coefficients. Specify 0 if all "
                        "rows are to be used. (default=%d)" % dflt_sample)
    

    return parser.parse_args()




def main():
    args = parse_options()
    inlist = [line.strip() for line in open(args.infile_list, "r")]

    for i in inlist:
        sys.stderr.write(i + "\n")
    sys.stderr.write("Loading count table\n")
  
    count_table, keep_list = load_data(inlist, args.seq,
                                       args.min_counts, args.skips)
    
    sys.stderr.write("Fitting coefficients\n")
    if args.fit_in_file:
        coefs_list = read_splines(args.fit_in_file)
    else:
        coefs_list = fit_splines(count_table)

    sys.stderr.write("Updating totals\n")
    if args.fit_out_file:
        write_splines(coefs_list, args.fit_out_file)
    else:
        outlist = [line.strip() for 
                   line in open(args.outfile_list,"r")]
        
        update_totals(inlist, outlist, count_table, coefs_list,
                      keep_list)

    sys.stderr.write("Finished!\n")

    
def open_files(file_list, r_w):
    files=[]
    
    for filename in file_list:
        if r_w.startswith("r"):
            if util.is_gzipped(filename):
                files.append(gzip.open(filename, r_w))
            else:
                files.append(open(filename, r_w))
        else:
            if filename.endswith(".gz"):
                files.append(gzip.open(filename, r_w))
            else:
                files.append(open(filename, r_w))
                
    
    return(files)


def get_at_gc_count(seq_h5, chrm, start, end):
    # seq HDF5 file contains ascii values for nucleotides
    # e.g. A = 65
    node = seq_h5.getNode("/%s" % chrm)
    vals = node[start-1:end]

    counts = np.bincount(vals)

    at_count = 0
    gc_count = 0
    
    if len(counts) >= ord("A"):
        at_count += counts[ord("A")]
    if len(counts) >= ord("T"):
        at_count += counts[ord("T")]
    if len(counts) >= ord("G"):
        gc_count += counts[ord("G")]
    if len(counts) >= ord("C"):
        gc_count += counts[ord("C")]
    
    return at_count, gc_count
        
        
    

    
def load_data(inlist, seq_h5_filename, min_counts, skips):
    infiles = open_files(inlist, "rt")

    seq_h5 = tables.openFile(seq_h5_filename, "r")
    
    end_of_file = False
    count_table = []
    keep_list = []

    info_list = []
    for infile in infiles:
        line = infile.readline()
        line = infile.readline()
        if not line:
            end_of_file = True
        else:
            info_list.append(line.strip().split())
            
    while not end_of_file:        
        snp_info = info_list[0]
        count_total = 0
        num_NAs = 0
        count_row = []
        chrm = snp_info[0]
        
        starts = [int(x) for x in snp_info[7].split(";")]
        ends = [int(x) for x in snp_info[8].split(";")]
        
        gc = 0
        at = 0
        for i in range(len(starts)):
            at_count, gc_count = get_at_gc_count(seq_h5, chrm,
                                                 starts[i],
                                                 ends[i])
            gc += gc_count
            at += at_count
                                
        for ind in range(len(infiles)):
            snp_info = info_list[ind]
            if snp_info[15] != "NA":
                new_count = int(snp_info[15])
                count_total += new_count
                count_row.append(new_count)
            else:
                num_NAs += 1
            for skip in range(skips):
                # skip lines
                line = infiles[ind].readline()
            line = infiles[ind].readline()

            if not line:
                end_of_file = True
            else:
                info_list[ind] = line.strip().split()

        if num_NAs == 0 and count_total >= min_counts and at+gc > 0:
            count_table.append([float(gc)/(at+gc), count_total] + 
                               count_row)
            keep_list.append(1)
        else:
            keep_list.append(0)
        for i in range(skips):
            keep_list.append(0)

    for infile in infiles:
        infile.close()
    count_table = np.array(count_table, dtype=np.float64)

    seq_h5.close()
    
    return count_table, keep_list



def read_splines(fit_in_file):
    spline_file = open(fit_in_file,"r")
    coefs = []
    for line in spline_file:
        coefs.append([float(x) for x in line.strip().split()])
    spline_file.close()
    return coefs


def write_splines(coefs,fit_out_file):
    spline_file = open(fit_out_file,"w")
    for ind in coefs:
        spline_file.write("\t".join([str(x) for x in ind]) + "\n")
    spline_file.close()
    return


def fit_splines(count_table):
    coefs=[]
    #iterate through each individual
    col_sums = np.sum(count_table,0)

    for col in range(2, count_table.shape[1]):
        slope_start = col_sums[col]/col_sums[1]
        sys.stderr.write(str(slope_start)+"\n")

        arg_list = (count_table[:,0], count_table[:,1], 
                    count_table[:,col])
        slope = fmin(splineone, [slope_start], args=arg_list, 
                     maxiter=500000, maxfun=500000)
        
        # Fit a slope only model (total read depth is different 
        # between individuls within an individual reads per region 
        # is linear
        coefs_slope = [0,0,0,0,0,0,slope[0],0,0,0]
        sys.stderr.write(str(coefs_slope)+"\n")
        
        # Read shift only (allowing more or fewer reads to fall 
        # in peaks)
        arg_list = ([0,0,0,0,0], count_table[:,0],
                    count_table[:,1], count_table[:,col])
        coefs_depth = fmin(splineread, [0, slope[0],0,0,0],
                           args=arg_list, maxiter=500000,
                           maxfun=500000)
        sys.stderr.write(str(coefs_depth) + "\n")

        # GC correct and total count only
        arg_list = ([0,slope[0],0,0,0], count_table[:,0],
                    count_table[:,1], count_table[:,col])
        
        coefs_gc = fmin(splinegc, [0,0,0,0,0],
                        args=arg_list, maxiter=500000, maxfun=500000)
        sys.stderr.write(str(coefs_gc)+"\n")
        
        # GC correct plus count shift
        arg_list = (count_table[:,0],count_table[:,1],
                    count_table[:,col])
        coefs_full = fmin(splinefit, 
                          list(coefs_gc) + list(coefs_depth),
                          args=arg_list, 
                          maxiter=500000, maxfun=500000)
        sys.stderr.write(str(coefs_full) + "\n")
        coefs.append(coefs_full)
    
    return(coefs)

    
def splineone(slope,x1,x2,y):
    return splinefit([0,0,0,0,0,0,slope[0],0,0,0],x1,x2,y)


def splinegc(gcargs,otherargs,x1,x2,y):
    a=np.hstack((gcargs,otherargs))
    return splinefit(a,x1,x2,y)


def splineread(otherargs,gcargs,x1,x2,y):
    a=np.hstack((gcargs,otherargs))
    return splinefit(a,x1,x2,y)


def splinefit(arg, x1, x2,  y):
    expecteds=[max(0.000001, x) 
               for x in calc_adjusted_totals(x1, x2, arg)]
  
    resids = y - expecteds


    sse=sum(resids**2)
    loglike=-sum(y*[math.log(x) for x in expecteds]-expecteds)

    return loglike #sse


def update_totals(inlist, outlist, count_table, coefs_table,
                  keep_list):
    infiles = open_files(inlist,"rt")
    outfiles = open_files(outlist,"w")
    row = 0
    count_row = 0

    # Write the header line
    for ind in range(len(infiles)):
        line = infiles[ind].readline()
        outfiles[ind].write(line)
        
    while row < len(keep_list):
        # calc expected counts for each individual
        adj_tot_list = []
        has_NAs = False

        if keep_list[row] == 1:
            for ind in range(len(infiles)):

                if count_table[count_row,2+ind]<0:
                    has_NAs=True
                    break

                adj_tot = calc_adjusted_totals(count_table[count_row, 0],
                                               count_table[count_row, 1],
                                               coefs_table[ind])
                adj_tot_list.append(max(adj_tot, -1000000))
                
            for ind in range(len(infiles)):
                line = infiles[ind].readline()
                snp_info = line.strip().split()

                if has_NAs:
                    outfiles[ind].write(line)
                else:
                    adj_tot_ind = adj_tot_list[ind] # / min(adj_tot_list)

                    outfiles[ind].write("\t".join(snp_info[:16]) +
                                        "\t"+str(adj_tot_ind) + "\n")
                outfiles[ind].flush()
            count_row+=1
        else:
            for ind in range(len(infiles)):
                line = infiles[ind].readline()
        row += 1
    for infile in infiles:
        infile.close()
    for outfile in outfiles:
        outfile.close()                           

        
def array_exp(x):
    if hasattr(x, "__len__"):
        return [math.exp(i) for i in x]
    else:
        return math.exp(x)

    
def calc_adjusted_totals(gc, tot, coefs):
    return array_exp(coefs[0] + gc*coefs[1] + gc**2*coefs[2]
                     + gc**3*coefs[3] + gc**4*coefs[4]) * \
                     (0 + tot*coefs[6] + tot**2*coefs[7] + 
                      tot**3*coefs[8] + tot**4*coefs[9])


main()
