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
usage: update_total_depth.py input_file_list output_directory

positional arguments:
  input_file_list       name of the .txt file containing a list of files for input
  out_directory         directory to output the updated files


optional arguments:
  None

This script is used to update the expected read depths of the input files for the
Combined Haplotype Test.  GC content biases (potential caused by the PCR step)
and peakiness of the data (potentially caused by differences in ChIP efficiency
or the expression of a few highly expressed genes) can vary greatly between
individuals, sequencing libraries, and lanes.  To correct for this, for each
individual, we fit quartic splines across regions for GC content and combined (
across all individuals) read depth.  This results in a fitted expected read depth
for each region x individual which is then used for the CHT.
"""


import numpy as np
import sys, pdb, gzip, argparse, genome.db, math
from scipy.optimize import *

def main():
    args=parse_options()
    inlist=[i for i in open(args.infile_list,"r")]

    for i in inlist:
        sys.stderr.write(i+"\n")
    sys.stderr.write("Loading count table\n")
  
    count_table,keep_list=load_data(inlist,args.min_counts, args.skips)
    
    #countfile=open("count_file.txt","w")
    #for i in range(count_table.shape[0]):
    #    countfile.write(" ".join([str(x) for x in count_table[i,]])+"\n")
    #exit()
    sys.stderr.write("Fitting coefficients\n")
    if args.fit_in_file:
        coefs_list=read_splines(args.fit_in_file)
    else:
        coefs_list=fit_splines(count_table)    

    sys.stderr.write("Updating totals\n")
    if args.fit_out_file:
        write_splines(coefs_list,args.fit_out_file)
    else:
        outlist=[i for i in open(args.outfile_list,"r")]
        update_totals(inlist,outlist,count_table,coefs_list,keep_list)

    sys.stderr.write("Finished!\n")
    
def open_files(file_list,r_w):
    files=[]
    for f in file_list:
        split_name=f.strip().split(".")
        if f.strip()[-3:]==".gz":
            files.append(gzip.open(f.strip(),r_w))
        else:
            files.append(open(f.strip(),r_w))
    return(files)

def load_data(inlist,min_counts,skips):
    infiles=open_files(inlist,"r")
    gdb=genome.db.GenomeDB(assembly="hg19")
    seq_track=gdb.open_track("seq")
    end_of_file=False
    count_table=[]
    keep_list=[]

    info_list=[]
    for infile in infiles:
        line=infile.readline()
        line=infile.readline()
        if not line:
            end_of_file=True
        else:
            info_list.append(line.strip().split())
        
    while not end_of_file:
        snp_info=info_list[0]
        count_total=0
        num_NAs=0
        count_row=[]
        chrm=snp_info[0]
        
        starts=[int(x) for x in snp_info[7].split(";")]
        ends=[int(x) for x in snp_info[8].split(";")]
        
        gc=0
        at=0
        for i in range(len(starts)):
            sequence=seq_track.get_seq_str(chrm,starts[i],ends[i])
            for base in list(sequence):
                if base == "C" or base == "G":
                    gc+=1
                if base == "A" or base == "T":
                    at+=1
        
        for ind in range(len(infiles)):
            snp_info=info_list[ind]
            if snp_info[15]!="NA":
                new_count=int(snp_info[15])
                count_total+=new_count
                count_row.append(new_count)
            else:
                num_NAs+=1
            for skip in range(skips):
                line=infiles[ind].readline()
            line=infiles[ind].readline()
            #sys.stderr.write(snp_info[1]+"\n")
            if not line:
                end_of_file=True
            else:
                info_list[ind]=line.strip().split()

        if num_NAs==0 and count_total>=min_counts and at+gc > 0:
            count_table.append([float(gc)/(at+gc),count_total]+count_row)
            keep_list.append(1)
        else:
            keep_list.append(0)
        for i in range(skips):
            keep_list.append(0)

    for infile in infiles:
        infile.close()
    count_table=np.array(count_table,dtype=np.float64)
    return count_table, keep_list

def read_splines(fit_in_file):
    spline_file=open(fit_in_file,"r")
    coefs=[]
    for line in spline_file:
        coefs.append([float(x) for x in line.strip().split()])
    spline_file.close()
    return coefs

def write_splines(coefs,fit_out_file):
    spline_file=open(fit_out_file,"w")
    for ind in coefs:
        spline_file.write("\t".join([str(x) for x in ind])+"\n")
    spline_file.close()
    return

def fit_splines(count_table):
    coefs=[]
    #iterate through each individual
    col_sums=np.sum(count_table,0)
    #pdb.set_trace()
    for col in range(2,count_table.shape[1]):
        slope_start=col_sums[col]/col_sums[1]
        sys.stderr.write(str(slope_start)+"\n")
        slope=fmin(splineone,[slope_start],args=(count_table[:,0],count_table[:,1],count_table[:,col]),maxiter=500000,maxfun=500000)
        # Fit a slope only model (total read depth is different between individuls
        # within an individual reads per region is linear
        coefs_slope=[0,0,0,0,0,0,slope[0],0,0,0]
        sys.stderr.write(str(coefs_slope)+"\n")
        # Read shift only (allowing more or fewer reads to fall in peaks)
        coefs_depth=fmin(splineread,[0,slope[0],0,0,0],args=([0,0,0,0,0],count_table[:,0],count_table[:,1],count_table[:,col]),maxiter=500000,maxfun=500000)
        sys.stderr.write(str(coefs_depth)+"\n")
        # GC correct and total count only
        coefs_gc=fmin(splinegc,[0,0,0,0,0],args=([0,slope[0],0,0,0],count_table[:,0],count_table[:,1],count_table[:,col]),maxiter=500000,maxfun=500000)
        sys.stderr.write(str(coefs_gc)+"\n")
        # GC correct plus count shift
        coefs_full=fmin(splinefit,list(coefs_gc)+list(coefs_depth),args=(count_table[:,0],count_table[:,1],count_table[:,col]),maxiter=500000,maxfun=500000)
        sys.stderr.write(str(coefs_full)+"\n")
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
    #pdb.set_trace()
    #if arg[5]<0:
    #    return 100000000000000000000
    
    expecteds=[max(0.000001,x) for x in calc_adjusted_totals(x1,x2,arg)]
  
    resids=y-expecteds

    #print arg
    sse=sum(resids**2)
    loglike=-sum(y*[math.log(x) for x in expecteds]-expecteds)
    #print sse
    #sys.stdout.flush()
    return loglike #sse

def update_totals(inlist,outlist,count_table,coefs_table,keep_list):
    infiles=open_files(inlist,"r")
    outfiles=open_files(outlist,"w")
    row=0
    count_row=0

    # Write the header line
    for ind in range(len(infiles)):
        line=infiles[ind].readline()
        outfiles[ind].write(line)
        
    while row < len(keep_list):
        # calc expected counts for each individual
        adj_tot_list=[]
        has_NAs=False
        #pdb.set_trace()
        if keep_list[row]==1:
            for ind in range(len(infiles)):
                #sys.stderr.write(str(count_table[row,2+ind])+"\n")
                if count_table[count_row,2+ind]<0:
                    has_NAs=True
                    #sys.stderr.write("Updating totals\n")
                    break
                adj_tot=calc_adjusted_totals(count_table[count_row,0],count_table[count_row,1],coefs_table[ind])
                adj_tot_list.append(max(adj_tot,-1000000))
            for ind in range(len(infiles)):
                line=infiles[ind].readline()
                snp_info=line.strip().split()

                if has_NAs:
                    outfiles[ind].write(line)
                else:
                    adj_tot_ind=adj_tot_list[ind]#/min(adj_tot_list)
                    #sys.stderr.write("Updating totals\n")
                    outfiles[ind].write("\t".join(snp_info[:16])+"\t"+str(adj_tot_ind)+"\n")
                outfiles[ind].flush()
            count_row+=1
        else:
            for ind in range(len(infiles)):
                line=infiles[ind].readline()
        row+=1
    for infile in infiles:
        infile.close()
    for outfile in outfiles:
        outfile.close()                           
        
def array_exp(x):
    if hasattr(x,"__len__"):
        return [math.exp(i) for i in x]
    else:
        return math.exp(x)

def calc_adjusted_totals(gc,tot,coefs):
    #pdb.set_trace()
    return array_exp(coefs[0]+gc*coefs[1]+gc**2*coefs[2]+gc**3*coefs[3]+gc**4*coefs[4]) * (0+tot*coefs[6]+tot**2*coefs[7]+tot**3*coefs[8]+tot**4*coefs[9])
    #return array_exp(coefs[0]+gc*coefs[1]+gc**2*coefs[2]+gc**3*coefs[3]+gc**4*coefs[4]) * (coefs[5]+tot*coefs[6]+tot**2*coefs[7]+tot**3*coefs[8]+tot**4*coefs[9])

def parse_options():
    parser=argparse.ArgumentParser()
    parser.add_argument("infile_list", action='store', default=None)
    parser.add_argument("outfile_list", action='store', default=None)
    parser.add_argument("-i", action='store',
                        dest='fit_in_file', default=None,
                        help="read coefficients from specified file  ")
    parser.add_argument("-o", action='store',
                        dest='fit_out_file', default=None,
                        help="only fit the model and write to specified file ")
    parser.add_argument("-m", action='store', type=int,
                        dest='min_counts', default=0,
                        help="minimum counts to use row ")
    parser.add_argument("--skip", action='store', type=int,
                        dest='skips', default=0,
                        help="lines to skip between each used for calculation ")
    return parser.parse_args()

main()
