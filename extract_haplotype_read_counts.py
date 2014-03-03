#!/bin/env python
#
# Copyright 2013 Graham McVicker and Bryce van de Geijn
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
usage: extract_haplotype_read_counts.py [-h] [--assembly ASSEMBLY]
                                        [--target_region_size TARGET_REGION_SIZE]
                                        [--sample_file SAMPLE_FILE]
                                        [--homozygous_as_counts {zero,rand_hap,rand_allele}]
                                        track_prefix pop individual input_file

positional arguments:
  track_prefix          prefix of tracks to extract reads from (e.g.
                        10_IND/PolII/read_counts/PolII_18505)
  pop                   population prefix for genotype tracks (YRI or CEU)
  individual            individual to extract read counts for (e.g. 18505)
  input_file            bed-like file to read coordinates of test SNP and
                        target region from

optional arguments:
  -h, --help            show this help message and exit
  --assembly ASSEMBLY   genome assembly that reads were mapped to (e.g. hg18)
  --target_region_size TARGET_REGION_SIZE
                        override target region size that is specified by input
                        file
  --sample_file SAMPLE_FILE
                        path to file containing ordered list of genotyped
                        individuals
  --homozygous_as_counts {zero,rand_hap,rand_allele}
                        how to report AS counts at linked het SNPs when test
                        SNP genotype is homozygous or unknown. zero (default):
                        set allele-specific counts to 0; rand_hap: randomly
                        choose one of the haplotypes to be 'reference';
                        rand_allele: choose random allele at each SNP to be
                        reference

This script is used to generate input files for the combined haplotype
test script.  It depends on a number of datafiles, which may make it
difficult for other people to use. More specifically this script reads
data from HDF5 files (a.k.a. tracks) and uses code from the 'genome'
library (https://github.com/gmcvicker/genome) to access them.

The script reads from the following HDF5 tracks. <PREFIX> and <POP> are specified by
positional command line arguments "track_prefix", "pop":
  <PREFIX>_AS_ref_count - number of allele-specific reads that match ref allele at each SNP
  <PREFIX>_AS_alt_count - number of allele-specific reads that match alt allele at each SNP
  <PREFIX>_AS_other_count - number of reads that match neither ref nor alt allele at each SNP
  <PREFIX>_read_start_count - number of aligned reads that start at each position

  <PREFIX>_impute2/snps - table with info about each SNP including alleles, and position
  <PREFIX>impute2/snp_index - mapping from genomic position to index in snps table
  impute2/<POP>_geno_probs - genotype probabilites for each individual
  impute2/<POP>_haplotypes - phasing information for alleles
  
"""

import argparse
import numpy as np
import sys
import gzip

import genome.db
import genome.coord

SNP_UNDEF = -1
HAP_UNDEF = -1


class SNP(object):
    def __init__(self, chrom, pos, name, ref_allele, alt_allele):
        self.chrom = chrom
        self.pos = pos
        self.name = name
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        
    


class DataTracks(object):
    def __init__(self, gdb, track_prefix, pop):

        # open tracks that read counts will be pulled from
        ref_as_name = "%s_AS_ref_count" % track_prefix

        alt_as_name = "%s_AS_alt_count" % track_prefix

        oth_as_name = "%s_AS_other_count" % track_prefix

        read_count_name = "%s_read_start_count" % track_prefix

        self.ref_count_track = gdb.open_track(ref_as_name)
        self.alt_count_track = gdb.open_track(alt_as_name)
        self.other_count_track = gdb.open_track(oth_as_name)
        self.read_count_track = gdb.open_track(read_count_name)

        # open tracks where SNP information can be extracted
        self.snp_track = gdb.open_track("impute2/snps")
        self.snp_index_track = gdb.open_track("impute2/snp_index")
        self.geno_track = gdb.open_track("impute2/%s_geno_probs" % 
                                         pop.lower())
        self.hap_track = gdb.open_track("impute2/%s_haplotypes" % 
                                        pop.lower())

    
    def close(self):
        """closes all of the data tracks"""

        self.ref_count_track.close()
        self.alt_count_track.close()
        self.other_count_track.close()
        self.read_count_track.close()

        self.snp_track.close()
        self.snp_index_track.close()
        self.geno_track.close()
        self.hap_track.close()
        


        
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--assembly", help="genome assembly that reads "
                        "were mapped to (e.g. hg18)", default=None)

    parser.add_argument('--target_region_size', 
                        help='override target region size that is '
                        'specified by input file',
                        type=int, default=None)

    parser.add_argument("--sample_file",
                        help="path to file containing ordered list of "
                        "genotyped individuals",
                        default=None)
    
    parser.add_argument("--homozygous_as_counts",
                        help="how to report AS counts at linked het SNPs when "
                        "test SNP genotype is homozygous or unknown. "
                        "zero (default): set allele-specific counts to 0; "
                        "rand_hap: randomly choose one of the haplotypes "
                        "to be 'reference'; "
                        "rand_allele: choose random allele at each SNP to "
                        "be reference", default="zero",
                        choices=("zero", "rand_hap", "rand_allele"))
                        

    parser.add_argument("track_prefix", 
                        help="prefix of tracks to extract reads from "
                        "(e.g. 10_IND/PolII/read_counts/PolII_18505)")

    parser.add_argument("pop",
                        help="population prefix for genotype tracks "
                        "(YRI or CEU)")
    
    parser.add_argument("individual", 
                        help="individual to extract read counts "
                        "for (e.g. 18505)")
    
    parser.add_argument("input_file", 
                        help="bed-like file to read coordinates of "
                        "test SNP and target region from")

    args = parser.parse_args()

    return args





def get_region_snps(dt, region_list, ind_idx):
    """Retrieves all of the SNPs in the requested regions. 
    The test SNP is also returned."""

    if len(region_list) == 0:
        raise genome.coord.CoordError("expected at least one coordinate, got 0")

    chrom = region_list[0].chrom
    
    snp_tab = dt.snp_track.h5f.getNode("/%s" % chrom.name)
    hap_tab = dt.hap_track.h5f.getNode("/%s" % chrom.name)
    geno_tab = dt.geno_track.h5f.getNode("/%s" % chrom.name)

    region_snps = []

    for region in region_list:
        if region.chrom.name != chrom.name:
            raise CoordError("only regions on same chromosome are supported")
        
        snp_idx = dt.snp_index_track.get_nparray(chrom, region.start, region.end)

        offsets = np.where(snp_idx != SNP_UNDEF)[0]

        test_snp = None

        for offset in offsets:
            i = snp_idx[offset]
            snp_row = snp_tab[i]

            # extract geno probs and haplotypes for this individual
            geno_probs = geno_tab[i, (ind_idx*3):(ind_idx*3 + 3)]
            haps = hap_tab[i, (ind_idx*2):(ind_idx*2 + 2)]

            snp = SNP(region.chrom, snp_row['pos'], 
                      snp_row['name'],
                      snp_row['allele1'],
                      snp_row['allele2'])

            # get heterozygote probability for SNP
            snp.het_prob = geno_probs[1]

            # linear combination of genotype probs:
            #     0*homo_ref + 1*het + 2*homo_alt
            snp.geno_sum = geno_probs[1] + 2.0*geno_probs[2]
            snp.haps = haps

            # TODO: set linkage probabilty properly
            snp.linkage_prob = 1.0

            region_snps.append(snp)
    
    return region_snps



def get_het_snps(snp_list):
    het_snp_list = []
    
    for snp in snp_list:
        if snp.haps[0] != snp.haps[1]:
            het_snp_list.append(snp)
    
    return het_snp_list
            
            
        

def lookup_individual_index(options, ind_name):
    """Gets the index of individual that is used 
    to lookup information in the genotype and haplotype tables"""

    if options.sample_file is None:
        sample_file = "/data/share/10_IND/IMPUTE/%s_samples.txt" % options.pop
    else:
        sample_file = options.sample_file

    sys.stderr.write("reading list of individuals from %s\n" % sample_file)
    f = open(sample_file)

    idx = 0
    for line in f:
        words = line.rstrip().split()

        name = words[0].replace("NA", "")
        if name == ind_name:
            f.close()
            return idx
        
        idx += 1

    raise ValueError("individual %s is not in samples file %s" %
                     (ind_name, options.sample_file))




def set_snp_counts(dt, region_list, snps, test_snp, options):
    """Sets counts of reference and alternate haplotype matching reads
    for each of the provided SNPs. Labeling of 'reference' or 'alternate'
    is with respect to the test SNP"""

    if test_snp and (test_snp.haps[0] != test_snp.haps[1]) and \
      (test_snp.haps[0] != HAP_UNDEF):      
        # test SNP is heterozygous: use this to phase counts that are
        # retrieved at linked het SNPs
        if test_snp.haps[0] == 0:
            # reference allele is first haplotype at test SNP
            ref_idx = 0
            alt_idx = 1
        else:
            # alt allele is first haplotype at test SNP
            ref_idx = 1
            alt_idx = 0
    else:
        # test SNP is homozygous or is undefined
        # so we have no way to tell which haplotype it is on
        if options.homozygous_as_counts == "rand_hap":
            # choose haplotype randomly
            if np.random.randint(2) == 0:
                ref_idx = 0
                alt_idx = 1
            else:
                ref_idx = 1
                alt_idx = 0
        else:
            ref_idx = None
            alt_idx = None
                

    for region in region_list:
        ref_counts = dt.ref_count_track.get_nparray(region.chrom, region.start,
                                                    region.end)
        
        alt_counts = dt.alt_count_track.get_nparray(region.chrom, region.start,
                                                    region.end)

        other_counts = dt.other_count_track.get_nparray(region.chrom, 
                                                        region.start,
                                                        region.end)
        for snp in snps:
            # we have het SNPs from several regions, but only want to consider
            # ones in current region
            if snp.pos >= region.start and snp.pos <= region.end:
                offset = snp.pos - region.start

                ref_count = ref_counts[offset]
                alt_count = alt_counts[offset]
                snp.other_count = other_counts[offset]

                if ref_idx is None:
                    if options.homozygous_as_counts == "zero":
                        snp.ref_hap_count = 0
                        snp.alt_hap_count = 0
                    elif options.homozygous_as_counts == "rand_allele":
                        # choose allele randomly to be reference
                        if np.random.randint(2) == 0:
                            snp.ref_hap_count = ref_count
                            snp.alt_hap_count = alt_count
                        else:
                            snp.ref_hap_count = alt_count
                            snp.alt_hap_count = ref_count
                    else:
                        raise ValueError("unknown homozygous_as_counts option %s" % 
                                         options.homozygous_as_counts)
                else:
                    if snp.haps[ref_idx] == 0:
                        # reference allele is on "reference" haplotype
                        snp.ref_hap_count = ref_count
                        snp.alt_hap_count = alt_count
                    elif snp.haps[ref_idx] == 1:
                        # reference allele is on "alternate" haplotype
                        snp.ref_hap_count = alt_count
                        snp.alt_hap_count = ref_count
                    else:
                        raise ValueError("expected haplotype to be defined")
    
        

def write_header(f):
    f.write("CHROM "
            "TEST.SNP.POS "
            "TEST.SNP.ID "
            "TEST.SNP.REF.ALLELE "
            "TEST.SNP.ALT.ALLELE "
            "TEST.SNP.GENOTYPE "
            "TEST.SNP.HAPLOTYPE "
            "REGION.START "
            "REGION.END "
            "REGION.SNP.POS "
            "REGION.SNP.HET.PROB "
            "REGION.SNP.LINKAGE.PROB "
            "REGION.SNP.REF.HAP.COUNT "
            "REGION.SNP.ALT.HAP.COUNT "
            "REGION.SNP.OTHER.HAP.COUNT "
            "REGION.READ.COUNT "
            "GENOMEWIDE.READ.COUNT\n")
    

def write_NA_line(f):
    f.write(" ".join(["NA"] * 15) + "\n")
    
    
def write_output(f, region_list, snps, test_snp, test_snp_pos, region_read_count, 
                 genomewide_read_count):


    chrom_name = region_list[0].chrom.name
    region_start_str = ";".join([str(r.start) for r in region_list])
    region_end_str = ";".join([str(r.end) for r in region_list])
    
    if test_snp is None:
        # the SNP did not exist, probably was removed between
        # 1000 genomes releases
        f.write("%s %d NA NA NA NA NA %s %s" %
                (chrom_name, test_snp_pos, region_start_str, region_end_str))
        
        f.write(" %s\n" % " ".join(["NA"] * 8))

        return
    
    f.write("%s %d %s %s %s %.2f %d|%d %s %s" %
            (test_snp.chrom.name, test_snp.pos, test_snp.name,
             test_snp.ref_allele, test_snp.alt_allele,
             test_snp.geno_sum, test_snp.haps[0], 
             test_snp.haps[1], region_start_str, region_end_str))
    
    # number of linked heterozygous SNPs that we can pull
    # haplotype-specific counts from
    n_het_snps = len(snps)

    if n_het_snps > 0:
        # write SNP positions
        f.write(" %s" % ";".join(["%d" % s.pos for s in snps]))

        # write SNP het probs
        f.write(" %s" % ";".join(["%.2f" % s.het_prob for s in snps]))

        # write SNP linkage probs
        f.write(" %s" % ";".join(["%.2f" % s.linkage_prob for s in snps]))

        # write SNP ref/alt/other haplotype counts
        f.write(" %s" % ";".join(["%d" % s.ref_hap_count for s in snps]))
        f.write(" %s" % ";".join(["%d" % s.alt_hap_count for s in snps]))
        f.write(" %s" % ";".join(["%d" % s.other_count for s in snps]))
    else:
        # no linked heterozygous SNPs
        f.write(" %s" % " ".join(["NA"] * 6))

    # write total read count for region and genome-wide read count
    f.write(" %d %d\n" % (region_read_count, genomewide_read_count))
    
    
    
    
def get_genomewide_count(gdb, track):
    stat = genome.trackstat.get_stats(gdb, track)
    return stat.sum


def get_region_read_counts(dt, region_list):
    total_count = 0

    for region in region_list:
        counts = dt.read_count_track.get_nparray(region.chrom, region.start, region.end)
        total_count += np.sum(counts)

    return total_count



def get_target_regions(args, chrom, words):
    """Parse start and end positions and return list of Coord objects representing
    target region(s)."""

    start_words = words[7].split(";")
    end_words = words[8].split(";")
    
    if len(start_words) != len(end_words):
        raise genome.coord.CoordError("number of start (%d) and end (%d) positions "
                                      "do not match" % (len(start_words), len(end_words)))
    
    n_coord = len(start_words)

    region_list = []
    for i in range(n_coord):
        start = int(start_words[i])
        end = int(end_words[i])        

        region = genome.coord.Coord(chrom, start, end)
        
        if args.target_region_size:
            if region.length() != args.target_region_size:
                # override the size of the target region
                # with size provided on command line
                mid = (region.start + region.end)/2
                region.start = mid - args.target_region_size/2
                region.end = mid + args.target_region_size/2
                if region.start < 1:
                    region.start = 1
                if region.end > chrom.length:
                    region.end = chrom.length

        region_list.append(region)

    return region_list



def main():
    args = parse_args()

    gdb = genome.db.GenomeDB(assembly=args.assembly)

    write_header(sys.stdout)
    
    individual = args.individual.replace("NA","").replace("GM", "")
    ind_idx = lookup_individual_index(args, individual)
    
    dt = DataTracks(gdb, args.track_prefix, args.pop)

    genomewide_read_counts = get_genomewide_count(gdb, dt.read_count_track)
    
    chrom_dict = gdb.get_chromosome_dict()

    if args.input_file.endswith(".gz"):
        f = gzip.open(args.input_file)
    else:
        f = open(args.input_file)

    line_count = 0

    if args.target_region_size:
        sys.stderr.write("setting target region size to %d\n" %
                         args.target_region_size)
    
    for line in f:
        line_count += 1
        if line_count % 1000 == 0:
            sys.stderr.write(".")

        if line.startswith("#"):
            continue
        
        words = line.rstrip().split()

        if words[1] == "NA":
            # no SNP defined on this line:
            write_NA_line(sys.stdout)
            continue

        
        chrom_name = words[0]
        chrom = chrom_dict[chrom_name]
        
        region_list = get_target_regions(args, chrom, words)

        snp_pos = int(words[1])
        snp_ref_base = words[3]
        snp_alt_base = words[4]
        # TODO: check that SNP ref/alt match?
                    
        snp_region = genome.coord.Coord(chrom, snp_pos, snp_pos)
        
        # pull out all of the SNPs in the target region(s)
        region_snps = get_region_snps(dt, region_list, ind_idx)

        # pull out test SNP
        test_snp_list = get_region_snps(dt, [snp_region], ind_idx)
        if len(test_snp_list) != 1:
            test_snp = None
            sys.stderr.write("WARNING: could not find test SNP at "
                             "position %s:%d\n" % (chrom.name, snp_pos))
            het_snps = []
        else:
            test_snp = test_snp_list[0]
                
            # pull out haplotype counts from linked heterozygous SNPs
            het_snps = get_het_snps(region_snps)
            set_snp_counts(dt, region_list, het_snps, test_snp, args)

        region_read_counts = get_region_read_counts(dt, region_list)

        write_output(sys.stdout, region_list, het_snps, test_snp, snp_pos,
                     region_read_counts, genomewide_read_counts)

    sys.stderr.write("\n")
    f.close()
    dt.close()

        

main()        

        
                             
