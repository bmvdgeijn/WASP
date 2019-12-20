import sys
import os

import numpy as np

import tables

import chromosome
import argparse

import gzip

import util

# TODO: change output dir--make command line arg?
OUTPUT_DIR = "."

SNP_UNDEF = -1


class SNP(object):
    def __init__(self, chrom, pos, name, ref_allele, alt_allele):
        self.chrom = chrom
        self.pos = pos
        self.name = name
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        


class SNPFiles(object):
    def __init__(self, args):
        # open tracks where SNP information can be extracted
        self.snp_tab_h5 = tables.open_file(args.snp_tab, "r")
        self.snp_index_h5 = tables.open_file(args.snp_index, "r")
        self.hap_h5 = tables.open_file(args.haplotype, "r")


    def close(self):
        self.snp_tab_h5.close()
        self.snp_index_h5.close()
        self.hap_h5.close()
        


class CombinedFiles(object):
    def __init__(self, output_dir, chrom_list):
        # combined allele-specific read counts
        as_count_filename = "%s/combined_as_count.h5" % output_dir
        self.as_count_h5 = tables.open_file(as_count_filename, "w")
        
        # combined mapped read counts
        read_count_filename = "%s/combined_read_count.h5" % output_dir
        self.read_count_h5 = tables.open_file(read_count_filename, "w")

        # counts of genotypes
        ref_count_filename = "%s/combined_ref_count.h5" % output_dir
        self.ref_count_h5 = tables.open_file(ref_count_filename, "w")
        
        alt_count_filename = "%s/combined_alt_count.h5" % output_dir
        self.alt_count_h5 = tables.open_file(alt_count_filename, "w")
        
        het_count_filename = "%s/combined_het_count.h5" % output_dir
        self.het_count_h5 = tables.open_file(het_count_filename, "w")
        
        self.filenames = [as_count_filename, read_count_filename,
                          ref_count_filename, alt_count_filename,
                          het_count_filename]

        self.h5_files = [self.as_count_h5, self.read_count_h5,
                         self.ref_count_h5, self.alt_count_h5, 
                         self.het_count_h5]

        # initialize all of these files
        atom = tables.UInt16Atom(dflt=0)
        
        for h5f in self.h5_files:
            for chrom in chrom_list:
                self.create_carray(h5f, chrom, atom)
                
    def create_carray(self, h5f, chrom, atom):
        zlib_filter = tables.Filters(complevel=1, complib="zlib")

        # create CArray for this chromosome
        shape = [chrom.length]
        carray = h5f.create_carray(h5f.root, chrom.name,
                                  atom, shape, filters=zlib_filter)

        return carray


    
    def add_counts(self, chrom_list, ind_files, snp_files, ind_idx):
        # add contribution from one individual to combined counts
        for chrom in chrom_list:
            node_name = "/%s" % chrom.name

            if node_name not in ind_files.ref_as_count_h5:
                continue
            if node_name not in ind_files.alt_as_count_h5:
                continue
            if node_name not in ind_files.read_count_h5:
                continue
            if node_name not in snp_files.hap_h5:
                continue
            if node_name not in snp_files.snp_index_h5:
                continue
            if node_name not in snp_files.snp_tab_h5:
                continue

            sys.stderr.write("  %s\n" % chrom.name)
            
            node = self.as_count_h5.get_node(node_name)
                
            ind_node = ind_files.ref_as_count_h5.get_node(node_name)
            node[:] += ind_node[:]

            ind_node = ind_files.alt_as_count_h5.get_node(node_name)
            node[:] += ind_node[:]

            node = self.read_count_h5.get_node(node_name)
            ind_node = ind_files.read_count_h5.get_node("/%s" % chrom.name)
            node[:] += ind_node[:]


            # get haplotypes for this individual
            hap_a_idx = ind_idx*2
            hap_b_idx = ind_idx*2+1            
            hap_tab = snp_files.hap_h5.get_node("/%s" % chrom.name)
            a_hap = hap_tab[:, hap_a_idx]
            b_hap = hap_tab[:, hap_b_idx]

            # determine genotype of SNPs for this individual
            is_homo_ref = (a_hap == 0) & (b_hap == 0)
            is_het      = ((a_hap == 0) & (b_hap == 1)) | ((a_hap == 1) & (b_hap == 0))
            is_homo_alt = (a_hap == 1) & (b_hap == 1)
            
            # get genomic location of SNPs
            i = snp_files.snp_index_h5.get_node("/%s" % chrom.name)[:]
            chrom_idx = np.where(i != SNP_UNDEF)[0]
            snp_idx = i[chrom_idx]

            # add to total genotype counts
            node = self.ref_count_h5.get_node("/%s" % chrom.name)
            node[chrom_idx] += is_homo_ref[snp_idx]
            
            node = self.het_count_h5.get_node("/%s" % chrom.name)
            node[chrom_idx] += is_het[snp_idx]
            
            node = self.alt_count_h5.get_node("/%s" % chrom.name)
            node[chrom_idx] += is_homo_alt[snp_idx]
            

    def close(self):
        for h5f in self.h5_files:
            h5f.close()

        # remove tempory combined files
        for filename in self.filenames:
            os.unlink(filename)
                         



class CountFiles(object):    
    def __init__(self, read_count_dir, individual):
        # open read count tracks for a single individual
        self.ref_as_count_h5 = tables.open_file("%s/ref_as_counts.%s.h5" % 
                                            (read_count_dir, individual), "r")
        self.alt_as_count_h5 = tables.open_file("%s/alt_as_counts.%s.h5" % 
                                            (read_count_dir, individual), "r")
        self.read_count_h5 = tables.open_file("%s/read_counts.%s.h5" %
                                             (read_count_dir, individual), "r")

        
    def close(self):
        """closes all of the data files"""
        self.ref_as_count_h5.close()
        self.alt_as_count_h5.close()
        self.read_count_h5.close()
        




        
def parse_args():
    parser = argparse.ArgumentParser(description="Uses read count and SNP information to generate "
                                     "a list of target regions and "
                                     "test SNPs that match specified "
                                     "criteria. The target regions"
                                     "are output to a BED-like file "
                                     "that can be used as input "
                                     "to the extract_haplotype_read_counts.py "
                                     "script.")
    
    parser.add_argument('--target_region_size', 
                        help='size of the target regions to output',
                        type=int, default=None)

    parser.add_argument('--min_as_count',
                        default=10, type=int,
                        help="minimum number of allele-specific reads "
                        "in target region (summed across individuals)")

    parser.add_argument('--min_read_count',
                        default=100, type=int,
                        help="minimum number of mapped reads in target region"
                        " (summed across individuals)")

    parser.add_argument('--min_het_count',
                        default=1, type=int,
                        help="minimum number of individuals that are heterozygous for "
                        "test SNP")

    parser.add_argument('--min_minor_allele_count',
                        default=1, type=int,
                        help='minimum number of test SNP minor alleles '
                        'present in individuals')
    
    parser.add_argument("--chrom",
                        help="Path to chromInfo.txt file (may be gzipped) "
                        "with list of chromosomes for the relevant genome "
                        "assembly. Each line in file should contain "
                        "tab-separated chromosome name and chromosome length "
                        "(in basepairs). chromInfo.txt files can be "
                        "downloaded from the UCSC genome browser. For "
                        "example, a chromInfo.txt.gz file for hg19 can "
                        "be downloaded from "
                        "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/",
                        metavar="CHROM_TXT_FILE",
                        required=True)
        
    parser.add_argument("--samples",
                        help="Path to text file containing identifiers for the "
                        "complete set of individuals with genotyping data. The "
                        "ordering of identifiers must be consistent with the "
                        "column order in the haplotype file. "
                        "The samples file is assumed to have one identifier "
                        "per line in the first column (other columns are "
                        "ignored).",
                        required=True,
                        metavar="SAMPLES_TXT_FILE",
                        default=None)
    
    parser.add_argument("--individuals",
                        help="A list of identifiers for the individuals that "
                        "will be used by the combined haplotype test. Read "
                        "counts for these individuals are read from HDF5 files "
                        "in the --read_count_dir. The argument can be a comma-"
                        "delimited list of identifiers (provided on command "
                        "line) or a path to a file containing a single "
                        "identifier per line. The individuals must be a subset "
                        "of those provided with the --samples argument.", 
                        metavar="INDIVIDUAL", required=True, default=None)
    
    parser.add_argument("--snp_tab",
                        help="Path to HDF5 file to read SNP information "
                        "from. Each row of SNP table contains SNP name "
                        "(rs_id), position, allele1, allele2.",
                        metavar="SNP_TABLE_H5_FILE",
                        required=True)
    
    parser.add_argument("--snp_index",
                        help="Path to HDF5 file containing SNP index. The "
                        "SNP index is used to convert the genomic position "
                        "of a SNP to its corresponding row in the haplotype "
                        "and snp_tab HDF5 files.",
                        metavar="SNP_INDEX_H5_FILE",
                        required=True)
    
    parser.add_argument("--haplotype",
                        help="Path to HDF5 file to read phased haplotypes "
                        "from.",
                        metavar="HAPLOTYPE_H5_FILE",
                        required=True,
                        default=None)

    parser.add_argument("--read_count_dir",
                        help="Path to directory containing HDF5 files with read counts "
                        "for each individual. These files are written by the bam2h5.py "
                        "program. The files must be named like: "
                        "ref_as_counts.<INDIVIDUAL>.h5, alt_as_counts.<INDIVIDUAL>.h5, "
                        "other_as_counts.<INDIVIDUAL>.h5, read_counts.<INDIVIDUAL>.h5",
                        metavar="READ_COUNT_DIR", required=True, default=None)

    parser.add_argument("--output_file",
                        help="Path to output file. If not specified output is written to stdout.",
                        default=None)

    args = parser.parse_args()

    return args




def get_samples_index(options):
    """Gets dictionary of sample_id => index mappings that is used 
    to lookup information in the genotype and haplotype tables"""
    sys.stderr.write("reading list of individuals from %s\n" % 
                     options.samples)
    
    f = open(options.samples, "rt")

    ind_dict = {}
    
    idx = 0
    for line in f:
        words = line.rstrip().split()

        # name = words[0].replace("NA", "")
        name = words[0]

        if name in ind_dict:
            raise ValueError("sample identifier '%s' appears multiple "
                             "times in file %s" % (name, options.samples))
        
        ind_dict[name] = idx
                
        idx += 1

    return ind_dict



def read_individuals(options, samp_idx):
    if os.path.exists(options.individuals):
        # read individuals from specified file
        ind_list = []
        f = open(options.individuals, "rt")

        for line in f:
            words = line.rstrip().split()
            identifier = words[0]
            
            ind_list.append(identifier)
    else:
        # no file specified, assume ','-delimited file
        # provided on command line
        ind_list = options.individuals.split(",")

    if len(ind_list) == 0:
        raise ValueError("could not read any individuals from '%s'" %
                         options.individuals)

    # check that all individuals were present in samples file
    for i in ind_list:
        if not i in samp_idx:
            raise ValueError("individual %s is not in samples file '%s'" % 
                             (i, options.samples))

    return ind_list



def main():

    sys.stderr.write("cmd: %s\n" % " ".join(sys.argv))
    
    args = parse_args()

    out_f = None
    if args.output_file:
        if args.output_file.endswith(".gz"):
            out_f = gzip.open(args.output_file, "wt")
        else:
            out_f = open(args.output_file, "wt")
    else:
        out_f = sys.stdout

    
    # make dictionary of identifier => index mapping
    samp_idx = get_samples_index(args)

    # read individuals
    individuals = read_individuals(args, samp_idx)
    
    chrom_list = chromosome.get_all_chromosomes(args.chrom)
    chrom_dict = chromosome.get_chromosome_dict(args.chrom)


    combined_files = CombinedFiles(OUTPUT_DIR, chrom_list)
    snp_files = SNPFiles(args)

    # STEP 1: make combined HDF5 files of AS counts, 
    # total mapped read counts, and genotype counts
    sys.stderr.write("summing genotypes and read counts across individuals\n")
    for ind in individuals:
        # open count files for this indivudal
        sys.stderr.write("individual: %s\n" % ind)
        count_files = CountFiles(args.read_count_dir, ind)

        ind_idx = samp_idx[ind]
        
        # add counts to combined totals
        combined_files.add_counts(chrom_list, count_files, snp_files, ind_idx)

        count_files.close()
        

    sys.stderr.write("generating list of target regions\n")
    
    # STEP 2: generate list of target regions centered on test SNPs:
    write_target_regions(out_f, args, chrom_list, combined_files, snp_files)

    combined_files.close()
    snp_files.close()




def write_target_regions(out_f, args, chrom_list, combined_files, snp_files):
    """Choose test SNPs and corresponding target regions and write to output file.
    Selected test SNPs:
       (a) are heterozygous in at least min_het_count individuals
       (b) have a minor allele count of at least min_minor_allele_count

    Selected target regions are centered on test SNPs and
         (a) have at least min_as_count allele-specific reads (across all individuals)
         (b) have at least min_read_count mapped reads (across all individuals)
    """

    for chrom in chrom_list:        
        node_name = "/%s" % chrom.name
        if node_name not in snp_files.snp_index_h5:
            continue
        if node_name not in snp_files.snp_tab_h5:
            continue
        
        sys.stderr.write("  %s\n" % chrom.name)

        sys.stderr.write("  getting genotype counts\n")
        ref_geno_count = combined_files.ref_count_h5.get_node(node_name)[:]
        het_geno_count = combined_files.het_count_h5.get_node(node_name)[:]

        ref_allele_count = ref_geno_count * 2 + het_geno_count
        # free memory as it is no longer needed
        del ref_geno_count

        alt_geno_count = combined_files.alt_count_h5.get_node(node_name)[:]
        alt_allele_count = alt_geno_count * 2 + het_geno_count

        del alt_geno_count

        sys.stderr.write("  getting minor allele counts\n")

        minor_count = np.amin(np.vstack([ref_allele_count, alt_allele_count]),
                              axis=0)
        
        idx = np.where((minor_count >= args.min_minor_allele_count) &
                       (het_geno_count >= args.min_het_count))[0]

        del het_geno_count
        del minor_count
        
        sys.stderr.write("  %d possible test SNPs\n" % idx.shape[0])

        read_counts = combined_files.read_count_h5.get_node(node_name)[:]
        as_read_counts = combined_files.as_count_h5.get_node(node_name)[:]

        snp_idx = snp_files.snp_index_h5.get_node(node_name)
        snp_tab = snp_files.snp_tab_h5.get_node(node_name)
        
        n_region = 0
                            
        for i in idx:
            start = int(max(1, i+1 - args.target_region_size/2))
            end = int(min(chrom.length, i+1 + args.target_region_size/2))

            n_reads = np.sum(read_counts[start-1:end])
            n_as_reads = np.sum(as_read_counts[start-1:end])

            snp_row = snp_tab[snp_idx[i]]
            
            if (n_reads >= args.min_read_count) and (n_as_reads >= args.min_as_count):
                # keep this target region

                # NOTE: currently this filter just uses total count of AS reads in region.
                # Would be better to take into account genotypes of each individual, since AS reads
                # are only useful for test in individuals that are heterozygous for test SNP
                out_f.write("%s %d %d %s %s + %s.%d %d %d\n" % 
                            (chrom.name, i+1, i+2, snp_row['allele1'],
                             snp_row['allele2'], chrom.name, start+1,
                             start, end))

                n_region += 1

        sys.stderr.write("  wrote %d test SNP / target region pairs\n" % n_region)





    
main()
