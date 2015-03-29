
import sys

import tables

import util


# VCF format from: http://www.1000genomes.org/node/101
# Fixed VCF fields are in this order:
# 0 - CHROM - name of chromosome
# 1 - POS - chromosome position (1st base numbered 1)
# 2 - ID - SNP identifier such as rsid 
# 3 - REF - reference allele
# 4 - ALT - alternate allele
# 5 - QUAL - phred-scaled prob that site is variable?
# 6 - FILTER - PASS or some other value to indicate site has not
#              not passed quality filters
# 7 - INFO - semicolon-separated series of short keys with 
#            optional values in the format: <key>=<data>[,data]
# 8 - FORMAT - colon-delimited format of following genotype columns
#              GT - genotype
#              DP - read depth
#              GL - genotype log10 likelihoods
#              GQ - genotype quality
#              HQ - haplotype quality    


MAX_ALLELE = 1024
MAX_INFO = 256

class SNPDesc(tables.IsDescription):
    name = tables.StringCol(16)
    pos = tables.Int32Col()

    # this will truncate the longest variants...
    # potentially want to have separate table for SNPs and indels
    # to save space, and allow for longer indels
    allele1 = tables.StringCol(MAX_ALLELE)
    allele2 = tables.StringCol(MAX_ALLELE)



def guess_chromosome_name(chrom_tab, filename):
    """Given table of chromosome names and a filename,
    attempts to guess what chromosome the file is for"""
    longest_match = ""
    
    for chrom in chrom_tab.get_all_chromosomes():
        # find longest chromosome name that is present in
        # filename. Use longest name because otherwise 'chr1' will
        # match filenames like 'chr10.vcf'
        if len(chrom.name) > len(longest_match):
            if chrom.name in filename:
                longest_match = chrom.name

    if longest_match == "":
        raise ValueError("could not determine name of "
                         "chromosome from file %s" % chrom.name)

    return longest_match



def read_vcf_header(f):
    """Skips over header lines for VCF"""
    for line in f:
        if line.startswith("##"):
            pass
        elif line.startswith("#CHROM"):
            words = line.rstrip()[1:].split()

            headers = words[0:9]
            sample_ids = words[9:]

            break
        else:
            raise ValueError("expected last line in VCF header to "
                             "start with #CHROM")

    return headers, sample_ids



def import_samples(vcf_filename, output_filename):
    """Extracts sample names from VCF header, and writes them to a
    text file."""
    f = util.check_open(vcf_filename)
    header, sample_ids = read_vcf_header(f)

    out_f = util.check_open(output_filename, mode="w")
    for samp in sample_ids:
        out_f.write("%s\n" % samp)

    out_f.close()
    f.close()

    return


def import_snps(chrom_tab, vcf_filenames, output_filename):
    """Creates an HDF5 file containing a single table per chromosome.
    Each table row contains information about a single SNP including
    identifier, position, reference allele, alternate allele.
    """
    # create an HDF5 file
    h5f = tables.openFile(output_filename, "w")

    chrom_dict = chrom_tab.get_chromosome_dict()
    
    for vcf_file in vcf_filenames:

        # guess name of chromosome from filename
        chrom_name = guess_chromosome_name(chrom_tab, vcf_file)
        chrom = chrom_dict[chrom_name]

        # open file, read header lines
        f = util.check_open(vcf_file)
        header, sample_ids = read_vcf_header(f)
        
        sys.stderr.write("importing SNPs for chromosome %s\n" % chrom.name)
        
        # create table to hold SNPs for this chromosome
        snp_tab = h5f.createTable("/", chrom_name, SNPDesc, "SNPs")
    
        row = snp_tab.row
        
        for line in f:
            # just care about first 9 fixed fields
            max_split = 9
            words = line.split(None, max_split)

            if len(words) < 9:
                raise ValueError("Expected at least 9 columns in "
                                 "each VCF row")

            pos = int(words[1])

            if pos > chrom.length or pos < 1:
                raise ValueError("position (%d) for SNP %s "
                                 "is outside of chromosome range 1-%d"
                                 % pos, words[2], chrom.length)
            
            row['pos'] = pos
            row['name'] = words[2]
            row['allele1'] = words[3]
            row['allele2'] = words[4]

            row.append()

        snp_tab.flush()
        f.close()
    
    h5f.close()
    

def check_samples(samples_filename, vcf_filename, sample_ids, vcf_sample_ids):
    """Checks that sample identifiers from sample file and VCF file match"""
    
    # make sure sample identifiers from headers in provided VCF
    # match those in existing samples file
    i = 0 
    if len(sample_ids) != len(vcf_sample_ids):
        raise ValueError("samples file (%s) has %d samples, but VCF file "
                         "(%s) has %d samples" % 
                         (samples_filename, len(sample_ids), 
                          vcf_file, len(vcf_sample_ids)))

    for s1, s2 in zip(sample_ids, vcf_sample_ids):
        i += 1            
        if s1 != s2:
            raise ValueError("sample identifier mismatch between samples "
                             "file (%s) and VCF file (%s): "
                             "[%d]: '%s' != '%s'" % 
                             (samples_filename, vcf_file, i, s1, s2))

    
    


def read_geno_probs(n_snp, n_samp, f):
    sys.stderr.write("allocating memory to hold genotype probabilities\n")
    cols_per_sample = 3

    data_matrix = np.empty((n_snp, n_samp * cols_per_sample),
                           dtype=np.float32)

    sys.stderr.write("reading genotype probabilities\n")
    warn_gl = False
    snp_num = 0

    for line in f:
        if snp_num >= n_snp:
            raise ValueError("VCF file contains more rows (%d) than SNP "
                             "file" % n_snp)

        words = line.rstrip().split()

        if len(words) < 9:
            raise ValueError("Expected at least 9 columns in each VCF row")

        pos = int(words[1])
        snp_id = words[2]
        ref_allele = words[3]
        alt_allele = words[4]
        gt_format = words[8].split(":")

        gt_fields = words[9:]

        if len(gt_fields) != n_samp:
            raise ValueError("number of genotype columns (%d) "
                             "does not match number of samples (%d)"
                             % (len(gt_fields), n_samp))

        # find index of GT and GL fields from ":"-delimited format string
        gt_idx = None
        gl_idx = None
        for i in range(len(gt_format)):
            if gt_format[i] == "GT":
                # genotype field
                gt_idx = i
            elif gl_format[i] == "GL":
                # genotype-likelihood field
                gl_idx = i

        if gl_idx is None and not warn_gq:
            sys.stderr.write("WARNING: Genotype likelihood (GL) "
                             "field is not present.\nUsing GT field "
                             "Assuming genotype probabilities of 0.99\n")
            warn_gl = True

        if gl_idx is None and gt_idx is None:
            raise ValueError("either GT or GL field is required and "
                             "neither are defined in format column")

        all_probs = []
        if gl_idx:
            # parse out genotype log-likelihoods, convert to probs
            all_probs = []
            for col in gt_fields:
                # should be 3 likelihoods per individual
                probs = [10.0**float(x) 
                         for x in col.split(":")[gl_idx].split(",")]

                if len(probs) != 3:
                    raise ValueError("expected 3 genotype likelihoods "
                                     "per sample")

                sum_probs = sum(probs)
                if sum_probs < 0.99 or sum_probs > 1.01:
                    raise ValueError("expected genotype likelihoods to "
                                     "sum to 1.0 (%s)" % 
                                     ",".join([str(x) for x in probs]))

                all_probs.extend(probs)
        else:
            # use genotypes, convert from 1.0 to 0.99 -- not ideal!
            for col in gt_fields:
                gt_val = col.split(":")[gt_idx]

                if gt_val == "0/0" or gt_val == "0|0":
                    # homozygous reference
                    all_probs.extend([0.99, 0.005, 0.005])
                elif(gt_val == "0/1" or gt_val == "0|1" or
                     gt_val == "1/0" or gt_val == "1|1"):
                    # heterozygous
                    all_probs.extend([0.005, 0.99, 0.005])
                elif gt_val == "1/1" or gt_val == "1|1":
                    # homozygous alternate
                    all_probs.extend([0.005, 0.005, 0.99])

        if len(all_probs) != n_samples * cols_per_sample:
            raise ValueError("expected %d likelihoods for %d samples "
                             "but got %d" % ((n_samples*cols_per_sample), 
                                             len(all_probs)))

        data_matrix[snp_num, :] = all_probs

        snp_num += 1

    if snp_num != n_snps:
        raise ValueError("number of SNP rows in VCF does not match "
                         "number of SNPs in %d"
                         (n_snp, snp_filename))

    return data_matrix




def write_carray(h5f, chrom, data_matrix):
    if data_matrix.dtype == "float32":
        atom = tables.Float32Atom(dflt=np.nan)
    elif data_matrix.dtype == "int8":
        atom = tables.Int8Atom(dflt=-1)
    else:
        raise ValueError("datatype not implemented %s" % 
                         str(data_matrix.dtype))

    zlib_filter = tables.Filters(complevel=1, complib="zlib")

    shape = data_matrix.shape
    carray = track.h5f.createCArray(h5f.root, chrom.name,
                                    atom, shape, filters=zlib_filter)

    carray[:,:] = data_matrix

    

def import_genos(chrom_tab, samples_filename, snp_filename, 
                 vcf_filenames, output_filename):
    """Creates an HDF5 file containing a single matrix of genotype 
    probabilities for each chromosome. Each row in the matrix corresponds 
    to a SNP in the snp_tab.  There are 3M columns, where M is number of 
    samples. The samples must be ordered in the same way as in the 
    samples.txt file. 

    Columns 1-3 are for the first sample, columns 4-6 are for the 
    second sample, etc.  The columns for each individual are 
    "homozygous reference probability", "heterozygous probability", 
    and "homozygous non-reference probability" and should sum to 1.0.
    """
    snp_h5 = tables.openFile(snp_filename)
    geno_h5 = tables.openFile(output_filename, "w")

    # read sample identifiers
    samp_f = open(samples_filename)
    sample_ids = [l.strip() for l in samp_f]
    samp_f.close()

    for vcf_file in vcf_filenames:
        f = util.check_open(vcf_filename)

        header, vcf_sample_ids = read_vcf_header(f)
        check_samples(samples_filename, vcf_file, sample_ids, vcf_sample_ids)

        # guess name of chromosome from filename
        chrom_name = guess_chromosome_name(chrom_tab, vcf_file)
        chrom = chrom_dict[chrom_name]

        # get SNP table for this chromosome
        if chrom_name not in snp_h5.root:
            raise ValueError("There is no SNP data for chromosome %s in "
                             "file %s" % (chrom_name, snp_filename))
        
        snp_tab = snp_h5.getNode("/%s" % chrom_name)
        n_snp = snp_tab.shape[0]
        n_samp = len(sample_ids)

        geno_probs = read_geno_probs(f, n_snp, n_samp)

        write_carray(geno_h5, chrom, geno_probs)
        
        f.close()

    snp_h5.close()
    geno_h5.close()
    
