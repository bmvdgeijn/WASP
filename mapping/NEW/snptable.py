import sys
import numpy as np
import gzip
import pysam


NUCLEOTIDES = set(['A', 'C', 'T', 'G'])
SNP_UNDEF = -1


# codes for CIGAR string
BAM_CMATCH     = 0   # M - match/mismatch to ref M
BAM_CINS       = 1   # I - insertion in read relative to ref
BAM_CDEL       = 2   # D - deletion in read relative to ref
BAM_CREF_SKIP  = 3   # N - skipped region from reference (e.g. intron)
BAM_CSOFT_CLIP = 4   # S - soft clipping (clipped sequence present in seq)
BAM_CHARD_CLIP = 5   # H - hard clipping (clipped sequence NOT present in seq)
BAM_CPAD       = 6   # P - padding (silent deletion from padded reference)
BAM_CEQUAL     = 7   # = - sequence match
BAM_CDIFF      = 8   # X - sequence mismatch

class SNPTable(object):
    def __init__(self):
        self.clear()

    def clear(self):
        # snp_index and indel_index are arrays of length
        # max(snp_pos, indel_pos) that provide lookup
        # into snp_pos, snp_allele1, etc. by chromosome position.
        # For example, if the first and second snps on the chromosome are
        # at positions 1234, 1455 then elements 1233 and 1444 of the
        # snp_index array will be 0 and 1 (and can be used to lookup
        # info for the SNP in snp_pos, snp_allele1, snp_allele2 arrays)
        self.snp_index = np.array([], dtype=np.int32)
        self.snp_pos = np.array([], dtype=np.int32)
        self.snp_allele1 = np.array([], dtype="|S10")
        self.snp_allele2 = np.array([], dtype="|S10")
        self.haplotypes = None
        self.n_snp = 0
        


    def read_h5(self, snp_tab_h5, snp_index_h5, hap_h5, chrom_name):
        """read in SNPs and indels from HDF5 input files"""

        node_name = "/%s" % chrom_name
        
        if node_name not in snp_tab_h5:
            sys.stderr.write("WARNING: chromosome %s is not "
                             "in snp_tab.h5 file, assuming no SNPs "
                             "for this chromosome\n" % chrom_name)
            self.clear()
            return

        # get numpy array of SNP idices
        node = snp_index_h5.getNode(node_name)
        self.snp_index = node[:]

        # get numpy array of SNP positions
        node = snp_tab_h5.getNode(node_name)
        self.snp_pos = node[:]['pos']
        self.snp_allele1 = node[:]['allele1']
        self.snp_allele2 = node[:]['allele2']
        self.n_snp = self.snp_pos.shape[0]

        self.haplotypes = hap_h5.getNode(node_name)
        
        
        

    def is_snp(self, allele1, allele2):
        """returns True if alleles appear to be 
        single-nucleotide polymorphism, returns false
        if appears to be an indel"""

        if (len(allele1) == 1) and (len(allele2) == 1):
            if allele1 in NUCLEOTIDES and allele2 in NUCLEOTIDES:
                # this is a SNP
                return True
            else:
                if ("-" in allele1) or ("-" in allele2):
                    # 1bp indel
                    return False
                else:
                    sys.stderr.write("WARNING: unexpected character "
                                     "in SNP alleles:\n%s/%s\n" %
                                     (allele1, allele2))
                    return False                
        
        return False
        


        
    def read_file(self, filename):
        """read in SNPs and indels from text input file"""
        try:
            if filename.endswith(".gz"):
                f = gzip.open(filename)
            else:
                f = open(filename, "r")
        except IOError:
            sys.stderr.write("WARNING: unable to read from file '%s', "
                             "assuming no SNPs for this chromosome\n" %
                             filename)
            self.clear()
            return
        
        snp_pos_list = []
        snp_allele1_list = []
        snp_allele2_list = []
        max_pos = 0

        for line in f:
            words = line.split()

            if(len(words) < 3):
                raise ValueError("expected at least 3 values per SNP "
              			 "file line but got %d:\n"
                                 "%s\n" % (len(words), line))

            pos = int(words[0])
            a1 = words[1].upper().replace("-", "")
            a2 = words[2].upper().replace("-", "")

            if pos <= 0:
                raise ValueError("expected SNP position to be >= 1:\n%s\n" %
                                 line)

            if pos > max_pos:
                max_pos = pos

            snp_pos_list.append(pos)
            snp_allele1_list.append(a1)
            snp_allele2_list.append(a2)

        f.close()

        # convert lists to numpy arrays, which allow for faster
        # lookups and use less memory
        self.snp_pos = np.array(snp_pos_list, dtype=np.int32)
        del snp_pos_list
        self.snp_allele1 = np.array(snp_allele1_list, dtype="|S10")
        del snp_allele1_list
        self.snp_allele2 = np.array(snp_allele2_list, dtype="|S10")
        del snp_allele2_list

        # make another array that makes it easy to lookup SNPs by their position
        # on the chromosome
        self.snp_index = np.empty(max_pos, dtype=np.int32)
        self.snp_index[:] = SNP_UNDEF
        self.snp_index[self.snp_pos-1] = np.arange(self.snp_pos.shape[0])

        self.n_snp = self.snp_pos.shape[0]

        # currently haplotypes can only be read from HDF5 file
        self.haplotypes = None

    
    def get_overlapping_snps(self, read):
        """Returns several lists: 
        [1] indices of SNPs that this read overlaps,
        [2] positions in read sequence that overlap SNPs, 
        [3] indices for indels that read overlaps, 
        [4] positions in read sequence that overlap indels. 
        First base of read is position 1."""
        
        # read.cigar is a list of tuples. Each tuple has two entries. The first
        # entry specifies the character in the cigar and the second entry
        # specifies the length of that character. The values are
        # M       BAM_CMATCH      0
        # I       BAM_CINS        1
        # D       BAM_CDEL        2
        # N       BAM_CREF_SKIP   3
        # S       BAM_CSOFT_CLIP  4
        # H       BAM_CHARD_CLIP  5
        # P       BAM_CPAD        6
        # =       BAM_CEQUAL      7
        # X       BAM_CDIFF       8
        # E.g. (0, 5) means 5 matches, and (4, 2) means a soft clip of 2bp
        read_start = 0
        read_end = 0
        genome_start = read.pos
        genome_end = read.pos

        # index into combined SNP/indel table for overlapping SNPs
        snp_idx = []
        # positions in read of overlapping SNPs
        snp_read_pos = []
        # index into combined SNP/indel table for overlapping indels
        indel_idx = []
        # positions in read of overlapping SNPs
        indel_read_pos = []
        
        for cigar in read.cigar:
            op = cigar[0] # CIGAR 'operation'
            op_len  = cigar[1] # length of operation
            
            if (op == BAM_CMATCH) or (op == BAM_CEQUAL) or (op == BAM_CDIFF):
                # match or mismatch to reference
                read_start = read_end + 1
                read_end = read_start + op_len - 1
                genome_start = genome_end + 1
                genome_end = genome_start + op_len - 1

                # check for SNP in this genome segment
                s = genome_start - 1
                e = min(genome_end, self.snp_index.shape[0])
                s_idx = self.snp_index[s:e]
                offsets = np.where(s_idx != SNP_UNDEF)[0]
                
                if offsets.shape[0] > 0:
                    # there are overlapping SNPs and/or indels
                    
                    for offset in offsets:
                        read_pos = offset + read_start
                        allele1 = self.snp_allele1[s_idx[offset]]
                        allele2 = self.snp_allele2[s_idx[offset]]
                        if self.is_snp(allele1, allele2):
                            snp_idx.append(s_idx[offset])
                            snp_read_pos.append(read_pos)
                        else:
                            indel_idx.append(s_idx[offset])
                            indel_read_pos.append(read_pos)

            elif op == BAM_CINS:
                # insert in read relative to reference
                read_start = read_end + 1
                read_end = read_start + op_len - 1

                # Genome sequence does not advance, no possibility
                # for read to overlap SNP, since these bases do
                # not exist in reference.
                # INDELs here should be picked up
                # by one of flanking match segments

            elif op == BAM_CDEL:
                # deletion in read relative to reference
                genome_start = genome_end + 1
                genome_end   = genome_start + op_len - 1

                # Read sequence does not advance, no possibility
                # for read to overlap SNP, since these bases do
                # not exist in read

                # in most cases deletion should be picked up
                # by flanking match segment, but there could be
                # nested indels

                s = genome_start - 1
                e = min(genome_end, self.snp_index.shape[0])
                
                # check for INDEL in this genome segment
                s_idx = self.snp_index[s:e]
                offsets = np.where(s_idx != SNP_UNDEF)[0]
                
                if offsets.shape[0] > 0:
                    # there are overlapping SNPs and/or indels
                    for offset in offsets:
                        read_pos = offset + read_start
                        allele1 = self.snp_allele1[s_idx[offset]]
                        allele2 = self.snp_allele2[s_idx[offset]]
                        if self.is_snp(allele1, allele2):
                            # ignore SNP
                            pass
                        else:
                            indel_idx.append(s_idx[offset])
                            # position in read is where we last left off
                            # in read sequence
                            indel_read_pos.append(read_end)
            elif op == BAM_CREF_SKIP:
                # section of skipped reference, such as intron
                genome_end = genome_end + op_len
                genome_start = genome_end

                # do nothing with SNPs/indels in this region
                # since they are skipped
                
            elif op == BAM_CSOFT_CLIP:
                # this part of read skipped
                read_start = read_end + 1
                read_end = read_start + op_len - 1

                # This is like insert, but at the beginning of the read.
                # TODO: handle indels? Sometimes a read can be softclipped
                # because it contains insert relative to reference, but
                # in these cases, presumably reference version of read
                # would map to same location (with higher score).

            elif seq_type == BAM_CHARD_CLIP:
                # these bases not included in read or genome
                pass

            elif seq_type == BAM_CPAD:
                # like an insert, likely only used in multiple-sequence
                # alignment where inserts may be of different lengths
                # in different seqs
                read_start += read_end + 1
                read_end = read_start + op_len - 1

            else:
                raise ValueError("unknown CIGAR code %d" % op)

        if read_end != len(read.seq):
            raise ValueError("length of read segments in CIGAR %d "
                             "does not add up to query length (%d)" %
                             (read_end, len(read.seq)))
        
        
        return snp_idx, snp_read_pos, indel_idx, indel_read_pos
