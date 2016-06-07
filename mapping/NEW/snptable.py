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
        self.snp_allele1 = np.array([], dtype="|S1")
        self.snp_allele2 = np.array([], dtype="|S1")
        self.indel_index = np.array([], dtype=np.int32)
        self.indel_pos = np.array([], dtype=np.int32)
        self.indel_allele1 = []
        self.indel_allele2 = []
        self.n_snp = 0
        
    
    def read_file(self, filename):
        """read in SNPs and indels from text input file"""
        sys.stderr.write("reading SNPs from file '%s'\n" % filename)

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
        indel_pos_list = []
        indel_allele1_list = []
        indel_allele2_list = []

        for line in f:
            words = line.split()

            if(len(words) < 3):
                raise ValueError("expected at least 3 values per SNP "
              			 "file line but got %d:\n"
                                 "%s\n" % (len(words), line))

            pos = int(words[0])
            a1 = words[1].upper()
            a2 = words[2].upper()

            if pos <= 0:
                raise ValueError("expected SNP position to be >= 1:\n%s\n" %
                                 line)

            if pos > max_pos:
                max_pos = pos

            if (len(a1) == 1) and (len(a2) == 1):
                if a1 in NUCLEOTIDES and a2 in NUCLEOTIDES:
                    # this is a SNP
                    snp_pos_list.append(pos)
                    snp_allele1_list.append(a1)
                    snp_allele2_list.append(a2)
                else:
                    if ("-" in a1) or ("-" in a2):
                        # 1bp indel
                        # reads overlapping indels are thrown out, although
                        # we will likely handle indels better soon
                        indel_pos_list.append(pos)
                        indel_allele1_list.append(a1.replace("-", ""))
                        indel_allele2_list.append(a2.replace("-", ""))
                    else:
                        sys.stderr.write("WARNING: unexpected character "
                                         "in SNP:\n%s\n" % line)
            else:
                # this is an indel
                indel_pos_list.append(pos)
                indel_allele1_list.append(a1)
                indel_allele2_list.append(a2)

        f.close()

        # convert lists to numpy arrays, which allow for faster
        # lookups and use less memory
        self.snp_pos = np.array(snp_pos_list, dtype=np.int32)
        del snp_pos_list
        self.snp_allele1 = np.array(snp_allele1_list, dtype="|S1")
        del snp_allele1_list
        self.snp_allele2 = np.array(snp_allele2_list, dtype="|S1")
        del snp_allele2_list

        self.indel_pos = np.array(indel_pos_list, dtype=np.int32)
        del indel_pos_list

        # make another array that makes it easy to lookup SNPs by their position
        # on the chromosome
        self.snp_index = np.empty(max_pos, dtype=np.int32)
        self.snp_index[:] = SNP_UNDEF
        self.snp_index[self.snp_pos-1] = np.arange(self.snp_pos.shape[0])

        self.indel_index = np.empty(max_pos, dtype=np.int32)
        self.indel_index[:] = SNP_UNDEF
        self.indel_index[self.indel_pos-1] = np.arange(self.indel_pos.shape[0])
        self.indel_allele1 = indel_allele1_list
        self.indel_allele2 = indel_allele2_list

        self.n_snp = self.snp_pos.shape[0]

    
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

        # index into SNP table for overlapping SNPs
        snp_idx = []
        # positions in read of overlapping SNPs
        snp_read_pos = []
        # index into indel table for overlapping indels
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
                genome_start = genome_start + 1
                genome_end = genome_start + op_len - 1

                # check for SNP in this genome segment
                s = genome_start - 1
                e = min(genome_end, self.snp_index.shape[0])
                s_idx = self.snp_index[s:e]
                offsets = np.where(s_idx != SNP_UNDEF)[0]
                
                if offsets.shape[0] > 0:
                    # there are overlapping SNPs
                    snp_idx.extend(s_idx[offsets])
                    # get the offset of the SNPs into the read
                    read_pos = offsets + read_start
                    snp_read_pos.extend(read_pos)

                # check for INDEL in this genome segment
                i_idx = self.indel_index[s:e]
                offsets = np.where(i_idx != SNP_UNDEF)[0]
                if offsets.shape[0] > 0:
                    indel_idx.extend(i_idx[offsets])
                    read_pos = offsets + read_start
                    indel_read_pos.extend(read_pos)

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
                genome_start = genome_start + 1
                genome_end   = genome_start + op_len - 1

                # Read sequence does not advance, no possibility
                # for read to overlap SNP, since these bases do
                # not exist in read

                # in most cases deletion should be picked up
                # by flanking match segment, but there could be
                # nested indels

                s = genome_start - 1
                e = min(genome_end, self.indel_index.shape[0])
                
                # check for INDEL in this genome segment
                i_idx = self.indel_index[s:e]
                offsets = np.where(i_idx != SNP_UNDEF)[0]
                if offsets.shape[0] > 0:
                    indel_idx.extend(i_idx[offsets])
                    # position in read is where we last left off
                    # in read sequence
                    indel_read_pos.extend(read_end)

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
