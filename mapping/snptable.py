import sys
import numpy as np
import gzip
import pysam
import operator

import util


NUCLEOTIDES = {b'A', b'C', b'T', b'G'}
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
        self.phase = None
        self.n_snp = 0
        self.samples = []
        


    def read_h5(self, snp_tab_h5, snp_index_h5, hap_h5, chrom_name,
                samples=None):
        """read in SNPs and indels from HDF5 input files"""

        node_name = "/%s" % chrom_name
        phase_node_name = "/phase_%s" % chrom_name
        
        if node_name not in snp_tab_h5:
            # try to remove chr prefix
            if chrom_name.startswith("chr"):
                node_name = "/%s" % chrom_name[3:]
                phase_node_name = "/phase_%s" % chrom_name[3:]
            # try to add chr prefix
            else:
                node_name = "/chr%s" % chrom_name
                phase_node_name = "/phase_chr%s" % chrom_name
            # if it still doesn't work, ignore this chromosome
            if node_name not in snp_tab_h5:
                sys.stderr.write("WARNING: chromosome %s is not "
                                 "in snp_tab.h5 file, assuming no SNPs "
                                 "for this chromosome\n" % chrom_name)
                self.clear()
                return

        # get numpy array of SNP idices
        node = snp_index_h5.get_node(node_name)
        self.snp_index = node[:]

        # get numpy array of SNP positions
        node = snp_tab_h5.get_node(node_name)
        self.snp_pos = node[:]['pos']
        self.snp_allele1 = node[:]['allele1']
        self.snp_allele2 = node[:]['allele2']
        self.n_snp = self.snp_pos.shape[0]
        self.samples = self.get_h5_samples(hap_h5, chrom_name)
        self.haplotypes = hap_h5.get_node(node_name)
        if phase_node_name in hap_h5:
            self.phase = hap_h5.get_node(phase_node_name)
        
        if samples:
            # reduce set of SNPs and indels to ones that are
            # polymorphic in provided list of samples
            samp_idx_dict, samp_idx = self.get_h5_sample_indices(hap_h5, chrom_name, samples)

            if len(samp_idx) == 0:
                # gracefully handle situation where there are no matching
                # samples on this chromosome
                sys.stderr.write("WARNING: chromosome %s haplotype file "
                                 "has no samples that match provided "
                                 "sample names, assuming no SNPs for "
                                 "this chromosome\n" % chrom_name)

                self.clear()
                return
                
            hap_idx = np.empty(samp_idx.shape[0]*2, dtype=np.int)
            hap_idx[0::2] = samp_idx*2
            hap_idx[1::2] = samp_idx*2 + 1
            haps = self.haplotypes[:, hap_idx]
            if self.phase:
                phase = self.phase[:, samp_idx]

            # count number of ref and non-ref alleles,
            # ignoring undefined (-1s)
            nonref_count = np.apply_along_axis(np.sum, 1, haps == 1)
            ref_count = np.apply_along_axis(np.sum, 1, haps == 0)
            total_count = nonref_count + ref_count
            is_polymorphic = (ref_count > 0) & (ref_count < total_count)

            # reduce to set of polymorphic positions
            sys.stderr.write("reducing %d SNPs on chromosome "
                             "%s to %d positions that are polymorphic in "
                             "sample of %d individuals\n" %
                             (haps.shape[0], chrom_name, 
                              np.sum(is_polymorphic), len(samples)))

            # make filtered and ordered samples for this chromosome
            # that corresponds to order of haplotypes
            sorted_samps = sorted(list(samp_idx_dict.items()),
                                  key=operator.itemgetter(1))
            self.samples = [x[0] for x in sorted_samps]
            
            self.haplotypes = haps[is_polymorphic,]
            if self.phase:
                self.phase = phase[is_polymorphic,]
            self.snp_pos = self.snp_pos[is_polymorphic]
            self.snp_allele1 = self.snp_allele1[is_polymorphic]
            self.snp_allele2 = self.snp_allele2[is_polymorphic]
            self.n_snp = self.snp_pos.shape[0]

            # regenerate index to point to reduced set of polymorphic SNPs
            self.snp_index[:] = -1                
            self.snp_index[self.snp_pos-1] = np.arange(self.n_snp,
                                                       dtype=np.int32)
                

    
    def get_h5_samples(self, h5f, chrom_name):
        """Reads list of samples that are present in 'samples' table 
        from haplotype HDF5 file"""
        samples = None

        node_name = "/samples_%s" % chrom_name
        
        if node_name not in h5f:
            # try to remove chr prefix
            if chrom_name.startswith("chr"):
                node_name = "/samples_%s" % chrom_name[3:]
            # try to add chr prefix
            else:
                node_name = "/samples_chr%s" % chrom_name
            # if it still doesn't work, raise an error
            if node_name not in h5f:
                raise ValueError("Cannot retrieve haplotypes for "
                                 "specified samples, because haplotype "
                                 "file %s does not contain '%s' table. "
                                 "May need to regenerate haplotype HDF5 file "
                                 "using snp2h5 with --samples option" %
                                 (h5f.filename, node_name))
                return samples

        node = h5f.get_node(node_name)
        samples = [row["name"].decode("utf-8") for row in node]

        # sys.stderr.write("SAMPLES: %s\n" % samples)
        
        return samples

    
    
    def get_h5_sample_indices(self, hap_h5, chrom_name, samples):
        """returns the indices of the the specified samples in the 
        HDF5 haplotype file. Indices are returned in a dictionary
        keyed on sample and as an array. Samples that are not 
        found in the haplotype HDF5 file for the specified chromosome 
        are not included in the dict or the array."""
        hap_samples = self.get_h5_samples(hap_h5, chrom_name)
        not_seen_samples = set(samples)
        seen_samples = set([])
        samp_idx = []
        samp_idx_dict = {}
        
        # get haplotype table indices of samples
        for i in range(len(hap_samples)):
            if hap_samples[i] in seen_samples:
                sys.stderr.write("WARNING: sample %s is present multiple "
                                 "times in haplotype table\n" % hap_samples[i])
            elif hap_samples[i] in not_seen_samples:
                # record index of this sample, add to set of samples
                # we have already observed
                samp_idx.append(i)
                samp_idx_dict[hap_samples[i]] = i
                not_seen_samples.remove(hap_samples[i])
                seen_samples.add(hap_samples[i])
            else:
                # this haplotype sample not in requested list
                pass

        if len(not_seen_samples) > 0:
            sys.stderr.write("WARNING: the following samples are not "
                             "present in haplotype table for chromosome "
                             "%s: %s\n" %
                             (chrom_name, ",".join(not_seen_samples)))
        
        return samp_idx_dict, np.array(samp_idx, dtype=np.int)

        

    def is_snp(self, allele1, allele2):
        """returns True if alleles appear to be 
        single-nucleotide polymorphism, returns false
        if appears to be an indel"""

        if (len(allele1) == 1) and (len(allele2) == 1):
            if allele1 in NUCLEOTIDES and allele2 in NUCLEOTIDES:
                # this is a SNP
                return True
            else:
                if (b"-" in allele1) or (b"-" in allele2):
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
            if util.is_gzipped(filename):
                f = gzip.open(filename, "rt")
            else:
                f = open(filename, "rt")
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

                # This is like insert, but at end of the read.
                # Sequence was not considered in alignment.
                # Usually this is because bases at end of read
                # were low quality. One option would be to 
                # pretend soft-clipped part of read was aligned
                # like match/mismatch and to consider SNPs in this
                # region. We have decided to not consider SNPs 
                # because this part of read is not actually aligned.

            elif op == BAM_CHARD_CLIP:
                # these bases not included in read or genome
                pass

            elif op == BAM_CPAD:
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
