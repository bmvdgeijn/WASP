import argparse
import gzip
import sys

import array
import pysam

class SNP:
    """SNP objects hold data for a single SNP"""
    def __init__(self, snp_line):
        """
        Initialize SNP object.
        
        Parameters:
        -----------

        snp_line : str
            Line from SNP file.

        Attributes:
        -----------

        pos : int
            Genomic position of SNP.

        alleles : list
            List of two alleles.

        ptype : str
            Type of polymorphism (snp or indel). If there are multiple alleles
            and one is an indel, ptype will be "indel". If the alleles are all
            single nucleotide variants, ptype will be "snp".

        max_len : int
            Maximum allele length. If greater than one, implies an insertion.

        """
        snp_split = snp_line.strip().split()
        self.pos = int(snp_split[0]) - 1
        self.alleles = [snp_split[1], snp_split[2]]
        self.ptype = "snp"
        self.max_len = 0
        
        for i in range(len(self.alleles)):
            if self.alleles[i] == "-":
                self.alleles[i] = ""
                self.ptype = "indel"
            elif len(self.alleles[i]) > self.max_len:
                self.max_len = len(self.alleles[i])
                
        if self.max_len > 1:
            self.ptype = "indel"
        
    def add_allele(self, new_alleles):
        """
        Add new alleles for a snp or indel.
        
        Parameters:
        -----------

        new_alleles : str
            List of alleles (each allele is a string like "A").

        """
        for new_allele in new_alleles:
            if new_allele == "-":
                self.ptype = "indel"
                new_allele = ""
            # Only add the new allele if it doesn't already exist.
            if not (new_allele in self.alleles):
                self.alleles.append(new_allele)
                if len(new_allele) > self.max_len:
                    self.max_len = len(new_allele)
                    
        if self.max_len > 1:
            self.ptype = "indel"

    def shift_indel(self):
        """
        Currently not used anywhere.
        """
        self.pos += 1
        self.max_len -= 1
        i = 0
        while i < len(self.alleles):
            if len(self.alleles) <= 1:
                self.alleles.pop(i)
            else:
                self.alleles[i] = self.alleles[i][1:]
                i += 1
        self.alleles.append("")

class BamScanner:
    """
    Class to keep track of all the information read in from the bamfile/snpfile.
    """

    def __init__(self, is_paired_end, max_window, file_name, keep_file_name,
                 remap_name, remap_num_name, fastq_names, snp_dir):
        """
        Constructor: opens files, creates initial table.
        
        Attributes:
        -----------

        is_paired_end : boolean
            Boolean indicating whether input data are paired end.

        snp_dir : str
            Path to directory that contains gzipped SNP files (one per
            chromosome).

        bamfile : pysam.Samfile
            Input bam file that we are reading.

        keep_bam : pysam.Samfile
            Output bam file of reads that do not overlap SNPs.

        remap_bam : pysam.Samfile
            Output bam file of reads that do overlap SNPs and need to be
            remapped.

        remap_num_file : gzip.GzipFile
            File to write XXX to.

        fastqs : list
            List of gzip.GzipFile objects for the different fastqs that will
            contain the reads to remap.

        read_table : list
            List of lists. Sublist i contains the reads whose positions are
            [real read position] % max_window.

        cur_read : pysam.XXX
            Current read from the bam file.

        end_of_file : boolean
            Boolean indicating whether we've reached the end of the bam file.

        remap_num : int
            A counter for the number of reads to be remapped. This starts at one
            and is incremented when a read (pair) is written to the fastq
            file(s).

        ref_match : int
            This is incremented everytime a read sequence matches a SNP
            genotype. Note that a particular read sequence can be looked at
            multiple times if it has multiple SNPs, so this is somewhat hard to
            interpret.

        alt_match : int
            This is initialized but not used anywhere.

        no_match : int
            This is incremented everytime a read sequence doesn't matche a SNP
            genotype. Note that a particular read sequence can be looked at
            multiple times if it has multiple SNPs, so this is somewhat hard to
            interpret.

        toss : int
            Number of reads tossed.

        nosnp : int
            Number of reads with no SNPs.

        remap : int
            Number of reads to remap.

        tot : int
            I think this is total number of reads, although it is only
            incremented in empty_slot_single(self) so it doesn't seem to be
            fully implemented right now.

        printstats : boolean
            Boolean for some print statements, currently not used.

        num_reads : int
            Number of reads for a given window that we have read but not yet
            written. This number is incremented when we read in a read and
            decremented when we pop a read out of the read table.

        cur_snp : SNP
            The current SNP to be or being parsed. 

        pos : int
            The current genomic position we are analyzing.

        chr_num : int
            Bam file ID number of the current chromosome we are analyzing.

        chr_name : str
            Name of the current chromosome we are analyzing.

        max_window : int
            Size of the window in base pairs to process reads. All of the reads
            and SNPs within max_window base pairs are processed at once. Any
            junction-spanning reads (i.e. with N in the cigar) that extend
            outside of the window are thrown out.

        """
        
        self.is_paired_end = is_paired_end
        
        # Read in all input files and create output files
        self.snp_dir = snp_dir
        self.bamfile = pysam.Samfile(file_name,"rb")
        self.keep_bam = pysam.Samfile(keep_file_name, "wb",
                                      template=self.bamfile)
        self.remap_bam = pysam.Samfile(remap_name, "wb", template=self.bamfile)
        self.remap_num_file = gzip.open(remap_num_name, "w")
        self.fastqs = [gzip.open(fqn,"w") for fqn in fastq_names]
        try:
            self.cur_read = self.bamfile.next()
        except:
            sys.stderr.write("No lines available for input")
            return()
        self.end_of_file = False

        self.remap_num = 1
        self.ref_match = 0
        self.alt_match = 0
        self.no_match = 0
        self.toss = 0
        self.nosnp = 0
        self.remap = 0
        self.tot = 0
        self.printstats = True
        
        self.num_reads = 0

        self.cur_snp = None
        
        self.pos = self.cur_read.pos
        self.chr_num = self.cur_read.tid
        self.chr_name = self.bamfile.getrname(self.cur_read.tid)
        self.max_window = max_window
                
        # Initialize the read tracking tables.
        self.read_table = [[] for x in range(self.max_window)]
        
        # Initialize the SNP and indel tracking tables.
        self.switch_chr()
        
        # Fill all tables.
        self.fill_table()
        
    def fill_table(self): 
        """
        Fills the table of reads starting from the current position and
        extending for the next <max_window> base pairs. The read table is a
        list of lists of length max_window. If the position of the current read
        is 100, the first sublist contains all of the reads at position 100, the
        next sublist contains all of the reads at position 101, etc.
        """
        if self.end_of_file:
            return()
            
        # For first read we need to set self.pos and initialize the SNP table.
        if self.num_reads == 0:
            self.pos = self.cur_read.pos
            self.init_snp_table()
            #self.num_reads+=1000
            
        while (self.cur_read.tid == self.chr_num) and \
          (self.cur_read.pos < (self.pos + self.max_window)):
            self.num_reads += 1
            self.read_table[self.cur_read.pos % self.max_window].append(self.cur_read)
           
            # Get a new read and check for the end of the file. 
            try:
                self.cur_read = self.bamfile.next()
            except:
                self.empty_table()
                self.end_of_file = True
                return()
        
        # Check to see if we've come across a new chromosome.
        if self.cur_read.tid != self.chr_num:
            self.empty_table()
            self.chr_num = self.cur_read.tid
            try:
                self.chr_name = self.bamfile.getrname(self.chr_num)
            except:
                sys.stderr.write("Problem with tid: " + str(self.chr_num) + "\n")
                self.skip_chr()
            self.pos = self.cur_read.pos
            self.switch_chr()
            self.fill_table()


    def switch_chr(self):
        """Switches to looking for SNPs on next chromosome."""
        chr_match = False
        while not chr_match and not self.end_of_file:
            try:
                self.snpfile = gzip.open("%s/%s.snps.txt.gz" 
                                         % (self.snp_dir,self.chr_name))
                sys.stderr.write("Starting on chromosome " + self.chr_name+"\n")
                chr_match = True
            except:
                sys.stderr.write("SNP file for chromosome " + 
                                 self.chr_name + " is not found. Skipping these reads.\n")
                self.skip_chr()
        
        self.end_of_snp_file = False
        self.get_next_snp()

    def init_snp_table(self):
        """
        Creates an empty SNP table starting from the position of the current
        and extending max_window base pairs. The SNP table is max_window long
        and has a zero if there are no variants overlapping a position or
        contains a SNP object if there is variant that overlaps a given
        position. 
        
        Also creates an indel table which is a list of lists of length
        max_window.

        Also creates an indel dict which is a dict whose keys are genomic
        positions and whose values are SNP objects whose ptype is indel.
        """
        # Number of SNPs in this table. I think this is total number of
        # different alleles across the whole table.
        self.num_snps = 0
        self.indel_dict = {}
        self.snp_table = [0 for x in range(self.max_window)]
        self.indel_table = [[] for x in range(self.max_window)]
        # Get SNPs in this window but skip SNPs that are upstream of the current
        # read.
        while not self.end_of_snp_file and self.cur_snp.pos < self.pos:
            self.get_next_snp()

        # Add SNPs downstream of the current read and within the current window.
        while not self.end_of_snp_file and (self.cur_snp.pos < self.pos + self.max_window):
            if self.cur_snp.ptype == "snp":
                self.add_snp()
            else:
                self.add_indel()
            self.get_next_snp()
        
        #sys.stderr.write(str(self.num_snps)+"\n")

    def add_snp(self):
        """
        Add a SNP to the SNP table. If the SNP table has a zero at this
        position, the SNP object will replace the zero. If the SNP table
        already has a SNP object at this position, the SNP will be added as new
        alleles. 
        """
        cur_pos = self.cur_snp.pos % self.max_window
        if self.snp_table[cur_pos] == 0:
            self.num_snps += 1
            self.snp_table[cur_pos] = self.cur_snp
        elif isinstance(self.snp_table[cur_pos], SNP):
            self.snp_table[cur_pos].add_allele(self.cur_snp.alleles)     
            
    def add_indel(self):
        """
        Add an indel to the indel table and indel dict. If there is already an
        indel in the indel dict at this position, add the alleles from cur_snp.
        TODO: Finish figuring out what's going on at the end of this method.
        """
        position = self.cur_snp.pos
        if self.indel_dict.has_key(position):
            start = self.indel_dict[position].max_len
            self.indel_dict[position].add_allele(self.cur_snp.alleles)
        else:
            self.indel_dict[position] = self.cur_snp
            start = 0
        end = self.indel_dict[position].max_len
        i = start
        while (i < end) and ((self.cur_snp.pos + i) < (self.pos + self.max_window)):
            self.indel_table[(self.cur_snp.pos + i) % self.max_window].append(position)
            i += 1

    def get_next_snp(self):
        """Read in next SNP or signal end of file."""
        snp_line = self.snpfile.readline()
        if snp_line:
            self.cur_snp = SNP(snp_line)
        else:
            self.end_of_snp_file = True

    def skip_chr(self):
        """Skips all of the reads from the chromosome of the current read and
        moves on to the next chromosome. Used if the SNP file can't be
        located."""
        while self.cur_read.tid == self.chr_num:
            try:
                self.cur_read = self.bamfile.next()
            except:
                self.empty_table()
                self.end_of_file = True
                return

        self.chr_num = self.cur_read.tid
        try:
            self.chr_name = self.bamfile.getrname(self.chr_num)
        except:
            sys.stderr.write("Problem with tid: " + str(self.chr_num) + "\n")
            self.skip_chr()

    def empty_slot_single(self):
        """Processes all reads that map to the current position and
        removes them from the read table Treats reads as single-end"""        
        cur_slot = self.pos % self.max_window
        while len(self.read_table[cur_slot]) > 0:
            self.tot += 1
            read = self.read_table[cur_slot].pop()
            self.num_reads -= 1
            seqs = self.check_for_snps(read, 0)
            num_seqs = len(seqs)
            if (num_seqs == 0) or (num_seqs > 10):
                continue
            if (num_seqs == 1):
                self.keep_bam.write(read)
            else:
                self.remap_num_file.write("%i\n" % (num_seqs-1))
                self.remap_num_file.flush()
                self.remap_bam.write(read)
                for seq in seqs[1:]:
                    loc_line = "%i:%s:%i:%i" % (self.remap_num, 
                                                self.chr_name, read.pos,num_seqs-1)
                    self.fastqs[0].write("@%s\n%s\n+%s\n%s\n" % (loc_line, seq, loc_line,read.qual))
                self.remap_num += 1
        #if self.printstats: 
        #    sys.stderr.write(str(self.tot)+" "+str(self.nosnp)+" "+str(self.remap)+" "+str(self.toss)+"\n")
        #    self.printstats = False
        #sys.stderr.write(str(self.ref_match)+" "+str(self.alt_match)+" "+str(self.no_match)+"\n")
        self.shift_SNP_table()

    def empty_slot_paired(self):
        """Processes all reads that map to the current position and
        removes them from the read table Treats reads as paired-end"""
        
        cur_slot = self.pos % self.max_window

        # While there are reads in this slot...
        while len(self.read_table[cur_slot]) > 0:
            # Pop the first read in the slot
            read = self.read_table[self.pos % self.max_window].pop()
            self.num_reads -= 1

            # Figure out the matching read position
            pair_chr_num = read.rnext 
            pair_pos = read.mpos 
            if (pair_chr_num != self.chr_num) or \
              ((pair_pos - self.pos) > self.max_window):
                continue

            # Find the slot the matching read is in
            pair_slot = pair_pos % self.max_window
            for indx in range(len(self.read_table[pair_slot])):

                #if self.read_table[pair_slot][indx].qname==read.qname:
                # for testing purposes???
                if self.read_table[pair_slot][indx].qname.split(":")[-1] == read.qname.split(":")[-1]:
                    pair_read = self.read_table[pair_slot].pop(indx)
                    self.num_reads -= 1
                    seq1s = self.check_for_snps(read, 0)
                    seq2s = self.check_for_snps(pair_read, read.mpos - read.pos)
                    num_seqs = len(seq1s)*len(seq2s)
                    if (num_seqs == 0) or (num_seqs > 32):
                        break
                    if (num_seqs == 1):
                        self.keep_bam.write(read)
                        self.keep_bam.write(pair_read)
                    else:
                        self.remap_bam.write(read)
                        self.remap_bam.write(pair_read)
                        self.remap_num_file.write("%i\n" % (2*(num_seqs-1)))
                        first = True
                        for seq1 in seq1s:
                            for seq2 in seq2s:
                                if first:
                                    left_pos = min(read.pos,pair_read.pos)
                                    right_pos = max(read.pos,pair_read.pos)
                                    loc_line = "%i:%s:%i:%i:%i" % (self.remap_num, self.chr_name, 
                                                                   left_pos, right_pos, num_seqs-1)
                                    self.fastqs[0].write("@%s\n%s\n+%s\n%s\n" % 
                                                         (loc_line, seq1, loc_line, read.qual))
                                    
                                    self.fastqs[1].write("@%s\n%s\n+%s\n%s\n" % 
                                                         (loc_line, self.reverse_complement(seq2),
                                                          loc_line, pair_read.qual))
                                first = False
                        self.remap_num += 1

                    # stop searching for the pair since it was found
                    break
        
        self.shift_SNP_table()

    def check_for_snps(self, read, start_dist):
        """
        Checks a single aligned read for overlapping SNPs and creates
        alternative sequences for remapping.

        Parameters
        ----------
        read : pysam.AlignedRead
            Read to check for SNPs in.

        start_dist : int
            I think this is the distance from the current position of the
            BamScanner to the start of the read.

        Returns
        -------
        seqs : list
            List of read sequences. This first entry is the read sequence from
            the bam file. Any subsequent sequences are the read sequence from
            the bam file except one base that overlapped a SNP is switched to
            the other allele. If the list is empty, the read overlaps an indel
            or has a CIGAR character besides N or M so we throw it out.
        """
        indx = read.pos % self.max_window
        # p keeps track of the number of read bases we've already analyzed. When
        # p = length of the read, we are done with this read.
        p = 0
        # num_snps is the number of SNPs in this read.
        num_snps = 0
        # I think seg_len is the distance from the current position of the
        # BamScanner to where we are 
        seg_len = start_dist
        seqs = [read.seq]
        if start_dist > 0:
            # has_junc indicates whether the read has an N in the CIGAR although
            # this doesn't seem to be used anywhere.
            has_junc = False
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
        # So a tuple (0, 5) means five matches and (4, 2) means a soft clip of
        # two.

        # We'll go through each cigar tuple one at a time.
        for cigar in read.cigar:
            seg_len += cigar[1]
            # Check whether this cigar segment is longer than the max window.
            # This generally happens if there is a junction read longer than the
            # max window.
            if seg_len > self.max_window:
                #TODO: count the occurrences of this error and just report at
                # the end.
                sys.stderr.write("Segment distance (from read pair and junction separation) "
                                 "is too large. A read has been thrown out. Consider increasing "
                                 "the max window size.\n")
                return([])

            if cigar[0] == 4:
                # CIGAR indicates a soft-clipping
                p = p + cigar[1]
            elif cigar[0] == 0:
                # CIGAR indicates a match alignment to the reference genome.
                # Since there is a match, let's go through each matched base and
                # see whether it contains a SNP.
                for i in range(cigar[1]):  
                    if len(self.indel_table[indx]) == 0:
                        snp = self.snp_table[indx]
                        if snp != 0:
                            num_snps += 1
                            if num_snps > 10:
                                # If there are more than 10 snps overlapping,
                                # throw out the read to prevent memory blow-up.
                                # TODO: should we increment self.toss here?
                                return([])
                            for seq in list(seqs):
                                matches = 0

                                for geno in snp.alleles:
                                    if seq[p] == geno:
                                        matches += 1
                                        for alt_geno in snp.alleles:
                                            if not alt_geno == geno:
                                                new_seq = (seq[:p] + alt_geno +
                                                           seq[p+1:])
                                                seqs.append(new_seq)
                                if matches == 0:
                                    self.no_match += 1
                                else:
                                    self.ref_match += 1
                    else:
                        # It's an indel, throw it out.
                        self.toss += 1
                        return([])
                    indx = (indx + 1) % self.max_window
                    p += 1
            elif cigar[0] == 3:
                # Skipped in the reference genome (splice junction).
                indx = (indx + cigar[1]) % self.max_window
                has_junc = True
            else:
                # There is a non-N/M in the read CIGAR--throw out the read.
                self.toss += 1
                return([])
                
        if len(seqs) == 1:
            self.nosnp += 1
        else:
            self.remap += 1
        return seqs
    
    def shift_SNP_table(self):             
        """Shifts the SNP table over one position and makes sure that
        indels are not lost."""
        
        self.pos += 1

        # Current slot to fill is the position + max_window - 1
        cur_slot=(self.pos-1) % self.max_window

        # Delete indels that are no longer used (if they ended at the previous position)
        for indel_pos in self.indel_table[cur_slot]:
            if (indel_pos + self.indel_dict[indel_pos].max_len-1) == (self.pos-1):
                del self.indel_dict[indel_pos]
        
        self.indel_table[cur_slot]=[]
        
        # Carry over indels from the previous slot
        for indel_pos in self.indel_table[cur_slot-1]:
            if (indel_pos + self.indel_dict[indel_pos].max_len-1) >= (self.pos+self.max_window-1):
                self.indel_table[cur_slot].append(indel_pos)

        if self.snp_table[cur_slot] != 0:
            self.num_snps -= 1
        self.snp_table[cur_slot] = 0

        # See if there is a SNP overlapping the current spot.
        while not self.end_of_snp_file and self.pos + self.max_window-1 > self.cur_snp.pos:
            sys.stderr.write(str(self.num_snps) + " " + str(self.pos) + " " + 
                             str(self.cur_snp.pos)+" !!!\n")
            sys.stderr.write("SNP out of order has been skipped\n")
            self.get_next_snp()

        while not self.end_of_snp_file and (self.cur_snp.pos == (self.pos + self.max_window - 1)):
            if self.cur_snp.ptype == "snp":
                self.add_snp()
            else:
                self.add_indel()
                if not self.cur_snp.pos in self.indel_table[cur_slot]:
                    self.indel_table[cur_slot].append(cur_snp.pos)
            self.get_next_snp()

    def empty_table(self):
        """Completely empties the read_table by repeatedly calling
        empty_slot function"""
        end_pos = self.pos + self.max_window
        while self.pos < end_pos:
            if self.is_paired_end:
                self.empty_slot_paired()
            else:
                self.empty_slot_single()

    def complement(self, letter):
        if letter == 'A':
            return('T')
        elif letter == 'T':
            return('A')
        elif letter == 'C':
            return('G')
        elif letter == 'G':
            return('C')
        else:
            return(letter)

    def reverse_complement(self,read):
        reverse = ""
        for letter in read:
            reverse = self.complement(letter)+reverse
        return reverse
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", action='store_true', 
                        dest='is_paired_end', default=False)
    parser.add_argument("-m", action='store', 
                        dest='max_window', type=int, default=100000)
    parser.add_argument("infile", action='store')
    parser.add_argument("snp_dir", action='store')
    
    options = parser.parse_args()
    infile = options.infile
    snp_dir = options.snp_dir
    name_split = infile.split(".")
    
    if len(name_split) > 1:
        pref = ".".join(name_split[:-1])
    else:
        pref = name_split[0]

    pysam.sort(infile, pref+".sort")

    sort_file_name = pref + ".sort.bam"
    keep_file_name = pref + ".keep.bam"
    remap_name = pref + ".to.remap.bam"
    remap_num_name = pref + ".to.remap.num.gz"

    if options.is_paired_end:
        fastq_names = [pref + ".remap.fq1.gz",
                       pref + ".remap.fq2.gz"]
    else:
        fastq_names = [pref + ".remap.fq.gz"]

    bam_data = BamScanner(options.is_paired_end, options.max_window, 
                          sort_file_name, keep_file_name, remap_name, 
                          remap_num_name, fastq_names, snp_dir)
    
    # TODO: Wrap stuff below in a function or method.
    # TODO: Necessary to call fill_table here? Called with bam_data is
    # initialized.
    bam_data.fill_table()
    
    #i=0
    while not bam_data.end_of_file:
        #i+=1
        #if i>50000:
            #sys.stderr.write(str(asizeof.asizeof(bam_data))+"\t"+str(asizeof.asizeof(bam_data.snp_table))+"\t"+str(asizeof.asizeof(bam_data.read_table))+"\t"+str(bam_data.num_reads)+"\t"+str(bam_data.num_snps)+"\n")
            #sys.stderr.write(str(asizeof.asizeof(bam_data))+"\t"+str(bam_data.num_reads)+"\t"+str(bam_data.num_snps)+"\t"+str(len(bam_data.indel_dict))+"\n")
            #i=0
        if options.is_paired_end:
            bam_data.empty_slot_paired()
        else:
            bam_data.empty_slot_single()
        bam_data.fill_table()
    
    sys.stderr.write("Finished!\n")

if __name__ == '__main__':
    main()
