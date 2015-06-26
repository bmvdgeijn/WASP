
from chromosome import Chromosome

class CoordError(Exception):
    """An exception indicating that something is wrong with a coordinate"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return str(self.value)


class Coord(object):
    def __init__(self, chrom, start, end, strand=0,
                 score=None, idnum=None, name=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.idnum = idnum
        self.score = score
        self.name = name

        if start > end:
            raise CoordError("start (%d) should be less than or "
                             "equal to end (%d)" % (start, end))

        if start < 1:
            raise CoordError("start (%d) should not be less "
                             "than 1" % start)

        if end > chrom.length:
            raise CoordError("end (%d) should not be greater than "
                             "length of chromosome "
                             "(%d)" % (end, chrom.length))
        
        if strand != 0 and strand != 1 and strand != -1:
            raise CoordError("strand should be one of (-1, 0, 1)")


    def key(self, use_strand=False):
        """Returns a tuple to be used for sorting of coordinates. If
        use_strand is True the tuple consists of (chrom.idnum, strand, start),
        otherwise the tuple is (chrom.idnum, start)"""
        if use_strand:
            return (self.chrom.idnum, self.strand, self.start)

        return (self.chrom.idnum, self.start)
        
        
    def __str__(self):
        """returns a string representation of this coordinate"""
        if self.idnum:
            id_str = str(self.idnum) + " "
        else:
            id_str = ""
            
        if self.strand == 1:
            strand_str = "(+)"
        elif self.strand == -1:
            strand_str = "(-)"
        else:
            strand_str = "(.)"
        
        return(id_str + str(self.chrom) + ":" + str(self.start) + "-" + \
               str(self.end) + strand_str)
    
    
    def length(self):
        """Returns the size of the region spanned by the
        coordinate in bases"""
        return self.end - self.start + 1
    
    def overlaps(self, other, use_strand=False):
        """Returns True if the coords overlap and False if they don't.
        If use_strand is True then the coordinates must be on the same
        strand to be considered overlapping."""
        if self.chrom.idnum != other.chrom.idnum:
            return False
        if use_strand and self.strand != other.strand:
            return False
        if self.start <= other.end and self.end >= other.start:
            return True
        return False

    def copy(self):
        """Creates a copy of this coordinate object and returns it.
        The copy is shallow in the sense that only the reference to
        the chromosome attribute is copied (i.e. a new Chromosome
        is not created)."""
        return Coord(self.chrom, self.start, self.end,
                     strand=self.strand, score=self.score,
                     idnum=self.idnum)


    def expand(self, n_bp):
        """Expands this coordinate by n_bp in both directions, but not
        exceeding the boundaries of the chromosome that this coordinate
        is on"""
        self.start = self.start - n_bp
        self.end = self.end + n_bp

        if self.chrom:
            # if there is an associated chromosome, don't go past the ends
            if self.start < 1:
                self.start = 1
            if self.end > self.chrom.length:
                self.end = self.chrom.length

    def within(self, other, use_strand=False):
        """checks whether one coordinate is completely contained
        within the other one. If use_strand is True then the
        this coordinate must also be on the same strand as the
        other one, in order to be considered 'contained'"""
        if self.chrom and other.chrom and (self.chrom.name != other.chrom.name):
            return False
        
        if use_strand and self.strand != other.strand:
            return False

        if self.start >= other.start and self.end <= other.end:
            return True

        return False
