import sys
import gzip
import re
import argparse


class Chromosome(object):
    """Represents a chromosome. Has a name, length and a few
    descriptive flags, such as whether the chromosome is a sex
    chromosome."""
    
    def __init__(self, idnum=None, name=None, length=None,
                 is_sex=False, is_rand=False, is_hap=False,
                 is_mito=False, is_x=False, is_y=False, is_auto=False):
        self.idnum = idnum
        self.name = name
        self.length = length

        # is this a "random" chromosome?
        self.is_rand = is_rand
        
        # is this a sex chromosome?
        self.is_sex = is_sex

        # is this an alternate haplotype chromosome?
        self.is_hap = is_hap

        # is this the mitochondrial chromosome?
        self.is_mito = is_mito

        # is this the X chromosome
        self.is_x = is_x

        # is this the Y chromosome
        self.is_y = is_y

        # is this an autosome
        self.is_auto = is_auto

        
    def copy(self):
        """Creates a new chromosome object with the same attributes
        as this one"""
        return Chromosome(idnum=self.idnum,
                          name=None, length=None,
                          is_sex=False, is_rand=False, is_hap=False,
                          is_mito=False, is_x=False, is_y=False,
                          is_auto=False)
                          
    def __str__(self):
        """returns a string representatin of this object"""
        return self.name

    def __cmp__(self, other):
        return cmp(self.idnum, other.idnum)


    def key(self):
        """Returns a key for sorting chromosomes based on their name"""
        m = re.match(r"^chr(\d+)", self.name)
        if m:
            # make sure autosomes are sorted numerically by padding with
            # leading 0s
            num = m.groups()[0]

            if len(num) < 3:
                name = ("0" * (3-len(num))) + num
            else:
                name = num
        else:
            # otherwise just sort lexigraphically
            name = self.name

        # first take non-haplo, non-rand, non-sex chromosomes, then
        # sort by name
        return (self.is_hap, self.is_rand, self.is_mito, self.is_sex, name)




def get_chromosome_dict(filename):
    """Returns a dictionary of all chromosomes in the chromInfo.txt 
    filename"""
    chrom_list = get_all_chromosomes(filename)
    chrom_dict = {}
    for chrom in chrom_list:
        chrom_dict[chrom.name] = chrom
        
    return chrom_dict




def get_chromosomes(filename, get_rand=False, get_auto=True, 
                    get_sex=True, get_x=True, get_y=False, 
                    get_hap=False, get_mito=False):
    """Returns a filtered list of Chromosomes read from the 
    specified chromInfo.txt file. Optional flags specify the 
    subset of Chromosomes that are returned. By default the 22 
    autosomes and chrX are retrieved (but chrY, the mitochondrial 
    chromosome, alternate haplotypes, and 'random' chromosomes are not)"""
    all_chrom = get_all_chromosomes(filename)
    chrom = None

    chrom_list = []
    for chrom in all_chrom:
        # use flags as negative filtering criteria only
        # e.g. if get_auto is False exclude all autosomes (including 
        # random and alternative haplotype chromosomes that are also
        # autosomes, regardless of get_rand and get_hap flags).
        if get_rand is False and chrom.is_rand:
            continue
        if get_auto is False and chrom.is_auto:
            continue
        if get_sex is False and chrom.is_sex:
            continue
        if get_x is False and chrom.is_x:
            continue
        if get_y is False and chrom.is_y:
            continue
        if get_hap is False and chrom.is_hap:
            continue
        if get_mito is false and chrom.is_mito:
            continue
    
        chrom_list.append(chrom)

    return chrom_list
    


def get_chromosomes_from_args(filename, args):
    """Convenience function, returns a list of chromosome objects
    in the provided file, that are obtained from parsing command line 
    args that may be in any one of the following forms 'chr2' or '3',
    or '1-22'"""
    if isinstance(args, str):
        # if we were given a string, put it in a list
        args = [args]

    if len(args) < 1:
        raise ValueError("expected at least one value, got 0")

    chrom_name_dict = get_chromosome_dict(filename)

    chrom_id_dict = dict([(str(c.idnum), c) for c \
                          in chrom_name_dict.values()])

    chrom_list = []
    for arg in args:
        words = arg.split("-")
        if len(words) == 2:
            # try parsing argument as a range of chromosome
            # ids like 1-22
            try:
                start = int(words[0])
                end = int(words[1])

                if start <= end:
                    vals = [str(x) for x in range(start, end+1)]
                else:
                    vals = [arg]
            except:
                vals = [arg]
        else:
            vals = [arg]

        for val in vals:
            if val in chrom_name_dict:
                chrom_list.append(chrom_name_dict[val])
            elif val in chrom_id_dict:
                chrom_list.append(chrom_id_dict[val])
            else:
                raise ValueError("unknown chromosome %s" % val)

    return chrom_list



def get_chromosome(filename, name):
    """Retrieves a single chromosome by name from the 
    provided chromInfo.txt file"""
    chrom_dict = get_chromosome_dict(filename)
    return chrom_dict[name]

        
        
def get_all_chromosomes(filename):
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    chrom_list = []
    
    for line in f:
        words = line.rstrip().split()

        if len(words) < 2:
            raise ValueError("expected at least two columns per line\n")
        
        chrom = Chromosome(name=words[0], length=int(words[1]))
        chrom_list.append(chrom)

        lc_name = chrom.name.lower()

        # determine whether this is autosome, sex or mitochondrial chrom
        if re.match('^chr(\d+)', lc_name):
            chrom.is_auto=True
        elif re.match("^chr[W-Zw-z]", lc_name):
            chrom.is_sex = True
        elif lc_name.startswith("chrm"):
            chrom.is_mito = True
        elif lc_name.startswith("chrun") or lc_name.startswith("chrur"):
            chrom.is_rand = True
        else:
            sys.stderr.write("WARNING: could not determine chromosome type "
                             "(autosome, sex, mitochondrial) from name "
                             "'%s'. Assuming 'random'\n" % chrom.name)
            chrom.is_rand = True

        if "rand" in chrom.name:
            # random chromosome
            chrom.is_rand = True

        if "hap" in chrom.name:
            # alt haplotype chromosome
            chrom.is_hap = True

    chrom_list.sort(key=Chromosome.key)

    idnum = 1
    for chrom in chrom_list:
        chrom.idnum = idnum
        idnum += 1

    f.close()

    return chrom_list




