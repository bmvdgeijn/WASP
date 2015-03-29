import sys
import gzip
import re
import argparse

import tables


class ChromDesc(tables.IsDescription):
    """This class defines the format of the chromosome table that is
    stored in HDF5 files"""
    idnum = tables.Int32Col()
    name = tables.StringCol(32)
    length = tables.Int32Col()
    is_sex = tables.BoolCol(dflt=False)
    is_auto = tables.BoolCol(dflt=True)
    is_rand = tables.BoolCol(dflt=False)
    is_hap = tables.BoolCol(dflt=False)
    is_mito = tables.BoolCol(dflt=False)
    is_y = tables.BoolCol(dflt=False)
    is_x = tables.BoolCol(dflt=False)



class ChromTab(object):
    """Represents a table of chromosomes that can be queried to get 
    specific subsets of chromosomes out."""


    def __init__(self, h5_filename, mode="r"):
        """Opens or creates a new chromosome table. If mode is "r"
        an existing HDF5 file containing chromosome data is opened. If mode
        is "w" a new, empty HDF5 file is created"""

        if mode == "w":
            # create a new file
            self.h5f = tables.openFile(h5_filename, "w")
            
            self.chrom_tab = self.h5f.createTable("/", "chromosome", 
                                                  ChromDesc, "chromosomes")
        else:
            self.h5f = tables.openFile(h5_filename, "r")
            self.chrom_tab = self.h5f.root.chromosome


    def load_from_file(self, input_filename):
        """Parses chromosomes from chromInfo.txt file (can be 
        downloaded from UCSC) and stores them in HDF5 file"""
        chrom_list = parse_chromosomes(input_filename)
        self.load_chromosomes(chrom_list)
        self.chrom_tab.flush()

        
    def load_chromosomes(self, chrom_list):
        """Stores a list of chromosomes in HDF5 file"""
        row = self.chrom_tab.row

        for chrom in chrom_list:
            row['idnum'] = chrom.idnum
            row['name'] = chrom.name
            row['length'] = chrom.length
            row['is_auto'] = chrom.is_auto
            row['is_sex'] = chrom.is_sex
            row['is_rand'] = chrom.is_rand
            row['is_hap'] = chrom.is_hap
            row['is_mito'] = chrom.is_mito
            row['is_y'] = chrom.is_y
            row['is_x'] = chrom.is_x                
            row.append()

    
    def get_chromosome_dict(self):
        """Returns a dictionary of all chromosomes in the HDF5 file,
        keyed on chromosome name"""
        chrom_list = self.get_all_chromosomes()
        chrom_dict = {}
        for chrom in chrom_list:
            chrom_dict[chrom.name] = chrom

        return chrom_dict


    def get_chromosomes(self, get_rand=False, get_auto=True, 
                        get_sex=True, get_x=True, get_y=False, 
                        get_hap=False, get_mito=False):
        """Returns a filtered list of Chromosomes from the
        HDF5 file. Optional flags specify the subset of Chromosomes
        that are returned. By default the 22 autosomes and chrX are
        retrieved (but chrY, the mitochondrial chromosome, alternate
        haplotypes, and 'random' chromosomes are not)"""
        chrom_list = []
        chrom = None

        flags = {'is_rand' : get_rand,
                 'is_auto' : get_auto,
                 'is_sex' : get_sex,
                 'is_x' : get_x,
                 'is_y' : get_y,
                 'is_hap' : get_hap,
                 'is_mito' :  get_mito}

        conditions = []
        for k, v in flags.items():
            if not v:
                # don't get this chromosome type
                conditions.append("(%s == False)" % k)

        row_iter = None
        if len(conditions) > 0:
            query = " & ".join(conditions)
            row_iter = self.chrom_tab.where(query)
        else:
            row_iter = self.chrom_tab

        for row in row_iter:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            chrom_list.append(chrom)

        return chrom_list
    

    def get_all_chromosomes(self):
        """Returns an unfiltered list of all of the chromosomes in the
        database"""
        chrom_list = []
        chrom = None

        for row in self.chrom_tab:
            chrom = Chromosome(idnum=row['idnum'],
                               name=row['name'],
                               length=row['length'],
                               is_auto=row['is_auto'],
                               is_hap=row['is_hap'],
                               is_mito=row['is_mito'],
                               is_sex=row['is_sex'],
                               is_x=row['is_x'],
                               is_y=row['is_y'])

            chrom_list.append(chrom)

        return chrom_list


    def get_chromosomes_from_args(self, args):
        """Convenience function, returns a list of chromosome objects
        that are obtained from parsing command line args that may be
        in any one of the following forms 'chr2' or '3', or '1-22'"""
        if isinstance(args, str):
            # if we were given a string, put it in a list
            args = [args]

        if len(args) < 1:
            raise ValueError("expected at least one value, got 0")

        chrom_name_dict = self.get_chromosome_dict()

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



    def get_chromosome(self, name):
        """Retrieves a single chromosome by name"""
        chrom_dict = self.get_chromosome_dict()
        return chrom_dict[name]


    def close(self):
        """Closes Chromosome table by closing underlying HDF5 file"""
        self.h5f.close()


        
        


    


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



        

def parse_chromosomes(filename):
    if filename.endswith(".gz"):
        f = gzip.open(filename)
    else:
        f = open(filename)

    chrom_list = []
    
    for line in f:
        words = line.rstrip().split()

        if len(words) < 2:
            raise ValueError("expected at least two columns per line\n")
        
        chrom = Chromosome(name=words[0], length=words[1])
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




