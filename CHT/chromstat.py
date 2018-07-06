
import sys
import numpy as np


class ChromStats(object):
    def __init__(self):
        self.n = 0
        self.n_nan = 0
        self.sum = 0
        self.min = None
        self.max = None


    def mean(self):
        """Calculates mean of sites that are not nan
        on this chromsome"""
        n = self.n - self.n_nan
        if n == 0:
            return np.inf
        
        return self.sum / float(n)


    def set_from_vals(self, vals):
        self.n = vals.size

        if str(vals.dtype).startswith('float'):
            nan_vals = np.isnan(vals)
            self.n_nan = np.sum(nan_vals)

            if self.n_nan < self.n:
                self.min = np.min(vals[~nan_vals])
                self.max = np.max(vals[~nan_vals])
                self.sum = np.sum(vals[~nan_vals])
        else:
            self.min = np.min(vals)
            self.max = np.max(vals)
            self.sum = np.sum(vals)

        

    def add(self, other):
        self.n += other.n
        self.n_nan += other.n_nan
        self.sum += other.sum
        
        if (self.min is None) or (other.min is not None and 
                                  self.min > other.min):
            self.min = other.min

        if (self.max is None) or (other.max is not None and
                                  self.max < other.max):
            self.max = other.max


    def __str__(self):
        return "n=%d n_nan=%s min=%s max=%s sum=%s" % \
            (self.n, str(self.n_nan), str(self.min), str(self.max), 
             str(self.sum))



def calc_stats(h5f, chrom_list, verbose=False):
    """Calculates stats for each chromosome in provided list as well
    as combined stats."""
    
    combined = ChromStats()

    for chrom in chrom_list:
        chrom_stat = ChromStats()
        node_name = "/%s" % chrom.name
        if node_name in h5f:
            node = h5f.get_node("/%s" % chrom.name)
            vals = node[:]
            chrom_stat.set_from_vals(vals)
            if verbose:
                sys.stderr.write("%s %s\n" % (str(chrom), str(chrom_stat)))
        else:
            sys.stderr.write("skipping chromosome %s because "
                             "not present in HDF5 file" % chrom.name)
        
        combined.add(chrom_stat)

    return combined


def set_stats(h5f, chrom_list, verbose=False):
    """Calculates stats for each chromosome and entire track and
    stores them as attributes on the chromosome nodes. The
    provided HDF5 file handle must have been opened in append mode"""
    combined = ChromStats()

    for chrom in chrom_list:
        node_name = "/%s" % chrom.name
        if node_name in h5f:
            chrom_stat = ChromStats()

            node = h5f.get_node(node_name)
            chrom_stat.set_from_vals(node[:])

            node.attrs.n = chrom_stat.n
            node.attrs.n_nan = chrom_stat.n_nan
            node.attrs.min = chrom_stat.min
            node.attrs.max = chrom_stat.max
            node.attrs.sum = chrom_stat.sum
            node.flush()

            if verbose:
                sys.stderr.write("%s %s\n" % (str(chrom), str(chrom_stat)))
            combined.add(chrom_stat)
        else:
            sys.stderr.write("skipping chromosome %s because "
                             "not present in HDF5 file\n" % chrom.name)
    
    return combined



def get_stats(h5f, chrom_list, verbose=False):
    """Retrieves stats that are stored as attributes for the specified
    set of chromosomes."""
    combined = ChromStats()
    chrom_stat = ChromStats()
    
    for chrom in chrom_list:
        node_name = "/%s" % chrom.name
        if node_name in h5f:
            node = h5f.get_node(node_name)
            if 'n' not in node.attrs:
                raise ValueError("Stat attributes are not set for track %s"
                                 % track.name)

            chrom_stat.n = node.attrs.n
            chrom_stat.n_nan = node.attrs.n_nan
            chrom_stat.min = node.attrs.min
            chrom_stat.max = node.attrs.max
            chrom_stat.sum = node.attrs.sum

            if verbose:
                sys.stderr.write("%s %s\n" % (str(chrom), str(chrom_stat)))
            combined.add(chrom_stat)
        else:
            sys.stderr.write("skipping chromosome %s because "
                             "not present in HDF5 file\n" % chrom.name)
            

    return combined
        
