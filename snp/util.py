import subprocess
import os
import gzip



def is_gzipped(filename):
    """Checks first two bytes of provided filename and looks for
    gzip magic number. Returns true if it is a gzipped file"""
    f = open(filename, "rb")

    # read first two bytes
    byte1 = f.read(1)
    byte2 = f.read(1)
    
    f.close()

    # check against gzip magic number 1f8b
    return (byte1 == chr(0x1f)) and (byte2 == chr(0x8b))



def check_open(filename, mode="r"):
    """Tries to open file and return filehandle. Takes into account
    that file may be gzipped. Raises exception if mode is write and 
    file already exists."""
    if mode.startswith("w") and os.path.exists(filename):
        raise IOError("file %s already exists" % filename)

    if mode == "w" or mode == "wb":
        if filename.endswith(".gz"):
            # create a gzipped file
            mode = "wb"
            return gzip.open(filename, mode)
    elif mode.startswith("r"):
        if is_gzipped(filename):
            # open a gzipped file
            mode = "rb"
            return gzip.open(filename, mode)

    return open(filename, mode)



def count_lines(filename):

    if not os.path.exists(filename):
        raise IOError("file '%s' does not exist" % filename)

    if not os.path.isfile(filename):
        raise IOError("'%s' is not a regular file" % filename)
    
    if is_gzipped(filename):
        p1 = subprocess.Popen(['zcat', filename],
                              stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['wc', '-l'], stdin=p1.stdout,
                               stdout=subprocess.PIPE)
        
        p1.stdout.close()

        out = p2.communicate()[0]
    else:
        out = subprocess.Popen(['wc', '-l', filename],
                               stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT).communicate()[0]

    return int(out.split()[0])
