
def check_pysam_version(min_pysam_ver="0.8.4"):
    """Checks that the imported version of pysam is greater than
    or equal to provided version. Returns 0 if version is high enough,
    raises ImportWarning otherwise."""
    import pysam

    min_ver = [int(x) for x in min_pysam_ver.split(".")]
    pysam_ver = [int(x) for x in pysam.__version__.split(".")]

    n_ver = min(len(pysam_ver), len(min_pysam_ver))
    
    for i in range(n_ver):
        if pysam_ver[i] < min_ver[i]:
            raise ImportWarning("pysam version is %s, but pysam version %s "
                                "or greater is required" % (pysam.__version__,
                                min_pysam_ver))
        if pysam_ver[i] > min_ver[i]:
            # version like 1.0 beats version like 0.8
            break
        
    return 0


def check_pytables_version():
    """Checks that PyTables version 3 is being used. PyTables version 3 
    changes the names of many functions and is not backwards compatible
    with PyTables 2. Previous versions of WASP used version 2, but switch
    to version 3 was made at same time as switch to python3."""
    import tables

    pytables_ver = [int(x) for x in tables.__version__.split(".")]

    if pytables_ver[0] < 3:
        raise ImportWarning("pytables version is %s, but pytables version "
                            ">=3 is required" % (tables.__version__))

    return 0


def is_gzipped(filename):
    """Checks first two bytes of provided filename and looks for
    gzip magic number. Returns true if it is a gzipped file"""
    f = open(filename, "rb")

    # read first two bytes
    byte1 = f.read(1)
    byte2 = f.read(1)
    
    f.close()

    # check against gzip magic number 1f8b
    return (byte1 == b'\x1f') and (byte2== b'\x8b')
