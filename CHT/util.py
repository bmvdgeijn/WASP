

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
