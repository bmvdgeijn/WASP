
import tables
import sys

import argparse


def main():
    parser = argparse.ArgumentParser(description="Writes names of samples "
                                     "contained in HDF5 file to stdout")
                                   
    parser.add_argument("h5file", help="HDF5 file containing /samples table")

    options = parser.parse_args()
    
    h5f = tables.openFile(options.h5file)

    for node in h5f.root:
        if node.name.startswith("samples"):
            _, chr_name = node.name.split("_", 1)

            sys.stdout.write("%s:\n" % chr_name)
            for row in node:
                sys.stdout.write("  %s\n" % row['name'])
            sys.stdout.write("\n")
    else:
        sys.stderr.write("%s does not contain samples table\n" % options.h5file)
        exit(2)

    h5f.close()


if __name__ == "__main__":
    main()
