# Edit the following variables as needed
HDF_INSTALL = $(HOME)/anaconda
EXTLIB = -L$(HDF_INSTALL)/lib
CC          = gcc
# CFLAGS      = -g -DH5_USE_16_API
CFLAGS      = -g -DH5_USE_16_API
LIB         = -lz -lm -lhdf5 -lhdf5_hl

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB)

objects=vcf.o impute.o util.o memutil.o err.o chrom.o snptab.o seq.o nuc.o

default: all 

all: snp2h5 fasta2h5 $(objects)

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDE)


snp2h5: snp2h5.o $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) snp2h5.c $(INCLUDE) $(LIBSHDF) $(LIB)


fasta2h5: fasta2h5.o $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) fasta2h5.c $(INCLUDE) $(LIBSHDF) $(LIB)

clean:
	rm -f $(objects) snp2h5 fasta2h5
