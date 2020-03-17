# Edit the following variables as needed
HDF_INSTALL=$(CONDA_PREFIX)
EXTLIB=-L$(HDF_INSTALL)/lib -Wl,-rpath,$(HDF_INSTALL)/lib
CC=gcc
LIB=-lz -lm -lhdf5 -lhdf5_hl

INCLUDE=-I$(HDF_INSTALL)/include
LIBSHDF=$(EXTLIB)
CFLAGS=-DH5_USE_16_API $(INCLUDE) -Wall -g

objects=vcf.o impute.o util.o memutil.o err.o chrom.o snptab.o seq.o nuc.o sampletab.o

default: all

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDE)

snp2h5: snp2h5.o $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) snp2h5.c $(LIBSHDF) $(LIB)

fasta2h5: fasta2h5.o $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) fasta2h5.c $(LIBSHDF) $(LIB)

all: snp2h5 fasta2h5 $(objects)

clean:
	rm -f $(objects) snp2h5.o fasta2h5.o snp2h5 fasta2h5
