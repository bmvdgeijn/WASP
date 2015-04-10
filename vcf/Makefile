# Edit the following variables as needed
HDF_INSTALL = $(HOME)/anaconda
EXTLIB = -L$(HDF_INSTALL)/lib
CC          = gcc
# CFLAGS      = -g -DH5_USE_16_API
CFLAGS      = -g -DH5_USE_16_API
LIB         = -lz -lm -lhdf5 -lhdf5_hl

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB)

objects=vcf.o util.o memutil.o err.o chrom.o snptab.o

default: all 

all: vcf2h5 $(objects)

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@ $(INCLUDE)


vcf2h5: vcf2h5.o $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) vcf2h5.c $(INCLUDE) $(LIBSHDF) $(LIB)


clean:
	rm -f $(objects) vcf2h5
