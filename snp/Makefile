# Edit the following variables as needed
HDF_INSTALL = $(HOME)/anaconda
EXTLIB = -L$(HDF_INSTALL)/lib
CC          = gcc
# CFLAGS      = -g -DH5_USE_16_API
CFLAGS      = -g -DH5_USE_16_API
LIB         = -lz -lm -lhdf5

INCLUDE   = -I$(HDF_INSTALL)/include
LIBSHDF   = $(EXTLIB)

objects=vcf.o util.o memutil.o err.o

default: all 

all: h5test vcf_to_h5genoprob $(objects)

$(objects): %.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@

h5test: h5test.c
	$(CC) $(CFLAGS) -o $@ h5test.c $(INCLUDE) $(LIBSHDF) $(LIB)


vcf_to_h5genoprob: vcf_to_h5genoprob.o $(objects)
	$(CC) $(CFLAGS) -o $@ $(objects) vcf_to_h5genoprob.c $(INCLUDE) $(LIBSHDF) $(LIB)



clean:
	rm $(objects) h5test vcf_to_h5genoprob
