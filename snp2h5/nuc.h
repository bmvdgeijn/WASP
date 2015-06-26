#ifndef __NUC_H__
#define __NUC_H__

#include <stdio.h>
#include <ctype.h>
#include "seq.h"

enum nucleotide {NUC_A=0, NUC_C, NUC_G, NUC_T, 
		 NUC_M, /* [M]ethyl: A or C */
		 NUC_K, /* [K]eto: G or T */
		 NUC_R, /* pu[R]ine: A or G */
		 NUC_Y, /* p[Y]rimidine: C or T */
		 NUC_W, /* [W]eak: A or T */
		 NUC_S, /* [S]trong: G or C */
		 NUC_GAP, NUC_N, NUM_NUCS};

/* the number of "true" nucleotides, i.e. not gaps or ambiguity chars */
#define NUM_REAL_NUCS 4

char nuc_id_to_char(const unsigned char id);

unsigned char nuc_char_to_id(const char nuc);

int nuc_is_real(const unsigned char id);

int nuc_is_ambi(const unsigned char id);

unsigned char nuc_comp(const unsigned char nuc_id);

char *nuc_ids_to_str(char *buf, const unsigned char *ids, const long len);

int nuc_ids_have_n(unsigned char *ids, const long len);

unsigned char *nuc_str_to_ids(unsigned char *buf, const char *str, 
			      const long len);

void nuc_ids_revcomp(unsigned char *nuc_ids, const long len);


#endif
