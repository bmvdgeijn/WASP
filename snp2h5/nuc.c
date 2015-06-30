#include <stdio.h>
#include <string.h>

#include "memutil.h"
#include "nuc.h"
#include "util.h"
#include "err.h"


static const char NUC_SYMBOL[NUM_NUCS] = {'A', 'C', 'G', 'T', 
					  'M', 'K', 'R', 'Y', 'W', 'S', 
					  '-', 'N'};


/**
 * Converts a nucleotide ID into a character.
 */
char nuc_id_to_char(const unsigned char id) {
  if(id > NUM_NUCS) {
    my_err("Invalid nucleotide identifier %d", id);
  }
  return NUC_SYMBOL[id];
}


/** 
 * returns TRUE if this is an id for an unambiguous nucleotide 
 * (one of A, C, T, G)
 */
int nuc_is_real(const unsigned char id) {
  return (id < NUC_M);
}

/*
 * returns TRUE if this is an id for a nucleotide ambiguity code
 * (one of M, K, R, Y, W, S).
 */
int nuc_is_ambi(const unsigned char id) {
  return (id > NUC_T) && (id < NUC_GAP);
}


/**
 * Converts a nucleotide to a unique nucleotide integer identifier.
 */
unsigned char nuc_char_to_id(const char nuc) {
  switch(nuc) {
  case('A'): case('a'): return NUC_A;
  case('C'): case('c'): return NUC_C;
  case('T'): case('t'): return NUC_T;
  case('G'): case('g'): return NUC_G;
  case('M'): case('m'): return NUC_M;
  case('K'): case('k'): return NUC_K;
  case('R'): case('r'): return NUC_R;
  case('Y'): case('y'): return NUC_Y;
  case('W'): case('w'): return NUC_W;
  case('S'): case('s'): return NUC_S;
  case('.'): case('-'): case('*'): return NUC_GAP;
  case('N'): case('n'): return NUC_N;
  }

  my_warn("%s:%d: unknown nucleotide character %c",
	  __FILE__, __LINE__, nuc);

  return NUC_N;
}


/**
 * returns the complement of the provided nucleotide ID
 */
unsigned char nuc_comp(const unsigned char nuc_id) {
  switch(nuc_id) {
  case(NUC_A): return NUC_T;
  case(NUC_C): return NUC_G;
  case(NUC_G): return NUC_C;
  case(NUC_T): return NUC_A;
  case(NUC_M): return NUC_K;
  case(NUC_K): return NUC_M;
  case(NUC_R): return NUC_Y;
  case(NUC_Y): return NUC_R;
  case(NUC_W): return NUC_W;
  case(NUC_S): return NUC_S;
  }  
  /* for N, GAP, complement is same */
  return nuc_id;
}


/**
 * Fills the provided buffer with a null-terminated string representation 
 * of the provided array of nucleotide arrays. The provided buffer
 * must be at least len+1 bytes long, and the provided nucleotide
 * id array must be at least len bytes long.
 *
 * If the provided buffer is NULL, a new buffer of length len+1 is
 * allocated and returned.
 */
char *nuc_ids_to_str(char *buf, const unsigned char *ids, const long len) {
  long i;
  
  if(buf == NULL) {
    buf = my_new(char, len+1);
  }

  for(i = 0; i < len; i++) {
    buf[i] = nuc_id_to_char(ids[i]);
  }
  buf[len] = '\0';
  
  return buf;
}



/**
 * Fills the provided buffer with nucleotide ids corresponding to the
 * string representation of DNA sequence provided. The returned array
 * is NOT null terminated. The len argument should correspond to the
 * length of the provided DNA character string. The provided buf must
 * be at least len bytes long. If the provided buf is NULL a new buffer
 * of length len is allocated and returned.
 *
 * The return value is just the ptr to the provided buf.
 */
unsigned char *nuc_str_to_ids(unsigned char *buf, const char *str, 
			      const long len) {
  long i;
  
  if(buf == NULL) {
    buf = my_new(unsigned char, len);
  }

  for(i = 0; i < len; i++) {
    buf[i] = nuc_char_to_id(str[i]);
  }
  return buf;
}




/**
 * Does an in-place reverse complement of an array of
 * nucleotide ids
 */
void nuc_ids_revcomp(unsigned char *nuc_ids, long len) {
  long i;

  util_breverse(nuc_ids, len);
  for(i = 0; i < len; i++) {
    nuc_ids[i] = nuc_comp(nuc_ids[i]);
  }
}


/**
 * Returns true if any nucleotides in provided array are NUC_N
 */
int nuc_ids_have_n(unsigned char *nucs, const long len) {
  long i;
  for(i = 0; i < len; i++) {
    if(nucs[i] == NUC_N) {
      return TRUE;
    }
  }
  return FALSE;
}

