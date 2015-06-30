#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>

#include "nuc.h"
#include "seq.h"
#include "util.h"
#include "memutil.h"
#include "err.h"



/**
 * Expands the current sequence buffer, by doubling it's size
 */
void seq_expand(Seq *seq) {
  seq->buf_sz = seq->buf_sz * 2;
  seq->sym = realloc(seq->sym, seq->buf_sz);

  if(!seq->sym) {
    my_err("%s:%d: failed to expand sequence buffer to %ld bytes", 
	  __FILE__, __LINE__, seq->buf_sz);
  }
}



/**
 * Initializes an empty sequence
 */
Seq *seq_new() {
  Seq *seq;

  seq = my_new(Seq, 1);
  seq->name[0] = '\0';
  seq->len = 0;
  seq->buf_sz = SEQ_DEFAULT_BUF_SZ;
  seq->sym = my_new(unsigned char, seq->buf_sz);

  return seq;
}




static int read_fasta_header(Seq *seq, gzFile f) {
  int hdr_idx = 0;
  int truncated = FALSE;
  int errnum;

  char c;

  /* expect first character to be '>' */
  c = gzgetc(f);
  if(c != '>') {
    if(gzeof(f)) {
      /* we've reached the end of the file */
      return 0;
    }
    if(c == -1) {
      const char *err = gzerror(f, &errnum);
      my_err("%s:%d: failed to read from fasta file: %s (errnum=%d)\n",
	     __FILE__, __LINE__, err, errnum);
    } else {
      my_err("%s:%d: expected fasta record to start with '>' not \\%d", 
	     __FILE__, __LINE__, c);
    }
  }

  /* fill in header */
  hdr_idx = 0;
  while((c = gzgetc(f)) != EOF) {
    if(c == '\n' || c == '\0' || c == '\r') {
      /* newline terminates header */
      break;
    }

    if(hdr_idx < SEQ_MAX_NAME_SZ) {
      seq->name[hdr_idx] = c;
      hdr_idx++;
    } else {
      truncated = TRUE;
    }
  }
  seq->name[hdr_idx] = '\0';

  if(truncated) {
    my_warn("truncated sequence name to max size of %d: '%s'",
	    SEQ_MAX_NAME_SZ, seq->name);
  }

  return 1;
}


/**
 * Sets the provided Seq's sym array with nucleotide symbols parsed
 * from the provided NULL-terminated sequence string.  The sequence
 * bufffer is expanded as necessary.
 */
void seq_read_str(Seq *seq, char *seq_str) {
  long i;

  /* parse sequence */
  i = 0;
  while(seq_str[i] != '\0') {
    if(i >= seq->buf_sz) {
      seq_expand(seq);
    }
    seq->sym[i] = nuc_char_to_id(seq_str[i]);
    i++;
  }

  seq->len = i;
}


/**
 * Reads a single FASTA record from the file with the provided name.
 * Returns length of read sequence, or -1 on failure.
 */
long seq_read_fasta_from_file(Seq *seq, const char *filename) {
  gzFile f;
  long seq_len;
  

  f = util_must_gzopen(filename, "rb");

  seq_len = seq_read_fasta_record(seq, f);

  gzclose(f);

  return seq_len;
}


/**
 * Reads the next fasta record from a file into the provided sequence,
 * expanding the sequence buffer as necessary. Returns number of
 * bases in sequence if sequence was read or -1 if at the end of the file.
 */
long seq_read_fasta_record(Seq *seq, gzFile f) {
  long seq_idx;
  char c;

  fprintf(stderr, "reading from fasta file\n");

  
  if(!read_fasta_header(seq, f)) {
    seq->len = 0;
    return -1;
  }

  /* read fasta sequence */
  seq_idx = 0;
  while((c = gzgetc(f)) != EOF) {
    if(isspace(c)) {
      /* skip whitespace */
      continue;
    }
    if(c == '>') {
      /* this is the end of the FASTA record */
      if(gzungetc(c, f) == EOF) {
	my_err("%s:%d: could not unget '>' character", __FILE__, __LINE__);
      }
      break;
    }

    if(seq_idx >= seq->buf_sz) {
      /* expand the seq buffer */
      seq_expand(seq);
    }
    seq->sym[seq_idx] = nuc_char_to_id(c);
    seq_idx += 1;
  }
  seq->len = seq_idx;


  fprintf(stderr, "read %ld bp of sequence\n", seq->len);
  
  return seq->len;
}


/**
 * Writes the sequence to the provided file in fasta format
 */
void seq_write_fasta_record(Seq *seq, gzFile f) {
  long i;
  int line_len;

  /* write header */
  if(seq->name == NULL) {
    gzprintf(f, ">\n");
  } else {
    gzprintf(f, ">%s\n", seq->name);
  }

  line_len = 0;

  /* write nucleotides */
  for(i = 0; i < seq->len; i++) {
    gzprintf(f, "%c", nuc_id_to_char(seq->sym[i]));
    line_len += 1;

    if(line_len >= SEQ_FASTA_LINE_LEN) {
      /* start a new line */
      gzprintf(f, "\n");
      line_len = 0;
    }
  }

  if(line_len != 0) {
    /* end last line */
    gzputc(f, '\n');
  }
}



/*
 * Returns a string representation of the nucleotides of this
 * sequence. The returned string should be freed when it is no longer
 * needed.
 */
char *seq_get_seqstr(Seq *seq) {
  long i;

  char *seqstr;
  seqstr = my_new(char, seq->len+1);

  for(i = 0; i < seq->len; i++) {
    seqstr[i] = nuc_id_to_char(seq->sym[i]);
  }
  seqstr[seq->len] = '\0';

  return seqstr;
}


/*
 * Fills the provided buffer with a null-terminated string
 * representation of the nucleotides in this sequence. The buffer must
 * be of length seq->len + 1 or greater.
 */
char *seq_get_seqstr_buf(Seq *seq, char *buf) {
  long i;

  for(i = 0; i < seq->len; i++) {
    buf[i] = nuc_id_to_char(seq->sym[i]);
  }
  buf[seq->len] = '\0';

  return buf;
}





/*
 * Does an in-place complement (but not reverse) of the provided sequence.
 * Does not alter the coordinates of the sequence.
 */
void seq_comp(Seq *seq) {
  long i;
  for(i = 0; i < seq->len; i++) {
    seq->sym[i] = nuc_comp(seq->sym[i]);
  }
}


/**
 * Performs an in-place of reversal (but not complement!) of the
 * provided sequence. Does not alter the coordinates of the sequence.
 */
void seq_rev(Seq *fwd) {
  util_breverse(fwd->sym, fwd->len);
}



/**
 * Performs an in-place reverse-compliment of the provided sequence.
 * and flips the strand of the sequence coordinate.
 */
void seq_revcomp(Seq *fwd_seq) {
  seq_rev(fwd_seq);
  seq_comp(fwd_seq);
}


/**
 * Creates a new copy of a sequence and returns it. The sequence
 * should be freed when it is no longer needed.
 */
Seq *seq_dup(Seq *seq) {
  Seq *new_seq;

  new_seq = my_new(Seq, 1);
  new_seq->len = seq->len;

  strncpy(new_seq->name, seq->name, SEQ_MAX_NAME_SZ);

  new_seq->sym = my_new(unsigned char, seq->len);
  new_seq->sym = memcpy(new_seq->sym, seq->sym, seq->len);
  new_seq->buf_sz = seq->len;

  return new_seq;
}


/**
 * Frees an allocated sequence structure. 
 */  
void seq_free(Seq *seq) {
  my_free(seq->sym);
  my_free(seq);
}


/**
 * Frees an array of sequence structures. 
 */  
void seq_array_free(Seq *seqs, long num) {
  long i;

  for(i = 0; i < num; i++) {
    my_free(seqs[i].sym);
  }
  my_free(seqs);
}

