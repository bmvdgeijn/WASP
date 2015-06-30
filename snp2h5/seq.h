#ifndef __SEQ_H__
#define __SEQ_H__

#include <stdio.h>
#include <zlib.h>


#define SEQ_MAX_NAME_SZ 127
#define SEQ_NAME_BUF_SZ 128
#define SEQ_DEFAULT_BUF_SZ 65535
#define SEQ_FASTA_LINE_LEN 70

/* 
 * A DNA, RNA, or amino acid sequence 
 */
typedef struct {
  char name[SEQ_NAME_BUF_SZ];      /* name of the sequence    */
  unsigned char *sym;  /* array of nucleotide ids that make up the sequence */
  long len;            /* length of the sequence */
  
  size_t buf_sz; /* size of the sym buffer */
} Seq;


Seq *seq_new(void);
void seq_expand(Seq *seq);
void seq_read_str(Seq *seq, char *seq_str);
long seq_read_fasta_from_file(Seq *seq, const char *filename);
long seq_read_fasta_record(Seq *seq, gzFile f);
void seq_write_fasta_record(Seq *seq, gzFile f);
void seq_rev(Seq *fwd_seq);
void seq_comp(Seq *seq);
void seq_revcomp(Seq *fwd_seq);
void seq_nucs_revcomp(unsigned char *nuc_ids, long len);

char *seq_get_seqstr_buf(Seq *seq, char *buf);
char *seq_get_seqstr(Seq *seq);



Seq *seq_dup(Seq *seq);
void seq_free(Seq *seq);
void seq_array_free(Seq *seqs, long num);


#endif
