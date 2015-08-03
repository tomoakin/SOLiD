/* basecache.h */
/* defines structures and functions associated with base and color sequence cache*/
#if !defined(seqmap_h_included)
#define seqmap_h_included
#include<stdio.h>
#include<stdint.h>
typedef struct seqmap_r_ {
  char * name;
  uint64_t offset;
  unsigned int length;
  int pseudo_chr_number;
} seqmap_r;
extern const int seqmap_nchr_max;
typedef struct refpos_t_ {
  char *name;
  int32_t pos;
} refpos_t;

seqmap_r find_target_fromfile(FILE* seqmapfile, const char*target_name);
seqmap_r scan_seqmap_record(const char*str);
int readseqmapfile(FILE* seqmapfile, char * linebuf, seqmap_r *seqlocarray, int *chr_sizes, int*seqlocarray_nelm, int *seqlocarray_size, int linebufsize);

refpos_t get_refpos(seqmap_r*, int, int);
void write_sq_lines(FILE* outfile, seqmap_r *seqlocarray, int seqlocarray_nelm, char*addstr);

#endif
