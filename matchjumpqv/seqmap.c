#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"seqmap.h"

const int seqmap_nchr_max = 64;

seqmap_r
scan_seqmap_record(const char*str)
{
  seqmap_r ret;
  int i;
  uint64_t t=0;
  for(i=0; str[i] != '\0'; i++){
    if(str[i] == '\t') break;
  }
  ret.name = (char*)malloc(i+1);
  strncpy(ret.name, str, i);
  ret.name[i]='\0';
  i++;
  while(isdigit(str[i])){
    t = t*10 + (str[i] - '0');
    i++;
  }
  ret.pseudo_chr_number = t;
  t = 0;
  while(isspace(str[i])){
    i++;
  }
  while(isdigit(str[i])){
    t = t*10 + (str[i] - '0');
    i++;
  }
  ret.offset = t;
  t=0;
  while(isspace(str[i])){
    i++;
  }
  while(isdigit(str[i])){
    t = t*10 + (str[i] - '0');
    i++;
  }
  ret.length = t;
  return ret;
}

seqmap_r find_target_fromfile(FILE* seqmapfile, const char*target_name)
{
  seqmap_r ret;
  const int bufsize = 1024;
  char stringbuffer[bufsize];
  char *p;
  int len = strlen(target_name);
  ret.name = NULL;
  ret.offset = 0;
  ret.length =0;
  ret.pseudo_chr_number =0;
  while((p = fgets(stringbuffer, bufsize, seqmapfile))){
    if(0 == strncmp(p, target_name, len)){
      return scan_seqmap_record(p);
    }
  }
  return ret;
}
int 
readseqmapfile(FILE* seqmapfile, char * linebuf, seqmap_r *seqlocarray, int *chr_sizes, int*seqlocarray_nelm, int * seqmap_list_size, int linebufsize)
{
  register int i = 0;
  register int j = 0;
  int nchr = 0;
  char *p;
  while((p = fgets(linebuf, linebufsize, seqmapfile))){
    if(i >= *seqmap_list_size){
      seqmap_r * nseqlocarray;
      nseqlocarray = (seqmap_r*)realloc(seqlocarray, (sizeof(seqmap_r)) * *seqmap_list_size * 2);
      if(nseqlocarray){
        seqlocarray = nseqlocarray;
      }else{
        fputs("memory allocation failure while reading seqmap file\n", stderr);
        exit(EXIT_FAILURE);
      }
    }
    seqlocarray[i] = scan_seqmap_record(p);
    j = seqlocarray[i].pseudo_chr_number;
    if(j>nchr) nchr=j;
    chr_sizes[j] = seqlocarray[i].offset + seqlocarray[i].length;
    i++;
  }
  *seqlocarray_nelm = i;
  return nchr;
}

/* returns sequence name and 1 based coordinate */
refpos_t 
get_refpos(seqmap_r * seqlocarray, int seqlocarray_nelm, int location)
{
  int sign;
  int abslocation;
  int i;
  refpos_t retval = {NULL, 0};
  if(location <0){
    abslocation = -location;
    sign = -1;
  }else{
    abslocation = location;
    sign = 1;
  }
  for(i = 0 ; i < seqlocarray_nelm; i++){
    int64_t offset = seqlocarray[i].offset;
    if(offset < abslocation + 30  && abslocation -30 < offset +seqlocarray[i].length){
       retval.name = seqlocarray[i].name;
       retval.pos = sign * (abslocation + 1 - offset);
       return retval;
    }
  }
  return retval;
}

void write_sq_lines(FILE* outfile, seqmap_r *seqlocarray, int seqlocarray_nelm, char*addstr)
{
  int i;
  for(i = 0; i < seqlocarray_nelm; i++){
    fprintf(outfile, "@SQ\tSN:%s\tLN:%u%s\n", seqlocarray[i].name, seqlocarray[i].length, addstr);
  }
}
