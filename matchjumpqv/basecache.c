/*
 * Copyright 2008, 2009 Tomoaki Nishiyama, Kanazawa University
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h> /* for open()*/

#include "basecache.h"

void
assigncache(cache_st*cache, int cachefd)
{
  int pagesize;
  struct stat sb;
  char *p;
  off_t mmapsize;

  fstat(cachefd,&sb);
  pagesize = getpagesize();
  mmapsize = (sb.st_size / pagesize +1) * pagesize;
  p = (char*)mmap(NULL, pagesize, PROT_READ, MAP_SHARED, cachefd, 0);/*preamble*/
  sscanf(p, "%i %" PRIu64" %" PRIu64, &pagesize, &cache->totalnuc, &cache->reverse_base);
  munmap(p, pagesize);
  cache->data = (chunk32*)mmap(NULL, mmapsize, PROT_READ, MAP_SHARED, cachefd, pagesize);
}

int
printseq(FILE*outfile, cache_st*cache, int start, int length, int zerobase)
{
  int i = 0;
  if(start < 0){
    start += cache -> reverse_base;
    if(zerobase) start -=2;
  }
  {
    uint64_t selector;
    int shiftcount;
    if(zerobase){
      selector = (start)/32;
      shiftcount = (start)%32 *2;
    }else{
      selector = (start -1)/32;
      shiftcount = (start-1)%32 *2;
    }
    chunk32* chunkp = cache -> data + selector;
    uint64_t curchunk_v = chunkp->validity >> shiftcount;
    uint64_t curchunk_n = chunkp->nucseq >> shiftcount;
    while(i < length){
      if(curchunk_v & 2){
        int curnuc = curchunk_n & 3;
        fputc(nuc2char[curnuc], outfile);
      }else{
        fputc('N', outfile);
      }
      shiftcount += 2;
      curchunk_v >>= 2;
      curchunk_n >>= 2;
      if(shiftcount >= 64){
        shiftcount = 0;
        chunkp ++;
        curchunk_v = chunkp -> validity;
        curchunk_n = chunkp -> nucseq;
      }
      i ++;
    }
    fputc('\n', outfile);
  }
  return 0;
}

int
printcolorseq(FILE*outfile, cache_st*cache, int start, int length, int zerobase)
{
  int i = 0;
  if(start < 0){
    start += cache -> reverse_base;
    if(zerobase) start -=2;
  }
  {
    uint64_t selector;
    int shiftcount;
    if(zerobase){
      selector = (start)/32;
      shiftcount = (start)%32 *2;
    }else{
      selector = (start -1)/32;
      shiftcount = (start-1)%32 *2;
    }
    chunk32* chunkp = cache->data+selector;
    uint64_t curchunk_v = chunkp->validity >> shiftcount;
    uint64_t curchunk_c = chunkp->colorseq >> shiftcount;
    while(i < length){
      if(curchunk_v & 1){
        int curcolor = curchunk_c & 3;
        fputc('0' + curcolor, outfile);
      }else{
        fputc('*', outfile);
      }
      shiftcount += 2;
      curchunk_v >>= 2;
      curchunk_c >>= 2;
      if(shiftcount >= 64){
        shiftcount = 0;
        chunkp ++;
        curchunk_v = chunkp -> validity;
        curchunk_c = chunkp -> colorseq;
      }
      i ++;
    }
    fputc('\n', outfile);
  }
  return 0;
}


typedef struct seqmap_r_ {
  char * name;
  uint64_t offset;
  unsigned int length;
  int pseudo_chr_number;
} seqmap_r; 

void
printtarget16(FILE* out, target16*target)
{
  int i;
/*  int length;
  uint64_t validity[2];
  uint64_t ns[2];
  uint64_t cs[2];*/
  int shiftcount = 0;
  uint64_t curchunk_v = target->validity;
  uint64_t curchunk_n = target->ns;
  uint64_t curchunk_c;
  for(i = 0; i < target->length; i++){
    if(curchunk_v & 2){
      int curnuc = curchunk_n & 3;
      fputc(nuc2char[curnuc], out);
    }else{
      fputc('N', out);
    }
    shiftcount += 2;
    curchunk_v >>= 2;
    curchunk_n >>= 2;
  }
  fputc('\n', out);
  shiftcount = 0;
  curchunk_v = target->validity;
  curchunk_c = target->cs;
  for(i = 0; i < target->length; i++){
    if(curchunk_v & 1){
      int curcolor = curchunk_c & 3;
      fputc('0' + curcolor, out);
    }else{
      fputc('N', out);
    }
    shiftcount += 2;
    curchunk_v >>= 2;
    curchunk_c >>= 2;
  }
  fputc('\n', out);
}

