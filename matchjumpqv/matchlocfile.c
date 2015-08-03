
/* copyright 2009, Tomoaki Nishiyama, Kanazawa University */
/* matchlocfile.c */
/* ioroutines for cachefile of gapped alined tags */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <sys/mman.h>
#include <inttypes.h>
#include "matchlocfile.h"
int
WriteTaglocToFile(int outfd, matchtagindex*mti, uint64_t totalalignment, uint64_t totalmappedtag, uint64_t totalpos)
{
  size_t pagesize;
  char* head;
  cursorcount *cursorindex;
  uint64_t genomesize;
  posdata *positionlevelarray;
  readtag *readtagarray;
  uint32_t position_cursor = 0;
  uint64_t readtag_cursor = 0;
  int64_t cursorindex_start;
  int64_t positionlevelarray_start; 
  int64_t positionlevelarray_size; 
  int64_t readtagarray_start;
  int64_t readtagarray_size;
  int64_t cursorindex_size;
  int i,j;
  char *byte_order;
  uint64_t b[2] = {0,0};
  b[0] = ((((((((((((((0LL + ('1' <<8)) + '2')<< 8) + '3')<<8)+'4')
         <<8)+'5')<<8)+'6')<<8)+'7')<<8)+'8'); 
  byte_order = (char*)b;

  genomesize = mti->genomesize;
  pagesize=getpagesize();
  head = mmap(NULL, pagesize, PROT_WRITE, MAP_SHARED, outfd, 0);
  if(head ==MAP_FAILED){
    perror("head");
    exit(EXIT_FAILURE);
  }
  /* construct fixed structure */
  cursorindex_start = pagesize;
  cursorindex_size =sizeof(cursorcount)*mti->arraysize;
  cursorindex = mmap(NULL, roundinps(cursorindex_size, pagesize), PROT_WRITE, MAP_SHARED, outfd, cursorindex_start);
  if(cursorindex ==MAP_FAILED){
    perror("cursorindex");
    exit(EXIT_FAILURE);
  }
  madvise(cursorindex, roundinps(cursorindex_size, pagesize), MADV_SEQUENTIAL);

  positionlevelarray_start = cursorindex_start + roundinps(cursorindex_size, pagesize);
  positionlevelarray_size = sizeof(posdata) * totalpos;
  positionlevelarray = mmap(NULL, roundinps(positionlevelarray_size, pagesize), PROT_WRITE, MAP_SHARED, outfd, positionlevelarray_start);
  if(positionlevelarray ==MAP_FAILED){
    perror("positionlevelarray");
    exit(EXIT_FAILURE);
  }
  madvise(positionlevelarray, roundinps(positionlevelarray_size, pagesize), MADV_SEQUENTIAL);

  readtagarray_size = sizeof(readtag) * totalalignment;
  readtagarray_start = positionlevelarray_start + roundinps(positionlevelarray_size, pagesize);
  readtagarray = mmap(NULL, roundinps(readtagarray_size, pagesize), PROT_WRITE, MAP_SHARED, outfd, readtagarray_start);
  if(readtagarray ==MAP_FAILED){
    perror("readtagarray");
    exit(EXIT_FAILURE);
  }
  madvise(readtagarray, roundinps(readtagarray_size, pagesize), MADV_SEQUENTIAL);

  ftruncate(outfd, readtagarray_start + roundinps(readtagarray_size, pagesize));
  sprintf(head, "Gapped alignment storage file. Format version 1.\n"
         "pagesize = %" PRId64 "\n"
         "genomesize = %" PRIu64 "\n"
         "arraysize = %" PRIu64 "\n"
         "upperbitmask = %" PRIX64 "\n"
         "lowerbitmask = %" PRIX64 "\n"
         "shiftbitcount = %i\n"
         "cursorindex_size = %" PRId64 "\n"
         "positinlevelarray_start = %" PRId64 "\n"
         "positionlevelarray_size = %" PRId64 "\n"
         "readtag_start = %" PRId64 "\n"
         "readtag_count = %" PRId64 "\n"
         "byte_order = %s\n",
         pagesize, mti->genomesize,mti->arraysize, mti->upperbitmask, mti->lowerbitmask, 
         mti->shiftbitcount,
         cursorindex_size, positionlevelarray_start, positionlevelarray_size,
         readtagarray_start, totalalignment, byte_order);
  /*  mti.each do */
  for(i = 0; i < mti->arraysize; i++){
    matchtaglist *curmatchtaglist = mti->pary[i];
    cursorindex[i].cursor = position_cursor;
    cursorindex[i].count = 0;
    while(curmatchtaglist){
      genomeloc_t loc = curmatchtaglist->location;
      positionlevelarray[position_cursor].loc = loc;
      positionlevelarray[position_cursor].count = curmatchtaglist->array_content_number;
      positionlevelarray[position_cursor].read_tag_head = readtag_cursor;
      for(j = 0; j < curmatchtaglist->array_content_number; j++){
//        fprintf(stderr, "readtag_cursor = %li, j = %d\n", readtag_cursor, j);
        readtagarray[readtag_cursor] = curmatchtaglist->alignedtag[j];
        readtag_cursor += 1;
        cursorindex[i].count += 1000.0 / curmatchtaglist->alignedtag[j].nhitpos;
      }
      curmatchtaglist = curmatchtaglist -> next;
      position_cursor += 1;
    }
  }
  cursorindex_start = pagesize;
  cursorindex_size =sizeof(cursorcount)*mti->arraysize;
  munmap(head, pagesize);
  munmap(cursorindex, roundinps(cursorindex_size, pagesize));
  munmap(positionlevelarray, roundinps(positionlevelarray_size, pagesize));
  munmap(readtagarray, roundinps(readtagarray_size, pagesize));
  return 0;
}
