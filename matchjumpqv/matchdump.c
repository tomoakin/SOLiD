/*
 * Copyright 2008, 2009 Tomoaki Nishiyama, Kanazawa University
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <float.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "matchlocfile.h"
#include "basecache.h"
#include "colorbase.h"
#include "seqmap.h"
#include "gd.h"
#include "gdfonts.h"
int printspace(FILE*out, int i)
{
  while(i>0){
    fputc(' ', out);
    i -=1;
  }
  return i;
}

static int seqmap_list_size = 2500;
static const int linebufsize = 1024000;
/*
extern char*optarg;
extern int optind;
*/
static void usage(void);

typedef struct match_ {
  int referencesequencenumber;
  int nmismatch;
  int64_t location;
} match;

#if 0
typedef struct tagdata50_{
  uint32_t threebitseq[5];/* color sequence with invalid state (.) */
} tagdata50;
typedef struct tag50align_{ /*static relocatable structure*/
  uint8_t alignedlen;
  uint8_t gappos;
  uint16_t score;
  uint16_t nhitpos; /* number of hit positions */
  uint16_t gaplen;
  genomeloc_t nextpos;
  tagdata50 seq;
} tag50align;
#endif


int
printtagdata(FILE* out, readtag* tag)
{
  fprintf(stdout, "aligned len: %i, score: %i nhitpos: %i, gappos: %i, gaplen: %i, nextpos:%d\n",
            tag->alignedlen, tag->score, tag->nhitpos, tag->gappos, tag->gaplen, tag-> nextpos);
  printtagdata50(stdout, tag->seq);
  fputc('\n',out);
  return 0;
}



static void usage()
{
  fputs("matchcachedump -b matchcache -c genomeseqcache -m seqmapfile -t target -s start -e end -w width [-o outfile] [-S strand] [-u unit]\n", stderr);
}

typedef struct  hit_density_pair_{
  double top;
  double reverse;
} hit_density_pair;


int Nmappablepositions = 453441362;
int
main(int argc, char **argv)
{
  int optc;
  int logtransform=0;
  int start_addr = 1, end_addr = 1000, unit=0, width=0;
  int ostart_addr = 1, oend_addr = 1000;
  int req_start_addr, req_end_addr;
  char * target_name = NULL;
  int scaff_length;
  char*linebuf;
  FILE *outfile = stdout;

  FILE *seqmapfile = NULL;
  seqmap_r * seqlocarray;
  seqmap_r  scaff_info; 
  int *chr_sizes;

  int strand = 3;
  double* out_density_array = NULL;
  char* cachefilename = NULL;
  int bincachefd = -1;
  uint64_t totaltagcount = 1000000; /* 1e6 */
  off_t scaff_start_offset;
  int flip_image=0;
  char*head;
  cursorcount *cursorindex;
  posdata *positionlevelarray;
  readtag *readtagarray;
  int64_t cursorindex_size;
  int64_t positionlevelarray_start; 
  int64_t positionlevelarray_size; 
  int64_t readtagarray_start;
  int64_t readtagarray_size;
  int64_t readtagmmapsize;
  uint64_t readtag_count = 0;
  uint64_t b[2] = {0,0};
  char*byte_order;
  int pagesize = getpagesize();
  int fpagesize;
  int cachefd = 0;
  uint64_t genomesize;
  uint64_t arraysize;
  uint64_t  upperbitmask, lowerbitmask;
  int shiftcount; 
  cache_st cache;

  byte_order = (char*)b;
  while((optc = getopt(argc,argv,"b:o:t:s:e:u:w:S:lp:c:m:")) != -1){
    switch(optc){
    case 'b':
        cachefilename = strdup(optarg);
        bincachefd = open(cachefilename, O_RDONLY ,644);
        if(bincachefd < 0){
          perror(argv[0]);
          fputs("bincachefile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'o':
        outfile = fopen(optarg, "w");
        if(outfile == NULL){
          perror(argv[0]);
          fputs("outfile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'c':
        cachefd = open(optarg, O_RDONLY, 0755);
        if(cachefd == -1){
          perror(argv[0]);
          fputs("cache file open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'm':
        seqmapfile = fopen(optarg, "r");
        if(seqmapfile == NULL){
          perror(argv[0]);
          fputs("seqmapfile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 't':  /*target name*/
      target_name = strdup(optarg);
      break;
    case 's': /*start*/
      req_start_addr = strtol(optarg, NULL, 0);
      start_addr = req_start_addr;
      break;
    case 'S': /*start*/
      strand = strtol(optarg, NULL, 0);
      break;
    case 'e': /*end*/
      req_end_addr = strtol(optarg, NULL, 0);
      end_addr = req_end_addr;
      break;
    case 'u': /*unit*/
      unit = strtol(optarg, NULL, 0);
      break;
    case 'w': /*width*/
      width = strtol(optarg, NULL, 0);
      break;
    case 'p': /*total number of locations 35-bp without N in the genome*/
      Nmappablepositions = strtol(optarg, NULL, 0);
      break;
    case 'l':
      logtransform = 1;
      break;
    default:
      usage();
      exit(EXIT_FAILURE);
    }
  }
  if(seqmapfile == NULL){
    fputs("seqmap file unspecified!\n", stderr);
    usage();
    exit(EXIT_FAILURE);
  }
  if(target_name == NULL){
    fputs("target name unspecified!\n", stderr);
    usage();
    exit(EXIT_FAILURE);
  }
  if(!cachefd){
    fputs("cache file unspecified!\n", stderr);
    usage();
    exit(EXIT_FAILURE);
  }
  assigncache(&cache, cachefd);
  initchar2nucmap();

  
  linebuf = malloc(linebufsize);
  if(linebuf == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }
  
  seqlocarray = malloc(sizeof(seqmap_r) * seqmap_list_size);
  if(seqlocarray == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }
  chr_sizes = malloc(sizeof(int)*seqmap_nchr_max);
  if(chr_sizes == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }

  scaff_info = find_target_fromfile(seqmapfile, target_name);
  if(req_end_addr>0){
    start_addr += scaff_info.offset;
    if(req_end_addr < scaff_info.length){
      end_addr += scaff_info.offset;
    }else{
      end_addr = scaff_info.offset + scaff_info.length;
    }
    ostart_addr = start_addr;
    oend_addr = end_addr;
  }else{
    start_addr -= scaff_info.offset;
    ostart_addr = start_addr + cache.reverse_base;
    end_addr -= scaff_info.offset;
    oend_addr = end_addr + cache.reverse_base;
  }
  
  head = mmap(NULL, pagesize, PROT_READ, MAP_PRIVATE, bincachefd, 0);
#if 0
 /* the file header format */
Gapped alignment storage file. Format version 1.
pagesize = 4096
genomesize = 480330409
arraysize = 16777216
upperbitmask = 3FFFFFC0
lowerbitmask = 3F
shiftbitcount = 6
cursorindex_size = 134217728
positinlevelarray_start = 134225920
positionlevelarray_size = 261128432
readtag_start = 395358208
readtag_count = 146618299
byte_order = 87654321
#endif
  sscanf(head, "Gapped alignment storage file. Format version 1.\n"
         "pagesize = %d" "\n"
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
         &fpagesize, &genomesize, &arraysize, &upperbitmask, &lowerbitmask, &shiftcount,  &cursorindex_size, &positionlevelarray_start, &positionlevelarray_size,
         &readtagarray_start, &readtag_count, byte_order);
  fprintf(stdout, "Gapped alignment storage file. Format version 1.\n"
         "pagesize = %d"  "\n"
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
         pagesize, genomesize, arraysize, upperbitmask, lowerbitmask, shiftcount, 
         cursorindex_size, positionlevelarray_start, positionlevelarray_size,
         readtagarray_start, readtag_count, byte_order);

  cursorindex = mmap(NULL, roundinps(cursorindex_size, pagesize), PROT_READ, MAP_PRIVATE, bincachefd, pagesize);
  positionlevelarray = mmap(NULL, roundinps(positionlevelarray_size, pagesize), PROT_READ, MAP_PRIVATE, bincachefd, positionlevelarray_start);
  readtagmmapsize = roundinps(readtag_count * sizeof(readtag), pagesize);
  readtagarray = mmap(NULL, readtagmmapsize, PROT_READ, MAP_PRIVATE, bincachefd, readtagarray_start);
  
  {
    int i;
    for(i = ((start_addr&upperbitmask)>>shiftcount); i < ((end_addr&upperbitmask)>>shiftcount); i++){
      posdata pos;
      int j,k;
      if(cursorindex[i].count == 0) continue;
      for(j = cursorindex[i].cursor; j <= cursorindex[i+1].cursor; j++){
        fprintf(stdout, "cursorindex[%i] pos: %d\n", i, j);
        pos = positionlevelarray[j];
        fprintf(stdout, "psary[%i]: loc: %d, count: %i, read_tag_head: %ld\n", j, pos.loc, pos.count, pos.read_tag_head);
        if(pos.loc < start_addr) continue;
        if(pos.loc > end_addr) break;
        for(k=0; k < pos.count; k++){
          readtag * t =readtagarray + pos.read_tag_head + k;
          if(t->gappos >0){
            printseq(stdout, &cache, pos.loc, t->gappos + 3, 1);
            printcolorseq(stdout, &cache, pos.loc, t->gappos, 1);
          }else{
            printseq(stdout, &cache, pos.loc, 50, 1);
            printcolorseq(stdout, &cache, pos.loc, 50, 1);
          }
          printtagdata(stdout, t);
          if (t->gappos > 0){
            printspace(stdout, t->gappos - 3);
            printseq(stdout, &cache, pos.loc + t->gaplen +t-> gappos - 3, 50 - t->gappos + 3, 1);
            printspace(stdout, t->gappos);
            printcolorseq(stdout, &cache, pos.loc + t->gaplen + t->gappos, 50 - t->gappos, 1);
          }
        }
      }
    }
  }
  return 0;
}

