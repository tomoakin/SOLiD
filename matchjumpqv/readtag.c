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
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/time.h>
#include <sys/resource.h>
#include "seqmap.h"
#include "matchlocations.h"
#include "matchlocfile.h"

static int max_match_points = 1100;
static int seqmap_list_size = 2500;
static int nchr_max = 64;
/*
extern char*optarg;
extern int optind;
*/
static void usage(void);

/*
static inline int printmatch(FILE*, match);
static inline match str2match(const char*);
*/

inline tag50alignlocation
str2tag50alignlocation(const char *str)
{
  tag50alignlocation ret;
  char *p = NULL;
  int64_t num;
  ret.location = 0;
  num = strtoll(str, &p, 10);
  ret.location = num;
  num = strtoll(p+1, &p, 10);
  ret.tag.score=num;
  num = strtoll(p+1, &p, 10);
  ret.tag.gappos=num;
  num = strtoll(p+1, &p, 10);
  ret.tag.gaplen=num;
  num = strtoll(p+1, &p, 10);
  ret.tag.alignedlen=num;
  num = strtoll(p+1, &p, 10);
  ret.tag.adapter_pos=num;
  ret.tag.nextpos=ret.location;
  ret.tag.nhitpos=0;
  ret.tag.seq.threebitseq[0]=0;
  ret.tag.seq.threebitseq[1]=1;
  ret.tag.seq.threebitseq[2]=2;
  ret.tag.seq.threebitseq[3]=3;
  ret.tag.seq.threebitseq[4]=4;
  return ret;
}

static inline int
printmatch(FILE*outfile, tag50alignlocation matchloc)
{
  const char *strand = "T";
  int64_t location = matchloc.location;
  if(matchloc.location < 0) {
    strand = "R";
    location = -matchloc.location;
  }
  return fprintf(outfile,"%s\t%li\n", strand, location);
}

static void usage()
{
  fputs("intali2icache -s seqmap [-o outfile] matchfiles ...\n", stderr);
}


int
main(int argc, char **argv)
{
  int optc;
  char*linebuf, *p;
  FILE *matchfile = stdin;
  int outfd = 0;
  FILE *cmapfile;
  FILE *seqmapfile = NULL;
  seqmap_r * seqlocarray;
  int seqlocarray_nelm;
  int nchr;
  int *chr_sizes;
  uint64_t total_seq_size;
  uint64_t totalmappedtag = 0;
  uint64_t totalmappedtag1 = 0;
  uint64_t totalnalign = 0;
  uint64_t totalnalign1 = 0;
  uint64_t totalmappedpos = 0;
  uint64_t linesprocessed = 0;
  struct rusage rused;
  const int linebuf_size= 1024000;
  matchtagindex*mti = NULL;
  while((optc = getopt(argc,argv,"c:o:s:z:")) != -1){
    switch(optc){
    case 'c':
        cmapfile = fopen(optarg, "r");
        if(cmapfile == NULL){
          perror(argv[0]);
          fputs("cmapfile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 's':
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
    case 'o':
        outfd = open(optarg, O_RDWR|O_CREAT|O_TRUNC, 0644);
        if(outfd < 0){
          perror(argv[0]);
          fputs("outfile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'r':
      break;
    case 'z':
      max_match_points = strtol(optarg, NULL, 0);
      break;
    default:
      usage();
      exit(EXIT_FAILURE);
    }
  }
  if(!outfd){
    usage();
    exit(EXIT_FAILURE);
  }
  if(!seqmapfile){
    usage();
    exit(EXIT_FAILURE);
  }
  
  linebuf = malloc(linebuf_size);
  if(linebuf == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }
  
  seqlocarray = malloc(sizeof(seqmap_r) * seqmap_list_size);
  if(seqlocarray == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }
  chr_sizes = malloc(sizeof(int)*nchr_max);
  if(chr_sizes == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }
  
  nchr = readseqmapfile(seqmapfile, linebuf, seqlocarray, chr_sizes, &seqlocarray_nelm, &seqmap_list_size, linebuf_size);

  {
    register int i = 0;
    register uint64_t t = 1;
    for(i=1; i<= nchr; i++){
      t += chr_sizes[i];
    }
    total_seq_size = t;
  }
  mti = initmatchtagindex(mti, total_seq_size, 6);

    fprintf(stderr,"%lu positions\n", total_seq_size);
  /* check the end of each chromosome */
  /* prepare hit density table; 2 double for each position */
  /* one is for top strand, the other for the reverse strand */
  /* ca. 8 G of memory for 500 Mb */
  /* float may be sufficient */

  /**/
  /* prepare chr_offset to scaff offset conversion table */
  /* this is not implemented right now, as may become unncessary */
  argv += optind;
  argc -= optind;
  {
    tag50alignlocation *matcharray = malloc(sizeof(tag50alignlocation)*max_match_points);
    int matchcount = 0;
    while(argc >0){
      getrusage(RUSAGE_SELF, & rused);
      fprintf(stderr, "user: %li s, sys: %li s\n", rused.ru_utime.tv_sec, rused.ru_stime.tv_sec);
#if 0
      fprintf(stderr, "maxrss: %li, ixrss: %li, idrss: %li, isrss: %li\n", 
               rused.ru_maxrss, rused.ru_ixrss,
               rused.ru_idrss, rused.ru_isrss);
#endif
      fprintf(stderr, "processing %s\n", argv[0]);
      matchfile = freopen(argv[0], "r", matchfile);
      argc--;
      argv++;
      while((p = fgets(linebuf, 1024000, matchfile))){
        linesprocessed += 1;
        if(*p == '>'){
          matchcount = 0;
          totalmappedtag1 += 1;
          char*matchstr;
          strtok(p,",");
          while((matchstr = strtok(NULL,","))){
            matcharray[matchcount] = str2tag50alignlocation(matchstr);
            matchcount ++;
            totalnalign1 += 1;
            if(matchcount > max_match_points) break;
          }
        }else if(*p == 'T'){
          if(matchcount == 0) continue;
          tagdata50 tag50 = str2tag50(p);
          totalmappedtag+=1;
          totalnalign += matchcount;
          {
            int i;
            for(i = 0; i < matchcount; i++){
              matcharray[i].tag.seq = tag50;
              matcharray[i].tag.nhitpos = matchcount;
              matcharray[i].tag.nextpos = matcharray[(i+1)%matchcount].location;
              inserttag(mti, matcharray[i].location, matcharray[i].tag, &totalmappedpos);
            }
          }
          matchcount = 0;
        }else{
#if 0
        fprintf(stderr, "unrecognized line at %" PRIu64 ": %s\n", linesprocessed, p);
#endif
        }
      }
    }
    fprintf(stdout, "%" PRIu64 " lines processed\n", linesprocessed);
    fprintf(stdout, "%" PRIu64 " tags were mapped to %" PRIu64 " alignments\n", totalmappedtag, totalnalign);
    fprintf(stdout, "%" PRIu64 " tags were mapped to %" PRIu64 " alignments\n", totalmappedtag1, totalnalign1);
    WriteTaglocToFile(outfd, mti, totalnalign, totalmappedtag, totalmappedpos);
  }
  return 0;
}

