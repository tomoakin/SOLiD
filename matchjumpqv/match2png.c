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
#include <math.h>
#include "matchlocfile.h"
#include "basecache.h"
#include "colorbase.h"
#include "seqmap.h"
#include "gd.h"
#include "gdfonts.h"
#include "gdfontt.h"
#include "matchlist.h"
#include "view.h"
 
int getimagexcoord(int64_t location, view* curview)
{
  int t;
  if(curview->flipimage == 0){
    if(location < curview->view_start_addr)
      return 0;
    t = (int)lround((location - curview->view_start_addr) * curview->xconversion_const); 
    if(t > curview->image_width)
      return curview->image_width;
    return t;
  }else{
    if(location > curview->view_end_addr)
      return 0;
    t = (int)lround((curview->view_end_addr - location + 1) * curview->xconversion_const); 
    if(t > curview->image_width)
      return curview->image_width;
    return t;
  }
}
void
drawgenomecolor(gdImagePtr gdimg, cache_st *cache, int nthline, gdFontPtr font, int fontcolor, int gapcolor, int *colorcolor, view *curview)
{
  char* colortable="0123....";
  int x1=0, x2;
  int y1, y2;
  int i;
  int color;
  int centershift = 0;
  centershift = lround((curview->xconversion_const - font->w)/2);
  y1 = lround(1 + nthline * curview->yconversion_const);
  y2 = lround((nthline + 1) * curview->yconversion_const);
#if 0
  fprintf(stderr, "curview_start_addr: %ld, view_end_addr: %ld", curview->view_start_addr, curview->view_end_addr);
#endif
  if(curview->flipimage == 0){
    for(i = curview->view_start_addr; i <=  curview->view_end_addr; i++){
      color = cachegetcolor(cache, i);
      x1 = getimagexcoord(i, curview);
      x2 = getimagexcoord(i + 1, curview);
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
    }
  }else{
    for(i = curview->view_end_addr; i >=  curview->view_start_addr; i--){
      color = cachegetcolor(cache, i);
      x1 = getimagexcoord(i + 1, curview);
      x2 = getimagexcoord(i, curview);
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
  }
}
void
drawflipelementswithcolor(gdImagePtr gdimg, int nthline, element* cur_el, gdFontPtr font, int fontcolor, int gapcolor, int *colorcolor, view* curview)
{
  int x1;
  int y1, y2, ym;
  int colorarray[50];
  char colorchararray[51];
  char* colortable="0123....";
  unsigned char linebuf[10240];
  int color;
  int centershift = 0;
  int offsetcorrection = 2;
  readtag *tagp;
  if(cur_el == NULL) return;
  if(cur_el->location < 0){
    offsetcorrection = 0;
  }
  tagp = cur_el->tagp;
  y1 = lround(1 + nthline * curview->yconversion_const);
  y2 = lround((nthline + 1) * curview->yconversion_const);
  tagdata50tointarray(colorarray, tagp->seq);
  tagdata50tostr(colorchararray, tagp->seq);
  colorchararray[50]='\0';
  sprintf(linebuf, ":%d:%i:%" PRId32, tagp->nhitpos, tagp->score, tagp->nextpos);
  centershift = lround((curview->xconversion_const - font->w)/2);
  if(centershift < 0) centershift = 0;
  {
    int x2;
    x2 = getimagexcoord(cur_el->location + tagp->adapter_pos - 1 +offsetcorrection + tagp->gaplen, curview);
    x1 = getimagexcoord(cur_el->location + 49 + offsetcorrection + tagp->gaplen, curview);
    gdImageFilledRectangle(gdimg, x1, y1, x2, y2, gapcolor);
  }
  x1 = getimagexcoord(cur_el->location + offsetcorrection, curview);

  /*write remaining colors*/
  if(tagp->gappos == 0){
    int i;
    int x2;
    for(i = 1; i <=  tagp->alignedlen; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + offsetcorrection + i , curview);
      gdImageFilledRectangle(gdimg, x2, y1, x1, y2, colorcolor[color]);
#if 0
      fprintf(stderr, "x1:%i, y1:%i, x2:%i, y2:%i\n", x2, y1, x1, y2);
#endif
      gdImageChar(gdimg, font, x2 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
    for(; i < 50; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + offsetcorrection + i, curview);
      gdImageChar(gdimg, font, x2 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
  }else{
    int i;
    int x2;
    ym = y1 + y2;
    ym /= 2;
    for(i = 1; i <  tagp->gappos; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + offsetcorrection + i, curview);
      gdImageFilledRectangle(gdimg, x2, y1, x1, y2, colorcolor[color]);
      gdImageChar(gdimg, font, x2 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
    x2=getimagexcoord(cur_el->location + offsetcorrection - 1 + tagp->gappos + tagp->gaplen, curview);
    gdImageLine(gdimg, x1, ym, x2, ym, gapcolor);
    x1 = x2;
    for(;i <= tagp->alignedlen; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + tagp->gaplen + offsetcorrection + i, curview);
      gdImageFilledRectangle(gdimg, x2, y1, x1, y2, colorcolor[color]);
#if 0
      fprintf(stderr, "i:%i, x1:%i, y1:%i, x2:%i, y2:%i\n", i, x2, y1, x1, y2);
#endif
      gdImageChar(gdimg, font, x2 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    } 
    for(; i < 50; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + tagp->gaplen + offsetcorrection + i, curview);
      gdImageChar(gdimg, font, x2 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
  }
  /*write the first *base* */
  {
    int firstbase = getnextbase(3, colorarray[0]);
    x1 = getimagexcoord(cur_el->location + offsetcorrection, curview);
    gdImageChar(gdimg, font, x1 + centershift, y1, nuc2char[firstbase], fontcolor);
  }
  x1 = getimagexcoord(cur_el->location + offsetcorrection - 1, curview);
  gdImageString(gdimg, font, x1 + centershift, y1, linebuf, fontcolor);
  drawflipelementswithcolor(gdimg, nthline, cur_el->next, font, fontcolor, gapcolor, colorcolor, curview);
}
void
drawelementswithcolor(gdImagePtr gdimg, int nthline, element* cur_el, gdFontPtr font, int fontcolor, int gapcolor, int *colorcolor, view* curview)
{
  int x1;
  int y1, y2, ym;
  int colorarray[50];
  char colorchararray[51];
  char* colortable="0123....";
  unsigned char linebuf[10240];
  int color;
  int centershift = 0;
  int offsetcorrection = 1;
  readtag *tagp;
  if(cur_el == NULL) return;
  drawelementswithcolor(gdimg, nthline, cur_el->next, font, fontcolor, gapcolor, colorcolor, curview);
  if(cur_el->location < 0){
    offsetcorrection = -1;
  }
  tagp = cur_el->tagp;
  y1 = lround(1 + nthline * curview->yconversion_const);
  y2 = lround((nthline + 1) * curview->yconversion_const);
  ym = y1 + y2;
  ym /= 2;
  tagdata50tointarray(colorarray, tagp->seq);
  tagdata50tostr(colorchararray, tagp->seq);
  colorchararray[50]='\0';
  sprintf(linebuf, ":%d:%i:%" PRId32 ":%d", tagp->nhitpos, tagp->score, tagp->nextpos, tagp->adapter_pos);
  x1 = getimagexcoord(cur_el->location + offsetcorrection, curview);
  centershift = lround((curview->xconversion_const - font->w)/2);
  if(centershift < 0) centershift = 0;
  /*write first *base* */
  {
    int firstbase = getnextbase(3, colorarray[0]);
    gdImageChar(gdimg, font, x1 + centershift, y1, nuc2char[firstbase], fontcolor);
  }

  {
    int x2;
    x1 = getimagexcoord(cur_el->location + tagp->adapter_pos + offsetcorrection + tagp->gaplen, curview);
    x2 = getimagexcoord(cur_el->location + 49 + 1 + offsetcorrection + tagp->gaplen, curview);
    gdImageFilledRectangle(gdimg, x1, y1, x2, y2, gapcolor);
  }

  x1 = getimagexcoord(cur_el->location + 1 + offsetcorrection, curview);
  /*write remaining colors*/
  if(tagp->gappos == 0){
    int x2;
    int i;
    x1 = getimagexcoord(cur_el->location + 1 + offsetcorrection, curview);
    for(i = 1; i <=  tagp->alignedlen; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + i + 1 + offsetcorrection, curview);
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
    for(; i < 50; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + i + 1 + offsetcorrection, curview);
#if 0
/*may add out-of-alignment color */
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
#endif
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
  }else{
    int i;
    int x2;
    for(i = 1; i <  tagp->gappos; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + i + 1 + offsetcorrection, curview);
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
    x2=getimagexcoord(cur_el->location + offsetcorrection + tagp->gappos + tagp->gaplen, curview);
    gdImageLine(gdimg, x1, ym, x2, ym, gapcolor);

    x1 = getimagexcoord(cur_el->location + offsetcorrection + tagp->gappos + tagp->gaplen, curview);
    i = tagp->gappos;
    for(;i <= tagp->alignedlen; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + tagp->gaplen + i + 1 + offsetcorrection, curview);
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    } 
    for(; i < 50; i++){
      color = colorarray[i];
      if(color >3) color = 4;
      x2 = getimagexcoord(cur_el->location + tagp->gaplen + i + 1 + offsetcorrection, curview);
#if 0
/*may add out-of-alignment color */
      gdImageFilledRectangle(gdimg, x1, y1, x2, y2, colorcolor[color]);
#endif
      gdImageChar(gdimg, font, x1 + centershift, y1, colortable[color], fontcolor);
      x1 = x2;
    }
  }
  gdImageString(gdimg, font, x1 + centershift, y1, linebuf, fontcolor);
}

void
drawelements(gdImagePtr gdimg, int hight, element* cur_el, int matchcolor, int gapcolor, int mismatchcolor, view* curview)
{
  int x1;
  int offsetcorrection = 0;
  readtag *tagp;
  if(cur_el == NULL) return;
  if(cur_el->location < 0){
    offsetcorrection = -2;
  }
  tagp = cur_el->tagp;
  x1 = getimagexcoord(cur_el->location + offsetcorrection, curview);
  if(tagp->gappos == 0){
    int x2;
    x2 = getimagexcoord(cur_el->location + tagp->alignedlen + 1 + offsetcorrection, curview);
    gdImageLine(gdimg, x1, hight, x2, hight, matchcolor);
  }else{
    int x2;
    x2=getimagexcoord(cur_el->location + tagp->gappos + 1 + offsetcorrection, curview);
    gdImageLine(gdimg, x1, hight, x2, hight, matchcolor);
    x1 = x2;
    x2=getimagexcoord(cur_el->location + tagp->gappos + tagp->gaplen + 1 + offsetcorrection, curview);
    gdImageLine(gdimg, x1, hight, x2, hight, gapcolor);
    x1 = x2;
    x2 = getimagexcoord(cur_el->location + tagp->alignedlen + tagp->gaplen + 1 + offsetcorrection, curview);
    gdImageLine(gdimg, x1, hight, x2, hight, matchcolor);
  }
  drawelements(gdimg, hight, cur_el->next, matchcolor, gapcolor, mismatchcolor, curview);
}

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


static void usage()
{
  fputs("matchcache2png -b matchcache -c genomeseqcache -m seqmapfile -t target -s start -e end -w width [-o outfile] [-S strand] [-u unit]\n", stderr);
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
  view curview;
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
  if(!(width >0)){
    fputs("width must be specified as a positive number!\n", stderr);
    usage();
    exit(EXIT_FAILURE);
  }
  if(!unit){
    unit = lround((fabs(req_end_addr - req_start_addr) + 1)/width);
    if(unit == 0) unit = 1;
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
  curview.scaffold=target_name;
  curview.view_start_scaffold_addr = req_start_addr;
  curview.view_end_scaffold_addr = req_end_addr;
  curview.image_width = width;

  if(req_end_addr>0){
    start_addr += scaff_info.offset;
    if(req_end_addr < scaff_info.length){
      end_addr += scaff_info.offset;
    }else{
      end_addr = scaff_info.offset + scaff_info.length;
    }
  }else{
    start_addr -= scaff_info.offset;
    end_addr -= scaff_info.offset;
  }
  if(start_addr < end_addr){
    curview.view_start_addr = start_addr;
    curview.view_end_addr = end_addr;
    curview.flipimage = 0;
    flip_image = 0;
  }else{
    curview.view_start_addr = end_addr;
    curview.view_end_addr = start_addr;
    curview.flipimage = 1;
    flip_image = 1;
    start_addr = curview.view_start_addr;
    end_addr = curview.view_end_addr;
  }
  if(!flip_image){
    curview.xconversion_const = width / (double)((req_end_addr - req_start_addr + 1 ));
  }else{
    curview.xconversion_const = width / (double)((req_start_addr - req_end_addr + 1 ));
  }
  
#if 0
  fprintf(stderr, "%s:%i --- %i\n", target_name, start_addr, end_addr);
#endif
  
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
         "readtag_count = %" PRIu64 "\n"
         "byte_order = %s\n",
         &fpagesize, &genomesize, &arraysize, &upperbitmask, &lowerbitmask, &shiftcount,  &cursorindex_size, &positionlevelarray_start, &positionlevelarray_size,
         &readtagarray_start, &readtag_count, byte_order);

  cursorindex = mmap(NULL, roundinps(cursorindex_size, pagesize), PROT_READ, MAP_PRIVATE, bincachefd, pagesize);
  positionlevelarray = mmap(NULL, roundinps(positionlevelarray_size, pagesize), PROT_READ, MAP_PRIVATE, bincachefd, positionlevelarray_start);
  readtagmmapsize = roundinps(readtag_count * sizeof(readtag), pagesize);
  readtagarray = mmap(NULL, readtagmmapsize, PROT_READ, MAP_PRIVATE, bincachefd, readtagarray_start);
  if(unit > 10){
    fprintf(stderr, "unit = %i\n", unit);
    out_density_array = calloc((end_addr - start_addr + 1)/unit + 1, sizeof(double));
    {
      int i;
      if(unit < 30){
        for(i = ((start_addr&upperbitmask)>>shiftcount); i < ((end_addr&upperbitmask)>>shiftcount); i++){
          posdata pos;
          int j,k;
          if(cursorindex[i].count == 0) continue;
          for(j = cursorindex[i].cursor; j <= cursorindex[i+1].cursor; j++){
            pos = positionlevelarray[j];
            if(pos.loc < start_addr) continue;
            if(pos.loc > end_addr) break;
            for(k=0; k < pos.count; k++){
              readtag * t =readtagarray + pos.read_tag_head + k;
              out_density_array[(pos.loc-start_addr+1)/unit] += 1.0/t->nhitpos;
            }
          }
        }
      }else{ /* more rough view; short cut to check for multihit correction*/
        for(i = ((start_addr&upperbitmask)>>shiftcount); i < ((end_addr&upperbitmask)>>shiftcount); i++){
          posdata pos;
          int j,k;
          if(cursorindex[i].count == 0) continue;
          for(j = cursorindex[i].cursor; j <= cursorindex[i+1].cursor; j++){
            pos = positionlevelarray[j];
            if(pos.loc < start_addr) continue;
            if(pos.loc > end_addr) break;
            out_density_array[(pos.loc-start_addr+1)/unit] += pos.count;
          }
        }
      }
    }
    {
      int i;
      gdImagePtr im;
      int black, gray, gray50, white, red;
      int imagehight = 100;
      double hconversionconst, wconversionconst;
      double min=DBL_MAX; 
      double max=0;
      for(i=0; i < (end_addr - start_addr + 1)/unit; i++){
        if(min > out_density_array[i]) min = out_density_array[i];
        if(max < out_density_array[i]) max = out_density_array[i];
       /*fprintf(stderr, "%i: %f\n", i, out_density_array[i]);*/
      }
      if(max < 1) logtransform =0;
      if(logtransform){
        min = -1;
        for(i=0; i < (end_addr - start_addr + 1)/unit; i++){
          if(out_density_array[i] < 1){ 
            out_density_array[i]=min;
          }else{
            out_density_array[i]=log10(out_density_array[i]);
          }
        }
        max = log10(max);
      }
      im = gdImageCreate(width, imagehight);
      black = gdImageColorAllocate(im, 0, 0, 0);
      gray = gdImageColorAllocate(im, 170, 170, 170);
      gray50 = gdImageColorAllocate(im, 127, 127, 127);
      white = gdImageColorAllocate(im, 255, 255, 255);
      red = gdImageColorAllocate(im, 255, 0, 0);
      gdImageFilledRectangle(im, 0,0, width, imagehight, white);

      if(max > 0){
        int realimagewidth;
        int imagestart = 0;
        double scale_step;
        hconversionconst = (imagehight-1) / (max - min);
        wconversionconst = width / (double)((fabs(req_end_addr - req_start_addr) + 1 ) / unit);
        realimagewidth = (double)(end_addr - start_addr+1)/(double)(req_end_addr-req_start_addr+1)*width;
        if(flip_image){
          realimagewidth = (double)(end_addr - start_addr+1)/(double)(req_start_addr-req_end_addr+1)*width;
          imagestart = width - realimagewidth;
        }

        scale_step = pow(10,floor(log10(max-min)));/* though pow(10,x) may be more expensive than exp10(x), it should be portable and no problem for one call*/
        if(max/scale_step < 2)
           scale_step /=2;
        if(max/scale_step > 5)
           scale_step *=2;
        {
          double di;
          gdFontPtr scalefont = gdFontGetSmall();
          di = 1.0 * unit * totaltagcount / Nmappablepositions; 
          if(max > di){
            int yi;
            yi= imagehight -1 - lround((di-min) * hconversionconst);
            gdImageLine(im, imagestart, yi, imagestart + realimagewidth, yi, red);
          }
          for(di=0.0; di < max; di += scale_step){
            int yi;
            yi= imagehight -1 - lround((di-min) * hconversionconst);
            gdImageLine(im, imagestart, yi, imagestart + realimagewidth, yi, gray);
          }
        }

        for(i=0; i < (end_addr - start_addr + 1)/unit; i++){
          int x1,x2,y2;
          x1 = floor(i * wconversionconst);
          x2 = ceil((i+1)* wconversionconst - 1 -0.000000001);
          if(flip_image){
            int tmp1, tmp2;
            tmp1 = width - x2;
            tmp2 = width - x1;
            x1 = tmp1;
            x2 = tmp2;
          }
          y2 = imagehight - 1 - lround((out_density_array[i] - min) * hconversionconst);
          gdImageFilledRectangle(im, x1, y2, x2, imagehight - 1, black);
        }
        {
          double di;
          gdFontPtr scalefont = gdFontGetSmall();
          unsigned char text[1000];
          for(di=0.0; di < max; di += scale_step){
            int yi;
            yi= imagehight -1 - lround((di-min) * hconversionconst);
            if(logtransform){
              sprintf(text, "%g", pow(10,di));
            }else{
              sprintf(text, "%g/%g", di, di/totaltagcount * Nmappablepositions / unit);
            }
            gdImageString(im, scalefont, imagestart, yi, text, gray50);
          }
        }
      }
      gdImagePng(im, outfile);
    }
  }else{
    int totaltagtobedisplayed = 0;
    int i;
    line *lines=NULL;
    int32_t scanstart = start_addr - 2000;
    if(start_addr > 0){
      if(scanstart< 0) scanstart = 1;
    }
    for(i = (scanstart&upperbitmask)>>shiftcount; i < ((end_addr&upperbitmask)>>shiftcount); i++){
      posdata pos;
      int j,k;
      if(cursorindex[i].count == 0) continue;
      for(j = cursorindex[i].cursor; j <= cursorindex[i+1].cursor; j++){
        pos = positionlevelarray[j];
        if(pos.loc > end_addr) break;
        for(k=0; k < pos.count; k++){
          readtag * t =readtagarray + pos.read_tag_head + k;
          int64_t endlocation = pos.loc;
          endlocation += t->alignedlen + t->gaplen;
          if(endlocation > start_addr){
            totaltagtobedisplayed +=1;
            lines = insertnewalignment(lines, t, pos.loc, 20);
          }
        }
      }
    }
    if(lines){
      gdImagePtr im;
      int nlines;
      int i = 0;
      int black, gray, gray50, white, red;
      int imagehight = 100;
      nlines = count_lines(lines);
#if 0
      fprintf(stderr, "%i tags in %i lines", totaltagtobedisplayed, nlines);
#endif
      if(curview.xconversion_const < 5 || nlines > 200){ /* unit == 1 */
        imagehight = nlines + 2;
        curview.image_hight = imagehight;
        curview.image_width = width;
        curview.unit = unit;
        curview.yconversion_const = 1.0;
        im = gdImageCreate(width, imagehight);
        black = gdImageColorAllocate(im, 0, 0, 0);
        gray = gdImageColorAllocate(im, 170, 170, 170);
        gray50 = gdImageColorAllocate(im, 127, 127, 127);/* this is darker than gray)*/
        white = gdImageColorAllocate(im, 255, 255, 255);
        red = gdImageColorAllocate(im, 255, 0, 0);
        gdImageFilledRectangle(im, 0,0, width, imagehight, white);
        while(lines){
          drawelements(im,i, lines->head, black, gray, red, &curview); 
          i++;
          lines = lines->next;
        }
        gdImagePng(im, outfile);
      }else{
        int linehight = 10;
        int blue, green, yellow, red;
        int color[5];
        gdFontPtr codefont = gdFontGetTiny();
        linehight = codefont->h;
        if(nlines *linehight < 100){
          linehight = 100/nlines;
        }
#if 0
        fprintf(stderr, "%i tags in %i lines", totaltagtobedisplayed, nlines);
#endif
        imagehight = linehight * (nlines+2) + 2;
        curview.image_hight = imagehight;
        curview.image_width = width;
        curview.unit = unit;
        curview.yconversion_const = linehight;
        im = gdImageCreate(width, imagehight);
        black = gdImageColorAllocate(im, 0, 0, 0);
        gray = gdImageColorAllocate(im, 170, 170, 170);

        blue = gdImageColorAllocate(im, 128, 255, 255);
        green = gdImageColorAllocate(im, 0, 255, 0);
        yellow =  gdImageColorAllocate(im, 255, 255, 0);
        red = gdImageColorAllocate(im, 255, 0, 0);
        color[0] = blue;
        color[1] = green;
        color[2] = yellow;
        color[3] = red;
        color[4] = gray;
        white = gdImageColorAllocate(im, 255, 255, 255);
        gdImageFilledRectangle(im, 0,0, width, imagehight, white);
        drawgenomecolor(im, &cache, i, codefont, black, gray, color, &curview);
        i++;
        if(curview.flipimage == 0){
          while(lines){
            drawelementswithcolor(im, i, lines->head, codefont, black, gray, color, &curview); 
            i++;
            lines = lines->next;
          }
        }else{
          while(lines){
            drawflipelementswithcolor(im, i, lines->head, codefont, black, gray, color, &curview); 
            i++;
            lines = lines->next;
          }
        }
        drawgenomecolor(im, &cache, i, codefont, black, gray, color, &curview);
        gdImagePng(im, outfile);
      }
    }else{ /* no line to draw */
      gdImagePtr im;
      gdFontPtr codefont = gdFontGetTiny();
      int linehight = codefont->h;
      int imagehight = linehight + 2;
      int blue, green, yellow, red;
      int color[5];
      int black, gray, white;
      curview.image_hight = imagehight;
      curview.image_width = width;
      curview.yconversion_const = linehight;
      im = gdImageCreate(width, imagehight);
      black = gdImageColorAllocate(im, 0, 0, 0);
      gray = gdImageColorAllocate(im, 170, 170, 170);
      blue = gdImageColorAllocate(im, 128, 255, 255);
      green = gdImageColorAllocate(im, 0, 255, 0);
      yellow =  gdImageColorAllocate(im, 255, 255, 0);
      red = gdImageColorAllocate(im, 255, 0, 0);
      color[0] = blue;
      color[1] = green;
      color[2] = yellow;
      color[3] = red;
      color[4] = gray;
      white = gdImageColorAllocate(im, 255, 255, 255);
      gdImageFilledRectangle(im, 0,0, width, imagehight, white);
      drawgenomecolor(im, &cache, 0, codefont, black, gray, color, &curview);
      gdImagePng(im, outfile);
    }
  }
  return 0;
}

