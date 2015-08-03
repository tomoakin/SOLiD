/*
 * matchextend: Seek to the position specified by the ma file
 * find heuristacally best exon-exon junction of the genome,
 * and show the alignment.
 * Copyright 2008, 2009 Tomoaki Nishiyama, Kanazawa University
 *
 *  The ma file may contain . as undetermined color.
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include <float.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <fcntl.h>

#include "basecache.h"
#include "colorbase.h"
#include "seqmap.h"
#include "adapter.h"

static int max_intron_size = 1500;
static int max_match_points = 1024;
static int seqmap_list_size = 2500;
static int nchr_max = 64;
static void usage(void);


typedef struct score_t_{
  float score; /* score is handled in floating point number for simplicity */
               /* this may be changed to integer type for speed optimization */
  int32_t intron_len; /* intron lengh */
  int8_t intron_pos; /* intron position */
  int8_t align_len;
  int8_t adapter_pos;
  int8_t adapter_align_pos;
} score_t;

target32*
gettarget32(cache_st*cache, int64_t location)
{
  static target32 ret_of_gettarget32;
  ret_of_gettarget32.length = 32;
  chunk32 *chunk = cache->data;
  
  if(location < 0){
    location += cache -> reverse_base - 2;
  }
  {
    int selector = (location)/32;
    int shiftcount = (location)%32 * 2;
    if(shiftcount){
      ret_of_gettarget32.validity[0] = chunk[selector].validity >> shiftcount;
      ret_of_gettarget32.ns[0] = chunk[selector].nucseq >> shiftcount;
      ret_of_gettarget32.cs[0] = chunk[selector].colorseq >> shiftcount;
      selector ++;
      ret_of_gettarget32.validity[0] |= chunk[selector].validity << (64-shiftcount);
      ret_of_gettarget32.ns[0] |= chunk[selector].nucseq << (64-shiftcount);
      ret_of_gettarget32.cs[0] |= chunk[selector].colorseq << (64-shiftcount);
    }else{
      ret_of_gettarget32.validity[0] = chunk[selector].validity;
      ret_of_gettarget32.ns[0] = chunk[selector].nucseq;
      ret_of_gettarget32.cs[0] = chunk[selector].colorseq;
      selector ++;
    }
    ret_of_gettarget32.cs[0] |= 3;
    ret_of_gettarget32.cs[0] ^= 3;
    ret_of_gettarget32.cs[0] |= ret_of_gettarget32.ns[0] & 3;
    return &ret_of_gettarget32;
  }
}
inline target16
gettarget16(cache_st*cache, int tagsize, int64_t location)
{
  target16 ret_of_gettarget16;
  ret_of_gettarget16.length = tagsize;
  chunk32 *chunk = cache->data;
  
  if(location < 0){
    location += cache -> reverse_base - 2;
  }
  {
    int selector = (location)/32;
    int shiftcount = (location)%32 * 2;
    if(shiftcount){
      ret_of_gettarget16.validity = chunk[selector].validity >> shiftcount;
      ret_of_gettarget16.ns = chunk[selector].nucseq >> shiftcount;
      ret_of_gettarget16.cs = chunk[selector].colorseq >> shiftcount;
      selector ++;
      ret_of_gettarget16.validity |= chunk[selector].validity << (64-shiftcount);
      ret_of_gettarget16.ns |= chunk[selector].nucseq << (64-shiftcount);
      ret_of_gettarget16.cs |= chunk[selector].colorseq << (64-shiftcount);
    }else{
      ret_of_gettarget16.validity = chunk[selector].validity;
      ret_of_gettarget16.ns = chunk[selector].nucseq;
      ret_of_gettarget16.cs = chunk[selector].colorseq;
      selector ++;
    }
    ret_of_gettarget16.cs |= 3;
    ret_of_gettarget16.cs ^= 3;
    ret_of_gettarget16.cs |= ret_of_gettarget16.ns & 3;
    return ret_of_gettarget16;
  }
}

inline target64
gettarget64(cache_st*cache, int tagsize, int64_t location)
{
  target64 ret_of_gettarget64;
  ret_of_gettarget64.length = tagsize;
  chunk32 *chunk = cache->data;
  
  if(location < 0){
    location += cache -> reverse_base - 2;
  }
  {
    int selector = (location)/32;
    int shiftcount = (location)%32 * 2;
    if(shiftcount){
      ret_of_gettarget64.validity[0] = chunk[selector].validity >> shiftcount;
      ret_of_gettarget64.ns[0] = chunk[selector].nucseq >> shiftcount;
      ret_of_gettarget64.cs[0] = chunk[selector].colorseq >> shiftcount;
      selector ++;
      ret_of_gettarget64.validity[0] |= chunk[selector].validity << (64-shiftcount);
      ret_of_gettarget64.ns[0] |= chunk[selector].nucseq << (64-shiftcount);
      ret_of_gettarget64.cs[0] |= chunk[selector].colorseq << (64-shiftcount);
      ret_of_gettarget64.validity[1] = chunk[selector].validity >> shiftcount;
      ret_of_gettarget64.ns[1] = chunk[selector].nucseq >> shiftcount;
      ret_of_gettarget64.cs[1] = chunk[selector].colorseq >> shiftcount;
      selector ++;
      ret_of_gettarget64.validity[1] |= chunk[selector].validity << (64-shiftcount);
      ret_of_gettarget64.ns[1] |= chunk[selector].nucseq << (64-shiftcount);
      ret_of_gettarget64.cs[1] |= chunk[selector].colorseq << (64-shiftcount);
    }else{
      ret_of_gettarget64.validity[0] = chunk[selector].validity;
      ret_of_gettarget64.ns[0] = chunk[selector].nucseq;
      ret_of_gettarget64.cs[0] = chunk[selector].colorseq;
      selector ++;
      ret_of_gettarget64.validity[1] = chunk[selector].validity;
      ret_of_gettarget64.ns[1] = chunk[selector].nucseq;
      ret_of_gettarget64.cs[1] = chunk[selector].colorseq;
    }
    ret_of_gettarget64.cs[0] |= 3;
    ret_of_gettarget64.cs[0] ^= 3;
    ret_of_gettarget64.cs[0] |= ret_of_gettarget64.ns[0] &3;
  }
  return ret_of_gettarget64;
}

static target64*
gettarget64_dest(target64* dest, cache_st*cache, int tagsize, int64_t location){
  target64 t;
  if(dest == NULL){
    dest = malloc(sizeof(target64));
  }
  t = gettarget64(cache, tagsize, location);
  dest->length = t.length;
  dest->validity[0] = t.validity[0];
  dest->validity[1] = t.validity[1];
  dest->ns[0] = t.ns[0];
  dest->ns[1] = t.ns[1];
  dest->cs[0] = t.cs[0];
  dest->cs[1] = t.cs[1];
  return dest;
}
static target64
recombinetarget64(target64 first, target64 second, int pos){
/* get pos bases from first, and remainder from the second */
  target64 ret;
  uint64_t mask = 0xFFFFFFFFFFFFFFFFULL;
  int shift, shiftr;
  ret.length=second.length;
  if(0 < pos && pos < 32){/* recombination in the first half */
    shiftr = 64 - pos * 2;
    ret.validity[0] = first.validity[0] & (mask >> shiftr);
    ret.ns[0] = first.ns[0] & (mask >> shiftr);
    ret.cs[0] = first.cs[0] & (mask >> shiftr);
    shift = pos * 2;
    ret.validity[0] |= second.validity[0] & (mask << shift);
    ret.ns[0] |= second.ns[0] & (mask << shift);
    if(pos < 31){
      ret.cs[0] |= second.cs[0] & (mask << (shift + 2));
    }
    ret.cs[0] |= ((uint64_t)getcolor(first.ns[0] >> (shift - 2), second.ns[0] >> shift)) << shift;
    /* cast to 64-bit type is necessary to enable shift more than 31 bits*/
/*    fprintf(stderr, "junctioncolor: %i\n",  getcolor(first.ns[0] >> (shift - 2), second.ns[0] >> shift));*/
    ret.validity[1] = second.validity[1];
    ret.ns[1] = second.ns[1];
    ret.cs[1] = second.cs[1];
  }else if(pos<64){/* recombination in the second half */
    ret.validity[0] = first.validity[0];
    ret.ns[0] = first.ns[0];
    ret.cs[0] = first.cs[0];

    shiftr = 64 - (pos - 32) * 2;
    shift = (pos - 32) * 2;
    if(shift){
      ret.validity[1] = first.validity[1] & (mask >> shiftr);
      ret.ns[1] = first.ns[1] & (mask >> shiftr);
      ret.cs[1] = first.cs[1] & (mask >> shiftr);
      ret.validity[1] |= second.validity[1] & (mask << shift);
      ret.ns[1] |= second.ns[1] & (mask << shift);
      ret.cs[1] |= second.cs[1] & (mask << (shift + 2));
      ret.cs[1] |= ((uint64_t)getcolor(first.ns[1] >> (shift - 2), second.ns[1] >> shift)) << shift;
    }else{
      ret.validity[1] = second.validity[1];
      ret.ns[1] = second.ns[1];
      ret.cs[1] = second.cs[1] & (mask <<  2);
      ret.cs[1] |= getcolor(first.ns[0] >> 62, second.ns[1]) ;
    }
  }else{
    fputs("invalid position parameter for recombinetarget64", stderr);
    exit(EXIT_FAILURE);
  }
  return ret;
}

static target64
getfusedtarget64(cache_st*cache, int64_t location, int tagsize, int intron_pos, int intron_length)
{
  target64 ret_getfusedtarget64;
  target64 first;
  target64 second;
  if(( intron_pos == 0)||(intron_pos >= tagsize)){
    return gettarget64(cache, tagsize, location);
  }
  first = gettarget64(cache, intron_pos, location);
  second = gettarget64(cache, tagsize, location + intron_length);
#if 0
  printasns(stderr, first.ns[0], 0, 32);
  printasns(stderr, first.ns[1], 0, tagsize - 32);
  fputc('\n', stderr);
  printascolor(stderr, first.cs[0], 0, 32);
  printascolor(stderr, first.cs[1], 0, tagsize - 32);
  fputc('\n', stderr);
  printasns(stderr, second.ns[0], 0, 32);
  printasns(stderr, second.ns[1], 0, tagsize - 32);
  fputc('\n', stderr);
  printascolor(stderr, second.cs[0], 0, 32);
  printascolor(stderr, second.cs[1], 0, tagsize - 32);
  fputc('\n', stderr);
#endif
  ret_getfusedtarget64 = recombinetarget64(first, second, intron_pos);
  return ret_getfusedtarget64;
}

static inline border16n
getbordertarget16(cache_st*cache, int64_t location, int intron_pos, int intron_length)
{
  border16n ret_getborder;
  target16 first;
  target16 second;
  first = gettarget16(cache, 1, location + intron_pos);
  second = gettarget16(cache, 1, location + intron_pos + intron_length - 2);
//  printasns(stderr, first.ns, 0, 2);
 // printasns(stderr, second.ns, 0, 5);
  ret_getborder.ns = (first.ns & 0x0F) | ((second.ns& 0x0f) << 4);
  ret_getborder.validity = first.validity | second.validity << 4;
  return ret_getborder;
}

typedef struct match_ {
  int referencesequencenumber;
  int nmismatch;
  int64_t location;
} match;
static inline match str2match(const char *);
typedef struct matchlistr_{
  int i;
  float score;
  int intron_len;
  int intron_pos;
  int align_len;
  int adapter_pos;
  int next;
} matchlistr;
typedef struct matchlistc_{
  int size;
  int first;
  int last;
  int capacity;
  matchlistr* matcharray;
}matchlistc;
static matchlistc* newmatchlist(int capacity)
{
  matchlistc* ret = malloc(sizeof(matchlistc));
  ret->size=1;
  ret->first=0;
  ret->last=0;
  ret->capacity = capacity;
  ret->matcharray = malloc(sizeof(matchlistr)*capacity);
  ret->matcharray[0].score=FLT_MAX;
  ret->matcharray[0].next=0;
  return ret;
}
  
static void
destroymatchlist(matchlistc* matchlist)
{
  free(matchlist->matcharray);
  free(matchlist);
}

static void
insertmatch(matchlistc*matchlist, int index, score_t score)
{
  matchlistr* array;
  int last;
  int current;
  int size = matchlist->size;
  if(matchlist->size == matchlist->capacity){
    array= realloc(matchlist->matcharray, sizeof(matchlistr) * matchlist->capacity * 2);
    if(array){
      matchlist->matcharray = array;
      matchlist->capacity *= 2;
    }else{
      fputs("memory allocation failed!",stderr);
      exit(EXIT_FAILURE);
    }
  }
  array = matchlist->matcharray;
  array[size].i = index;
  array[size].score = score.score;
  array[size].intron_len = score.intron_len;
  array[size].intron_pos = score.intron_pos;
  array[size].align_len = score.align_len;
  array[size].adapter_pos = score.adapter_pos;
 
  last = matchlist->last;
  current = matchlist->first;
  while(1){
    if(array[current].score > score.score)break;
    last = current;
    if(array[current].next == current) break;
    current = array[current].next;
  }
  array[last].next=size;
  array[size].next=current;
  matchlist->size++;
  if(array[size].score < array[matchlist->first].score)
    matchlist->first = size;
}

typedef struct matchdata_{
  char* tag_id;
  char* colorstring_first;
  char* colorstring_last;
  int length_first;
  int length_last;
  int length;
  uint64_t cs_first;
  uint64_t cs_last[2];
  uint64_t cs[2];
  int* qv;
  int nmatchpoints_first;
  match* matchpoints_first;
  int nmatchpoints_last;
  match* matchpoints_last;
} matchdata;
static inline match* nthmatch(matchdata* entry, int n)
{
  return entry->matchpoints_first + n;
}
static inline int64_t
getlast32(matchdata *entry)
{
  int shiftcount;
  int selector;
  int64_t retval = 0;
  shiftcount = 64 - (entry->length % 32) * 2;
  selector = entry->length/32;
  /* shift left with shiftcount */
  if(selector){
    retval = entry->cs[0] >> (64 - shiftcount);
    retval |= entry->cs[1] << shiftcount;
  }else{
    retval = entry->cs[0] << shiftcount;
  }
  return retval;
}

static void
printentryhead(FILE* out, matchdata* entry)
{
  fprintf(out, ">%s", entry->tag_id);
}
static void
printentrycolor(FILE* out, matchdata* entry)
{
  fprintf(out, "%s%s\n", entry->colorstring_first, entry->colorstring_last + 2);
}
static void
printreadtag(FILE* out, matchdata* entry)
{
  int i;
  fprintf(out, ">%s\n%s %s\n", entry->tag_id, entry->colorstring_first, entry->colorstring_last);
  for(i=0; i < entry->length; i++)
    fprintf(out, "%i ", entry->qv[i]);
  fputs("\n",out);
}
static void
printtarget64n(FILE* out, target64*target, int length)
{
  int i;
/*  int length;
  uint64_t validity[2];
  uint64_t ns[2];
  uint64_t cs[2];*/
  int shiftcount = 0;
  uint64_t curchunk_v = target->validity[0];
  uint64_t curchunk_n = target->ns[0];
  for(i = 0; i < length; i++){
    if(curchunk_v & 2){
      int curnuc = curchunk_n & 3;
      fputc(nuc2char[curnuc], out);
    }else{
      fputc('N', out);
    }
    shiftcount += 2;
    curchunk_v >>= 2;
    curchunk_n >>= 2;
    if(shiftcount >= 64){
      shiftcount = 0;
      curchunk_v = target->validity[1];
      curchunk_n = target->ns[1];
    }
  }
}

static void
printtarget64(FILE* out, target64*target)
{
  int i;
/*  int length;
  uint64_t validity[2];
  uint64_t ns[2];
  uint64_t cs[2];*/
  int shiftcount = 0;
  uint64_t curchunk_v = target->validity[0];
  uint64_t curchunk_n = target->ns[0];
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
    if(shiftcount >= 64){
      shiftcount = 0;
      curchunk_v = target->validity[1];
      curchunk_n = target->ns[1];
    }
  }
  fputc('\n', out);
  shiftcount = 0;
  curchunk_v = target->validity[0];
  curchunk_c = target->cs[0];
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
    if(shiftcount >= 64){
      shiftcount = 0;
      curchunk_v = target->validity[1];
      curchunk_c = target->cs[1];
    }
  }
  fputc('\n', out);
}
static void
printmatchedpattern(FILE* outfile, matchdata* entry, matchlistr* match, cache_st *cache)
{
#if 0
  target64 target = getfusedtarget64(cache, entry->matchpoints_first[match->i].location, entry->length, match->intron_pos, match->intron_len);
  fprintf(outfile, "targetsequence at %lli\n", entry->matchpoints[match->i].location);
  printtarget64(outfile, target);
  fprintf(outfile, "junction: %d\n", match->bestj);
#endif
  fprintf(outfile, ",%" PRId64 ".%d.%d.%d.%d.%d",entry->matchpoints_first[match->i].location, (int)(match->score*-100), match->intron_pos, match->intron_len, match->align_len, match->adapter_pos);
#if 0
  printtarget64n(outfile, &target, entry->length);
  fprintf(outfile, "target+adapter\n");
  printtarget64(outfile, target);
  printtarget64n(outfile, target, match->bestj);
  fprintf(outfile, "mismatch: %i\tqvscore: %i\n", match->score, match->qvscore);
#endif
}
void writematchlist(FILE*, matchdata*, matchlistc*, cache_st *, double);
static int
getbestscore(matchlistc*matchlist){
  int i;
  i = matchlist -> first;
  return matchlist->matcharray[i].score;
}
void
writematchlist(FILE* outfile, matchdata*entry, matchlistc* matchlist, cache_st *cache, double minscore)
{
  int i;
  double score;
  printentryhead(outfile, entry);
#if 0
  printreadtag(outfile, entry);
  fputs("\n", outfile);
#endif
  i = matchlist->first;
  score = matchlist->matcharray[i].score;
  while(1){
    if (i == matchlist->last) break;
    if (matchlist->matcharray[i].score > score + 4.0){
#if 0
      fprintf(outfile, ";%f;", matchlist->matcharray[i].score);
#endif
      break;
    }
    if (matchlist->matcharray[i].score > minscore){
#if 0
      fputs(";", outfile);
#endif
      break;
    }
    printmatchedpattern(outfile, entry, matchlist->matcharray+i, cache);
    i = matchlist->matcharray[i].next;
  }
  fputs("\n", outfile);
  printentrycolor(outfile, entry);
}
static inline int
nmatchpoint(matchdata* match)
{
  return match->nmatchpoints_first;
}
static matchdata* 
getfastaentrywithqv(FILE* matchfile_first, FILE* matchfile_last, FILE* qualfile, char* linebuf)
{
  int intheentry;
  int shiftcount, selector, i;
  matchdata *ret;
  char *p, *q;
  ret = malloc(sizeof(matchdata));
  ret->matchpoints_first= calloc(max_match_points, sizeof(match));
  ret->matchpoints_last= calloc(max_match_points, sizeof(match));
  ret->tag_id = NULL;
  ret->nmatchpoints_first = 0;
  ret->nmatchpoints_last = 0;
  ret->length_first = 0;
  ret->length_last = 0;
  while((p = fgets(linebuf, 1024000, matchfile_first))){
    if(*p == '>'){
      char *matchstr;
      ret->tag_id = strdup(strtok(p + 1, ",\n"));
      while((matchstr = strtok(NULL, ",\n"))){
        ret->matchpoints_first[ret->nmatchpoints_first] = str2match(matchstr);
        ret->nmatchpoints_first++;
      }
    }else if(*p=='#'){
      continue;
    }else{
      /* process color string */
      size_t length = strlen(linebuf);
      ret->colorstring_first = q = malloc(length);
      if(is_validnuc(*p))
        (*(q++)) = (*(p++));
      while(*p){
        if((*p >= '0' && *p <= '3') || *p == '.'){
          *(q++) = *(p++);
          ret->length_first ++;
        }else if(isspace(*p)){
           p++;
        }else{
           break;
        }
      }
      *q = '\0';
      ret->cs_first = 0;
      ret->cs_first |=  getnextbase(char2nucmap[(unsigned)ret->colorstring_first[0]], char2color(ret->colorstring_first[1]));
      shiftcount = 2;
      for(i=2; i <= ret->length_first; i++){
        ret->cs_first |= (uint64_t)(3&char2color(ret->colorstring_first[i])) << shiftcount;
        shiftcount += 2;
        if(shiftcount == 64){
          shiftcount = 0;
          fprintf(stderr, "First part too long!\n");
          printascolor(stderr, ret->cs[0], 0, 31);
          printascolor(stderr, ret->cs[1], 0, 31);
          fprintf(stderr, "i: %i, ret->length: %i\n", i, ret->length_first);
          fprintf(stderr, "tag: %s %s\n", ret->colorstring_first, ret->colorstring_last);
          exit(EXIT_FAILURE);
        }
      }
      break;
    }
  }
  if(!ret->tag_id){
    free(ret);
    return NULL;
  }
  intheentry = 0;
  while((p = fgets(linebuf, 1024000, matchfile_last))){
    if(*p == '>'){
      /* check if the tag_id match the firstpart */
      p++;
      q = ret->tag_id;
      while(*q){
        if (*q != *p) break;
        q++; p++;
      }
      if(*q == '\0'){
        /* the right entry is found */
        intheentry = 1;
      }else if (strcmp(q, "_FILT") == 0){
//        fprintf(stdout, "first half was filtered: %s\n", linebuf);
        *q = '\0';
        intheentry = 1;
      }else if (strcmp(p, "_FILT") == 0){
 //       fprintf(stdout, "second half was filtered: %s\n", linebuf);
        intheentry = 1;
      }else {
//        fprintf(stdout, "second file specific tag: %s\n", linebuf);
        exit(EXIT_FAILURE);
      }
    }else if(*p=='#'){
      continue;
    }else{
      size_t length;
      if(!intheentry){
//        fprintf(stdout, "Entry id did not match\n");
 //       fprintf(stdout, "first tag id: %s\n", ret->tag_id);
        exit(EXIT_FAILURE);
      }
      /* process color string */
      length = strlen(linebuf);
      ret->colorstring_last = q = malloc(length);
      if(is_validnuc(*p))
        (*(q++)) = (*(p++));
      while(*p){
        if((*p >= '0' && *p <= '3') || *p == '.'){
          *(q++) = *(p++);
          ret->length_last ++;
        }else if(isspace(*p)){
           p++;
        }else{
           break;
        }
      }
      *q = '\0';
      ret->cs_last[0] = ret->cs_last[1] = 0;
      shiftcount = 0;selector = 0;
      /* the first two charactor T0 is meaningless for the last half*/
      for(i=2; i <= ret->length_last; i++){
        ret->cs_last[selector] |= (uint64_t)(3&char2color(ret->colorstring_last[i])) << shiftcount;
        shiftcount += 2;
        if(shiftcount == 64){
          shiftcount = 0;
          selector += 1;
          if(selector >1){
            fprintf(stderr, "Too long tag sequence!\n");
            printascolor(stderr, ret->cs[0], 0, 31);
            printascolor(stderr, ret->cs[1], 0, 31);
            fprintf(stderr, "i: %i, ret->length: %i\n", i, ret->length_last);
            fprintf(stderr, "tag: %s\n", ret->colorstring_last);
            exit(EXIT_FAILURE);
          }
        }
      }
      break;
    }
  }
  ret->cs[0] = ret->cs_first | (ret->cs_last[0]  << (ret->length_first * 2));
  ret->cs[1] = (ret->cs_last[0] >> (64 - ret->length_first * 2)) |(ret->cs_last[1]  << (ret->length_first * 2));
  ret->length = ret->length_first + ret->length_last - 1;
  ret->qv = calloc(ret->length, sizeof(int));
  intheentry = 0;
  while((p = fgets(linebuf, 1024000, qualfile))){
    if(*p == '>'){
      /* check if the tag_id match the csfasta */
      p++;
      q = ret->tag_id;
      while(*q){
        if (*q != *p) break;
        q++; p++;
      }
      if(*q == '\0'){
        /* the right entry is found */
        intheentry = 1;
      }else if (strcmp(q, "_FILT") == 0){
//        fprintf(stdout, "this tag was filtered: %s\n", linebuf);
        *q = '\0';
        intheentry = 1;
      }else {
//        fprintf(stdout, "QV file specific tag: %s\n", linebuf);
      }
    }else if(*p=='#'){
      continue;
    }else{
      if(intheentry){
        for(i = 0; i< ret->length; i++){
          ret->qv[i] = strtol(p, &q, 10);
          p=q;
        }
        break;
      }
    }
  }
  return ret;
}

static void
destroymatchdata(matchdata*entry)
{
  if(!entry)
    return;
  free(entry->tag_id);
  free(entry->colorstring_first);
  free(entry->colorstring_last);
  free(entry->qv);
  free(entry->matchpoints_first);
  free(entry->matchpoints_last);
  free(entry);
}

static inline match
str2match(const char *str)
{
  match ret;
  char *p = NULL;
  int64_t num;
  errno = 0;
  ret.referencesequencenumber = -1;
  ret.nmismatch = -1;
  ret.location = 0;
  num = strtoll(str, &p, 10);
  if (*p == '_'){
    ret.referencesequencenumber = num;
    num = strtoll(p+1, &p, 10);
  }
  if(*p == '.'){
    ret.location = num;
    num = strtoll(p+1, &p, 10);
    ret.nmismatch=num;
  }
  return ret;
}

static inline int
printmatch(FILE*outfile,match matchloc)
{
  const char *strand = "T";
  int64_t location = matchloc.location;
  if(matchloc.location < 0) {
    strand = "R";
    location = -matchloc.location;
  }
  if(matchloc.referencesequencenumber < 0){
    return fprintf(outfile,"-\t%s\t%" PRId64 "\n", strand, location);
  }else{
    return fprintf(outfile,"%i\t%s\t%" PRId64 "\n", matchloc.referencesequencenumber, strand, location);
  }
}


static void usage()
{
  fputs("matchjumpqv -m matchfile -q quality -s seqmapfile -c sequencecachefile [-o outfile]\n", stderr);
}

static int64_t
find_bestmatchlocation(cache_st *cache, int64_t firstlocation, uint64_t last32cs, int maxgapsize, int comparebases)
{
  int64_t nuccount;
  int64_t bestlocation = firstlocation;
  int minmismatch = 32;
  uint64_t selector;
  int shiftcount;
  uint64_t color32;
  uint64_t n32;
  uint64_t tcolor32;
  uint64_t tn32;
  uint64_t comp32;
  uint64_t mask = 0x5555555555555555LL;
  int i;
  if(!cache){
    fputs("find_bestmatchlocation was called with NULL pointer", stderr);
    exit(EXIT_FAILURE);
  }
  if(firstlocation >=0){
    nuccount = firstlocation;
  }else{
    nuccount = cache->reverse_base + firstlocation - 2;
  }
  selector = (nuccount)/32;
  shiftcount = (nuccount) % 32 * 2;
  tcolor32 = cache->data[selector].colorseq;
  tn32 = cache->data[selector].nucseq;
  color32 = tcolor32 >> shiftcount;
  n32 = tn32 >> shiftcount;
  selector += 1;
  tcolor32 = cache->data[selector].colorseq;
  tn32 = cache->data[selector].nucseq;
  color32 |= tcolor32 << (64-shiftcount);
  n32 |= tn32 << (64-shiftcount);
  mask <<= (64-comparebases*2);
  for(i=0; i<maxgapsize; i++){
    comp32 = color32 ^ last32cs;
    comp32 |=  comp32 >> 1;
    comp32 &= mask;
    comp32 += comp32 >> 2;
    comp32 &= 0x3333333333333333LL;
    comp32 += comp32 >> 4;
    comp32 += comp32 >> 8;
    comp32 &= 0x000F000F000F000FLL;
    comp32 += comp32 >> 16;
    comp32 += comp32 >> 32;
    comp32 &= 0x3FLL;

#if 0 
      fprintf(stderr, "mismatch: %i, i=%d\n", comp32, i);
      printasns(stderr, n32, 0, 32);
      fputc('\n', stderr);
      printascolor(stderr, color32, 0, 32);
      fputc('\n', stderr);
      printascolor(stderr, last32cs, 0, 32);
      fputc('\n', stderr);
      printascolor(stderr, color32^last32cs, 0, 32);
      fputc('\n', stderr);
      fputc('\n', stderr);
#endif
    if(comp32 < minmismatch){
      minmismatch = comp32;
      bestlocation = firstlocation + i;
    }

    color32 >>= 2;
    color32 |= (3LL << 62) & (tcolor32 << (62 - shiftcount));
    n32 >>= 2;
    n32 |= (3LL << 62) & (tn32 << (62 - shiftcount));
    shiftcount += 2;
    if (shiftcount >= 64){
      shiftcount = 0;
      selector += 1;
      tcolor32 = cache->data[selector].colorseq;
      tn32 = cache->data[selector].nucseq;
    }
  }
#if 0
 fprintf(stderr, "minmismatch: %i, i = %d\n", minmismatch, bestlocation - firstlocation);
#endif 
  return bestlocation;
}

static inline double mismatchpenalty(int qv)
{
  const double log075 = log(0.75)/log(2.0);
  const double log10 = log(10.0)/log(2.0);
  return log10*qv / 10.0 - log075; /* this formula should be reconsidered */
}
static inline double matchreward(int qv)
{
  const double rlog2 = 1/log(2.0);
  return 2.0 + log(1-pow(10, -qv/10.0))*rlog2; 
}

static inline double gappenalty(int gapsize, border16n acdo)
{
  const double rlog2 = 1/log(2.0);
  double score = 0.0;
  score += log(gapsize)*rlog2 + 4;
  if((gapsize > 16) && (acdo.ns & 0xFF) == (2|(3<<2)|(0<<4)|(2<<6))){
    /* intron accepted with canonical acceptor donor pair */
    /*GT -- AG*/
    score -= 8;
  }
  return score;
}

static inline score_t
reevalwithadapter(cache_st *cache, matchdata* entry, int64_t firstlocation, score_t score, tag *adapter)
{
  int i,j;
  score_t ret = score;
  double positionalscore[50];
  double opositionalscore[50];
  int best_junction = 50;
  double  best_terminal_score = FLT_MAX;
  uint64_t cd;
  target64 t;
  t = getfusedtarget64(cache, firstlocation, entry->length, score.intron_pos, score.intron_len);
  cd = t.cs[0] ^ entry->cs[0];
  for(i = 0; i < 50;){
    if(entry->qv[i] > 0){
      if(cd & 3){ /*mismatch*/
        positionalscore[i] = mismatchpenalty(entry->qv[i]);
      }else{ /*match*/
        positionalscore[i] = - matchreward(entry->qv[i]);
      }
    }else{ 
      positionalscore[i] = 0;
    }
    cd >>= 2;
    i++;
    if(i == 32) cd = t.cs[1] ^ entry->cs[1];
  }
  for(i = 1; i < 50 ; i++){
    positionalscore[i] += positionalscore[i-1];
  }
//  fputs("original positinal score: ", stdout);
  for(i= 0; i < 50; i++){
    opositionalscore[i] = positionalscore[i];
//    fprintf(stdout, "%i: qv %i, %f; ", i, entry->qv[i], positionalscore[i]);
  }
//  fputs("\n", stdout);
  for(i = 50; i >=32; i--){
    mergeadapter(&t, i, adapter);
    cd = t.cs[1] ^ entry->cs[1];
    cd >>= 2 * (i-32);
    for(j = i; j < 50; j++){
      if(entry->qv[j] > 0){
        if(cd & 3){
          positionalscore[j] = mismatchpenalty(entry->qv[j]);
        }else{
          positionalscore[j] = - matchreward(entry->qv[j]);
        }
        positionalscore[j] += positionalscore[j-1];
      }else{
        positionalscore[j] = positionalscore[j-1];
      }
      cd >>= 2;
    }
//    fprintf(stdout, "%i: %f; ", i, positionalscore[49]);
    if(positionalscore[49] < best_terminal_score){
      best_junction = i;
      best_terminal_score = positionalscore[49];
    }
  }
  for(; i >= 0; i--){
    mergeadapter(&t, i, adapter);
    cd = t.cs[0] ^ entry->cs[0];
    cd >>= 2 * (i-32);
    for(j = i; j < 50;){
      if(entry -> qv[j] > 0){
        if(cd & 3){
          positionalscore[j] = mismatchpenalty(entry->qv[j]);
        }else{
          positionalscore[j] = - matchreward(entry->qv[j]);
        }
        positionalscore[j] += positionalscore[j-1];
      }else{
        positionalscore[j] = positionalscore[j-1];
      }
      cd >>= 2;
      j++;
      if(j == 32) cd = t.cs[1] ^ entry->cs[1];
    }
//    fprintf(stdout, "%i: %f; ", i, positionalscore[49]);
    if(positionalscore[49]< best_terminal_score){
      best_junction = i;
      best_terminal_score = positionalscore[49];
    }
  }
  if(score.intron_pos > 0 && score.intron_pos < best_junction){
    double gapp;
    border16n acdo = getbordertarget16(cache, firstlocation, score.intron_pos, score.intron_len);
    gapp = gappenalty(score.intron_len, acdo);
    for(i = score.intron_pos; i < best_junction; i ++){
      opositionalscore[i] += gapp;
    }
  }
  {
    double bestscore = FLT_MAX;
    int bestscorepos =0;
    for(i = 0; i< best_junction; i ++){
//      fprintf(stdout, "%i: %f, %f; ", i, opositionalscore[i], bestscore);
      if(opositionalscore[i] < bestscore){
        bestscore = opositionalscore[i];
        bestscorepos = i;
      }  
    }
    ret.score = bestscore;
    ret.align_len= bestscorepos;
    ret.adapter_pos = best_junction;
    ret.adapter_align_pos = 50;
  }
  return ret;
}

static inline  score_t 
evalalignment(cache_st *cache,  matchdata* entry, int64_t firstlocation, int gappos, int gapsize, target64 t)
{
  int i;
  uint64_t cd;
  score_t retv;
  double score = 0.0;
  double minscore = FLT_MAX;
  double minscore_pregap = FLT_MAX;
  int minscore_pos = 0;
  int minscore_pos_pregap = 0;
  const double rlog2 = 1/log(2.0);
 
  cd = t.cs[0] ^ entry->cs[0];
  i = 0;
  while(i < gappos){
    if(entry->qv[i] > 0){
      if(cd & 3){ /*mismatch*/
        score += mismatchpenalty(entry->qv[i]);
      }else{ /*match*/
        score -= matchreward(entry->qv[i]);
        if(score < minscore){
          minscore = score;
          minscore_pos = i;
        }
      }
    }
    cd >>= 2;
    i++;
    if(i == 32) cd = t.cs[1] ^ entry->cs[1];
  }
  minscore_pregap = minscore;
  minscore_pos_pregap = minscore_pos;
  if(gapsize){
    border16n acdo = getbordertarget16(cache, firstlocation, gappos, gapsize);
    score += gappenalty(gapsize, acdo);
  }
  while(i < entry->length){
    if(entry->qv[i] > 0){
      if(cd & 3){ /*mismatch*/
        score += mismatchpenalty(entry->qv[i]);
      }else{ /*match*/
        score -= matchreward(entry->qv[i]); 
        if(score < minscore){
          minscore = score;
          minscore_pos = i;
        }
      }
    }
    cd >>= 2;
    i++;
    if(i == 32) cd = t.cs[1] ^ entry->cs[1];
  }
  retv.score = (float) minscore;
  retv.intron_len = gapsize;
  retv.intron_pos = (int8_t) gappos;
  retv.align_len = (int8_t) minscore_pos;
  if(minscore_pos < gappos){
    retv.intron_len = 0;
    retv.intron_pos = 0;
  }
  return retv;
}

static int64_t
find_best_gapposition(cache_st *cache, matchdata* entry, int64_t firstlocation, int gapsize)
{
  int i;
  int best_i = -1;
  score_t tscore;
  score_t score;
  target64 first, second;
  score.score = FLT_MAX;
  first = gettarget64(cache, entry->length, firstlocation);
  second = gettarget64(cache, entry->length, firstlocation + gapsize);
  for(i = 6; i < entry->length; i++){
    target64 t;
    t = recombinetarget64(first, second, i);
    tscore = evalalignment(cache, entry, firstlocation, i, gapsize, t); 
    if(tscore.score < score.score){
      score.score = tscore.score;
      best_i = i;
    }
//    score
  }
  return best_i;

}

int
main(int argc, char **argv)
{
  int optc;
  char *linebuf, *p;
  FILE *matchfile_first = stdin;
  FILE *matchfile_last = NULL;
  FILE *qualfile = NULL;
  FILE *outfile = stdout;
  FILE *statfile;
  char *statfilename;
  int cachefile;
  FILE *seqmapfile = NULL;
  seqmap_r *seqlocarray;
  int nchr;
  int *chr_sizes;
  uint64_t *chr_start_addr;
  cache_st cache;
  int maxgapsize = max_intron_size;
  double threshold = 40.0;
  tag* adapter = NULL;
  initchar2nucmap();
  while((optc = getopt(argc,argv,"a:c:f:m:l:o:s:t:q:z:")) != -1){
    switch(optc){
    case 'a':
      adapter = newadapter(optarg);
      break;
    case 'f':
    case 'm':
        matchfile_first = fopen(optarg, "r");
        if(matchfile_first == NULL){
          perror(argv[0]);
          fputs("matchfile_first open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'l':
        matchfile_last = fopen(optarg, "r");
        if(matchfile_last == NULL){
          perror(argv[0]);
          fputs("matchfile_last open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;

    case 'q':
        qualfile = fopen(optarg, "r");
        if(qualfile == NULL){
          perror(argv[0]);
          fputs("qualfile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'c':
        cachefile = open(optarg, O_RDONLY, 0644);
        if(cachefile < 0){
          perror(argv[0]);
          fputs("cachefile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
        assigncache(&cache, cachefile);
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
        outfile = fopen(optarg, "w");
        if(outfile == NULL){
          perror(argv[0]);
          fputs("outfile open failed: ", stderr);
          fputs(optarg, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
        statfilename = malloc(strlen(optarg) + 6);
        strcpy(statfilename,optarg);
        strcat(statfilename, ".stat");
        statfile = fopen(statfilename, "w");
        if(statfile == NULL){
          perror(argv[0]);
          fputs("outfile open failed: ", stderr);
          fputs(statfilename, stderr);
          fputc('\n', stderr);
          usage();
          exit(EXIT_FAILURE);
        }
      break;
    case 'r':
      break;
    case 't':
      threshold = strtod(optarg, NULL);
      break;
    case 'z':
      max_match_points = strtol(optarg, NULL, 0);
      break;
    default:
      usage();
      exit(EXIT_FAILURE);
    }
  }
  if(!seqmapfile){
    usage();
    exit(EXIT_FAILURE);
  }
  
  linebuf = malloc(1024000);
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
  chr_start_addr = malloc(sizeof(uint64_t)*nchr_max);
  if(chr_start_addr == NULL){
    perror(argv[0]);
    fputs("memory allocation failed: ", stderr);
  }
  
  {
    register int i = 0;
    register int j = 0;
    nchr = 0;
    while((p = fgets(linebuf, 1024000, seqmapfile))){
      seqlocarray[i] = scan_seqmap_record(p);
      j = seqlocarray[i].pseudo_chr_number;
      if(j>nchr) nchr=j;
      chr_sizes[j] = seqlocarray[i].offset + seqlocarray[i].length;
      i++;
    }
  }

  {
    register int i = 0;
    register uint64_t t = 1;
    chr_start_addr[0] = 0;
    for(i=1; i<= nchr; i++){
      chr_start_addr[i] = t;
      t += chr_sizes[i];
    }
  }

  {
    score_t score;
    matchdata *entry;
    while((entry = getfastaentrywithqv(matchfile_first, matchfile_last, qualfile, linebuf))){
      int i;
      matchlistc *matchlist;
      uint64_t last32cs;
      int64_t firstlocation, location, gappos;
      int64_t gapsize;
      matchlist = newmatchlist(max_match_points);
      last32cs = getlast32(entry);
#if 0
      printreadtag(outfile, entry);
      printascolor(outfile, last32cs, 0, 32);
      fputc('\n', outfile);
#endif
      for(i=0; i < nmatchpoint(entry); i++){
        target64 t;
        /* at first evaluate without gap*/
        firstlocation = entry->matchpoints_first[i].location;
        t = gettarget64(&cache, entry->length, firstlocation);
/*        printtarget64(stdout, &t);*/
        score = evalalignment(&cache, entry, firstlocation, 0, 0, t); 
        if(score.align_len < 45){
        /* if there is posibility to find a gapped alignment */
          int comparebases = 32;
          if(comparebases > 50 - score.align_len + 4){
            comparebases = 50 - score.align_len + 4;
          }
          location = find_bestmatchlocation(&cache, firstlocation + entry->length - 32, last32cs, maxgapsize, comparebases);
          gapsize = location - (firstlocation + entry->length - 32);
          if(gapsize > 0){
            gappos = find_best_gapposition(&cache, entry, firstlocation, gapsize);
            t = getfusedtarget64(&cache, firstlocation, entry->length, gappos, gapsize);
            score = evalalignment(&cache, entry, firstlocation, gappos, gapsize, t); 
          }
        }
#if 0
        fprintf(outfile, "firstlocation: %" PRId64", gapsize: %" PRId64 ", gappos: %" PRId64 "\n", firstlocation, gapsize, gappos);
        fprintf(outfile, "score: %f, intron_len: %i, intron_pos: %i, align_len: %i \n", score.score, score.intron_len, score.intron_pos,score.align_len); 
#endif
        if(score.score < - threshold){
          score = reevalwithadapter(&cache, entry, firstlocation, score, adapter);
          insertmatch(matchlist, i, score);
        }
#if 0
        fprintf(outfile, "firstlocation: %" PRId64", gapsize: %" PRId64 ", gappos: %" PRId64 "\n", firstlocation, gapsize, gappos);
        fprintf(outfile, "score: %f, intron_len: %i, intron_pos: %i, align_len: %i, adapter_pos: %i \n", score.score, score.intron_len, score.intron_pos,score.align_len, score.adapter_pos); 
#endif

      }
      writematchlist(outfile, entry, matchlist, &cache, -threshold);
      destroymatchlist(matchlist);
      destroymatchdata(entry);
    }
  }
  return 0;
}

