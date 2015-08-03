/* Copyright 2009, Tomoaki Nishiyama, Kanazawa University */
/* implementation of matchlocations*/
#include<stdint.h>
#include<stdlib.h>
#include<limits.h>
#include "matchlocations.h"
tagdata50
str2tag50(const char *str)
{
  tagdata50 ret;
  int selector = 0;
  int shift = 0;
  int color;
  const char * p = str;
  ret.threebitseq[0]=0;
  ret.threebitseq[1]=0;
  ret.threebitseq[2]=0;
  ret.threebitseq[3]=0;
  ret.threebitseq[4]=0;
  while(*p){
    switch(*p){
      case '0':
      case '1':
      case '2':
      case '3':
        color = *p - '0';
        ret.threebitseq[selector] |= color << shift;
        break;
      case '.':
        color = 4;
        ret.threebitseq[selector] |= color << shift;
        break;
      default:
        p++;
        continue;
    }
    shift += 3;
    if(shift > 27){
      shift = 0;
      selector += 1;
    }
    p++;
  }
  return ret;
}


matchtagindex*
initmatchtagindex(matchtagindex* this, uint64_t genomesize, int shiftbits)
{
  uint64_t fm =1;
  int maxshiftbit = 0;
  uint64_t as;
  if(this == NULL) this= malloc(sizeof(matchtagindex));
  this->genomesize = genomesize;
  this->shiftbitcount = shiftbits;
  while(fm<genomesize){
    maxshiftbit += 1;
    fm <<= 1;
  }
  fm >>= shiftbits - 1;
  as = fm;
  fm -= 1;
  fm <<= shiftbits;
  this->upperbitmask = fm;
  this->lowerbitmask = UINT64_MAX ^ (UINT64_MAX << shiftbits);
  this-> pary = calloc(sizeof(matchtaglist), as);
  if(this->pary){
    this -> arraysize = as;
    return this;
  }
  return NULL;
}
matchtaglist*
newmatchtaglist()
{
  matchtaglist * ret;
  ret = calloc(sizeof(matchtaglist),1);
  return ret;
}
matchtaglist*
inserttaginlist(matchtaglist* list, int location, tag50align tag, uint64_t* totalpos)
{
  if(list == NULL){
    list = newmatchtaglist();
  }
  if(list-> array_content_number  == 0){
    list -> location = location;
    list -> alignedtag = malloc(sizeof(tag50align));
    list -> allocated_array_size = 1;
    list -> array_content_number = 1;
    *totalpos += 1;
    *(list -> alignedtag) = tag;
  return list;
  }
  /* list-> array_content_nubmer != 0*/
  if(list -> location == location){
    if( list -> array_content_number >= list -> allocated_array_size){
      list->allocated_array_size *= 2;
      list->alignedtag = realloc(list->alignedtag, sizeof(tag50align) * list-> allocated_array_size);
    }
    list -> alignedtag[list->array_content_number] = tag;
    list -> array_content_number += 1;
    return list;
  }
  if(list -> location < location){
    list->next = inserttaginlist(list -> next, location, tag, totalpos);
    return list;
  }
  if(list -> location > location){
    matchtaglist * nc;
    nc = inserttaginlist(NULL, location, tag, totalpos);
    nc -> next = list;
    return nc;
  }
  return list;
}
int
inserttag(matchtagindex* mti, int location, tag50align tag, uint64_t *totalpos)
{
  unsigned int indexcursor= (location & mti->upperbitmask) >> mti->shiftbitcount;
  mti->pary[indexcursor] = inserttaginlist(mti->pary[indexcursor], location, tag, totalpos);
  return 0;
}
int
printtagdata50(FILE*out, tagdata50 tag)
{
  int i;
  int shift;
  int mask = 7;
  char table[]="0123....";
  for(i =0; i < 5; i++){
    for(shift=0; shift < 30; shift += 3){
      fputc(table[(tag.threebitseq[i] >> shift )& mask], out);
    }
  }
  return 0;
}

char*
tagdata50tostr(char*out, tagdata50 tag)
{
  int i;
  int shift;
  int mask = 7;
  char table[]="0123....";
  int j = 0;
  for(i =0; i < 5; i++){
    for(shift=0; shift < 30; shift += 3){
      out[j++]=table[(tag.threebitseq[i] >> shift )& mask];
    }
  }
  return 0;
}
int*
tagdata50tointarray(int*out, tagdata50 tag)
{
  int i;
  int shift;
  int mask = 7;
  int j = 0;
  for(i =0; i < 5; i++){
    for(shift=0; shift < 30; shift += 3){
      out[j++]=(tag.threebitseq[i] >> shift )& mask;
    }
  }
  return 0;
}
