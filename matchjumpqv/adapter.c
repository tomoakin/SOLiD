/*
 * adapter.c
 * Copyright 2008, 2009 Tomoaki Nishiyama, Kanazawa University
 */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <assert.h>
#include "colorbase.h"
#include "basecache.h"
#include "adapter.h"
tag* newadapter(char*tagstring)
{
  int selector;
  int length;
  int shiftcount;
  int lastnuc, curnuc, curcolor;
  tag *ret;
  char *p = tagstring;
  ret = (tag*)malloc(sizeof(tag));
  if(ret == NULL){
    fprintf(stderr, "memory allocation failed for tag!\n");
    exit(EXIT_FAILURE);
  }
  ret->cs[0] =0; ret->cs[1] =0; ret->cs[2] =0;
  ret->ns[0] =0; ret->ns[1] =0; ret->ns[2] =0;
  if(!(char2nucmap[0xff&*p]&4)){
    fprintf(stderr, "adapter not a nucleotide sequence!\n");
    fprintf(stderr, "adapter: %s\n", tagstring);
    exit(EXIT_FAILURE);
  }
  lastnuc = ret->head_nuc = char2nucmap[0xff&*p] & 3;
  length = 0;
  selector = 0;
  shiftcount = 0;
  while(char2nucmap[0xff&*p]&4){
    curnuc = char2nucmap[0xff&*p];
    ret->ns[selector] |= (uint64_t)(curnuc&3) << shiftcount;
    curcolor = getcolor(lastnuc, curnuc);
    ret->cs[selector] |= (uint64_t)curcolor << shiftcount;

    length ++;
    shiftcount +=2;
    if(shiftcount == 64){
      shiftcount = 0;
      selector +=1;
      if (selector >2){
        fprintf(stderr, "Too long adapter sequence!\n");
        fprintf(stderr, "adapter: %s\n", tagstring);
        exit(EXIT_FAILURE);
      }
    }
    p++; lastnuc = curnuc;
  }
    
  ret->length = length;
  return ret;
}

void
mergeadapter(target64*target, int j, tag*adapter)
{
  /*joint target seq with the adpter */
  uint64_t cs[2], ns[2];
  int selector, shiftcount;
  uint64_t targetmask[2];
  targetmask[0]= 0xFFFFFFFFFFFFFFFFULL;
  targetmask[1]= 0xFFFFFFFFFFFFFFFFULL;
  selector = 1;
  shiftcount = 128 - j * 2;
  if(shiftcount >= 64){
    targetmask[selector] = 0;
    selector -= 1;
    shiftcount -= 64;
  }
  if(shiftcount>0){
    targetmask[selector] >>= shiftcount;
  }
  cs[0] = target->cs[0] & targetmask[0];
  cs[1] = target->cs[1] & targetmask[1];
  ns[0] = target->ns[0] & targetmask[0];
  ns[1] = target->ns[1] & targetmask[1];
  if(j<32){
    cs[0] |= (uint64_t)getcolor((target->ns[0] >> (j-1)*2) & 3, adapter-> head_nuc) << (j)*2;
    cs[0] |= adapter->cs[0] << j*2;
    ns[0] |= adapter->ns[0] << j*2;
    cs[1] |= adapter->cs[0] >> (64-j*2);
    ns[1] |= adapter->ns[0] >> (64-j*2);
    cs[1] |= adapter->cs[1] << j*2;
    ns[1] |= adapter->ns[1] << j*2;
  }
  if(j == 32){
    cs[1] |= getcolor((target->ns[0] >> 31*2) & 3, adapter-> head_nuc) ;
    cs[1] |= adapter->cs[0];
    ns[1] |= adapter->ns[0];
  }
  if(j > 32){
    cs[1] |= (uint64_t)getcolor((target->ns[1] >> (j-32-1)*2) & 3, adapter-> head_nuc) << (j-32)*2;
    cs[1] |= adapter->cs[0] << (j-32)*2;
    ns[1] |= adapter->ns[0] << (j-32)*2;
  }
  target->cs[0] = cs[0];
  target->cs[1] = cs[1];
  target->ns[0] = ns[0];
  target->ns[1] = ns[1];
}
