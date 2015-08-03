/* matchlist.c*/
/*
 * Copyright 2008, 2009 Tomoaki Nishiyama, Kanazawa University
 */
#include <stdlib.h>
#include "matchlist.h"

int count_lines(line* lp)
{
  if(lp == NULL) return 0;
  return 1 + count_lines(lp->next);
}

line* newline()
{
  line* lp;
  lp = malloc(sizeof(line));
  lp->next = NULL;
  lp->head = NULL;
  lp->tail = NULL;
  lp->end_location = 0;
  return lp;
}
element* newelement(readtag* newtag, int64_t location)
{
  element* ep;
  ep = malloc(sizeof(element));
  ep->location = location;
  ep->next = NULL;
  ep->tagp = newtag;
  return ep;
}

line*
insertnewalignment(line *lines, readtag* newtag, int64_t location, int minspace)
{
  if(lines == NULL){
    lines = newline();
    lines->head = newelement(newtag, location);
    lines->tail = lines->head;
    lines->end_location = location + newtag->gaplen + 50; /* tag length insteadof align length */
    return lines;
  }
  if(lines->end_location + minspace < location){
    lines->tail->next = newelement(newtag, location);
    lines->tail = lines->tail->next;
    lines->end_location = location + newtag->gaplen + 50;
  }else{
    lines->next = insertnewalignment(lines->next, newtag, location, minspace);
  }
  return lines;
}
