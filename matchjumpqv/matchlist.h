/* matchlist.h*/
/*
 * Copyright 2008, 2009 Tomoaki Nishiyama, Kanazawa University
 */
#if !defined(matchlist_h_included)
#define matchlist_h_included
#include<inttypes.h>
#include "matchlocfile.h"
typedef struct element_{
  struct element_* next;
  int64_t location;
  readtag *tagp;
} element;

typedef struct line_{
  struct line_* next;
  element *head;
  element *tail;
  int64_t end_location; /* can be calculated from tagp but keep it outside*/
} line;
line* insertnewalignment(line *lines, readtag* newtag, int64_t location, int minspace);
int count_lines(line*);

#endif
