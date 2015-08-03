
#if !defined(matchlocfile_h_included)
#define matchlocfile_h_included

#include <stdio.h>
#include <stdint.h>
#include "matchlocations.h"

typedef struct cursorcount_{
  uint32_t cursor;
  uint32_t count;
} cursorcount;
typedef struct pos_{
  genomeloc_t loc;
  uint32_t count;
  uint64_t read_tag_head;
} posdata;
typedef tag50align readtag;
int WriteTaglocToFile(int, matchtagindex*, uint64_t, uint64_t, uint64_t);

static inline int64_t roundinps(int64_t ds, size_t pagesize)
{
  return (ds/pagesize +1)* pagesize;
}

#endif
