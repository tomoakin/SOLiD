/* basecache.h */
/* defines structures and functions associated with base and color sequence cache*/
#if !defined(basecache_h_included)
#define basecache_h_included
#include<inttypes.h>
#include<sys/types.h>
#include<sys/stat.h>
#include<sys/mman.h>
#include<limits.h>
#include "colorbase.h"
typedef struct chunk32_{
  uint64_t validity;
  uint64_t nucseq;
  uint64_t colorseq;
} chunk32;
typedef struct cache_st_{
  uint64_t totalnuc;
  uint64_t reverse_base;
  chunk32* data;
} cache_st;

void assigncache(cache_st*cache, int cachefd);


static inline void
setbase(chunk32* seqcachearray, int nuccount, int nuc, int color, int curnucvalid, int curcolorvalid)
{
  uint64_t selector = (nuccount - 1) /32;
  int shiftcount= (nuccount-1) % 32 * 2;
  seqcachearray[selector].validity |= ((curnucvalid?2:0) | (curcolorvalid?1:0)) << shiftcount;
  seqcachearray[selector].nucseq |= (nuc & 3) << shiftcount;
  seqcachearray[selector].colorseq |= (color & 3) << shiftcount;
}

static inline int
cachegetcolor(const cache_st * cache, int64_t location)
{
  const chunk32* seqcachearray;
  int64_t nuccount;
  uint64_t selector;
  int shiftcount;
  if(!cache) return 4;
  seqcachearray = cache->data;
  if((location < 0) && (-location > cache->totalnuc))
      return 4;
  if((location > 0) && (location > cache->totalnuc)){
    return 4;
  }
  if(location < 0){
    nuccount = location + cache->reverse_base;
  }else{
    nuccount = location;
  }
  selector = (nuccount - 1) / 32;
  shiftcount = (nuccount - 1) % 32 * 2;
  if((seqcachearray[selector].validity >> shiftcount) &1)
    return (seqcachearray[selector].colorseq >> shiftcount) & 3;
  else return 4;
}

typedef struct border16n_{
  uint32_t validity;
  uint32_t ns;
} border16n;

typedef struct target16_{
  int length;
  uint32_t validity;
  uint32_t ns;
  uint32_t cs;
} target16;

typedef struct target32_{
  int length;
  uint64_t validity[2];
  uint64_t ns[2];
  uint64_t cs[2];
} target32;

typedef struct target64_{
  int length;
  uint64_t validity[2];
  uint64_t ns[2];
  uint64_t cs[2];
} target64;

typedef struct tag_{
  int length;
  int head_nuc;
  uint64_t cs[3];
  uint64_t ns[3];
} tag;
target32* gettarget32(cache_st*, int64_t);
target16 gettarget16(cache_st*, int, int64_t);
target64 gettarget64(cache_st*, int, int64_t);
void printtarget16(FILE* out, target16*target);

int printseq(FILE*outfile, cache_st*cache, int start, int length, int zerobase);
int printcolorseq(FILE*outfile, cache_st*cache, int start, int length, int zerobase);
#endif
