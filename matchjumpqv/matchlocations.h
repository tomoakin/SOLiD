
#if !defined(matchlocations_h_included)
#define matchlocations_h_included

#include <stdio.h>
#include <stdint.h>
typedef int32_t genomeloc_t;

typedef struct tagdata50_{
  uint32_t threebitseq[5];/* color sequence with invalid state (.) */
} tagdata50;
typedef struct tag50align_{ /*static relocatable structure*/
  uint8_t alignedlen;
  uint8_t gappos;
  uint8_t adapter_pos;
  uint8_t nhitpos; /* number of hit positions */
  uint16_t score;
  uint16_t gaplen;
  genomeloc_t nextpos;
  tagdata50 seq;
} tag50align;
typedef struct tag50alignloc_{
  tag50align tag;
  genomeloc_t location;
}tag50alignlocation;
typedef struct matchtaglist_{
  struct matchtaglist_* next;
  genomeloc_t location; /* for genomes no more than 2G */
  int allocated_array_size; /*capacity of the allocated array*/
  int array_content_number; /*actual number of aligned tag*/
  tag50align *alignedtag; /* the content is relocatable*/
} matchtaglist;

typedef struct matchtagindex_{
  uint64_t genomesize;
  uint64_t arraysize;
  uint64_t upperbitmask;
  uint64_t lowerbitmask;
  int shiftbitcount;
  matchtaglist** pary;
} matchtagindex;

matchtagindex* initmatchtagindex(matchtagindex*, uint64_t, genomeloc_t);
int inserttag(matchtagindex*, int, tag50align, uint64_t*totalpos);
inline tag50alignlocation str2tag50alignlocation(const char*);
tagdata50 str2tag50(const char*);
int printtagdata50(FILE*out, tagdata50 tag);
char* tagdata50tostr(char*out, tagdata50 tag);
int* tagdata50tointarray(int*out, tagdata50 tag);

#endif
