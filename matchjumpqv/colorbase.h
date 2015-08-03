/* colorbase.h */
/* defines structures and functions associated with base and color sequence*/
#if !defined(colorbase_h_included)
#define colorbase_h_included
#include<stdio.h>
#include<stdint.h>
#include<limits.h>

extern char char2nucmap[CHAR_MAX+1];
extern const char* nuc2char;
extern const int dibasecolortable[16];
/*
 * 'aa' => '0', 'ac' => '1', 'ag' => '2', 'at' => '3',
 * 'ca' => '1', 'cc' => '0', 'cg' => '3', 'ct' => '2',
 * 'ga' => '2', 'gc' => '3', 'gg' => '0', 'gt' => '1',
 * 'ta' => '3', 'tc' => '2', 'tg' => '1', 'tt' => '0'
 */
/* 'a0' => 'a', 'a1' => 'c', 'a2' => 'g', 'a3' => 't',
 * 'c0' => 'c', 'c1' => 'a', 'c2' => 't', 'c3' => 'g',
 * 'g0' => 'g', 'g1' => 't', 'g2' => 'a', 'g3' => 'c',
 * 't0' => 't', 't1' => 'g', 't2' => 'c', 't3' => 'a'
 */
extern void initchar2nucmap();
static inline int
getcolor(int prevbase, int curbase)
{
  return dibasecolortable[((prevbase & 3) << 2 ) | (curbase & 3)];
}
static inline int
getnextbase(int prevbase, int color)
{
  return dibasecolortable[((prevbase & 3) << 2 ) | (color & 3)];
}

static inline int
is_validnuc(int c)
{
  if((c == 'T') || (c == 'G') || (c == 'C') ||(c == 'A'))
    return 1;
  else
    return 0;
}
static inline int
char2color(int c)
{
  if(c >= '0' && c <= '3')
    return c-'0';
  else
    return 0;
}
static inline void
printasns(FILE* out, uint64_t ns, int start, int end)
{
  int i;
  ns >>= start*2;
  for(i = start; i < end; i++){
    int curnuc = ns & 3;
    fputc(nuc2char[curnuc], out);
    ns >>=2;
  }
/*  fputc('\n', out);*/
}
static inline void
printascolor(FILE* out, uint64_t cs, int start, int end)
{
  int i;
  cs >>= start*2;
  for(i = start; i < end; i++){
    int curcolor = cs & 3;
    fputc('0' + curcolor, out);
    cs >>=2;
  }
/*  fputc('\n', out);*/
}
static inline uint64_t
reversecolor(uint64_t cs, int n)
{
  uint64_t retval = 0;
  uint64_t inval = cs;
  int i;
  for(i = n; i > 0; i--){
    retval = retval << 2 | (inval & 3);
    inval >>= 2;
  }
  return retval;
}


#endif

