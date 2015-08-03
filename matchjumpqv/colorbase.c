/* colorbase.c */
/* defines structures and functions associated with base and color sequence*/
#include "colorbase.h"
const char* nuc2char="ACGT";
const int dibasecolortable[16] = {0, 1, 2, 3, 1, 0, 3, 2, 2, 3, 0, 1, 3, 2, 1, 0};
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

char char2nucmap[CHAR_MAX+1];

void
initchar2nucmap()
{
  int i;
  for(i=0; i <= CHAR_MAX; i++) char2nucmap[i]=0;
  char2nucmap['a'] = 4 | 0;
  char2nucmap['A'] = 4 | 0;
  char2nucmap['c'] = 4 | 1;
  char2nucmap['C'] = 4 | 1;
  char2nucmap['g'] = 4 | 2;
  char2nucmap['G'] = 4 | 2;
  char2nucmap['t'] = 4 | 3;
  char2nucmap['T'] = 4 | 3;
  char2nucmap['u'] = 4 | 3;
  char2nucmap['U'] = 4 | 3;
  char2nucmap['N'] = 8;
  char2nucmap['n'] = 8;
  char2nucmap['Y'] = 8;
  char2nucmap['y'] = 8;
  char2nucmap['r'] = 8;
  char2nucmap['R'] = 8;
  char2nucmap['s'] = 8;
  char2nucmap['S'] = 8;
  char2nucmap['w'] = 8;
  char2nucmap['W'] = 8;
  char2nucmap['b'] = 8;
  char2nucmap['B'] = 8;
  char2nucmap['m'] = 8;
  char2nucmap['M'] = 8;
  char2nucmap['h'] = 8;
  char2nucmap['H'] = 8;
  char2nucmap['v'] = 8;
  char2nucmap['V'] = 8;
}

