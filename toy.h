#include <stdio.h>
#include <stdlib.h>
#include <time.h>


#define TK_K 3
#define TK_N 4
#define TK_Q 97
#define psize 4


void toy_polmul_naive(short *dst, short *a, short *b, int add_to_dst);
void toy_poladd(short *d, short *t, short *e);
void toy_poldiff(short *d, short *t, short *e);
void mat_vct_mul(short *p, short *a, short *b);
void toy_gen(short *A, short *t, short *s);
void swapm(short *mat, int z1, int z2);
void transpose(short *matrix);
void toy_enc( short *A,  short *t, int plain, short *u, short *v);
int toy_dec(short *s, short *u,  short *v);
