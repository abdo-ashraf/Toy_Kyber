#include "toy.h"


void toy_polmul_naive(short *dst, short *a, short *b, int add_to_dst){
    short tempo[4];

    tempo[0] = (a[0]*b[0] - a[3]*b[1] - a[2]*b[2] - a[1]*b[3]) % TK_Q;
    tempo[1] = (a[1]*b[0] + a[0]*b[1] - a[3]*b[2] - a[2]*b[3]) % TK_Q;
    tempo[2] = (a[2]*b[0] + a[1]*b[1] + a[0]*b[2] - a[3]*b[3]) % TK_Q;
    tempo[3] = (a[3]*b[0] + a[2]*b[1] + a[1]*b[2] + a[0]*b[3]) % TK_Q;

    if(add_to_dst==1){
        dst[0] += tempo[0];
        dst[1] += tempo[1];
        dst[2] += tempo[2];
        dst[3] += tempo[3];
    }
    else{
        dst[0] = tempo[0];
        dst[1] = tempo[1];
        dst[2] = tempo[2];
        dst[3] = tempo[3];
    }
}


void toy_poladd(short *d, short *t, short *e){
    d[0] = (t[0] + e[0]) % TK_Q;
    d[1] = (t[1] + e[1]) % TK_Q;
    d[2] = (t[2] + e[2]) % TK_Q;
    d[3] = (t[3] + e[3]) % TK_Q;
}


void toy_poldiff(short *d, short *t, short *e){
    d[0] = (t[0] - e[0]) % TK_Q;
    d[1] = (t[1] - e[1]) % TK_Q;
    d[2] = (t[2] - e[2]) % TK_Q;
    d[3] = (t[3] - e[3]) % TK_Q;
}


// p = a.b
void mat_vct_mul(short *p, short *a, short *b){
    short r_conv[TK_K*TK_K*TK_N];
    // convolute vertor on matrix
    for (int i=0; i < TK_K; i++) {
            toy_polmul_naive(&r_conv[(i*TK_K*TK_N) + (0*TK_N)], &a[(i*TK_K*TK_N) + (0*TK_N)], &b[0*TK_N], 0);
            toy_polmul_naive(&r_conv[(i*TK_K*TK_N) + (1*TK_N)], &a[(i*TK_K*TK_N) + (1*TK_N)], &b[1*TK_N], 0);
            toy_polmul_naive(&r_conv[(i*TK_K*TK_N) + (2*TK_N)], &a[(i*TK_K*TK_N) + (2*TK_N)], &b[2*TK_N], 0);
    }
    // sum r_conv on axis 0
        toy_poladd(&p[0*TK_N], &r_conv[0*TK_N], &r_conv[1*TK_N]);
        toy_poladd(&p[0*TK_N], &p[0*TK_N], &r_conv[2*TK_N]);

        toy_poladd(&p[1*TK_N], &r_conv[3*TK_N], &r_conv[4*TK_N]);
        toy_poladd(&p[1*TK_N], &p[1*TK_N], &r_conv[5*TK_N]);

        toy_poladd(&p[2*TK_N], &r_conv[6*TK_N], &r_conv[7*TK_N]);
        toy_poladd(&p[2*TK_N], &p[2*TK_N], &r_conv[8*TK_N]);
}


static void vector_fill(short *buf, int n){
    for (int k = 0; k < n; ++k) {
        short val = rand();
        val=(val>>1&1)-(val&1);
        if(val<0)
            val+=TK_Q;

        buf[k]=val;
    }
}


void toy_gen(short *A, short *t, short *s)
{
    short e1[TK_K*TK_N];
    for (int k = 0; k < TK_K*TK_K*TK_N; ++k) {
        A[k]=rand()%TK_Q;
    }
    vector_fill(s, TK_K*TK_N);
    vector_fill(e1, TK_K*TK_N);

    // t = A.s + e1
    mat_vct_mul(t, A, s); // t = A.s
    for(int i=0; i<TK_K*TK_N; i+=TK_N){
        toy_poladd(&t[i], &t[i], &e1[i]);
    }
}


void swapm(short *mat, int z1, int z2){
    z1 = z1*TK_N;
    z2 = z2*TK_N;
    short temp[TK_N];

    temp[0] = mat[z1+0];
    temp[1] = mat[z1+1];
    temp[2] = mat[z1+2];
    temp[3] = mat[z1+3];

    mat[z1+0] = mat[z2+0];
    mat[z1+1] = mat[z2+1];
    mat[z1+2] = mat[z2+2];
    mat[z1+3] = mat[z2+3];

    mat[z2+0] = temp[0];
    mat[z2+1] = temp[1];
    mat[z2+2] = temp[2];
    mat[z2+3] = temp[3];
}


void transpose(short *matrix) {
    swapm(matrix, 1,3);
    swapm(matrix, 5,7);
    swapm(matrix, 2,6);
}


void toy_enc( short *A,  short *t, int plain, short *u, short *v){

    short r[TK_K*TK_N];
    short e1[TK_K*TK_N];
    short e2[TK_N];
    for (int i = 0; i < TK_K*TK_N; i++) {
        int val = rand() & 3; // Uniform distribution
        r[i] = (val & 1) - ((val >> 1) & 1); // Binomial distribution
        r[i] %= TK_Q;
    }
    for (int i = 0; i < TK_K*TK_N; i++) {
        int val = rand() & 3; // Uniform distribution
        e1[i] = (val & 1) - ((val >> 1) & 1); // Binomial distribution
        e1[i] %= TK_Q;
    }
    for (int i = 0; i < TK_N; i++) {
        int val = rand() & 3; // Uniform distribution
        e2[i] = (val & 1) - ((val >> 1) & 1); // Binomial distribution
        e2[i] %= TK_Q;
    }

    transpose(A);
    mat_vct_mul(u, A,r);
    toy_poladd(u,u,e1);
    
    short msg_bits[psize];
    for(int i=0; i<psize; i++){
        msg_bits[i] = (plain>>i&1) * (TK_Q/2);
    }

    short temp[TK_K*TK_N];
    short vv[TK_N];
    toy_polmul_naive(&temp[0], &t[0], &r[0],0);
    toy_polmul_naive(&temp[4], &t[4], &r[4],0);
    toy_polmul_naive(&temp[8], &t[8], &r[8],0);

    toy_poladd(&vv[0*TK_N], &temp[0*TK_N], &temp[1*TK_N]);
    toy_poladd(&vv[0*TK_N], &vv[0*TK_N], &temp[2*TK_N]);

    toy_poladd(e2, vv, e2);
    toy_poladd(v, msg_bits, e2);
}


int toy_dec(short *s, short *u,  short *v){
    short temp[TK_K*TK_N];
    short vv[TK_N];
    toy_polmul_naive(&temp[0], &s[0], &u[0],0);
    toy_polmul_naive(&temp[4], &s[4], &u[4],0);
    toy_polmul_naive(&temp[8], &s[8], &u[8],0);

    toy_poladd(&vv[0*TK_N], &temp[0*TK_N], &temp[1*TK_N]);
    toy_poladd(&vv[0*TK_N], &vv[0*TK_N], &temp[2*TK_N]);

    short p[TK_N];

    toy_poldiff(p, v, vv);

    int plain=0;
    int val;
    int bit;
    for (int i=0;i<psize; i++){
        val = p[i];
        if(val>TK_Q/2)
            val -= TK_Q;
        bit = abs(val)>TK_Q/4;
        plain |= bit<<i;
    }
    return plain;
}