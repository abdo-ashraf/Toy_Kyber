#include "toy.h"

int main(int argc, char *argv[]){
    if(argc != 2){
        printf("plain text needed");
        return 0;
    }

    int int_plain = atoi(argv[1]);
    printf("in_plain: %d\n", int_plain);
    srand (time(NULL));
    short A[TK_K*TK_K*TK_N];
    short t[TK_K*TK_N];
    short s[TK_K*TK_N];
    toy_gen(A, t, s);
    printf("A: %d t: %d s: %d\n", A, t, s);

    short u[TK_K*TK_N];
    short v[TK_N];
    toy_enc(A, t, int_plain, u, v);
    printf("u: %d v: %d \n", u, v);

    int plain = toy_dec(s, u, v);
    printf("out_plain: %d", plain);
    return 1;
}