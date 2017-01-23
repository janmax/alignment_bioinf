#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdbool.h>

#define info(M, ...) fprintf(stdout, "[INFO] " M "\n", ##__VA_ARGS__)
#define max(a,b,c) (((a)>(b))?((a)>(c)?(a):(c)):((b)>(c)?(b):(c)))
#define s(x,y) blosum_mat[(x)][(y)]

#define HIGHEST_CHAR 'Z'
#define MAX_SEQSIZE 4096
#define NUCLIDS 20

enum {DIAGONAL, UP, LEFT, STOP};

struct sequence {
    size_t length;
    unsigned char sequence[MAX_SEQSIZE];
    char seq_name[32];
};

struct alignment {
    size_t length;
    int32_t score;
    char first[MAX_SEQSIZE * 2];
    char second[MAX_SEQSIZE * 2];
};

char blosum_mat[HIGHEST_CHAR][HIGHEST_CHAR];
extern const unsigned char nuclids[];

#endif // ALIGNMENT_H
