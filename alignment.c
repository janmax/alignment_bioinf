#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>

#define info(message, ...) printf(message, "\n", __VA_ARGS__);

static inline void fill_blosum62(char blosum_mat[128][128])
{
    FILE *fp = fopen("BLOSUM62", "r");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        exit(errno);
    }

    unsigned char p1, p2;
    int penalty;
    while(fscanf(fp, "%c %c %d\n", &p1, &p2, &penalty) != EOF) {
        blosum_mat[p2][p1] = penalty;
    }
}

int main(int argc, char *const *argv)
{
    int c;
    int gap_penalty = -1;
    FILE *output = stdout;
    char *input_value;
    while((c = getopt (argc, argv, "g:o:i:")) != -1) {
        switch (c) {
        case 'g':
            gap_penalty = strtol(optarg, NULL, 10);
            break;
        case 'o':
            output = fopen(optarg, "w");
            break;
        case 'i':
            input_value = optarg;
            break;
        }
    }

    if (input_value == NULL)
        printf("[INFO] You need to specify an input filename in fasta format.\n");

    printf("Gap %d\n", gap_penalty);
    printf("In %s\n", input_value);
    fprintf(output, "Out %s\n", "DEMO"); // where the code goes

    char blosum_mat[128][128];
    fill_blosum62(blosum_mat);

    printf("penalty %d\n", blosum_mat['P']['E']);
    printf("penalty %d\n", blosum_mat['P']['E']);
    return 0;
}