#include "print_functions.h"

inline void print_blosum()
{
    printf("%3s", "");
    for (unsigned i = 0; i < NUCLIDS; ++i)
        printf("%3c", amino_acids[i]);
    printf("\n");
    for (unsigned i = 0; i < NUCLIDS; ++i) {
        printf("%3c", amino_acids[i]);
        for (unsigned j = 0; j < NUCLIDS; ++j) {
            printf("%3d", blosum_mat[amino_acids[i]][amino_acids[j]]);
        }
        printf("\n");
    }
}

inline void print_matrix(int32_t **F,
                         int32_t **Z,
                         uint32_t n,
                         uint32_t m)
{
    for (unsigned i = 0; i <= n; ++i) {
        for (unsigned j = 0; j <= m; ++j) {
            if (j < m)
                printf("%4d%4s", F[i][j] ? F[i][j] : 0, Z[i][j + 1] == LEFT ? "<" : " " );
            else
                printf("%4d", F[i][j] ? F[i][j] : 0);
        }
        printf("\n");
        for (unsigned j = 0; j <= m && i < n; ++j) {
            if (j < m)
                printf("%4s%4s", Z[i + 1][j] == UP ? "^" : " ", Z[i + 1][j + 1] == DIAGONAL ? "\\" : " ");
            else
                printf("%4s", Z[i + 1][j] == UP ? "^" : " ");
        }
        printf("\n");
    }
}


inline void print_alignment(struct alignment * a)
{
    printf("Found an alignment of length %ld:\n", a->length);
    printf("%s\n", a->upper);
    printf("%s\n", a->lower);
    printf("Alignment score --> %d\n", a->score);
}

inline void print_help()
{
    printf("usage: ./main [-h]\n");
    printf("    [-g gap_penalty]\n");
    printf("    [-i fasta_file]\n");
    printf("    [-o output (default: stdout)]\n");
    printf("    [-b BLOSUM_file (default: BLOSUM62)]\n");
    printf("    [-p] -- print BLOSUM matrix.\n");
}
