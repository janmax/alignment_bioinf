#include "print_functions.h"

inline void print_blosum()
{
    printf("%3s", "");
    for (unsigned i = 0; i < NUCLIDS; ++i)
        printf("%3c", nuclids[i]);
    printf("\n");
    for (unsigned i = 0; i < NUCLIDS; ++i) {
        printf("%3c", nuclids[i]);
        for (unsigned j = 0; j < NUCLIDS; ++j) {
            printf("%3d", blosum_mat[nuclids[i]][nuclids[j]]);
        }
        printf("\n");
    }
}

inline void print_matrix(int32_t **F,
                         int32_t **Z,
                         uint32_t n,
                         uint32_t m)
{
    for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < m; ++j) {
            if (j < m - 1)
                printf("%4d%4s", F[i][j] ? F[i][j] : 0, Z[i][j + 1] == LEFT ? "<" : " " );
            else
                printf("%4d", F[i][j] ? F[i][j] : 0);
        }
        printf("\n");
        for (unsigned j = 0; j < m && i < n - 1; ++j) {
            if (j < m - 1)
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
    printf("%s\n", a->first);
    printf("%s\n", a->second);
    printf("Alignment score --> %d\n", a->score);
}
