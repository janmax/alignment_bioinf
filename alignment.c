#include "alignment.h"
#include "print_functions.h"

const unsigned char amino_acids[] = "ARNDCQEGHILKMFPSTWYV";

/*** static functions begin ***************************************************/
static void strrev(char *p)
{
    char *q = p;
    while (q && *q) ++q;
    for (--q; p < q; ++p, --q)
        *p = *p ^ *q,
         *q = *p ^ *q,
          *p = *p ^ *q;
}

static int32_t * allocate_2darray(int32_t ***arr, int32_t n, int32_t m)
{
    *arr = (int32_t**) calloc(1, n * sizeof(int32_t*));
    int32_t *arr_data = calloc(1, n * m * sizeof(int32_t));
    for (int32_t i = 0; i < n; i++)
        (*arr)[i] = arr_data + i * m;
    return arr_data; // free point
}

static void deallocate_2darray(int32_t*** arr, int32_t* arr_data)
{
    free(arr_data);
    free(*arr);
}
/*** static functions end *****************************************************/


/*** public member functions **************************************************/
void fill_blosum(char * filename)
{
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        exit(errno);
    }

    unsigned char p1, p2;
    int32_t penalty;
    while (fscanf(fp, " %c %c %d\n", &p1, &p2, &penalty) != EOF) {
        blosum_mat[p2][p1] = penalty;
    }
    fclose(fp);
}

void read_sequence(FILE * file_pos, struct sequence * s)
{
    fscanf(file_pos, ">%s\n", s->seq_name);
    fscanf(file_pos, "%s\n", s->sequence);
    s->length = strlen((char *) s->sequence);
}

void needleman_wunsch(struct alignment *a,
                      int32_t d, // gap_penalty
                      const struct sequence * s1,
                      const struct sequence * s2)
{
    const uint32_t n = s1->length;
    const uint32_t m = s2->length;

    int32_t **F, **Z;
    int32_t *F_free, *Z_free;
    F_free = allocate_2darray(&F, n + 1, m + 1);
    Z_free = allocate_2darray(&Z, n + 1, m + 1);

    // initialize borders
    for (uint32_t i = 1; i <= n; ++i) {
        F[i][0] = F[i - 1][0] - d;
        Z[i][0] = UP;
    }
    for (uint32_t j = 1; j <= m; ++j) {
        F[0][j] = F[0][j - 1] - d;
        Z[0][j] = LEFT;
    }

    // fill the 2d matrix and prepare backtracking
    for (unsigned i = 1; i <= n; ++i) {
        for (unsigned j = 1; j <= m; ++j) {
            int32_t f1 = F[i - 1][j - 1] + s(s1->sequence[i - 1], s2->sequence[j - 1]);
            int32_t f2 = F[i - 1][j] - d;
            int32_t f3 = F[i][j - 1] - d;
            int32_t max = F[i][j] = max(f1, f2, f3);
            if (max == f1)
                Z[i][j] = DIAGONAL;
            else if (max == f2)
                Z[i][j] = UP; // has to be swaped
            else
                Z[i][j] = LEFT;
        }
    }

    // backtrack the best match
    unsigned index = 0;
    int32_t i = n, j = m;
    while (i >= 0 && j >= 0) {
        if (Z[i][j] == DIAGONAL) {
            a->upper[index] = s1->sequence[i - 1];
            a->lower[index] = s2->sequence[j - 1];
            i--; j--;
        } else if (Z[i][j] == UP) {
            a->upper[index] = s1->sequence[i - 1];
            a->lower[index] = '-';
            i--;
        } else {
            a->upper[index] = '-';
            a->lower[index] = s2->sequence[j - 1];
            j--;
        }
        index++;
    }

    strrev(a->upper);
    strrev(a->lower);
    a->score = F[n][m];
    a->length = index;

    if (a->length < 20)
        print_matrix(F, Z, n, m);

    // free the used arrays
    deallocate_2darray(&F, F_free);
    deallocate_2darray(&Z, Z_free);
}

void smith_waterman(struct alignment *a,
                    int32_t d, // gap_penalty
                    const struct sequence * s1,
                    const struct sequence * s2)
{
    uint32_t n = s1->length;
    uint32_t m = s2->length;

    // arrays are allocated with calloc and therefore zeroed
    int32_t **F, **Z;
    int32_t *F_free, *Z_free;
    F_free = allocate_2darray(&F, n + 1, m + 1);
    Z_free = allocate_2darray(&Z, n + 1, m + 1);

    // final index
    int32_t max_i = 0;
    int32_t max_j = 0;

    // fill the 2d matrix and prepare backtracking
    int32_t f1, f2, f3, * max;
    for (unsigned i = 1; i <= n; ++i) {
        for (unsigned j = 1; j <= m; ++j) {
            max = &F[i][j];
            f1 = F[i - 1][j - 1] + s(s1->sequence[i - 1], s2->sequence[j - 1]);
            f2 = F[i - 1][j] - d;
            f3 = F[i][j - 1] - d;
            *max = max(f1, f2, f3);
            if (*max <= 0) {
                *max = 0;
                Z[i][j] = STOP;
            } else if (*max == f1)
                Z[i][j] = DIAGONAL;
            else if (*max == f2)
                Z[i][j] = UP; // has to be swaped
            else
                Z[i][j] = LEFT;

            // track maximum
            if (F[i][j] >= F[max_i][max_j]) {
                max_i = i; max_j = j;
            }
        }
    }

    // backtrack the best match
    unsigned index = 0;
    int32_t i = max_i;
    int32_t j = max_j;
    while (F[i][j] && i >= 0 && j >= 0) {
        if (Z[i][j] == DIAGONAL) {
            a->upper[index] = s1->sequence[i - 1];
            a->lower[index] = s2->sequence[j - 1];
            i--; j--;
        } else if (Z[i][j] == UP) {
            a->upper[index] = s1->sequence[i - 1];
            a->lower[index] = '-';
            i--;
        } else {
            a->upper[index] = '-';
            a->lower[index] = s2->sequence[j - 1];
            j--;
        }
        index++;
    }

    // write the alignment struct
    strrev(a->upper);
    strrev(a->lower);
    a->score = F[n][m];
    a->length = index;

    if (a->length < 20)
        print_matrix(F, Z, n, m);

    // free the used arrays
    deallocate_2darray(&F, F_free);
    deallocate_2darray(&Z, Z_free);
}
/*** public member functions end **********************************************/
