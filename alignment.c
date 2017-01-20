#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>

#define info(M, ...) fprintf(stdout, "[INFO] " M "\n", ##__VA_ARGS__)
#define max(a,b,c) (((a)>(b))?((a)>(c)?(a):(c)):((b)>(c)?(b):(c)))
#define HIGHEST_CHAR 'Z'
#define MAX_SEQSIZE 4096
#define s(x,y) blosum_mat[x][y]

enum {DIAGONAL, UP, DOWN};

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

static inline void print_alignment(struct alignment * a)
{
    printf("Found an alignment of length %ld:\n", a->length);
    printf("%s\n", a->first);
    printf("%s\n", a->second);
    printf("Alignment score --> %d\n", a->score);
}

static inline void fill_blosum62(char blosum_mat[HIGHEST_CHAR][HIGHEST_CHAR])
{
    FILE *fp = fopen("BLOSUM62", "r");
    if (fp == NULL) {
        printf ("File not created okay, errno = %d\n", errno);
        exit(errno);
    }

    unsigned char p1, p2;
    int penalty;
    while (fscanf(fp, "%c %c %d\n", &p1, &p2, &penalty) != EOF) {
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

static inline int32_t * allocate_2darray(int32_t ***arr, int32_t n, int32_t m)
{
    *arr = (int32_t**) calloc(1, n * sizeof(int32_t*));
    int32_t *arr_data = calloc(1, n * m * sizeof(int32_t));
    for (int32_t i = 0; i < n; i++)
        (*arr)[i] = arr_data + i * m;
    return arr_data;
}

static inline void deallocate_2darray(int32_t*** arr, int32_t* arr_data)
{
    free(arr_data);
    free(*arr);
}

struct alignment * needleman_wunsch(int32_t d, // gap_penalty
                                    const struct sequence s1,
                                    const struct sequence s2)
{
    struct alignment *a = malloc(sizeof(struct alignment));
    uint32_t n = s1.length;
    uint32_t m = s2.length;

    int32_t **F, **Z;
    int32_t *F_free, *Z_free;
    F_free = allocate_2darray(&F, n, m);
    Z_free = allocate_2darray(&Z, n, m);

    // initialize borders
    for (uint32_t i = 1; i < n; ++i)
        F[i][0] = F[i - 1][0] - d;
    for (uint32_t j = 1; j < m; ++j)
        F[0][j] = F[0][j - 1] - d;

    // fill the 2d matrix and prepare backtracking
    for (unsigned i = 1; i < n; ++i) {
        for (unsigned j = 1; j < m; ++j) {
            int32_t f1 = F[i - 1][j - 1] + s(s1.sequence[i], s2.sequence[j]);
            int32_t f2 = F[i - 1][j] - d;
            int32_t f3 = F[i][j - 1] - d;
            int32_t max = F[i][j] = max(f1, f2, f3);
            if (max == f1)
                Z[i][j] = DIAGONAL;
            else if (max == f2)
                Z[i][j] = DOWN; // has to be swaped
            else
                Z[i][j] = UP;
        }
    }

    // backtrack the best match
    char r1[n + m];
    char r2[n + m];
    unsigned index = 0;
    unsigned i = n - 1, j = m - 1;
    while (i + j != 0) {
        if (Z[i][j] == DIAGONAL) {
            r1[index] = s1.sequence[i];
            r2[index] = s2.sequence[j];
            i--; j--;
        } else if (Z[i][j] == UP) {
            r1[index] = '-';
            r2[index] = s2.sequence[j];
            j--;
        } else {
            r1[index] = s1.sequence[i];
            r2[index] = '-';
            i--;
        }
        index++;
    }

    // write the alignment struct
    for (unsigned i = 0; i < index; ++i) {
        a->first[i] = r1[index - i - 1];
        a->second[i] = r2[index - i - 1];
    }
    a->first[index] = '\0';
    a->second[index] = '\0';
    a->score = F[n - 1][m - 1];
    a->length = index;

    // free the used arrays
    deallocate_2darray(&F, F_free);
    deallocate_2darray(&Z, Z_free);
    return a;
}


int main(int argc, char *const *argv)
{
    int c;
    int gap_penalty = 1;
    FILE *input;
    FILE *output = stdout;
    char *input_name = NULL;

    /*** begin arg parsing ****************************************************/
    while ((c = getopt (argc, argv, "g:o:i:")) != -1) {
        switch (c) {
        case 'g':
            gap_penalty = strtol(optarg, NULL, 10);
            break;
        case 'o':
            output = fopen(optarg, "w");
            break;
        case 'i':
            input_name = optarg;
            break;
        }
    }

    if (input_name == NULL) {
        info("You need to specify an input filename in fasta format.");
        return 0;
    } else {
        input = fopen(input_name, "r");
    }
    /*** finished arg parsing *************************************************/

    // read gap penalties from BLOSUM file
    fill_blosum62(blosum_mat);

    // read two sequences
    struct sequence s1;
    struct sequence s2;
    read_sequence(input, &s1);
    read_sequence(input, &s2);
    fclose(input);

    // get global alignment with needleman wunsch
    struct alignment * a = needleman_wunsch(gap_penalty, s1, s2);

    print_alignment(a);
    free(a);
    return 0;
}
