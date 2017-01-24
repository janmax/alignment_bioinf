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

char blosum_mat[HIGHEST_CHAR][HIGHEST_CHAR];
extern const unsigned char amino_acids[];

struct alignment {
    size_t length;
    int32_t score;
    char upper[MAX_SEQSIZE * 2];
    char lower[MAX_SEQSIZE * 2];
};

/**
 * This function implements the Needleman-Wunsch algorithm as suggested in
 * Durbin et. al. -- Biological sequence analysis.
 *
 * @param a  The address of an empty struct alignment
 *           where the solution will be written to
 * @param d  The chosen gap penalty
 * @param s1 A pointer to a struct sequence that contains a sequence
 * @param s2 A pointer to a struct sequence that contains a sequence
 */
void needleman_wunsch(struct alignment * a,
                      int32_t d, // gap_penalty
                      const struct sequence * s1,
                      const struct sequence * s2);

/**
 * This function implements the Smith-Waterman algorithm as suggested in
 * Durbin et. al. -- Biological sequence analysis.
 *
 * @param a  The address of an empty struct alignment
 *           where the solution will be written to
 * @param d  The chosen gap penalty
 * @param s1 A pointer to a struct sequence that contains a sequence
 * @param s2 A pointer to a struct sequence that contains a sequence
 */
void smith_waterman(struct alignment * a,
                    int32_t d, // gap_penalty
                    const struct sequence * s1,
                    const struct sequence * s2);

/**
 * Fills the BLOSUM matrix that is defined above.
 *
 * @param filename A path to a file that contains a BLOSUM matrix as letter
 *                 tuples with an asociated score value.
 */
void fill_blosum(char * filename);

/**
 * Reading a sequence from a fasta file.
 * @param file_pos A vaild FILE pointer that points to the position of
 *                 a line in a fasta file beginning with > NAME
 * @param s        The struct sequence to fill with data.
 */
void read_sequence(FILE * file_pos, struct sequence * s);

#endif // ALIGNMENT_H
