#ifndef PRINT_FUNCTIONS_H
#define PRINT_FUNCTIONS_H

#include <stdlib.h>
#include "alignment.h"

void print_blosum();
void print_matrix(int32_t **F,
                                int32_t **Z,
                                uint32_t n,
                                uint32_t m);
void print_alignment(struct alignment * a);

#endif // PRINT_FUNCTIONS_H
