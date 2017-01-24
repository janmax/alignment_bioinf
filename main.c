#include "alignment.h"
#include "print_functions.h"

int main(int argc, char *const *argv)
{
    // default values
    int gap_penalty = 1;
    bool blosum_flag = false;
    char *BLOSUM_path = "BLOSUM62";
    struct sequence s[2];
    struct alignment nw, sm;

    // necessary for arg parsing
    int arg;
    FILE *input = NULL;
    FILE *output = stdout;

    /*** begin arg parsing ****************************************************/
    while ((arg = getopt (argc, argv, "g:o:i:b:ph")) != -1) {
        switch (arg) {
        case 'g':
            gap_penalty = strtol(optarg, NULL, 10);
            break;
        case 'o':
            output = fopen(optarg, "w");
            break;
        case 'i':
            input = fopen(optarg, "r");
            break;
        case 'b':
            BLOSUM_path = optarg;
            break;
        case 'p':
            blosum_flag = true;
            break;
        default:
        case 'h':
            print_help();
            return 0;
        }
    }

    // check for input arg
    if (input == NULL) {
        info("You need to specify an input filename in fasta format.");
        return 0;
    }
    /*** finished arg parsing *************************************************/

    // read gap penalties from BLOSUM file and print if requested
    fill_blosum(BLOSUM_path);
    if (blosum_flag)
        print_blosum();

    // read two sequences
    for (unsigned i = 0; i < 2; ++i)
        read_sequence(input, &s[i]);
    fclose(input);

    // get global alignment with needleman wunsch
    printf("Needleman-Wunsch output:\n");
    needleman_wunsch(&nw, gap_penalty, &s[0], &s[1]);
    print_alignment(&nw);

    // get local alignment with smith waterman
    printf("\nSmith-Waterman output:\n");
    smith_waterman(&sm, gap_penalty, &s[0], &s[1]);
    print_alignment(&sm);

    // all done -- exit
    return 0;
}
