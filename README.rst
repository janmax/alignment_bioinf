****************************************************************
Implementation of Needleman-Wnusch and Smith-Waterman algorithms
****************************************************************

Created by Jan Maximilian Michal as part of the lecture Algortithmen der
Bioinformatik by Prof. Dr Morgenstern of the University of Göttingen.

How to use the code?
====================

There is not much to see. Start by compiling the code. There are few
dependencies that should all be included in standard UNIX systems.

.. code-block:: bash

    make
    ./main -i fasta/demo.fasta

The program will take a BLOSUM matrix a fasta file with two alignments and find
a global and a local alignment for these sequences. The dynamic programming
matrix will be printed along with the optimal alignment and some information
about it.

You can use the following options to customize its behavior:

.. code-block:: text

    usage: ./main [-h]
        [-g gap_penalty]
        [-i fasta_file]
        [-o output (default: stdout)]
        [-b BLOSUM_file (default: BLOSUM62)]
        [-p] # prints the BLOSUM matrix.

The examples in Durbin et. al. can be reproduced by invoking

.. code-block:: bash

    ./main b BLOSUM50 -i fasta/demo.fasta -g 8

How is the code organized?
==========================

The structure is the following:

- **alignment.h** declares and documents the logic for the implemented
algorithms.

- **print_functions.h** just declares some helper functions that print out
some structs nicely. Very helpful for debugging as well.

- **main.c** glues everything together, reads data and passes it to the
algorithms, which therefore could be reused in another project easily.

Static functions are usually not commented, since they just encapsulate behavior
that would distract from the relevant code parts.

