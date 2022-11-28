These files generate all elliptic curves of given conductor(s) by solving a
large list of corresponding Thue--Mahler equations in parallel. This program
takes as input a single conductor N, a range of conductors [N1,...,N2], or, with
the flag -l, an arbitrary finite list of conductors:

Parameters
    N1 [N2]
        A single conductor, N1, or the range [N1,...,N2].
    [-l N1] [N2...]: [OPTIONAL]
        A finite, arbitrary list of conductors to generate [N1,N2,...].

These functions print all output, logfiles, and errors to a directory called
Data/${name}, created as the program proceeds. If such a directory alread exists
in Data/, these functions generate a directory Data/${name}i, where
Data/${name}j exists for all j < i.

To run this program using GNU parallel, running at most 40 jobs in parallel, run

$ cd TMEC
$ nohup Code/computeEllipticCurves.sh N1 [N2] &
