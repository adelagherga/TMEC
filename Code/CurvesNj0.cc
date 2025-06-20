// Program to list curves of given conductor with j-invariant 0

// Usage:
// ./CurvesNj0 <N>  # for one conductor N
// ./CurvesNj0 <N1> <N2>  # for all conductors from N1 to N2 inclusive

// NB the number of parameters on the command line is 1 or 2

// Output: one curve per line, conductor N then a-invariants [a1,a2,a3,a4,a6]

// Example:
//
// $ ./CurvesNj0 1 100
// 27 [0,0,1,0,0]
// 27 [0,0,1,0,-7]
// 36 [0,0,0,0,1]
// 36 [0,0,0,0,-27]

#include <eclib/egros.h>
#include "egros_cache.h"

int main (int argc, char *argv[])
{
  if ( (argc < 2) || (argc > 3) )
    {
      cerr << "Usage: CurvesNj0 N or CurvesNj0 N1 N2" << endl;
      return 0;
    }
  long n1, n2;
  char* t;
  n1 = strtol(argv[1], &t, 10); // 10 is the base
  if (*t)
    {
      cerr << "Usage: CurvesNj0 N or CurvesNj0 N1 N2" << endl;
      return 0;
    }
  if (argc==3)
    {
      n2 = strtol(argv[2], &t, 10); // 10 is the base
      if (*t)
        {
          cerr << "Usage: CurvesNj0 N or CurvesNj0 N1 N2" << endl;
          return 0;
        }
    }
  else
    n2 = n1;

  initprimes("PRIMES");

  // skip to next multiple of 3:
  while (n1%3!=0)
    n1++;
  for (long n=n1; n<=n2; n+=3)
    {
      bigint N(n);
      vector<bigint> suppN = pdivs(N);
      if (is_N_possible_j_0(N, suppN))
        {
          vector<CurveRed> Elist = get_egros_from_j_0(suppN);
          for (auto E: Elist)
            if (E.conductor() == N)
              cout << N << " " << (Curve)E <<endl;
        }
    }
}
