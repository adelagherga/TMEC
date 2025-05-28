// Program to list curves of given conductor and j-invariant

// Usage:
// ./CurvesNj <j> <N>  # for one conductor N
// ./CurvesNj <j> <N1> <N2>  # for all conductors from N1 to N2 inclusive

// NB the number of parameters on the command line is 2 or 3.  The
// j-invariant should entered as a single integer or num/den with no
// spaces.

// Output: one curve per line, conductor N then a-invariants [a1,a2,a3,a4,a6]

// Example:
//
// $ ./CurvesNj -4096/11 1 100
// 11 [0,-1,1,0,0]
// 99 [0,0,1,-3,-5]
// $ ./CurvesNj 10976 1 200
// 128 [0,-1,0,-9,-7]
// 128 [0,1,0,-2,-2]
// 128 [0,1,0,-9,7]
// 128 [0,-1,0,-2,2]

#include "egros.h"

#define VERBOSE 0

int main (int argc, char *argv[])
{
  if ( (argc < 3) || (argc > 4) )
    {
      cerr << "Usage: CurvesNj j N or CurvesNj j N1 N2" << endl;
      return 0;
    }
  // parse the rational j
  bigrational j;
  string js(argv[1]);
  istringstream jss(js);
  jss >> j;
#if VERBOSE
  cout << "j-invariant "<<j<<endl;
#endif
  long n1, n2;
  char* t;
  n1 = strtol(argv[2], &t, 10); // 10 is the base
  if (*t)
    {
      cerr << "Usage: CurvesNj0 N or CurvesNj0 N1 N2" << endl;
      return 0;
    }
  if (argc==4)
    {
      n2 = strtol(argv[3], &t, 10); // 10 is the base
      if (*t)
        {
          cerr << "Usage: CurvesNj0 N or CurvesNj0 N1 N2" << endl;
          return 0;
        }
#if VERBOSE
      cout << "conductors from "<<n1<<" to "<<n2<<endl;
#endif
    }
  else
    {
      n2 = n1;
#if VERBOSE
      cout << "conductor "<<n1<<endl;
#endif
    }

  initprimes("PRIMES");

  for (long n=n1; n<=n2; n++)
    {
      if (!is_valid_conductor(n))
        {
          continue;
        }
      bigint N(n);
      vector<CurveRed> Elist = egros_from_j(j, pdivs(N));
      for (auto Ei = Elist.begin(); Ei!=Elist.end(); ++Ei)
        if (getconductor(*Ei) == N)
          cout << N << " " << (Curve)(*Ei) <<endl;
    }
}
