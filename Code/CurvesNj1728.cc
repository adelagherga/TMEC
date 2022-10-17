// Program to list curves of given conductor with j-invariant 1728

// Usage:
// ./CurvesNj0 <N>  # for one conductor N
// ./CurvesNj0 <N1> <N2>  # for all conductors from N1 to N2 inclusive

// NB the number of parameters on the command line is 1 or 2

#include "egros.h"

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

  for (long n=n1; n<=n2; n++)
    {
      if (!is_valid_conductor(n))
        {
          continue;
        }
      bigint N(n);
      //vector<CurveRed> Elist = get_egros_from_j_1728(N);
      vector<CurveRed> Elist = get_egros_from_j_1728(pdivs(N));
      for (auto Ei = Elist.begin(); Ei!=Elist.end(); ++Ei)
        if (getconductor(*Ei) == N)
          cout << N << " " << (Curve)(*Ei) <<endl;
    }
}
