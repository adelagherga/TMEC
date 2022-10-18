// Program to list curves of given conductor and j-invariant

// Usage:
// ./CurvesNj0 <N>  # for one conductor N
// ./CurvesNj0 <N1> <N2>  # for all conductors from N1 to N2 inclusive

// NB the number of parameters on the command line is 1 or 2

#include "egros.h"

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
  cout << "j-invariant "<<j<<endl;
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
      cout << "conductors from "<<n1<<" to "<<n2<<endl;
    }
  else
    {
      n2 = n1;
      cout << "conductor "<<n1<<endl;
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
