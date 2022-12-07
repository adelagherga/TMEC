//////////////////////////////////////////////////////////////////////////////////////////////
//
// Program to list elliptic curves with 2-torsion and given conductor or conductor support
//
// Conditional on an explicit form of Szpiro's conjecture
//
// Author: John Cremona,  based on joint work with Ariel Pacetti
//
/////////////////////////////////////////////////////////////////////////////////////////////

#include <eclib/marith.h>
#include <eclib/curve.h>
#include <assert.h>
#include "curve_utils.h"

// Usage:
// ./twotorsioncurves <N> 0 (or just <N>) # for exact conductor N, default output directory
// ./twotorsioncurves <N1> 1  # for conductors N' with supp(N')=supp(N), default output directory
// ./twotorsioncurves <N1> 2  # for conductors N' from with supp(N')\subseteq supp(N), default output directory

#define VERBOSE 0

// The default directory to hold output files:
const string default_output_directory = "../Data/Curves";
// For each N the output file is output_directory/curves_{N}.txt

string output_filename(long n, const string& dir)
{
  ostringstream s;
  s << dir << "/curves_" << n << ".txt";
  return s.str();
}

int main (int argc, char *argv[])
{
  if ( (argc < 2) || (argc > 4) )
    {
      cerr << "Usage: twotorsioncurves N or twotorsioncurves N 1 or twotorsioncurves N 2" <<endl;
      return 0;
    }
  char* t;
  long n = strtol(argv[1], &t, 10); // 10 is the base
  int supp=0; // default: exact conductor
  if (argc==3)
    supp = strtol(argv[2], &t, 10); // 10 is the base

  initprimes("PRIMES");
  bigint N(n);
  vector<bigint> PP = pdivs(N);

#if VERBOSE > 1
  cout << "N = " << N <<endl;
  cout << "primes " << PP <<endl;
  cout << "support flag = " << supp << endl;
#endif
#if VERBOSE > 0
  if (supp==0)
    cout << "Curves with 2-torsion and conductor " << n << endl;
  if (supp==1)
    cout << "Curves with 2-torsion and bad reduction at " << PP << endl;
  if (supp==2)
    cout << "Curves with 2-torsion and good reduction outside " << PP << endl;
#endif


  vector<CurveRed> curves = CurvesWith2Torsion(N, supp);

  for (auto Ci = curves.begin(); Ci!=curves.end(); ++Ci)
    {
      bigint Ni = getconductor(*Ci);
      cout << Ni << " ";
      // cout << ": " << pdivs(Ni)<< " : ";
      cout << (Curve)(*Ci) <<endl;
    }
#if VERBOSE > 0
  cout << curves.size() << " curve(s)" <<endl;
#endif
}

