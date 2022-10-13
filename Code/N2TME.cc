// Program to list irreducible cubic forms associated to elliptic
// curves of given conductor with no 2-torsion

// Usage:
// ./N2TME <N>  # for one conductor N, default output directory
// ./N2TME <N1> <N2>  # for all conductors from N1 to N2 inclusive, default output directory
// ./N2TME <N> <dir> # for one conductor N, output directory dir
// ./N2TME <N1> <N2> <dir> # for all conductors from N1 to N2 inclusive, output directory dir

// NB the number of parameters on the command line is 1, 2 or 3, and
// if 2 we need to check whether the 2nd is an integer or not

#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>
#include <assert.h>
#include "TME.h"

#define VERBOSE 0

// The default directory to hold output files:
const string default_output_directory = "../Data/Forms";
// For each N the output file is output_directory/{N}Forms.csv

string output_filename(long n, const string& dir)
{
  ostringstream s;
  s << dir << "/" << n << "Forms.csv";
  return s.str();
}

// return 0 for an integer (put into n), 1 for a string (put into s)
int parse_int_or_string(const char* arg, long& n, string& s)
{
  // cout<<"parsing "<<arg<<" ..."<<flush;
  char* t;
  n = strtol(arg, &t, 10); // 10 is the base
  if (!*t)
    {
      // cout<<"long int n = "<<n<<", returning 0" <<  endl;
      return 0;
    }
  s = string(arg);
  // cout<<"string s = "<<s<<", returning 1" <<  endl;
  return 1;
}

int main (int argc, char *argv[])
{
  if ( (argc < 2) || (argc > 4) )
    {
      cerr << "Usage: N2TME N or N2TME N1 N2 or N2TME N dir or N2TME N1 N2 dir" <<endl;
      return 0;
    }
  long n1, n2, n;
  string dir, dummy;
  int t = parse_int_or_string(argv[1], n1, dir);
  if (t) // string dir is set, now parse one or two ints
    {
      t = parse_int_or_string(argv[2], n1, dummy);
      assert(t==0);
      if (argc == 4) // parse another int
        {
          t = parse_int_or_string(argv[3], n2, dummy);
          assert(t==0);
        }
      else
        n2 = n1;
    }
  else // int n1 is set
    {
      t = parse_int_or_string(argv[2], n2, dir);
      if (t) // string dir is set
        {
          if (argc == 4) // parse another int
            {
              t = parse_int_or_string(argv[3], n2, dummy);
              assert(t==0);
            }
          else
            n2 = n1;
        }
      else // we now have n1, n2
        {
          if (argc == 4) // parse a string
            {
              t = parse_int_or_string(argv[3], n, dir);
            }
          else
            dir = default_output_directory;
        }
    }

#if VERBOSE
  if (dir==default_output_directory)
    {
      if (argc==2)
        cerr << "TM equations for conductor " << n1 << endl;
      else
        cerr << "TM equations for conductors from " << n1 << " to " << n2 << endl;
    }
  else
    {
      if (argc==2)
        cerr << "TM equations for conductor " << n1 << " (output to "<<dir<<")"<<endl;
      else
        cerr << "TM equations for conductors from " << n1 << " to " << n2 << " (output to "<<dir<<")"<< endl;
    }

#endif
  initprimes("PRIMES");

  for (long n=n1; n<=n2; n++)
    {
#if VERBOSE
      cout << "N = " << n << ": " <<flush;
#endif
      if (!is_valid_conductor(n))
        {
#if VERBOSE
          cout << " -- not a valid conductor" << endl;
#endif
          continue;
        }
      string ofname = output_filename(n, dir);
      vector<TM_eqn> TM_eqns = get_TMeqnsN(n);

#if VERBOSE
      int neqns = TM_eqns.size();
      if (neqns)
        cout << neqns;
      else
        cout <<"No";
      cout << " TM equations for conductor " << n << endl;
      cout << "Writing output to " << ofname << endl;
#endif
      ofstream fout;
      fout.open(ofname.c_str());
      for (auto T = TM_eqns.begin(); T!=TM_eqns.end(); ++T)
        {
          // output to screen includes N and D:
          cout << T->as_string(3) << endl;
          // output to file includes N but not D:
          fout << T->as_string(1) << endl;
        }
      fout.close();
#if VERBOSE
      cout << neqns << " lines written to " << ofname << endl;
#endif
    }
}
