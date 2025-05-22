//////////////////////////////////////////////////////////////////////////////////////////////
//
// Program to test divisor_iterator class
//
/////////////////////////////////////////////////////////////////////////////////////////////

#include <eclib/marith.h>
#include "arith_utils.h"

int main (int argc, char *argv[])
{
  if ( (argc < 2) || (argc > 3) )
    {
      cerr << "Usage: divtest N" <<endl;
      return 0;
    }
  char* t;
  long n = strtol(argv[1], &t, 10); // 10 is the base

  initprimes("PRIMES");
  bigint N(n);
  int ok = test_divisor_iterator(N);
  if (ok)
    cout << "N = " << N << ": test passes" << endl;
  else
    cout << "N = " << N << ": test fails" << endl;
}

