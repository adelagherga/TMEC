// Program to compare two lists of TM equations read from files (to
// compare the C++ and Magma implementations)

#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>

#include "TME.h"

//#define VERBOSE

int main ()
{
  initprimes("PRIMES");
  int verbose=1;
  string file1, file2;
  cout << "Enter first filename: "; cin>>file1;
  cout << "Enter second filename: "; cin>>file2;

  vector<TM_eqn> list1 = read_TMeqns(file1);
  cout << "Read "<<list1.size()<<" equations from "<<file1<<endl;
  vector<TM_eqn> list2 = read_TMeqns(file2);
  cout << "Read "<<list2.size()<<" equations from "<<file2<<endl;

  int ok = compare_TM_eqn_lists(list1, list2, verbose);
  if (ok)
    cout<<"Lists agree!"<<endl;
  else
    cout<<"Lists do not agree"<<endl;
}
