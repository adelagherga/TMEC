// class for holding a conductor and associated factorization data:

#include "cubic_utils.h"

class Ndata {
public:
  bigint N;
  bigint N0; // prime-to-6 part of N
  int alpha; // valuation of 2
  int beta;  // valuation of 3
  vector<bigint> support; // all prime factors
  vector<bigint> Mprimes; // prime factors > 3 with valuation 2
  vector<bigint> Aprimes; // prime factors > 3 with valuation 1

  // constructors
  Ndata(int conductor)
    :N(conductor)
  {
    init();
  }
  Ndata(const bigint& conductor)
    :N(conductor)
  {
    init();
  }
  void init();
  Ndata() {;}

  // equality test:
  int operator== (const Ndata& N2) const {return (N==N2.N);}

};

// class for holding a discriminant and associated factorization data

class Ddata {
public:
  Ndata NN; // pointer to conductor data
  bigint D;
  bigint D0; // prime-to-6 part of |D|
  int alpha; // valuation of 2
  int beta;  // valuation of 3
  int s;     // sign
  // constructors
  Ddata(const Ndata& Ndat, const bigint& D23, int al, int be, int sg);
  Ddata(const Ndata& Ndat, const bigint& d);
  Ddata() {;}

  // equality test:
  int operator== (const Ddata& D2) const {return (NN==D2.NN && D==D2.D);}
};

// class for holding RHS data for a TM-equation

class TM_RHS {
public:
  bigint a;
  vector<bigint> plist;
  // constructor
  TM_RHS(const bigint& aa, const vector<bigint>& pl)
    :a(aa), plist(pl)
  {;}
  TM_RHS()
    :a(1), plist({})
  {;}
  // for output:
  operator string() const;

  // equality test:
  int operator== (const TM_RHS& R2) const {return (a==R2.a && plist==R2.plist);}
};

class TM_eqn {
public:
  Ddata DD; // discriminant data (which also contains conductor data)
  cubic F;
  TM_RHS RHS;

  TM_eqn(const Ddata& dd, const cubic& f, const TM_RHS& rhs)
    :DD(dd), F(f), RHS(rhs)
  {;}

  // Construct from a string "N,D,[a,b,c,d],A,plist"
  TM_eqn(const string& s);

  // local test: return 0 if impossible, else 1 (and the RHS may have changed)
  int local_test();

  // return an output string
  // (ND=0) "[a,b,c,d],A,plist"
  // (ND=1) "N,[a,b,c,d],A,plist"
  // (ND=2) "D,[a,b,c,d],A,plist"
  // (ND=3) "N,D,[a,b,c,d],A,plist"
  string as_string(int ND=3) const;

  // for long output (converts to string containing "N,D,[a,b,c,d],A,plist")
  operator string() const
  {
    return as_string(1);
  }

  // for input (parsing of a string of the same form as output)

  friend istream& operator>>(istream& s, TM_eqn& tme);

  // equality test:
  int operator== (const TM_eqn& T2) const {return (DD==T2.DD && RHS==T2.RHS && F==T2.F);}

  // equivalence test:
  int is_gl2_equivalent(const TM_eqn& T2) const {return (DD==T2.DD && RHS==T2.RHS && F.gl2_equivalent(T2.F));}
};

// Return a list of discriminants for one conductor
vector<Ddata> get_discriminants(const Ndata& NN);

// Return a list of RHSs (a, primes) for one discriminant
vector<TM_RHS> get_RHS(const Ddata& D);

// Return a list of irreducible cubic forms up to GL(2,Z)-equivalence) for one discriminant
vector<cubic> get_cubics(const Ddata& DD);

// for p||N and p not dividing D=disc(F) we require that F(u,v)=0 (mod
// p) has a nontrivial solution:
int local_test(const cubic& F, const Ddata& DD, const bigint& p);

// Return all TM equations for one discriminant
vector<TM_eqn> get_TMeqnsD(const Ddata& DD);

// Return all TM equations for one conductor
vector<TM_eqn> get_TMeqnsN(const Ndata& NN);

// Read TM equations from a file
vector<TM_eqn> read_TMeqns(const string& filename);

// Write TM equations to a file or stdout
void write_TMeqns(vector<TM_eqn> TMEs, const string& filename="stdout");

int compare_TM_eqn_lists(const vector<TM_eqn>& list1, const vector<TM_eqn>& list2, int verbose=1);
