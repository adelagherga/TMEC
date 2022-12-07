// class to iterate through divisors of a factored positive integer

class divisor_iterator {
public:
  int ok;            // flags that iteration not finished
  int np;            // number of primes
  vector<bigint> PP; // list of np primes
  vector<long> EE;   // list of np maximum exponents
  vector<long> ee;   // list of np current exponents
  vector<bigint> NN; // list of np+1 partial products
  // constructors
  divisor_iterator(const vector<bigint>& P, const vector<long>& E)
    :ok(1), PP(P), EE(E)
  {
    np = PP.size();
    ee.resize(np, 0);
    NN.resize(np+1, BIGINT(1));
  }
  divisor_iterator(const bigint& N)
    :ok(1)
  {
    PP = pdivs(N);
    for (auto pi=PP.begin(); pi!=PP.end(); ++pi)
      EE.push_back(val(*pi,N));
    np = PP.size();
    ee.resize(np, 0);
    NN.resize(np+1, BIGINT(1));
  }
  divisor_iterator()
    :ok(1), np(0)
  {
    NN.resize(1, BIGINT(1));
  }
  // increment if possible
  void increment();
  // check for end
  int is_ok() {return ok;}
  // deliver current value
  bigint value() {return NN[0];}
};

// test function for divisor_itersor
int test_divisor_iterator(const bigint& N);

// Compute N from its factorization (lists of primes and exponents) --
// (name taken from gp)
bigint factorback(const vector<bigint>&PP, const vector<int>& EE);

// Squarefree part of N
bigint sqf(const bigint& N);

// Maximum conductor for a given list of primes
bigint MaxN(const vector<bigint>&S);
