// class to iterate through divisors of a factored positive integer

class divisor_iterator {
protected:
  int ok;            // flags that iteration not finished
  int np;            // number of primes
  int nd;            // number of divisors
  vector<bigint> PP; // list of np primes
  vector<long> EE;   // list of np maximum exponents
  vector<long> ee;   // list of np current exponents
  vector<bigint> NN; // list of np+1 partial products

public:
  // constructors
  divisor_iterator(const vector<bigint>& P, const vector<long>& E);
  divisor_iterator(const bigint& N);
  divisor_iterator();

  // increment if possible
  void increment();
  // check for end
  int is_ok() {return ok;}
  // reset
  void rewind()
    {
      ee.resize(np, 0);
      NN.resize(np+1, bigint(1));
      ok=1;
    }
  // deliver current value
  bigint value() {return NN[0];}
  // total number of divisors
  long ndivs() {return nd;}
  // report on current status
  void report();
};

// test function for divisor_itersor
int test_divisor_iterator(const bigint& N);

// Compute N from its factorization (lists of primes and exponents) --
// (name taken from gp)
bigint factorback(const vector<bigint>&PP, const vector<int>& EE);

// Radical of N
bigint radical(const bigint& N);

// Maximum conductor for a given list of primes
bigint MaxN(const vector<bigint>&S);
