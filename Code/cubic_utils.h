// multiply all integers in a list by a constant:
vector<bigint> multiply_list(const bigint& a, const vector<bigint>& L);

// multiply all integers in L1 by all in L2:
vector<bigint> multiply_lists(const vector<bigint>& L1, const vector<bigint>& L2);

// multiply all integers in L by p^e for e in exponents:
vector<bigint> multiply_list_by_powers(const bigint& p, const vector<int>& exponents, const vector<bigint>& L);

// for q prime > 3, returns a list of representatives of the values of
// F(u,v) mod q modulo cubes
vector<bigint> image_mod_cubes(const cubic& F, const bigint& q);

// Similar to AG's magma function of the same name. Return 1 iff there
// exist (u,v) not (0,0) mod q and exponents e such that
// F(u,v)=a*prod(p^e) mod q.
int modpCheck(const cubic& F, const bigint& a, const vector<bigint>& primes, const bigint& q);

// similar to AG's modpcheckDivRHS. Return 1 iff there exists
// primitive (u,v) such that F(u,v)=0 (mod a).
int modaCheck(const cubic& F, const bigint& a);
