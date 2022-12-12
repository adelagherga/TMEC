// *conjectural* bound on the discriminant (up to quadratic twist) of
// curves with 2-torsion and good reduction outside the primes
// dividing N, based on the maximum possible condictor such a curve
// can have, together with an explaicit form of Szpiro's conjecture.
bigint SzpiroBound(const bigint& N);

// Minimal model of y^2 = x^3+a*x^2+b*x+c
CurveRed Eab(const bigint& a, const bigint& b, const bigint& c=BIGINT(0));

// Quadratic twist of an elliptic curve
CurveRed TwistD(const CurveRed& E, const bigint& D);

// Given a list of elliptic curves E, and one discriminant D, return the
// list of twists of the curves by D
vector<CurveRed> TwistsD(const vector<CurveRed>& EE, const bigint& D);

// Given a list of elliptic curves E, and one prime p, return the
// list of twists of the curves by:
// +p if p=1 (mod 4)
// -p if p=3 (mod 4)
// -4, 8 and -8 if p=2

vector<CurveRed> TwistsP(const vector<CurveRed>& EE, const bigint& p);

// Given a vector of elliptic curves, and a vector of primes, return a
// list of all twists of the curves by discriminants supported on
// those primes (including the original curves)

vector<CurveRed> TwistsPP(const   vector<CurveRed>& EE, const vector<bigint>& PP);

// All curves with 2-torsion and conductor N' with
// if support==0, N' = N
// if support==1, supp(N') = supp(N)
// if support==2, supp(N') contained in N

vector<CurveRed> CurvesWith2Torsion(const bigint& N, int support=0);

