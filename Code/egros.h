// egros.h - declarations of functions for curves with good reduction outside S

#include <eclib/curve.h>
#include <set>

int is_nth_power(const bigint& x, int n);
bigint prime_to_S_part(const bigint& x,  const vector<bigint>& S);
int is_S_unit(const bigint& x,  const vector<bigint>& S);
int has_good_reduction_outside_S(const CurveRed& C, const vector<bigint>& S);

int is_S_integral(const bigrational& j, const vector<bigint>& S);
int is_j_possible(const bigrational& j, const vector<bigint>& S);
Curve Curve_from_j(const bigrational& j); // one curve with this j-invariant
vector<bigint> twist_factors(const vector<bigint>& S, int n); // only intended for n=2,4,6
vector<CurveRed> egros_from_j_1728(const vector<bigint>& S);
vector<CurveRed> egros_from_j_0(const vector<bigint>& S);
vector<CurveRed> egros_from_j(const bigrational& j, const vector<bigint>& S);

// cached version of egros_from_j_0(S)
vector<CurveRed> get_egros_from_j_0(const vector<bigint>& S);
vector<CurveRed> get_egros_from_j_0(const bigint& N);

int test_conductor_j_0(const bigint& N, const vector<bigint>& support);
int test_conductor_j_1728(const bigint& N, const vector<bigint>& support);

// cached version of egros_from_j_1728(S)
vector<CurveRed> get_egros_from_j_1728(const vector<bigint>& S);
vector<CurveRed> get_egros_from_j_1728(const bigint& N);

// cached version of egros_from_j(S)
vector<CurveRed> get_egros_from_j(const bigrational& j, const vector<bigint>& S);
vector<CurveRed> get_egros_from_j(const bigrational& j, const bigint& N);
