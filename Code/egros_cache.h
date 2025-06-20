// egros_cache.h - declarations of cached versions of functions for curves with good reduction outside S
//////////////////////////////////////////////////////////////////////////

#include <eclib/curve.h>
#include <set>

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
