// egros_cache.cc:   implementations of cached versions of functions for curves with good reduction outside S
//////////////////////////////////////////////////////////////////////////

#include <eclib/egros.h>
#include "egros_cache.h"

// These four maps form a cache so that we don't compute more than
// once the list curves curves with given S or given conductor and j=0
// or j=1728

// Map from S to a list of curves with EGR outside S and j=0
map<vector<bigint>, vector<CurveRed>> Elists_0_by_S;

// Map from N to a list of curves with conductor N and j=0
map<bigint, vector<CurveRed>> Elists_0_by_N;

// Map from S to a list of curves with EGR outside S and j=1728
map<vector<bigint>, vector<CurveRed>> Elists_1728_by_S;

// Map from N to a list of curves with conductor N and j=1728
map<bigint, vector<CurveRed>> Elists_1728_by_N;

// These two maps are a cache of lists of curves with good reduction
// outside S, or conductor N, and one fixed j-invariant:

// Map from {j,S} to a list of curves with EGR outside S and j=j
map<pair<bigrational,vector<bigint>>, vector<CurveRed>> Elists_j_by_S;

// Map from {j,N} to a list of curves with conductor N and j=j (unused)
// map<pair<bigrational,bigint>, vector<CurveRed>> Elists_j_by_N;

// Return a list of curves with good reduction outside S and j=0,
// either retrieving from the cache or computing and caching first:
vector<CurveRed> get_egros_from_j_0(const vector<bigint>& S)
{
  auto Elist = Elists_0_by_S.find(S);
  if (Elist == Elists_0_by_S.end())
    {
      vector<CurveRed> Es = egros_from_j_0(S);
      Elists_0_by_S[S] = Es;
      return Es;
    }
  else
    {
      return Elist->second;
    }
}

// Return a list of curves with conductor N and j=0,
// either retrieving from the cache or computing and caching first:
vector<CurveRed> get_egros_from_j_0(const bigint& N)
{
  auto Elist = Elists_0_by_N.find(N);
  if (Elist == Elists_0_by_N.end())
    {
      vector<CurveRed> Es = get_egros_from_j_0(pdivs(N));
      for (auto E: Es)
        {
          bigint N1 = E.conductor();
          if (Elists_0_by_N.find(N1) == Elists_0_by_N.end()) // new conductor
            Elists_0_by_N[N1] = {E};
          else // existing conductor
            if (std::find(Elists_0_by_N[N1].begin(), Elists_0_by_N[N1].end(), E) == Elists_0_by_N[N1].end()) // new curve
              Elists_0_by_N[N1].push_back(E);
        }
      return Elists_0_by_N[N];
    }
  else
    {
      return Elist->second;
    }
}

// Return a list of curves with good reduction outside S and j=1728,
// either retrieving from the cache or computing and caching first:
vector<CurveRed> get_egros_from_j_1728(const vector<bigint>& S)
{
  auto Elist = Elists_1728_by_S.find(S);
  if (Elist == Elists_1728_by_S.end())
    {
      vector<CurveRed> Es = egros_from_j_1728(S);
      Elists_1728_by_S[S] = Es;
      return Es;
    }
  else
    {
      return Elist->second;
    }
}

// Return a list of curves with conductor N and j=1728,
// either retrieving from the cache or computing and caching first:
vector<CurveRed> get_egros_from_j_1728(const bigint& N)
{
  auto Elist = Elists_1728_by_N.find(N);
  if (Elist == Elists_1728_by_N.end())
    {
      vector<CurveRed> Es = get_egros_from_j_1728(pdivs(N));
      for (auto Ei = Es.begin(); Ei!=Es.end(); ++Ei)
        {
          bigint N1 = getconductor(*Ei);
          if (Elists_1728_by_N.find(N1) == Elists_1728_by_N.end()) // new conductor
            Elists_1728_by_N[N1] = {*Ei};
          else // existing conductor
            if (std::find(Elists_1728_by_N[N1].begin(), Elists_1728_by_N[N1].end(), *Ei) == Elists_1728_by_N[N1].end()) // new curve
              Elists_1728_by_N[N1].push_back(*Ei);
        }
      return Elists_1728_by_N[N];
    }
  else
    {
      return Elist->second;
    }
}

// cached version of egros_from_j(S)
vector<CurveRed> get_egros_from_j(const bigrational& j, const vector<bigint>& S)
{
  auto Elist = Elists_j_by_S.find({j,S});
  if (Elist == Elists_j_by_S.end())
    {
      vector<CurveRed> Es = egros_from_j(j, S);
      Elists_j_by_S[{j,S}] = Es;
      return Es;
    }
  else
    {
      return Elist->second;
    }
}
