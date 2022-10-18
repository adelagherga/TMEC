// egros.cc:   Implementation of functions for elliptic curves with good reduction outside S
//////////////////////////////////////////////////////////////////////////
//
// Copyright 2022 John Cremona
// 
// This file is part of the eclib package.
// 
// eclib is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// eclib is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with eclib; if not, write to the Free Software Foundation,
// Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA
// 
//////////////////////////////////////////////////////////////////////////

#include "egros.h"

const bigint zero(0), one(1), two(2), three(3);

int is_nth_power(const bigint& x, int n)
{
  if (is_zero(x))
    return 1;
  if ((x<0) && (n%2==0))
    return 0;
  vector<bigint> plist = pdivs(x);
  for (auto pi=plist.begin(); pi!=plist.end(); ++pi)
    if (val(*pi, x)%n!=0)
      return 0;
  return 1;
}

bigint prime_to_S_part(const bigint& x,  const vector<bigint>& S)
{
  if (is_zero(x))
    return x;
  bigint y(abs(x));
  for (auto pi=S.begin(); pi!=S.end() && y>1; ++pi)
    divide_out(y, *pi);
  return y;
}

int is_S_unit(const bigint& x,  const vector<bigint>& S)
{
  return prime_to_S_part(x, S)==1;
}

int is_S_integral(const bigrational& j, const vector<bigint>& S)
{
  return is_S_unit(den(j), S);
}

int has_good_reduction_outside_S(const CurveRed& C, const vector<bigint>& S)
{
  return is_S_unit(getconductor(C), S);
}

Curve Curve_from_j(const bigrational& j) // one curve with this j-invariant
{
  bigint a1, a2, a3, a4, a6;
  bigint n = num(j);
  bigint m = n-1728*den(j);
  if (is_zero(n)) // j=0
    {
      a3=1;
      a4=0;
      a6=0; // 27a3
    }
  else
    if (is_zero(m)) // j=1728
      {
        a3=0;
        a4=-1;
        a6=0; // 32a2
      }
    else // j!=0,1728
    {
      a3=0;
      a4 = -3*n*m;
      a6 = -2*n*m*m;
    }
  return Curve(zero, zero, a3, a4, a6);
}


int is_j_possible(const bigrational& j, const vector<bigint>& S)
{
  bigint n = num(j);
  bigint m = n-1728*den(j);
  if (is_zero(m)) // j=1728
    return 1;
  if (is_zero(n)) // j=0
    return std::find(S.begin(), S.end(), BIGINT(3)) != S.end();
  if (!is_S_integral(j, S))
    return 0;
  return
    is_nth_power(prime_to_S_part(n, S), 3)
    &&
    is_nth_power(prime_to_S_part(m, S), 2);
}

vector<bigint> twist_factors(const vector<bigint>& S, int n)
// only intended for n=2,4,6
{
  vector<bigint> wlist = {one,-one};
  for (auto pi=S.begin(); pi!=S.end(); ++pi)
    {
      bigint p = *pi;
      vector<bigint> ppowers = {one};
      for (int i=1; i<n; i++)
        ppowers.push_back(ppowers[i-1]*p);
      vector<bigint> pwlist;
      pwlist.reserve(wlist.size()*n);
      for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
        for (auto pp = ppowers.begin(); pp!=ppowers.end(); ++pp)
          pwlist.push_back((*pp)*(*wi));
      wlist = pwlist;
    }
  return wlist;
}

vector<CurveRed> egros_from_j_1728(const vector<bigint>& S)
{
  vector<CurveRed> Elist;
  int no2 = std::find(S.begin(), S.end(), two) == S.end();
  vector<bigint> wlist = twist_factors(S, 4);
  for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
    {
      bigint w = *wi;
      if (no2) w *= 4;
      Curve E(zero,zero,zero,w,zero);
      Curvedata Emin(E, 1);
      CurveRed Ered(Emin);
      if (has_good_reduction_outside_S(Ered, S))
        Elist.push_back(Ered);
    }
  return Elist;
}

vector<CurveRed> egros_from_j_0(const vector<bigint>& S)
{
  vector<CurveRed> Elist;
  int no3 = std::find(S.begin(), S.end(), three) == S.end();
  if (no3)
    return Elist;
  int no2 = std::find(S.begin(), S.end(), two) == S.end();
  bigint zero(0);
  vector<bigint> wlist = twist_factors(S, 6);
  for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
    {
      bigint w = *wi;
      if (no2) w *= 16;
      Curve E(zero,zero,zero,zero,w);
      Curvedata Emin(E, 1);
      CurveRed Ered(Emin);
      if (has_good_reduction_outside_S(Ered, S))
        Elist.push_back(Ered);
    }
  return Elist;
}

vector<CurveRed> egros_from_j(const bigrational& j, const vector<bigint>& S)
{
  vector<CurveRed> Elist;
  if (!is_j_possible(j, S))
    {
      // cout << "j="<<j<<" impossible for S="<<S<<endl;
      return Elist;
    }
  bigint n = num(j);
  if (is_zero(n))
    return egros_from_j_0(S);
  bigint m = n-1728*den(j);
  if (is_zero(m))
    return egros_from_j_1728(S);

  bigint a4 = -3*n*m;
  bigint a6 = -2*n*m*m;
  // cout<<"Base curve [0,0,0,"<<a4<<","<<a6<<"]\n";
  vector<bigint> Sx = S;
  vector<bigint> Sy = pdivs(n*m*(n-m));
  vector<bigint> extra_primes;
  for (auto pi=Sy.begin(); pi!=Sy.end(); ++pi)
    {
      if (std::find(Sx.begin(), Sx.end(), *pi) == Sx.end())
        {
          Sx.push_back(*pi);
          extra_primes.push_back(*pi);
        }
    }
  vector<bigint> wlist = twist_factors(Sx, 2);
  // cout << wlist.size() << " twist factors"<<endl;


  // We'll test twists of [0,0,0,a4,a6], whose discriminant is
  // 1728n^2m^3(n-m). For primes p>3 not in S we already have
  // ord_p(n)=0(3), ord_p(m)=0(2) and ord_p(n-m)=ord_p(denom(j))=0, so
  // ord_p(disc)=0(6).  For there to be any good twists we want
  // ord_p(disc)=0(12).  The twist by w is [0,0,0,w^2*a4,w^3*a6] which
  // has disc w^6 times the that of the base curve.  So for the
  // 'extra' primes p (not in S) with ord_p(n)=3(6) we must have ord_p(w) odd,
  // while for p with ord_p(m)=2(4) we must have ord_p(w) odd.

  vector<bigint> a4a6primes;
  for (auto pi=extra_primes.begin(); pi!=extra_primes.end(); ++pi)
    {
      bigint p = *pi;
      if ((p==two) || (p==three))
        continue;
      if ((val(p,n)%6==3) || (val(p,m)%4==2))
        a4a6primes.push_back(p);
    }
  // cout << "extra_primes = "<<a4a6primes<<endl;
  // cout << "a4a6primes =   "<<a4a6primes<<endl;
  bigint zero(0);
  int no2 = std::find(S.begin(), S.end(), BIGINT(2)) == S.end();

  for (auto wi=wlist.begin(); wi!=wlist.end(); ++wi)
    {
      bigint w = *wi;
      int ok=1;
      for (auto pi=a4a6primes.begin(); pi!=a4a6primes.end() && ok; ++pi)
        {
          ok = (val(*pi,w)%2==1);
        }
      if (!ok)
        {
          // cout << "Skipping w = " << w <<endl;
          continue;
        }
      // cout << "Using w = " << w <<endl;
      bigint w2 = w*w;
      bigint w3 = w*w2;
      if (no2)
        w *= 16;
      Curve E(zero,zero,zero,w2*a4,w3*a6);
      Curvedata Emin(E, 1);
      CurveRed Ered(Emin);
      if (has_good_reduction_outside_S(Ered, S))
        Elist.push_back(Ered);
    }
  return Elist;
}

// Map from S to a list of curves with EGR outside S and j=0
map<vector<bigint>, vector<CurveRed>> Elists_0_by_S;

// Map from N to a list of curves with conductor N and j=0
map<bigint, vector<CurveRed>> Elists_0_by_N;

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

vector<CurveRed> get_egros_from_j_0(const bigint& N)
{
  auto Elist = Elists_0_by_N.find(N);
  if (Elist == Elists_0_by_N.end())
    {
      vector<CurveRed> Es = get_egros_from_j_0(pdivs(N));
      for (auto Ei = Es.begin(); Ei!=Es.end(); ++Ei)
        {
          bigint N1 = getconductor(*Ei);
          if (Elists_0_by_N.find(N1) == Elists_0_by_N.end()) // new conductor
            Elists_0_by_N[N1] = {*Ei};
          else // existing conductor
            if (std::find(Elists_0_by_N[N1].begin(), Elists_0_by_N[N1].end(), *Ei) == Elists_0_by_N[N1].end()) // new curve
              Elists_0_by_N[N1].push_back(*Ei);
        }
      return Elists_0_by_N[N];
    }
  else
    {
      return Elist->second;
    }
}

// Test whether N is a possible conductor for j=0: 3|N, no p||N and
// usual bounds on ord_p(N)
int test_conductor_j_0(const bigint& N, const vector<bigint>& support)
{
  if (!div(three,N))
    return 0;
  int ok=1;
  for (auto pi=support.begin(); pi!=support.end() && ok; ++pi)
    {
      bigint p = *pi;
      int np = val(p,N);
      int okp = (np>=2) && (np<= (p==two? 8 : (p==three? 5 : 2)));
      ok = ok && okp;
    }
  return ok;
}

// Test whether N is a possible conductor for j=1728: 2|N, no p||N and
// usual bounds on ord_p(N)
int test_conductor_j_1728(const bigint& N, const vector<bigint>& support)
{
  if (!div(two,N))
    return 0;
  int ok=1;
  for (auto pi=support.begin(); pi!=support.end() && ok; ++pi)
    {
      bigint p = *pi;
      int np = val(p,N);
      int okp = (np>=2) && (np<= (p==two? 8 : (p==three? 5 : 2)));
      ok = ok && okp;
    }
  return ok;
}

// Map from S to a list of curves with EGR outside S and j=1728
map<vector<bigint>, vector<CurveRed>> Elists_1728_by_S;

// Map from N to a list of curves with conductor N and j=1728
map<bigint, vector<CurveRed>> Elists_1728_by_N;

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


// Map from {j,S} to a list of curves with EGR outside S and j=j
map<pair<bigrational,vector<bigint>>, vector<CurveRed>> Elists_j_by_S;

// Map from {j,N} to a list of curves with conductor N and j=j
map<pair<bigrational,bigint>, vector<CurveRed>> Elists_j_by_N;

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
