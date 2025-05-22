#include <eclib/marith.h>
#include <eclib/curve.h>
#include <eclib/isogs.h>
#include <eclib/parifact.h>
#include <assert.h>
#include "arith_utils.h"
#include "curve_utils.h"
#include <set>

bigint SzpiroBound(const bigint& N)
{
  vector<bigint> S = pdivs(N);
  bigint Nmax = MaxN(S);
  double logNmax = log(Nmax);
  double logSFN=log(sqf(N));
  // cout<<"Szpiro bound:"<<endl;
  // cout<<"N = "<<N<<endl;
  // cout<<"Nmax = "<<Nmax<<endl;
  // cout<<"log(Nmax) = "<<logNmax<<endl;
  // cout<<"logSFN = "<<logSFN<<endl;
  bigint M(1);
  for (auto pi=S.begin(); pi!=S.end(); ++pi)
    {
      bigint p=*pi;
      double logp = log(p);
      double factor(p==2? 6.1: 3);
      int e = floor(factor*(logNmax-logSFN+logp)/logp);
      // cout<<"(p,e) = ("<<p<<","<<e<<")\n";
      M *= pow(p,e);
    }
  if (ndiv(2,N))
    M *= 256;
  return M;
}

// Minimal model of y^2 = x^3+a*x^2+b*x+c
CurveRed Eab(const bigint& a, const bigint& b, const bigint& c)
{
  bigint zero(0);
  Curvedata C(zero,a,zero,b,c, 1);
  return CurveRed(C);
}

// Quadratic twist of an elliptic curve

CurveRed TwistD(const CurveRed& E, const bigint& D)
{
  bigint c4, c6;
  E.getci(c4, c6);
  CurveRed ED = Eab(bigint(0),-27*D*D*c4,-54*D*D*D*c6);
  // cout<<"Twisting E="<<(Curve)E<<", N_E="<<getconductor(E)<<" by "<<D<<": "<<(Curve)ED<<" N_ED="<<getconductor(ED)<<endl;
  return ED;
}

// Given a list of elliptic curves E, and one discriminant D, return the
// list of twists of the curves by D

vector<CurveRed> TwistsD(const vector<CurveRed>& EE, const bigint& D)
{
  vector<CurveRed> ans;
  for(auto Ei = EE.begin(); Ei!=EE.end(); ++Ei)
    ans.push_back(TwistD(*Ei,D));
  return ans;
}

// Given a list of elliptic curves E, and one prime p, return the
// list of twists of the curves by:
// +p if p=1 (mod 4)
// -p if p=3 (mod 4)
// -4, 8 and -8 if p=2

vector<CurveRed> TwistsP(const vector<CurveRed>& EE, const bigint& p)
{
  long p4 = posmod(p,4);
  if (p4==1)
    return TwistsD(EE,p);
  if (p4==3)
    return TwistsD(EE,-p);
  vector<CurveRed>
    ans1 = TwistsD(EE,bigint(-4)),
    ans2 = TwistsD(EE,bigint(-8)),
    ans3 = TwistsD(EE,bigint(8));
  ans1.insert(ans1.end(), ans2.begin(), ans2.end());
  ans1.insert(ans1.end(), ans3.begin(), ans3.end());
  return ans1;
}

// Given a vector of elliptic curves, and a vector of primes, return a
// list of all twists of the curves by discriminants supported on
// those primes (including the original curves)

vector<bigint> a_invariants(const CurveRed& C)
{
  vector<bigint> ai(5);
  C.getai(ai[0],ai[1],ai[2],ai[3],ai[4]);
  return ai;
}

vector<long> primes10 = {2,3,5,7,11,13,17,19,23,29};

vector<bigint> sort_key(const CurveRed& C)
{
  vector<bigint> key;
  vector<bigint> ainvs(5);
  bigint NC = getconductor(C);
  key.push_back(NC);
  C.getai(ainvs[0],ainvs[1],ainvs[2],ainvs[3],ainvs[4]);
  for(auto pi=primes10.begin(); pi!=primes10.end(); ++pi)
    {
      long p = *pi;
      if (posmod(NC,p)==0)
        continue;
      long ap = ellap(posmod(ainvs[0],p), posmod(ainvs[1],p), posmod(ainvs[2],p),
                      posmod(ainvs[3],p), posmod(ainvs[4],p), p);
      key.push_back(bigint(ap));
    }
  key.insert(key.end(), ainvs.begin(), ainvs.end());
  return key;
}

int operator==(const CurveRed& C1, const CurveRed& C2)
{
  return a_invariants(C1) == a_invariants(C2);
}

int operator<(const CurveRed& C1, const CurveRed& C2)
{
  return sort_key(C1) < sort_key(C2);
}

vector<CurveRed> TwistsPP(const vector<CurveRed>& EE, const vector<bigint>& PP)
{
  vector<CurveRed> ans = EE;
  for (auto pi=PP.begin(); pi!=PP.end(); ++pi)
    {
      vector<CurveRed> p_twists = TwistsP(ans, *pi);
      ans.insert(ans.end(), p_twists.begin(), p_twists.end());
    }
  // now remove any duplicates
  std::set<CurveRed> s;
  s.insert(ans.begin(), ans.end());
  ans.assign(s.begin(), s.end());
  return ans;
}

// return list of curves linked to C by 2-power isogenies (including C
// itself)
vector<CurveRed> twopowerisog(const CurveRed& C)
{
  vector<CurveRed> curves;
  curves.push_back(C);
  int more=1;
  while (more)
    {
      auto nc = curves.size();
      std::set<CurveRed> newcurves;
      for (auto Ci=curves.begin(); Ci!=curves.end(); ++Ci)
        {
          newcurves.insert(*Ci);
          vector<CurveRed> curves2 = twoisog(*Ci,0);
          newcurves.insert(curves2.begin(), curves2.end());
        }
      curves.assign(newcurves.begin(), newcurves.end());
      auto new_nc = curves.size();
      more = (new_nc > nc);
    }
  return curves;
}

// utility function: d1 will run over all divisors of D1=maxD1(b,PP)
// omitting those whose 2-valuation is 2
bigint maxD1(const bigint& b, const vector<bigint>& PP)
{
  vector<int> EE;
  for(auto pi = PP.begin(); pi!=PP.end(); ++pi)
    {
      bigint p = *pi;
      EE.push_back((val(p,b)==1?(p==2?3:1):0));
    }
  return factorback(PP,EE);
}

// utility function: d2 will run over all divisors of D2=maxD2(b,PP, Dmax)
bigint maxD2(const bigint& b, const vector<bigint>& PP, const bigint& Dmax)
{
  vector<int> EE;
  for(auto pi = PP.begin(); pi!=PP.end(); ++pi)
    {
      bigint p = *pi;
      EE.push_back((div(p,b)?0:val(p,Dmax)));
    }
  return factorback(PP,EE);
}

// All curves with 2-torsion and conductor N' with
// if support==0, N' = N
// if support==1, supp(N') = supp(N)
// if support==2, supp(N') contained in N

#define DEBUG 2

vector<CurveRed> CurvesWith2Torsion(const bigint& N, int support)
{
  bigint M = sqf(N);
  bigint N2 = N;
  int oddN = !div(2,N);
  if (oddN) N2 *=2;
  bigint Dmax = SzpiroBound(N);
  vector<bigint> PP = pdivs(Dmax);
  vector<long> EE, EE2;
  for (auto pi=PP.begin(); pi!=PP.end(); ++pi)
    {
      long e = val(*pi,Dmax);
      EE.push_back(e);
      EE2.push_back(e/2);
    }
#if DEBUG > 0
  cerr<<"Dmax = "<<Dmax<<endl;
  cerr<<"primes: "<<PP<<endl;
  cerr<<"exponents: "<<EE<<endl;
  cerr<<"exponents/2: "<<EE2<<endl;
  int nb=1;
  for(auto e=EE2.begin(); e!=EE2.end(); ++e)
    nb *= (1+*e);
  cerr<<nb<<" b values"<<endl;
#endif
  vector<vector<bigint>> ab_pairs;
  bigint a, asq, b, d1, d2, d1d2;
  // iterates through b s.t. b^2 | Dmax
  for (divisor_iterator b_iter(PP,EE2); b_iter.is_ok(); b_iter.increment())
    {
      b = b_iter.value();
      bigint mD1 = maxD1(b, PP);
      bigint mD2 = maxD2(b, PP, Dmax);
#if DEBUG > 1
      cout<<"b = "<<b<<endl;
      cout<<"maxD1 = "<<mD1<<endl;
      cout<<"maxD2 = "<<mD2<<endl;
#endif
      // iterate through d1 | maxD1
#if DEBUG > 1
      divisor_iterator d1_iter(mD1);
      cout<<d1_iter.ndivs()<< " d1 values"<<endl;
#endif
      for (divisor_iterator d1_iter(mD1); d1_iter.is_ok(); d1_iter.increment())
        {
          d1 = d1_iter.value();
#if DEBUG > 2
          cout<<" d1 = "<<d1<<endl;
#endif
          assert (!(d1==0));
          if (val(2,d1)==2)
            continue;
#if DEBUG > 2
          divisor_iterator d2_iter(mD2);
          cout<<d2_iter.ndivs()<< " d2 values"<<endl;
#endif
          // iterate through d2 | maxD2
          for (divisor_iterator d2_iter(mD2); d2_iter.is_ok(); d2_iter.increment())
            {
              d2 = d2_iter.value();
#if DEBUG > 3
              cout<<"  d2 = "<<d2<<endl;
#endif
              assert (!(d2==0));
              d1d2 = d1*d2_iter.value();
              if (isqrt(d1d2+4*b,a))
                ab_pairs.push_back({a,b});
              if (isqrt(-d1d2+4*b,a))
                ab_pairs.push_back({a,b});
              if (isqrt(d1d2-4*b,a))
                ab_pairs.push_back({a,-b});
            }
#if DEBUG > 2
          cout<<"  last d2 = " << d2<<" #(a,b) pairs now "<<ab_pairs.size()<<endl;
#endif
        }
    }
  std::set<vector<bigint>> sab;
  sab.insert(ab_pairs.begin(), ab_pairs.end());
  ab_pairs.assign(sab.begin(), sab.end());
#ifdef DEBUG2
  cout<<"Found "<<ab_pairs.size()<<" (a,b) pairs"<<endl;
#endif
  // convert (a,b) to curves y^2=x(x^2+ax+b):
  vector<CurveRed> curves;
  for(auto ab=ab_pairs.begin(); ab!=ab_pairs.end(); ++ab)
    curves.push_back(Eab((*ab)[0],(*ab)[1]));

#if DEBUG > 1
  for(auto Ci=curves.begin(); Ci!=curves.end(); ++Ci)
  cout<<getconductor(*Ci)<<" : "<<(Curve)(*Ci)<<endl;
#endif

  // include all twists:
  curves = TwistsPP(curves, PP);
#if DEBUG > 0
  cout<<"With all twists: "<<curves.size()<<" curves"<<endl;
#if DEBUG > 1
  for(auto Ci=curves.begin(); Ci!=curves.end(); ++Ci)
  cout<<getconductor(*Ci)<<" : "<<(Curve)(*Ci)<<endl;
#endif
#endif

  // If N is odd, discard curves with even conductor. Then:
  // if support==2 keep all;
  // if support==1 keep only those with bad reduction at all p|N;
  // if support==0 keep only those with conductor N.
  vector<CurveRed> curves1;
  for (auto Ci=curves.begin(); Ci!=curves.end(); ++Ci)
    {
      bigint NC = getconductor(*Ci);
      int keep;
      switch (support)
        {
        case 0:
          {
            keep = (NC==N);
            break;
          }
        case 1:
          {
            keep = (sqf(NC)==M);
            break;
          }
        default:
          {
            keep = !(oddN && div(2,NC));
          }
        }
      if (keep)
        curves1.push_back(*Ci);
    }
  curves = curves1;
#if DEBUG > 0
  cout<<"With support constraints: "<<curves.size()<<" curves"<<endl;
#endif

  // close under 2-isogenies:
  std::set<CurveRed> curves2;
  for (auto Ci=curves.begin(); Ci!=curves.end(); ++Ci)
    {
      vector<CurveRed> C2C = twopowerisog(*Ci);
      curves2.insert(C2C.begin(), C2C.end());
    }
  curves.assign(curves2.begin(), curves2.end());
#if DEBUG > 0
  cout<<"With 2-power-isogenous curves: "<<curves.size()<<" curves"<<endl;
#endif
  std::sort(curves.begin(), curves.end());
  return curves;
}
