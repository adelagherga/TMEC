//////////////////////////////////////////////////////////////////////////////////////////////
//
// Program to list elliptic curves with 2-torsion and given conductor or conductor support
//
// Conditional on an explicit form of Szpiro's conjecture
//
// Author: John Cremona,  based on joint work with Ariel Pacetti
//
/////////////////////////////////////////////////////////////////////////////////////////////

#include <eclib/marith.h>
#include <eclib/curve.h>
#include <eclib/isogs.h>
#include <assert.h>

// Usage:
// ./twotorsioncurves <N> 0 (or just <N>) # for exact conductor N, default output directory
// ./twotorsioncurves <N> 1  # for conductors N' with supp(N')=supp(N), default output directory
// ./twotorsioncurves <N> 2  # for conductors N' from with supp(N')\subseteq supp(N), default output directory

#define VERBOSE 0

// The default directory to hold output files:
const string default_output_directory = "../Data/Curves";
// For each N the output file is output_directory/curves_{N}.txt

// *conjectural* bound on the discriminant (up to quadratic twist) of
// curves with 2-torsion and good reduction outside the primes
// dividing N, based on the maximum possible conductor such a curve
// can have, together with an explicit form of Szpiro's conjecture.
bigint SzpiroBound(const bigint& N)
{
  vector<bigint> S = pdivs(N);  // support of N
  bigint Nmax = MaxN(S);        // largest possible conductor with this support
  double logNmax = log(Nmax);   // log of this
  double logSFN=log(radical(N));    // log of radical of N
  // cout<<"Szpiro bound:"<<endl;
  // cout<<"N = "<<N<<endl;
  // cout<<"Nmax = "<<Nmax<<endl;
  // cout<<"log(Nmax) = "<<logNmax<<endl;
  // cout<<"logSFN = "<<logSFN<<endl;
  bigint M(1);
  for (auto p: S)
    {
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

// Minimal model of y^2 = x^3+a*x^2+b*x
CurveRed Eab(const vector<bigint>& ab) // ab=[a,b]
{
  bigint zero(0);
  Curvedata C(zero, ab[0], zero, ab[1], zero, 1); // 1 means make model minimal
  return CurveRed(C);
}

// return list of curves linked to C by 2-power isogenies (including C itself)
vector<CurveRed> twopowerisog(const CurveRed& C)
{
  vector<CurveRed> curves;
  curves.push_back(C);
  while (1)
    {
      // newcurves will contain all E in curves and those 2-isogenous to them.
      // We use  new list as we'll be looping through the old one.
      vector<CurveRed> newcurves;
      for (auto Ci: curves)
        {
          vector<CurveRed> curves2 = twoisog(Ci,0);
          curves2.push_back(Ci);
          for (auto Cj: curves2)
            if (std::find(newcurves.begin(), newcurves.end(), Cj) == newcurves.end())
              newcurves.push_back(Cj);
        }
      // Now copy newcurves back into curves if it has more
      if (newcurves.size() > curves.size())
        curves.assign(newcurves.begin(), newcurves.end());
      else
        break;
    }
  return curves;
}

// utility function: d1 will run over all divisors of D1=maxD1(b,PP)
// omitting those whose 2-valuation is 2
bigint maxD1(const bigint& b, const vector<bigint>& PP)
{
  vector<int> EE;
  for(auto p: PP)
    {
      EE.push_back((val(p,b)==1?(p==2?3:1):0));
    }
  return factorback(PP,EE);
}

// utility function: d2 will run over all divisors of D2=maxD2(b,PP, Dmax)
bigint maxD2(const bigint& b, const vector<bigint>& PP, const bigint& Dmax)
{
  vector<int> EE;
  for(auto p: PP)
    {
      EE.push_back((div(p,b)?0:val(p,Dmax)));
    }
  return factorback(PP,EE);
}

// All curves with 2-torsion and conductor N' with
// if support==0, N' = N
// if support==1, supp(N') = supp(N)
// if support==2, supp(N') contained in N

#define DEBUG 0

vector<CurveRed> CurvesWith2Torsion(const bigint& N, int support)
{
  bigint M = radical(N);
  bigint N2 = N;
  int oddN = !div(2,N);
  if (oddN) N2 *=2;
  bigint Dmax = SzpiroBound(N);
  vector<bigint> PP = pdivs(Dmax);
  vector<long> EE, EE2;
  for (auto p: PP)
    {
      long e = val(p,Dmax);
      EE.push_back(e);
      EE2.push_back(e/2);
    }
#if DEBUG > 1
  cerr<<"N: "<<N<<endl;
  //cerr<<"Dmax = "<<Dmax<<endl;
  cerr<<"Dmax primes: "<<PP<<endl;
#if DEBUG > 2
  cerr<<"Dmax exponents: "<<EE<<endl;
  cerr<<"Dmax exponents/2: "<<EE2<<endl;
#endif
  int nb=1;
  for(auto e: EE2)
    nb *= (1+e);
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
#if DEBUG > 2
      cout<<"b = "<<b<<endl;
      cout<<"maxD1 = "<<mD1<<endl;
      cout<<"maxD2 = "<<mD2<<endl;
#endif
      // iterate through d1 | maxD1
#if DEBUG > 2
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
#if DEBUG > 1
  cout<<"Found "<<ab_pairs.size()<<" (a,b) pairs"<<endl;
#endif
  // convert (a,b) to curves y^2=x(x^2+ax+b):
  vector<CurveRed> curves;
  for(auto ab: ab_pairs)
    curves.push_back(Eab(ab));

#if DEBUG > 1
  for(auto C: curves)
    cout << C.conductor()<<" : "<<(Curve)C<<endl;
#endif

  // include all twists:
  curves = AllTwists(curves, PP);
#if DEBUG > 1
  cout<<"With all twists: "<<curves.size()<<" curves"<<endl;
#if DEBUG > 1
  for(auto C: curves)
    cout << C.conductor()<<" : "<<(Curve)C<<endl;
#endif
#endif

  // If N is odd, discard curves with even conductor. Then:
  // if support==2 keep all;
  // if support==1 keep only those with bad reduction at all p|N;
  // if support==0 keep only those with conductor N.
  vector<CurveRed> curves1;
  for (auto C: curves)
    {
      bigint NC = C.conductor();
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
            keep = (radical(NC)==M);
            break;
          }
        default:
          {
            keep = !(oddN && div(2,NC));
          }
        }
      if (keep)
        curves1.push_back(C);
    }
  curves.assign(curves1.begin(), curves1.end());
#if DEBUG > 0
  cout<<"With support constraints: "<<curves.size()<<" curves"<<endl;
#endif

  // close under 2-isogenies:
  vector<CurveRed> curves2;
  for (auto C: curves)
    {
      vector<CurveRed> C2C = twopowerisog(C);
      for (auto C2: C2C)
        if (std::find(curves2.begin(), curves2.end(), C2) == curves2.end())
          curves2.push_back(C2);
    }
  curves.assign(curves2.begin(), curves2.end());
#if DEBUG > 0
  cout<<"With 2-power-isogenous curves: "<<curves.size()<<" curves"<<endl;
#endif
  std::sort(curves.begin(), curves.end());
  return curves;
}

string output_filename(long n, const string& dir)
{
  ostringstream s;
  s << dir << "/curves_" << n << ".txt";
  return s.str();
}

int main (int argc, char *argv[])
{
  if ( (argc < 2) || (argc > 4) )
    {
      cerr << "Usage: twotorsioncurves N or twotorsioncurves N 1 or twotorsioncurves N 2" <<endl;
      return 0;
    }
  char* t;
  long n = strtol(argv[1], &t, 10); // 10 is the base
  int supp=0; // default: exact conductor
  if (argc==3)
    supp = strtol(argv[2], &t, 10); // 10 is the base

  initprimes("PRIMES");
  if (supp==0 && !is_valid_conductor(n))
    {
#if VERBOSE > 0
      cout << "No elliptic curves have conductor " << n << endl;
#endif
      exit(0);
    }
  bigint N(n);
  vector<bigint> PP = pdivs(N);

#if VERBOSE > 1
  cerr << "\nN = "<<N<<endl;
  cout << "N = " << N <<endl;
  cout << "primes " << PP <<endl;
  cout << "support flag = " << supp << endl;
#endif

#if(0) // uncomment this to skip N with too many prime factors
  if (PP.size()>5)
    {
      cout << "# prime factors = "<< PP.size()<<" is > 5 ("<<PP<<"), skipping"<<endl;
      exit(0);
    }
#endif
#if(0) // uncomment this to only run N with a specific number of prime factors
  if (PP.size()!=6)
    {
      cout << "# prime factors = "<< PP.size()<<" is not 6 ("<<PP<<"), skipping"<<endl;
      exit(0);
    }
#endif
#if VERBOSE > 0
  if (supp==0)
    cout << "Curves with 2-torsion and conductor " << n << endl;
  if (supp==1)
    cout << "Curves with 2-torsion and bad reduction at " << PP << endl;
  if (supp==2)
    cout << "Curves with 2-torsion and good reduction outside " << PP << endl;
#endif


  vector<CurveRed> curves = CurvesWith2Torsion(N, supp);

#if VERBOSE > 0
  cout << "---------------------------------------------------"<<endl;
  cerr << curves.size() << " curves found" <<endl;
  cout << curves.size() << " curves found" <<endl;
  cout << "---------------------------------------------------"<<endl;
#endif
  for (auto Ci = curves.begin(); Ci!=curves.end(); ++Ci)
    {
      bigint Ni = getconductor(*Ci);
      cout << Ni << " ";
      // cout << ": " << pdivs(Ni)<< " : ";
      cout << (Curve)(*Ci) <<endl;
    }
#if VERBOSE > 0
  cout << curves.size() << " curves found" <<endl;
#endif
}

