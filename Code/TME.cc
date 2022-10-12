#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>
#include <assert.h>

#include "TME.h"

const vector<int> signs = {1,-1};
const vector<vector<int>> alpha0list = {{2}, {2,3}, {2,4}, {2,3,4}, {2,3,4}, {2,3}, {2,3,4}, {3,4}, {3}};
const vector<vector<int>> beta0list = {{0}, {0,1}, {0,1,3}, {3}, {4}, {5}};
const vector<int> powersof2 = {1,2,4,8,16,32,64,128,256};
const vector<int> powersof3 = {1,3,9,27,81,243};
const bigint one(1), two(2), three(3), six(6);

void Ndata::init()
{
  support = pdivs(N);
  N0 = N;
  alpha = divide_out(N0, two);
  beta  = divide_out(N0, three);
  for (auto pi = support.begin(); pi!=support.end(); ++pi)
    {
      bigint p = *pi;
      if (p>3)
        {
          if (val(p,N)==1)
            Mprimes.push_back(p);
          else
            Aprimes.push_back(p);
        }
    }
}

Ddata::Ddata(const Ndata& Ndat, const bigint& D23, int al, int be, int sg)
  :NN(Ndat), D0(D23), alpha(al), beta(be), s(sg)
{
  D = D0 * powersof2[alpha] * powersof3[beta] * s;
}

Ddata::Ddata(const Ndata& Ndat, const bigint& d)
  :NN(Ndat), D(d), D0(d), s(sign(d))
{
  alpha = divide_out(D0, 2);
  beta  = divide_out(D0, 3);
}

vector<Ddata> get_discriminants(const Ndata& NN)
{
  vector<int> alpha0s = alpha0list[NN.alpha];
  vector<int> beta0s = beta0list[NN.beta];
  vector<bigint> N1s = posdivs(NN.N0);

  vector<Ddata> Dlist;
  for (auto N1i  = N1s.begin(); N1i!=N1s.end(); ++N1i)
    for (auto alpha0  = alpha0s.begin(); alpha0!=alpha0s.end(); ++alpha0)
      for (auto beta0  = beta0s.begin(); beta0!=beta0s.end(); ++beta0)
        for (auto s = signs.begin(); s!=signs.end(); s++)
          Dlist.push_back(Ddata(NN, *N1i, *alpha0, *beta0, *s));
  return Dlist;
}

// these two maps are not yet used in the code which follows

map<pair<int,int>, vector<int>> alpha_map = { {{0,2}, {0,3}},
                                              {{2,2}, {1}},
                                              {{2,4}, {0,1}},
                                              {{3,2}, {1,2}},
                                              {{3,3}, {2}},
                                              {{3,4}, {0,1}},
                                              {{4,4}, {0,1}},
                                              {{5,2}, {0}},
                                              {{5,3}, {1}},
                                              {{6,4}, {0,1}},
                                              {{7,3}, {0}},
                                              {{7,3}, {0}},
                                              {{8,3}, {1}}};
map<pair<int,int>, vector<int>> beta_map = { {{0,0}, {0}},
                                             {{2,3}, {0}},
                                             {{3,3}, {0,1}},
                                             {{4,4}, {0,1}},
                                             {{5,5}, {0,1}}};

vector<TM_RHS> get_RHS(const Ddata& D)
{
  vector<bigint> alist;
  alist.push_back(one);
  vector<bigint> plist;

  // p=2
  int aN = D.NN.alpha, aD = D.alpha;
  switch(aN*10+aD) {
  case 2:
    alist = multiply_list_by_powers(two, {0,3}, alist);
    break;
  case 22:
  case 53:
  case 83:
    alist = multiply_list_by_powers(two, {1}, alist);
    break;
  case 24:
  case 34:
  case 44:
  case 64:
    alist = multiply_list_by_powers(two, {0,1}, alist);
    break;
  case 32:
    alist = multiply_list_by_powers(two, {1,2}, alist);
    break;
  case 33:
    alist = multiply_list_by_powers(two, {2}, alist);
    break;
  case 52:
  case 73:
  case 74:
    // alist = multiply_list_by_powers(two, {0}, alist);
    break;
  default:
    plist.push_back(two);
  }
  // p=3
  int bN = D.NN.beta, bD = D.beta;
  switch(bN*10+bD) {
  case 0:
  case 23:
    // alist = multiply_list_by_powers(3, {0}, alist);
    break;
  default:
    if (bN>=3)
      alist = multiply_list_by_powers(three, {0,1}, alist);
    else
      plist.push_back(three);
  }

  // p>3, multiplicative
  for (auto pi=D.NN.Mprimes.begin(); pi!=D.NN.Mprimes.end(); ++pi)
    plist.push_back(*pi);

  // p>3, additive
  for (auto pi=D.NN.Aprimes.begin(); pi!=D.NN.Aprimes.end(); ++pi)
    {
      bigint p = *pi;
      if (val(p,D.D)==2)
        alist = multiply_list_by_powers(p, {0,1}, alist);
      else
        plist.push_back(p);
    }

  std::sort(plist.begin(), plist.end());
  vector<TM_RHS> RHSs;
  for (auto ai=alist.begin(); ai!=alist.end(); ++ai)
    RHSs.push_back(TM_RHS(*ai, plist));
  return RHSs;
}

TM_RHS::operator string() const
{
  ostringstream s;
  s << a << ",[";
  for (auto pi=plist.begin(); pi!=plist.end(); ++pi)
    {
      if (pi!=plist.begin())
        s << ",";
      s << *pi;
    }
  s << "]";
  return s.str();
}

//#define DEBUG12

vector<cubic> get_cubics(const Ddata& DD)
{
  // (12) says that if val(3,disc(F))>=3 then restrict to b=c=0 (mod 3):
  int eqn12 = (DD.beta < 3);
  vector<cubic> Flist = reduced_cubics(DD.D,
                                       0, // 0 to exclude reducibles
                                       1, // 1 for GL2-equivalence
                                       0); // verbosity level
#ifdef DEBUG12
  cout<<"D="<<DD.D<<": "<<Flist.size()<<" cubics before applying (12): "<<Flist<<endl;
#endif
  vector<cubic> Flist2;
  for(auto Fi = Flist.begin(); Fi!=Flist.end(); ++Fi)
    {
      cubic F = *Fi;
      if (eqn12 || (div(three,F.b()) && div(three,F.c())))
        {
          F.normalise(); // 1 for gl2
          Flist2.push_back(F);
        }
#ifdef DEBUG12
      else
        {
          cout << "omitting F="<<F<<" as (12) not satisfied"<<endl;
        }
#endif
    }
#ifdef DEBUG12
  cout<<"D="<<DD.D<<": "<<Flist2.size()<<" cubics after applying (12): "<<Flist2<<endl;
#endif
  return Flist2;
}

//#define DEBUG_LOCAL_TEST

int TM_eqn::local_test()
{
#ifdef DEBUG_LOCAL_TEST
  cout<<"Local testing of "<<(string)(*this)<<endl;
#endif
  if (!modaCheck(F, RHS.a))
    {
#ifdef DEBUG_LOCAL_TEST
      cout<<" - fails modaCheck"<<endl;
#endif
        return 0;
    }

  vector<bigint> newprimes;

  for (auto pi=RHS.plist.begin(); pi!=RHS.plist.end(); ++pi)
    {
      bigint p = *pi;
      if (F.has_roots_mod(p)) // then F mod p roots, so we keep p as an RHS prime
        {
#ifdef DEBUG_LOCAL_TEST
          cout<<" - p = "<<p<<" is kept"<<endl;
#endif
          newprimes.push_back(p);
        }
      else // F mod p has no roots...
        {
          // check (13) is satisfied, F(u,v) is not 0 mod p
          if ((p>=3) && (val(p,DD.NN.N)==1) && !div(p,DD.D))
            {
#ifdef DEBUG_LOCAL_TEST
              cout<<" - p = "<<p<<" fails condition (13)"<<endl;
#endif
              return 0;
            }
          // Now the equation has no solutions with a positive power
          // of p on RHS, and we check whether it is satisfiable mod p
          // with p not dividing the RHS:
          if (!modpCheck(F, RHS.a, RHS.plist, p))
            {
#ifdef DEBUG_LOCAL_TEST
              cout<<" - p = "<<p<<" fails"<<endl;
#endif
              return 0; // no it is not, discard this RHS
            }
          // otherwise it is, we do *not* keep this prime and hence
          // only look for TME solutions with no p dividing the RHS
#ifdef DEBUG_LOCAL_TEST
          cout<<" - p = "<<p<<" passes and is removed"<<endl;
#endif
        }
    }
  RHS.plist = newprimes;
  return 1;
}

// for output (puts either N,D,[a,b,c,d],A,plist or [a,b,c,d],A,plist into string)
string TM_eqn::as_string(int ND) const
{
  ostringstream s;
  if (ND)
    s << DD.NN.N << "," << DD.D << ",";
  F.output(s);
  s << "," << (string)RHS;
  return s.str();
}

//#define DEBUG_IMPRIMITIVE

vector<TM_eqn> get_TMeqnsD(const Ddata& DD)
{
  vector<TM_eqn> TMeqns;
  vector<cubic> Flist = get_cubics(DD);
  vector<TM_RHS> RHSlist = get_RHS(DD);
#ifdef DEBUG12
  cout<<"D = "<<DD.D<<": "<<endl;
  cout<<"  "<<Flist.size()<<" cubics: "<<Flist<<endl;
  cout<<"  "<<RHSlist.size()<<" RHSs"<<endl;
#endif
  for (auto Fi=Flist.begin(); Fi!=Flist.end(); ++Fi)
    {
      cubic F = *Fi;
      bigint g = F.content();
      assert (div(g,six));
      int imprimitive = (g!=one);
      Ddata DD1 = DD;
// if F is imprimitive with content g, we check that g divides a (skip
// this equation of not), and fivide F and a by g, adjusting D:
      if (imprimitive)
        {
#ifdef DEBUG_IMPRIMITIVE
          cout << " - F="<<F<<" is not primitive, content is "<<g<<endl;
#endif
          F = F / g;
          DD1 = Ddata(DD.NN, DD.D/pow(g,4));
#ifdef DEBUG_IMPRIMITIVE
          cout << " - primitive F="<<F<<", with discriminant "<<DD1.D<<endl;
#endif
        }
      for (auto RHSi=RHSlist.begin(); RHSi!=RHSlist.end(); ++RHSi)
        {
          TM_RHS RHS = *RHSi;
          if (imprimitive && !div(g, RHS.a))
            {
#ifdef DEBUG_IMPRIMITIVE
              cout << " - skipping RHS value "<<RHS.a<<" as not a multiple of the content "<<g<<endl;
#endif
              continue; // skip this (F,RHS) pair
            }
          else
            {
#ifdef DEBUG_IMPRIMITIVE
              cout << " - dividing RHS value "<<RHS.a<<" by the content "<<g<<endl;
#endif
              RHS.a /= g;
            }
          TM_eqn T(DD1, F, RHS);
          if (T.local_test())
            TMeqns.push_back(T);
        }
    }
  return TMeqns;
}

vector<TM_eqn> get_TMeqnsN(const Ndata& NN)
{
  vector<TM_eqn> TMeqns;
  vector<Ddata> Dlist = get_discriminants(NN);
  for (auto Di=Dlist.begin(); Di!=Dlist.end(); ++Di)
    {
      vector<TM_eqn> TMeqnsD = get_TMeqnsD(*Di);
      TMeqns.insert(TMeqns.end(), TMeqnsD.begin(), TMeqnsD.end());
    }
  return TMeqns;
}

// for p||N and p not dividing D=disc(F) we require that F(u,v)=0 (mod
// p) has a nontrivial solution:
int local_test(const cubic& F, const Ddata& DD, const bigint& p)
{
  return (div(p,DD.D) || F.has_roots_mod(p));
}

istream& operator>>(istream& s, TM_eqn& tme)
{
  string tme_string;
  s >> tme_string;
  tme = TM_eqn(tme_string);
  return s;
}

// Construct from a string N,D,[a,b,c,d],A,plist

//#define DEBUG_PARSER

const string ignore_chars = "[](),";
const string space = " ";

TM_eqn::TM_eqn(const string& s)
{
#ifdef DEBUG_PARSER
  cout << "parsing string: "<<s<<endl;
#endif
  // Make a copy of s
  string ss(s);
  // replace commas and brackets by spaces
  for (auto ci=ss.begin(); ci!=ss.end(); ++ci)
    {
      char c = *ci;
      if (std::find(ignore_chars.begin(), ignore_chars.end(), c) != ignore_chars.end())
        ss.replace(ci, ci+1, space);
    }
#ifdef DEBUG_PARSER
  cout << " - after replacing "<<ignore_chars<<" by spaces: "<<ss<<endl;
#endif
  // read the first 7 integers:
  istringstream is(ss);
  bigint N, D, a, b, c, d, A;
  is >> N >> D >> a >> b >> c >> d >> A >> ws;
  Ndata NN(N);
  DD = Ddata(NN, D);
  F = cubic(a,b,c,d);
  unimod m;
  F.sl2_reduce(m);
  F.normalise();
  // now read the primes (there may be none)
  vector<bigint> plist;
  bigint p;
  while (!is.eof())
    {
      is >> p;
      plist.push_back(p);
      is >> ws;
    }
  // Sort the primes in the list:
  std::sort(plist.begin(), plist.end());
#ifdef DEBUG_PARSER
  cout << " N="<<N<<", D="<<D<<", F=["<<a<<","<<b<<","<<c<<","<<d<<"], A="<<A<<", plist = "<<plist<<endl;
#endif
  RHS = TM_RHS(A, plist);
}

// Read TM equations from a file
vector<TM_eqn> read_TMeqns(const string& filename)
{
  vector<TM_eqn> TM_eqns;
  string line;
  ifstream data(filename.c_str());
  if (data)
    {
      while (data>>line)
        {
          TM_eqn T(line);
          TM_eqns.push_back(T);
        }
      data.close();
    }
  return TM_eqns;
}

// Write TM equations to a file or stdout
void write_TMeqns(vector<TM_eqn> TMEs, const string& filename)
{
  if (filename.compare("stdout"))
    {
      ofstream out(filename.c_str());
      for (auto Ti = TMEs.begin(); Ti!=TMEs.end(); ++Ti)
        out << (string)(*Ti) << endl;
      out.close();
    }
  else
    {
      for (auto Ti = TMEs.begin(); Ti!=TMEs.end(); ++Ti)
        cout << (string)(*Ti) << endl;
    }
}

int compare_TM_eqn_lists(const vector<TM_eqn>& list1, const vector<TM_eqn>& list2, int verbose)
{
  int n1 = list1.size();
  int n2 = list2.size();
  if (n1!=n2)
    {
      if (verbose)
        cout<<"List 1 has size "<<n1<<" but list 2 has size "<<n2<<endl;
      else
        return 0;
    }
  int nmatches=0;
  vector<int> check1(n1,0), check2(n2,0);
  vector<int>::iterator c1=check1.begin();
  for (auto T1i=list1.begin(); T1i!=list1.end(); ++T1i)
    {
      TM_eqn T1 = *T1i;
      int found = 0;
      vector<int>::iterator c2=check2.begin();
      for (auto T2i=list2.begin(); T2i!=list2.end() && (!found); ++T2i)
        {
          if (*c2==0)
            {
              if(verbose>1)
                cout<<"comparing "<<(string)T1<<" and "<<(string)(*T2i)<<"..."<<flush;
              found = T1.is_gl2_equivalent(*T2i);
              if(found)
                {
                  *c1=1;
                  *c2=1;
                  if(verbose>1)
                    cout<<"equivalent"<<endl;
                }
              else
                {
                  if(verbose>1)
                    cout<<"different"<<endl;
                }
            }
          c2++;
        }
      if (found)
        {
          nmatches++;
        }
      else
        {
          if (verbose)
            cout<<"No match for "<<(string)T1<<" in list 2!"<<endl;
        }
      c1++;
      if(verbose>1)
        {
          cout<<"After testing "<<(string)T1<<", check1 = "<<check1<<" and check2 = "<<check2<<endl;
        }
    }
  int ok = (nmatches==n1) && (nmatches==n2);
  if (verbose && !ok)
    {
      if (nmatches<n1)
        cout<<(n1-nmatches)<<" equations in list 1 not in list 2"<<endl;
      if (nmatches<n2)
        {
          cout<<(n2-nmatches)<<" equations in list 2 not in list 1:"<<endl;
          vector<int>::iterator c2=check2.begin();
          for (auto T2i=list2.begin(); T2i!=list2.end(); ++T2i, ++c2)
            if (*c2==0)
              cout<<(string)(*T2i)<<endl;
        }
    }
  return ok;
}
