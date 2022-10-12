#include <eclib/marith.h>
#include <eclib/unimod.h>
#include <eclib/polys.h>
#include <eclib/cubic.h>
#include <assert.h>

#include "cubic_utils.h"

// multiply all integers in a list by a constant:
vector<bigint> multiply_list(const bigint& a, const vector<bigint>& L)
{
  vector<bigint> aL = L;
  for (auto x=aL.begin(); x!=aL.end(); ++x)
    (*x) *= a;
  return aL;
}

// multiply all integers in a list by all in a second list:
vector<bigint> multiply_lists(const vector<bigint>& L1, const vector<bigint>& L2)
{
  vector<bigint> L3;
  L3.reserve(L1.size()*L2.size());
  for (auto x=L1.begin(); x!=L1.end(); ++x)
    for (auto y=L2.begin(); y!=L2.end(); ++y)
      L3.push_back((*x)*(*y));
  return L3;
}

// multiply all integers in L by p^e for e in exponents:
vector<bigint> multiply_list_by_powers(const bigint& p, const vector<int>& exponents, const vector<bigint>& L)
{
  vector<bigint> peL;
  peL.reserve(L.size()*exponents.size());
  for (auto e=exponents.begin(); e!=exponents.end(); ++e)
    {
      bigint x = pow(p, *e);
      for (auto y=L.begin(); y!=L.end(); ++y)
        peL.push_back(x*(*y));
    }
  return peL;
}

int is_cube(const bigint& a, const bigint& q)
{
  if (div(q,a) || div(q,a-1) || div(3,q+1))
    return 1;
  bigint b;
  power_mod(b, a, (q-1)/3, q); // b = a^(q-1)/3 mod q
  return div(q,b-1);
}

// for q prime, returns a list of representatives of the values of
// F(u,v) mod q modulo cubes
vector<bigint> image_mod_cubes(const cubic& F, const bigint& q)
{
  vector<bigint> images;

  // first see if 0 is a value:

  //  vector<bigint> coeffs = {F.a(), F.b(), F.c(), F.d()};
  if (F.has_roots_mod(q))
    images.push_back(BIGINT(0));

  // if q=2 (mod 3) or q=3, then all nozero values occur and all are cubes:

  if (div(3,q+1) || q==3)
    {
      images.push_back(BIGINT(1));
      return images;
    }

  // Now q=1 (mod 3) and we must see which of the three nonzero cosets mod cubes are hit:

  // We keep track of the inverses of the images found to simplify the coset check
  vector<bigint> inverses;

  // check F(1,0):
  bigint v = F.a();
  if (v!=0)
    {
      images.push_back(v);
      inverses.push_back(invmod(v, q));
    }

  // check F(0,x) for all x mod q:
  for (bigint x(0); x<q; x++)
    {
      v = F.eval(x);
      if (v==0)
        continue;
      int repeat=0;
      for (auto w = inverses.begin(); !repeat && w!=inverses.end(); ++w)
        repeat = is_cube((*w) * v, q);
      if (!repeat) // we have a new coset
        {
          images.push_back(v);
          inverses.push_back(invmod(v, q));
        }
      if (inverses.size()==3)
        break;
    }
  return images;
}

// (function name as in AG's magma code) Return 1 iff there exist
// (u,v) not (0,0) mod q and exponents e such that F(u,v)=a*prod(p^e)
// mod q.
int modpCheck(const cubic& F, const bigint& a, const vector<bigint>& primes, const bigint& q)
{
  // So we don't have to construct copies of primes with q removed:
  if (std::find(primes.begin(), primes.end(), q)!=primes.end())
    return 1;

  if (div(2,q+1)) // then F takes one and hence all nonzero values since all are cubes
    return 1;

  for (auto pi = primes.begin(); pi!=primes.end(); ++pi)
    if (!is_cube(*pi,q)) // then powers of p cover all cosets mod cubes
      return 1;

  // Now all p are cubes mod q so can be ignored, we just check if F
  // takes the value a (mod cubes)
  vector<bigint> images = image_mod_cubes(F, q);
  bigint b = invmod(a,q);
  for (auto ci = images.begin(); ci!=images.end(); ++ci)
    {
      if (is_zero(*ci)) // if 0 is a value, ignore it
        continue;
      if (is_cube(b*(*ci),q)) // i.e. a=c mod cubes
        return 1;
    }
  return 0;
}

// similar to AG's modpCheckdivrhs. Return 1 iff there exists
// primitive (u,v) such that F(u,v)=0 (mod a).
int modaCheck(const cubic& F, const bigint& a)
{
  if (a==1)
    return 1;

  bigint g = F.content();
  if (!div(g,a))
    return 0;

  vector<bigint> plist = pdivs(a);

  if (plist.size()>1) // use CRT
    {
      for (auto pi = plist.begin(); pi!=plist.end(); ++pi)
        {
          bigint q = pow(*pi, val(*pi, a));
          if (!modaCheck(F, q))
            return 0;
        }
      return 1;
    }

  // Now a is a prime power p^e
  bigint p = plist[0];
  if (!F.has_roots_mod(p))
    return 0;
  if (val(p,a)==1) // a=p, nothing more to do
    return 1;

  // Now we must see if any projective root mod p lifts to p^e. Note
  // that in our application we expect p to divide the discriminant of
  // F, so roots mod p will not be simple, otherwise Hensel would make
  // this redundant.

  bigint b = a/p; // = p^{e-1}
  bigint one(1);
  // Test for roots above (1:0) mod p:
  if (div(p,F.a()))
    {
      // Test (u,v) = (1, p*w) for w mod p^{e-1}
      for (bigint w(0); w<b; w++)
        if (div(a,F.eval(one, p*w)))
          return 1;
    }
  // find affine roots r mod p:
  vector<bigint> roots = F.roots_mod(p);
  for (auto r = roots.begin(); r!=roots.end(); ++r)
    {
      // Test (u,v) = (r+p*w, 1) for w mod p^{e-1}
      for (bigint w(0); w<b; w++)
        if (div(a,F.eval(p*w + *r, one)))
          return 1;
    }
  // None of the roots did lift, so we fail
  return 0;
}
