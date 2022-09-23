/*
multGroup.m

Given an ideal fb in OK and elements alpha_1,...,alpha_n of K^* whose support
is disjoint from fb, these functions determine an abstract abelian group
isomorphic to (OK/fb)^* and the images of alpha_1,...,alpha_n in this group.
This is to replace Magma's UnitGroup(quo<OK|fb>), which is too slow when fb is
divisible by a very large power of a prime ideal, and is based on Section 4.2 of
H. Cohen, "Advanced Topics in Computational Algebraic Number Theory", Springer.

Authors
    Adela Gherga <adelagherga@gmail.com>
    Samir Siksek <samir.siksek@gmail.com>
Created
    18 August 2022
*/

residueClassImage:=function(delta,modfa,fa)
    /*
      Determine the image of delta, a non-zero element of K, in OK/fa.
      Here, the denominator ideal of delta is coprime to fa.

      Parameters
          delta: FldOrdElt
              A non-zero element of K.
          modfa: Map
              The residue class map OK --> OK/fa.
          fa: RngOrdIdl
              An ideal of OK.
      Returns
          u: FldFinElt
              The image of delta in OK/fa.
   */
    K:=Universe([delta] cat Basis(fa));
    OK:=MaximalOrder(K);
    if delta in OK then
	return modfa(delta);
    end if;
    facts:=Factorisation(Denominator(delta)*OK);
    facts:=[fact[1] : fact in facts];
    facts:=[fact : fact in facts | Valuation(delta*OK,fact) lt 0];
    ppow:=&*[fact^(-Valuation(delta*OK,fact)) : fact in facts];
    assert ppow+fa eq 1*OK;
    beta:=CRT(fa,ppow,OK!1,OK!0);
    // beta is congruent to 1 mod fa
    // and 0 mod ppow (which the denominator ideal of delta).
    assert beta in OK;
    alpha:=OK!(delta*beta);
    return modfa(alpha)/modfa(beta);
end function;

multQuo:=function(fp,a,b,alphaList)
    /*
      Given a prime ideal fp in OK, integers a, b such that b>a, and elements
      alpha_1,...,alpha_n of (1+fp^a), determine an abstract abelian group
      isomorphic to (1+fp^a)/(1+fp^b) and the images of alpha_1,...,alpha_n in
      this group.

      Parameters
          fp: RngOrdIdl
              A prime ideal of OK.
          a: RndIntElt
          b: RndIntElt
              An integer such that b>a.
          alphaList: SeqEnum
              A list of elements of (1+fp^a).
      Returns
          D: GrpAb
	      An abstract abelian group isomorphic to (1+fp^a)/(1+fp^b).
	  alphaListImage: SeqEnum
	      A list of the images of the elements in alphaList in D.
   */
    assert b gt a;
    K:=NumberField(Universe(Basis(fp)));
    n:=Degree(K);
    OK:=MaximalOrder(K);
    assert &and[Valuation((alpha-1)*OK,fp) ge a : alpha in alphaList];
    fpb:=fp^b;
    _,modfpb:=quo<OK | fpb>;
    p:=Characteristic(ResidueClassField(fp));
    e:=Valuation(p*OK,fp);
    assert e ge 1;
    k0:=1+Floor(e/(p-1));
    if (b gt 2*a) and (a lt k0) then
	D2,l2,D2gens:=$$(fp,a,2*a,alphaList);
	betaList:=[];
	assert #l2 eq #alphaList;
	for i in [1..#alphaList] do
	    alpha:=alphaList[i];
	    coords:=Eltseq(l2[i]);
	    assert #coords eq #D2gens;
	    beta:=(modfpb(alpha)/
		   &*[modfpb(D2gens[j])^coords[j] : j in [1..#coords]])@@modfpb;
	    assert Valuation(beta-1,fp) ge 2*a;
	    Append(~betaList,beta);
	end for;
	s:=#OrderedGenerators(D2);
	invsD2:=Invariants(D2);
	betaList:=[D2gens[j]^invsD2[j] : j in [1..s]] cat betaList;
	D1,l1:=$$(fp,2*a,b,betaList);
	r:=#OrderedGenerators(D1);
	A:=FreeAbelianGroup(r+s);
	invsD1:=Invariants(D1);
	rels:=[invsD1[i]*A.i : i in [1..r]];
	for j in [1..s] do
	    Append(~rels,A!(Eltseq(l1[j]) cat
			    [0 : i in [1..s]])-invsD2[j]*A.(r+j));
	end for;
	D,pi:=quo<A | rels>;
	l1:=l1[s+1..#l1];
	assert #l1 eq #alphaList;
	assert #l2 eq #alphaList;
	l:=[];
	for j in [1..#alphaList] do
	    Append(~l,pi(A!(Eltseq(l1[j]) cat Eltseq(l2[j]))));
	end for;
	return D,l;
    end if;
    if b le 2*a then
	fpb:=fp^b;
	fpa:=fp^a;
	GOK:=FreeAbelianGroup(n);
	A:=sub<GOK | [GOK!Eltseq(OK!l) : l in Basis(fpa)]>;
	B:=sub<GOK | [GOK!Eltseq(OK!l) : l in Basis(fpb)]>;
	C,modB:=quo<GOK|B>;
	D:=sub<C | [modB(h) : h in OrderedGenerators(A)]>;
	assert #C eq Norm(fp)^b;
	assert #D eq Norm(fp)^(b-a);
	alphaList:=[GOK!(Eltseq(OK!(alpha-1))) : alpha in alphaList];
	assert &and[alpha in A : alpha in alphaList];
	alphaList:=[D!(modB(alpha)) : alpha in alphaList];
	Dgens:=[OK!(Eltseq(d@@modB)) : d in OrderedGenerators(D)];
	assert &and[d in fpa : d in Dgens];
	Dgens:=[1+d : d in Dgens];
	return D,alphaList,Dgens;
    end if;
    if a ge k0 then
	fpb:=fp^b;
	fpa:=fp^a;
	GOK:=FreeAbelianGroup(n);
	A:=sub<GOK | [GOK!Eltseq(OK!l) : l in Basis(fpa)]>;
	B:=sub<GOK | [GOK!Eltseq(OK!l) : l in Basis(fpb)]>;
	C,modB:=quo<GOK|B>;
	D:=sub<C | [modB(h) : h in OrderedGenerators(A)]>;
	assert #C eq Norm(fp)^b;
	assert #D eq Norm(fp)^(b-a);
	OKfp,mapToComp:=Completion(OK,fp : Precision:=b+1);
	alphaList:=[(Log(mapToComp(alpha)))@@mapToComp : alpha in alphaList];
	assert &and[alpha in fpa : alpha in alphaList];
	alphaList:=[GOK!Eltseq(OK!alpha) : alpha in alphaList];
	assert &and[alpha in A : alpha in alphaList];
	alphaList:=[D!(modB(alpha)) : alpha in alphaList];
	return D,alphaList;
    end if;
end function;

multGroupPrimePower:=function(fp,k,alphaList)
    /*
      Given a prime ideal fp in OK and elements alpha_1,...,alpha_n of OK whose
      support is disjoint from fp, determine an abstract abelian group
      isomorphic to (OK/fp)^* and the images of alpha_1,...,alpha_n in this
      group. This algorithm is based on Proposition 4.2.4 of H. Cohen,
      "Advanced Topics in Computational Algebraic Number Theory", Springer.

      Parameters
          fp: RngOrdIdl
              A prime ideal of OK.
          k: RndIntElt
              A positive integer.
          alphaList: SeqEnum
              A list of elements of OK whose support is coprime to fp.
      Returns
          D: GrpAb
	      An abstract abelian group isomorphic to (OK/fp)^*.
	  alphaListImage: SeqEnum
	      A list of the images of the elements in alphaList in D.
   */
    K:=NumberField(Universe(Basis(fp)));
    n:=Degree(K);
    OK:=MaximalOrder(K);
    assert &and[Valuation(alpha*OK,fp) eq 0 : alpha in alphaList];
    assert k ge 1;
    if k eq 1 then
	R,modfp:=ResidueClassField(fp);
	D1,phi1:=MultiplicativeGroup(R);
	return D1,[(modfp(alpha))@@phi1 : alpha in alphaList];
    end if;
    q:=Norm(fp);
    fpk:=fp^k;
    _,modfpk:=quo<OK | fpk>;
    prec:=Ceiling(Log(k)/Log(2));
    xList:=[];
    yList:=[];
    for alpha in alphaList do
	z:=modfpk(alpha);
	x:=z;
	for i in [1..prec] do
	    x:=x - (x^(q-1)-1)/(modfpk(q-1)*x^(q-2));
	end for;
	assert x^(q-1) eq 1;
	y:=(z/x)@@modfpk;
	x:=x@@modfpk;
	assert Valuation(x-alpha,fp) ge 1;
	assert Valuation(y-1,fp) ge 1;
	Append(~xList,x);
	Append(~yList,y);
	// We've written each alpha as alpha = x y mod fp^k where
	// x is the Teichmueller lift of alpha, and y = 1 mod fp^k.
    end for;
    R,modfp:=ResidueClassField(fp);
    D1,phi1:=MultiplicativeGroup(R);
    xList:=[(modfp(x))@@phi1 : x in xList];
    D2,yList:=multQuo(fp,1,k,yList);
    D,phi1,phi2:=DirectSum(D1,D2);
    assert #xList eq #yList;
    return D,[phi1(xList[i])+phi2(yList[i]) : i in [1..#xList]];
end function;

multGroup:=function(fb,alphaList)
    /*
      Given an ideal fb in OK and elements alpha_1,...,alpha_n of K^* whose
      support is disjoint from fb, determine an abstract abelian group
      isomorphic to (OK/fb)^* and the images of alpha_1,...,alpha_n in this
      group. This is to replace Magma's UnitGroup(quo<OK|fb>), which is too
      slow when fb is divisible by a very large power of a prime ideal, and is
      based on Section 4.2 of H. Cohen, "Advanced Topics in Computational
      Algebraic Number Theory", Springer.

      Parameters
          fb: RngOrdIdl
              An ideal of OK.
          alphaList: SeqEnum
              A list of elements of K^* whose support is coprime to fb.
      Returns
          D: GrpAb
	      An abstract abelian group isomorphic to (OK/fb)^*.
	  alphaListImage: SeqEnum
	      A list of the images of the elements in alphaList in D.
   */
    K:=NumberField(Universe(Basis(fb)));
    n:=Degree(K);
    OK:=MaximalOrder(K);
    facts:=Factorisation(fb);
    D:=AbelianGroup([]);
    alphaListImage:=[D!0 : i in [1..#alphaList]];
    if #facts eq 0 then
	return D,alphaListImage;
    end if;
    _,modfb:=quo<OK | fb>;
    alphaList:=[(residueClassImage(alpha,modfb,fb))@@modfb
		: alpha in alphaList];
    // Replaced each alpha by an element of OK representing the same class
    // modulo fb.
    for fact in facts do
	fp:=fact[1];
	k:=fact[2];
	D1,alphaListImage1:=multGroupPrimePower(fp,k,alphaList);
	D,phi,psi:=DirectSum(D,D1);
	assert #alphaListImage eq #alphaListImage1;
	alphaListImage:=[phi(alphaListImage[i])+psi(alphaListImage1[i])
			 : i in [1..#alphaListImage]];
    end for;
    return D,alphaListImage;
end function;
