/*
equationsInK.m

These functions determine tau, delta_1,...,delta_r in K, giving all possible
equations a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}
(as in equation (15)).

Authors
    Adela Gherga <adelagherga@gmail.com>
    Samir Siksek <samir.siksek@gmail.com>
Created
    9 August 2022
*/

adequateLM:=function(alpha,beta,p)
    /*
      For a rational prime p, determine sets Lp and Mp, adequate for
      (alpha,beta), as in Algorithm 1.

      Parameters
          alpha: FldOrdElt
          beta: RngIntElt
          p: RndIntElt
              A rational prime.
      Returns
          Lp: SetEnum
              A list of ideals fb.
          Mp: SetEnum
              A list of pairs (fb,fp), where fp is a prime ideal with
	      ramification index = inertial degree = 1.
   */
    OK:=MaximalOrder(Universe([alpha,beta]));
    facts:=Factorisation(p*OK);
    facts:=[fact[1] : fact in facts]; // The prime ideals above p.
    cB:=[];
    cU:={};
    fb:=1*OK;
    for fp in facts do
        vbeta:=Valuation(beta*OK,fp);
        valpha:=Valuation(alpha*OK,fp);
        fb:=fb*fp^(vbeta+Min(valpha,0));
        if valpha ge 0 then
            R,modfp:=ResidueClassField(fp);
            Rp:=PrimeField(R);
            u:=residueClassImage(-alpha,modfp,fp);
            if u in Rp then
                Append(~cB,fp);
                u:=Integers()!(Rp!u);
                assert (u ge 0) and (u le p-1);
                Include(~cU,u);
            end if;
        end if;
    end for;
    if #cB eq 0 then
        return {fb},{};
    end if;
    if #cB eq 1 then
        P:=cB[1];
        if InertiaDegree(P) eq 1 and RamificationDegree(P) eq 1 then
            return {},{<fb,P>};
        end if;
    end if;
    M:={};
    if #cU eq p then
        L:={};
    else
        L:={fb};
    end if;
    for u in cU do
        add1,add2:=$$((u+alpha)/p,p*beta,p);
        L:=L join add1;
        M:=M join add2;
    end for;
    return L join {1*OK},M;
end function;

adequateLMRefined:=function(a0,theta,p)
    /*
      For a rational prime p, determine sets Lp and Mp, adequate
      for (-theta/a_0,a_0), including refinements.

      Parameters
          a0: RngIntElt
          theta: RngOrdElt
              A generator of K.
          p: RndIntElt
              A rational prime.
      Returns
          Lp: SetEnum
              A list of ideals fb.
          Mp: SetEnum
              A list of pairs (fb,fp), where fp is a prime ideal with
	      ramification index = inertial degree = 1.
   */
    OK:=MaximalOrder(Parent(theta));
    Lp,Mp:=adequateLM(-theta/a0,a0,p);

    // Apply refinements.
    Mp:={<pr[1]/pr[2]^Valuation(pr[1],pr[2]), pr[2]> : pr in Mp};
    for pr in Mp do
	fbd:=pr[1];
	fp:=pr[2];
	Lp:={fb : fb in Lp | IsIntegral(fb/fbd) eq false or
			     fb/fbd ne fp^Valuation(fb/fbd,fp)};
    end for;
    if IsDivisibleBy(a0,p) then
	d:=Degree(NumberField(OK));
	Mp1:={};
	Lp:={fb : fb in Lp | Valuation(Norm(fb),p) ge (d-1)*Valuation(a0,p)};
	for pr in Mp do
	    fbd:=pr[1];
	    fp:=pr[2];
	    while Valuation(Norm(fbd),p) le (d-1)*Valuation(a0,p) do
		if Valuation(Norm(fbd),p) eq (d-1)*Valuation(a0,p) then
		    Mp1:=Mp1 join {<fbd,fp>};
		    break;
		else
		    fbd:=fbd*fp;
		end if;
	    end while;
	end for;
	Mp:=Mp1;
    end if;
    return Lp,Mp;
end function;

normInv:=function(a0,theta,R)
    /*
      Determine all ideals of OK with norm R for some (a,S) in cZ,
      including refinements.

      Parameters
          a0: RngIntElt
          theta: RngOrdElt
              A generator of K.
          R: RndIntElt
              A positive integer.
      Returns
          cR: SetEnum
              All ideals of OK with norm R.
   */
    OK:=MaximalOrder(Parent(theta));
    assert R in Integers();
    R:=Integers()!R;
    assert R ge 1;
    if R eq 1 then
	return {1*OK};
    end if;
    p:=Max(PrimeDivisors(R));
    cR:=[];

    Lp,Mp:=adequateLM(-theta/a0,a0,p);
    for fb in Lp do
	if Valuation(Norm(fb),p) eq Valuation(R,p) then
	    Append(~cR,fb);
	end if;
    end for;
    for pr in Mp do
	fb:=pr[1];
	fp:=pr[2];
	while Valuation(Norm(fb),p) le Valuation(R,p) do
	    if Valuation(Norm(fb),p) eq Valuation(R,p) then
		Append(~cR,fb);
		break;
	    else
		fb:=fb*fp;
	    end if;
	end while;
    end for;

    if IsEmpty(cR) then
	return {};
    else
	return &join{{fp*fa : fa in $$(a0,theta, R div Norm(fp))} : fp in cR};
    end if;
end function;

idealEquations:=function(alist,a,primelist)
    /*
      Determine all possible sets (fa,S). Here, fa is an ideal of K and S is a
      list of prime ideals of K, fp_1,...,fp_s, such that
      (a_0 X - theta Y)OK = fa fp_1^{n_1} ... fp_s^{n_s} (as in equation (13)).

      Parameters
          alist: SeqEnum
              A list of integer coefficients a_0, a_1,...,a_d.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
      Returns
          afplist: List
              A list of pairs (fa,S), where S is a list of prime ideals of K,
	      fp_1,...,fp_s.
   */
    assert &and[IsPrime(p) : p in primelist];
    assert &and[Valuation(a,p) eq 0 : p in primelist];
    assert &and[a_i in Integers() : a_i in alist];
    a0:=Integers()!alist[1];
    assert a0 ne 0;
    d:=#alist-1;
    assert d ge 3;
    QUV<U,V>:=PolynomialRing(Rationals(),2);
    Qx<x>:=PolynomialRing(Rationals());
    F:=&+[alist[i+1]*U^(d-i)*V^i : i in [0..d]];
    assert IsHomogeneous(F);
    f:=a0^(d-1)*Evaluate(F,[x/a0,1]);
    assert IsMonic(f);
    assert Degree(f) eq d;
    assert IsIrreducible(f);
    assert &and[a_i in Integers() : a_i in Coefficients(f)];
    K<theta>:=NumberField(f);
    OK:=MaximalOrder(K);
    theta:=OK!theta;

    afplist:=[* [* 1*OK, [] *] *];
    for p in primelist do
	afplistNew:=[* *];
	Lp,Mp:=adequateLMRefined(a0,theta,p);
	afplistNew1:=[* [* pr[1]*fb,pr[2] *] : fb in Lp,pr in afplist *];
	afplistNew2:=[* [* pr[1]*qr[1],pr[2] cat [qr[2]] *]
		      : qr in Mp, pr in afplist *];
	afplist:=afplistNew1 cat afplistNew2;
    end for;

    R:=Integers()!(a*a0^(d-1));
    afplistNew:=[* *];
    for pr in afplist do
	af:=pr[1];
	R0:=AbsoluteValue(R div GCD(Norm(af),R));
	cR:=normInv(a0,theta,R0);
	afplistNew:=afplistNew cat [* [* af*I, pr[2] *] : I in cR *];
    end for;
    afplist:=afplistNew;
    for pr in afplist do // Running sanity checks.
	af:=pr[1];
	fplist:=pr[2];
	assert &and[InertiaDegree(fq)*RamificationDegree(fq) eq 1
		    : fq in fplist];
	normlist:=[Norm(fq) : fq in fplist];
	assert #Set(normlist) eq #normlist;
	assert Set(normlist) subset Set(primelist);
	Naf:=Norm(af);
	assert IsDivisibleBy(Naf,a*a0^(d-1));
	assert Set(PrimeDivisors(Naf div (a*a0^(d-1)))) subset Set(primelist);
    end for;
    return afplist;
end function;

principalize:=function(fa,S)
    /*
      From the pair (fa,S), where S = (fp_1,...,fp_r) and
      (a_0 X - theta Y)OK = fa fp_1^{n_1} ... fp_s^{n_s} is an equation of
      ideals in K, determine delta_1,...,delta_r and all possible tau such that
      a_0 X - Y theta = tau \delta_1^{b_1} ... \delta_r^{b_r} is an equation
      in K (as in equation (15)).
      If the class [fa] is not contained in the subgroup of the class group
      generated by [fp_1,...,fp_r] then then is is impossible; the function
      returns false,[],[].

      Parameters
          fa: RngOrdIdl
              An ideal of OK.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
      Returns
          tf: BoolElt
              A true/false value. If the class [fa] is not contained in the
	      subgroup of the class group generated by [fp_1,...,fp_r] then
	      then is is impossible; return false, otherwise return true.
          tauList: SeqEnum
              The set of possibilities for tau in K. This is equal to 1/2 of the
              number of roots of unity in K (since we can absorb \pm 1
	      in (X,Y)).
          deltaList: SeqEnum
              A list of elements delta_1,...,delta_r in K which form the basis
              for the S-unit group OS^* modulo roots of unity.
   */
    K:=NumberField(Parent(Basis(fa)[1]));
    OK:=MaximalOrder(K);
    s:=#S;
    Zs:=FreeAbelianGroup(s);
    ClK,eta:=ClassGroup(K);
    // Given an ideal fa, we obtain its class in ClK by fa@@eta.
    phi:=hom<Zs -> ClK | [ClK | fp@@eta : fp in S ]>;
    img:=Image(phi);
    if -fa@@eta in img eq false then
	return false,[K | ],[K | ];
    end if;
    if s eq 0 then
	fa2:=fa;
    else
	rr:=(-fa@@eta)@@phi;
	rr:=Eltseq(Zs!rr);
	fa2:=fa*&*[S[i]^rr[i] : i in [1..s]];
    end if;
    tf,alpha:=IsPrincipal(fa2);
    assert tf; // The ideal really is principal, and equal to alpha.
    if #S eq 0 then
	U,mu:=UnitGroup(K);
    else
	U,mu:=SUnitGroup(S);
    end if;
    deltaList:=[K!(mu(u)) : u in OrderedGenerators(U)];
    rtun:=deltaList[1];
    assert IsRootOfUnity(rtun);
    deltaList:=deltaList[2..#deltaList];

    // Reorder deltaList so that the units apear first.
    deltaList:=[delta : delta in deltaList | Abs(Norm(delta)) eq 1] cat
	       [delta : delta in deltaList | Abs(Norm(delta)) ne 1];
    r:=Order(U.1) div 2;
    assert rtun^r eq OK!(-1);
    tauList:=[alpha*rtun^i : i in [0..(r-1)]];
    return true,tauList,deltaList;
end function;

equationsInK:=function(alist,a,primelist)
    /*
      Determine a list of all possible elements in K, tau, delta_1,...,delta_r
      so that a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r} in K
      (as in equation (15)).

      Parameters
          alist: SeqEnum
              A list of coefficients a_0, a_1,...,a_d.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
      Returns
          tauDeltaList: List
              A list of pairs (tau,deltaList), where deltaList is
              the list delta_1,...,delta_r.
   */
    afplist:=idealEquations(alist,a,primelist);
    if (#afplist eq 1) then
	printf "There is %o ideal equations to solve.\n",#afplist;
    else
	printf "There are %o ideal equations to solve.\n",#afplist;
    end if;
    tauDeltaList:=[* *];
    for pr in afplist do
	af:=pr[1];
	fplist:=pr[2];
	tf,tauList,deltaList:=principalize(af,fplist);
	if tf then
	    for tau in tauList do
		tauDeltaList:=tauDeltaList cat [* [* tau,deltaList *] *];
	    end for;
	end if;
    end for;
    ranks:=[#pr[2] : pr in tauDeltaList];
    // Sort the elements in order of rank for SmallSieveInfo.
    Sort(~ranks,~permut);
    return [* tauDeltaList[i^permut] : i in [1..#tauDeltaList] *];
end function;
