/*
sift.m

For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, these functions
sieve for all solutions (X,Y,b_1,...,b_r). Here, the solution vectors
(b_1,...,b_r) are known to live in the intersection of certain translated
lattices, having norm <= cBf; the sieve recovers all such vectors.

Authors
    Adela Gherga <adelagherga@gmail.com>
    Samir Siksek <samir.siksek@gmail.com>
Created
    12 August 2022
*/

detSVP:=function(Zr,L)
    /*
      Compute the determinant of L regarded as a lattice of Z^r and the norm
      of the shortest vector contained in L.

      Parameters
          Zr: GrpAb
              The abelian group isomorphic to Z^r.
          L: GrpAb
              An abelian group defining a lattice of Z^r.
      Returns
	  det: RngIntElt
              The determinant of the lattice L.
          lambdaL: RngIntElt
              The squared L^2 norm of the shortest vector of the lattice L.
   */
    assert L subset Zr;
    r:=#Invariants(Zr);
    ZZr:=StandardLattice(r);
    LL:=sub<ZZr | [ Eltseq(Zr!l) : l in OrderedGenerators(L)]>;
    lambdaL:=Norm(ShortestVectors(LL)[1]);
    return Determinant(LL),lambdaL;
end function;

cVP:=function(L,w,distSq)
    /*
      Determine all vectors v of the lattice L such that the squared L^2 norm
      of (v - w) is <= distSq.

      Parameters
          L: Lat
	  w: LatElt
          distSq: RngIntElt
      Returns
	  cvp: SeqEnum
              A list of vectors v of the lattice L such that the squared L^2
	      norm of (v - w) is <= distSq.
   */
    assert distSq ge 0;
    if distSq gt 0 then
	cvp:=[];
	P:=CloseVectorsProcess(L,w,distSq);
	while (IsEmpty(P) eq false) do
	    Append(~cvp,NextVector(P));
	end while;
    else
	if w in L then
	    cvp:=[w];
	else
	    cvp:=[];
	end if;
    end if;
    return cvp;
end function;

vectorsInCoset:=function(Zr,L,w,distSq)
    /*
      Determine all vectors v of the translated lattice w + L having squared
      L^2 norm <= distSq.

      Parameters
          Zr: GrpAb
              The abelian group isomorphic to Z^r.
          L: GrpAb
              An abelian group defining a lattice of Z^r.
	  w: GrpAbElt
              An element of Z^r.
          distSq: RngIntElt
      Returns
	  vecs: SeqEnum
              A list of vectors v of the translated lattice w + L having squared
	      L^2 norm <= distSq.
   */
    assert L subset Zr;
    r:=#Invariants(Zr);
    ZZr:=StandardLattice(r);
    ww:=ZZr!(Eltseq(w));
    LL:=sub<ZZr | [ Eltseq(Zr!l) : l in OrderedGenerators(L)] >;
    vecs:=cVP(LL,-ww,distSq);
    vecs:=[ww+vv : vv in vecs];
    vecs:=[Eltseq(ZZr!vv) : vv in vecs];
    return vecs;
end function;

cosetIntersect:=function(Zr,w1,L1,w2,L2)
    /*
      Given lattices L1, L2 and vectors w1, w2, determine the lattice L3 and
      vector w3 of Z^r such that  w3 + L3 = (w1 + L1) intersect (w2 + L2),
      where possible.

      Parameters
          Zr: GrpAb
              The abelian group isomorphic to Z^r.
	  w1: GrpAbElt
              An element of Z^r.
          L1: GrpAb
              An abelian group defining a lattice of Z^r.
	  w2: GrpAbElt
              An element of Z^r.
          L2: GrpAb
              An abelian group defining a lattice of Z^r.

      Returns
          tf: BoolElt
              A true/false value. This value is true if
	      w3 + L3 = (w1 + L1) intersect (w2 + L2), and false if
              emptyset = (w1 + L1) intersect (w2 + L2).
          w3: GrpAbElt
          L3: GrpAb
   */
    if w1-w2 in L1+L2 eq false then
	return false,Zr!0,L1;
    end if;
    D,i1,i2,p1,p2:=DirectSum(L1,L2);
    kappa:=hom<D->Zr | [ p2(d)-p1(d) : d in OrderedGenerators(D)]>;
    L3:=L1 meet L2;
    dv:=(w1-w2)@@kappa;
    l1:=p1(dv);
    l2:=p2(dv);
    assert w1+l1 eq w2+l2;
    return true,w1+l1,L3;
end function;

deeperSift:=function(Zr,Lc,wc,cBfsq,bigInfs,depth)
    /*
      Determine lattice vectors b in the translated lattice wc + Lc such that
      the squared L^2 norm of b is <= cBfsq, sieved with prime ideals outside
      of S (as in Procedure 1).

      Parameters
          Zr: GrpAb
              The abelian group isomorphic to Z^r.
          Lc: GrpAb
              An abelian group defining a lattice of Z^r.
          wc: GrpAbElt
              An element of Z^r.
          cBfsq: FldReElt
              The final, squared upper bound for the L^2 norm of the vector
	      (b_1,...,b_r) after reduction.
	  bigInfs: List
	      For each entry (q,Bq,Piq,SqShifted) of smallInfs, this is the
	      pair (phiq,Tq), where phiq defines the map Z^r --> Bq, and Tq
	      is the set Sq intersected with piq(Z^r).
          depth: RngIntElt
              An integer denoting the depth of the recusion.
      Returns
          vecs: SeqEnum
              A list of vectors b of the translated lattice wc + Lc having
	      squared L^2 norm <= cBfsq.
   */
    rk:=#Invariants(Lc);
    det,lambdaL:=detSVP(Zr,Lc);
    if (lambdaL gt 4*cBfsq) or (#bigInfs eq 0) then
	vecs:=vectorsInCoset(Zr,Lc,wc,cBfsq);
	for I in bigInfs do
	    phiq:=I[1];
	    Tq:=I[2];
	    vecs:=[v : v in vecs | phiq(Domain(phiq)!v) in Tq];
	end for;
	return vecs;
    end if;
    imLList:=[#I[2]/Index(Zr,Kernel(I[1])+Lc) : I in bigInfs];
    Sort(~imLList,~permut);
    bigInfs:=[* bigInfs[i^permut] : i in [1..#bigInfs] *];
    bigInf:=bigInfs[1];
    bigInfs:=bigInfs[2..#bigInfs];
    phiq:=bigInf[1];
    Tq:=bigInf[2];
    phiqLc:=phiq(Lc);
    phiwc:=phiq(wc);
    Tqcut:=[t : t in Tq | (t-phiwc) in phiqLc];
    K:=Kernel(phiq);
    vecs:=[];
    for t in Tqcut do
	s:=t@@phiq;
	tf,wcNew,LcNew:=cosetIntersect(Zr,wc,Lc,s,K);
	if tf then
	    vecs:=vecs cat $$(Zr,LcNew,wcNew,cBfsq,bigInfs,depth+1);
	end if;
    end for;
    return vecs;
end function;

sift:=function(tau,deltaList,Zr,Lc,wc,SLeft,rangeLeft,cBfsq,bigInfs,depth)
    /*
      Determine lattice vectors b in the translated lattice wc + Lc such that
      the squared L^2 norm of b is <= cBfsq, sieved with prime ideals of S and
      those outside of S (as in Procedure 1).

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          Zr: GrpAb
              The abelian group isomorphic to Z^r.
          Lc: GrpAb
              An abelian group defining a lattice of Z^r.
          wc: GrpAbElt
              An element of Z^r.
          SLeft: SeqEnum
              A subset of the list of prime ideals fp_1,...,fp_s of K.
          rangeLeft: SeqEnum
              The mode up-to-date lower and upper bounds for the fp-valuation
              of eps, for fp in SLeft.
          cBfsq: FldReElt
              The final, squared upper bound for the L^2 norm of the vector
	      (b_1,...,b_r) after reduction.
	  bigInfs: List
	      For each entry (q,Bq,Piq,SqShifted) of smallInfs, this is the
	      pair (phiq,Tq), where phiq defines the map Z^r --> Bq, and Tq
	      is the set Sq intersected with piq(Z^r).
          depth: RngIntElt
              An integer denoting the depth of the recusion.
      Returns
          vecs: SeqEnum
              A list of vectors b of the translated lattice wc + Lc having
	      squared L^2 norm <= cBfsq.
   */
    assert #rangeLeft eq #SLeft;
    rk:=#Invariants(Lc);
    det,lambdaL:=detSVP(Zr,Lc);
    if lambdaL gt 4*cBfsq then // This is lambda(Lc)^2 > (2Bf)^2.
	vecs:=vectorsInCoset(Zr,Lc,wc,cBfsq);
	for I in bigInfs do
	    phiq:=I[1];
	    Tq:=I[2];
	    vecs:=[v : v in vecs | phiq(Domain(phiq)!v) in Tq];
	end for;
	return vecs;
    end if;
    if #SLeft eq 0 then
	return deeperSift(Zr,Lc,wc,cBfsq,bigInfs,depth+1);
    end if;
    K:=NumberField(Universe([tau] cat deltaList));
    theta:=K.1;
    d:=Degree(K);
    assert d ge 3;
    OK:=MaximalOrder(K);
    r:=#deltaList;
    fp:=SLeft[1];
    SLeft:=SLeft[2..#SLeft];
    krange:=rangeLeft[1];
    rangeLeft:=rangeLeft[2..#rangeLeft];
    assert #SLeft eq #rangeLeft;
    kmin:=krange[1];
    kmax:=krange[2]; // k is the valution of (a_0 X - theta Y) at fp.
    assert kmin ge 0;
    assert kmax ge kmin;
    p:=Characteristic(ResidueClassField(fp));
    assert Degree(fp) eq 1;
    assert Valuation(p*OK,fp) eq 1;
    fa:=p*OK/fp;
    for fact in Factorisation(fa) do
	assert &and[Valuation(delta*OK,fact[1]) eq 0 : delta in deltaList];
    end for;
    mpol:=MinimalPolynomial(theta);
    lcf:=LeadingCoefficient(mpol);
    assert lcf eq 1; // Surely theta is an algebraic integer!
    henPrec:=2*(Valuation(Discriminant(mpol),p))+kmax+1;
    rts:=Roots(mpol,pAdicRing(p,henPrec));
    rts:=[Integers()!r[1] : r in rts];
    rts:=[r : r in rts | Valuation((r-theta)*OK,fp) ge kmax+1];
    assert #rts ge 1;
    theta0:=rts[1]; // Thus ord_{fp}(theta - theta_0) >= kmax+1.
    taufacts:=Factorisation(tau*OK);
    taunumfacts:=[fact : fact in taufacts | fact[2] gt 0];
    // The factorization of the numerator ideal of tau*OK.
    if #taunumfacts eq 0 then
	taunum:=1*OK;
    else
	taunum:=&*[fact[1]^(fact[2]) : fact in taunumfacts];
    end if;
    // Now taunum is the numerator ideal of tau.
    Z1:=FreeAbelianGroup(1);
    eta:=hom<Zr->Z1 | [Valuation(delta,fp)*Z1.1 : delta in deltaList]>;
    etaLc:=eta(Lc);
    assert #OrderedGenerators(etaLc) eq 1;
    modulus:=Eltseq(Z1!(etaLc.1))[1];
    class:=Valuation(tau,fp)+Eltseq(Z1!(eta(wc)))[1];
    krange:=[kmin..kmax];
    krange:=[k : k in krange | IsDivisibleBy(k-class,modulus)];
    if #krange eq 0 then
	return [];
    end if;
    H:=Kernel(eta);
    vecs:=[];
    kmin:=Minimum(krange);
    assert kmin ge 0;

    if kmin eq 0 then
	vp:=(-Valuation(tau,fp))*Z1.1;
	tf,wcNew,LcNew:=cosetIntersect(Zr,vp@@eta,H,wc,Lc);
	if tf then
	    invs:=Invariants(LcNew);
	    assert #invs eq rk-1; // Checking that the rank has gone done by 1.
	    vecs:=vecs cat $$(tau,deltaList,Zr,LcNew,wcNew,SLeft,
			      rangeLeft,cBfsq,bigInfs,depth+1);
	end if;
    end if;

    krange:=[k : k in krange | k ge 1];
    krange:=[k : k in krange | (taunum+fa^k) eq ((theta-theta0)*OK+fa^k)];
    for k in krange do
	fb:=fa^k/(taunum+fa^k);
	tau0:=(theta0-theta)/tau;
	for fact in Factorisation(fb) do
	    assert Valuation(tau0*OK,fact[1]) eq 0; // Sanity check!
	end for;
	if (fb eq ideal<OK|1>) then
	    kprime:=0;
	else
	    kprime:=Max([Ceiling(fact[2]/Valuation(p*OK,fact[1])) :
			 fact in Factorisation(fb)]);
	end if;
	pkprime:=p^kprime;
	assert Denominator(pkprime*OK/fb) eq 1;
	if kprime ge 1 then
	    assert Denominator((p^(kprime-1)*OK)/fb) ne 1;
	end if;
	// kprime really is the smallest positive integer such that
	// fb divides p^kprime.
	Zmod:=Integers(pkprime);
	U,eps:=UnitGroup(Zmod);
	gens:=[Integers()!(u@eps) : u in OrderedGenerators(U)];
	if IsOdd(p) and kprime ge 1 then
	    assert #gens eq 1; // Sanity check: there is a primitive root.
	end if;
	alphaList:=[K!tau0] cat [K!g : g in gens] cat
		   [K!delta: delta in deltaList];
	D,alphaList:=multGroup(fb,alphaList);
	// D is (isomorphic to) the group (OK/fb)^* and alphaList is the image
	// in D of the original alphaList defined above, i.e. the images in D
	// of tau0, the generators of Zp'^*, and the delta_i.
	tau0im:=alphaList[1];
	alphaList:=alphaList[2..#alphaList];
	gensim:=alphaList[1..#gens];
	alphaList:=alphaList[(#gens+1)..#alphaList];
	assert #alphaList eq #deltaList;
	G,pi:=quo<D | sub<D | gensim>>;
	tau0im:=pi(tau0im);
	alphaList:=[pi(alpha) : alpha in alphaList];
	assert #alphaList eq r;
	phi:=hom<Zr->G | alphaList>;
	// This the map \phi in the paper.
	if tau0im in Image(phi) then
	    w:=tau0im@@phi;
	    L:=Kernel(phi);
	    tf,wcNew,LcNew:=cosetIntersect(Zr,w,L,wc,Lc);
	    if tf then
		invs:=Invariants(LcNew);
		assert #invs eq rk; // Checking that the rank has is same.
		vp:=(k-Valuation(tau,fp))*Z1.1;
		v:=vp@@eta;
		tf,wcNew,LcNew:=cosetIntersect(Zr,v,H,wcNew,LcNew);
		if tf then
		    invs:=Invariants(LcNew);
		    assert #invs eq rk-1;
		    // Checking that the rank has is reduced by 1.
		    vecs:=vecs cat $$(tau,deltaList,Zr,LcNew,wcNew,SLeft,
				      rangeLeft,cBfsq,bigInfs,depth+1);
		end if;
	    end if;
	end if;
    end for;
    return vecs;
end function;
