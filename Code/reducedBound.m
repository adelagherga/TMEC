/*
reducedBound.m

These functions determine a very large bound for B = max(|b_1|,...,|b_r|), where
a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}. This bound is then
successively reduced using lattice approximation lattices until no additional
improvement is possible.

Authors
    Adela Gherga <adelagherga@gmail.com>
    Samir Siksek <samir.siksek@gmail.com>
Created
    10 August 2022
*/

degL:=function(K)
    /*
      Given a number field K = Q(theta) having degree at least 3, determine the
      smallest possible degree for L = Q(theta_1,theta_2,theta_3), where
      theta_1, theta_2, theta_3 are three distinct conjugates of theta.

      Parameters
          K: FldNum
              A number field K = Q(theta) of degree at least 3.
      Returns
          D: RngIntElt
              The smallest possible degree of L = Q(theta_1,theta_2,theta_3).
   */
    d:=Degree(K);
    assert d ge 3;
    G:=GaloisGroup(K);
    S:=[[i,j,k] : i,j,k in [1..d] | i lt j and j lt k];
    inds:=[Index(G,Stabilizer(G,l)) : l in S];
    D:=Min(inds);
    if IsTransitive(G,3) then
	assert D eq d*(d-1)*(d-2); // Sanity check!
    end if;
    return D;
end function;

initialBound:=function(tau,deltaList,S)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      a bound for B = max(|b_1|,...,|b_r|) using the bounds of Matveev and Yu.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
      Returns
          c17: FldReElt
              A positive real constant used to derive c20.
          c20: FldReElt
              The upper bound for B = max(|b_1|,...,|b_r|) derived using the
              bounds of Matveev and Yu.
   */
    K:=Universe([tau] cat deltaList);
    K:=NumberField(K);
    theta:=K.1;
    OK:=MaximalOrder(K);
    assert &and[ IsPrime(P) : P in S ];
    assert #Set(S) eq #S;
    assert &and[ RamificationIndex(P) eq 1 : P in S];
    assert &and[ InertiaDegree(P) eq 1 : P in S];
    assert &join[Support(delta*OK) : delta in deltaList] eq Set(S);
    u,v:=Signature(K);
    assert #deltaList eq #S+u+v-1; // This is a sanity check.
    // We want deltaList to be a basis for the S-unit group OS^* modulo torsion,
    // and the S-unit group has rank #S+u+v-1.
    T:=[Characteristic(ResidueClassField(P)) : P in S];
    assert #Set(T) eq #T;
    d:=Degree(K);
    D:=degL(K);
    assert IsDivisibleBy(D,d);
    e:=Exp(1);
    pi:=Pi(RealField());
    htheta:=AbsoluteLogarithmicHeight(theta);
    htau:=AbsoluteLogarithmicHeight(tau);
    c7:=Log(2) + 2*htheta + htau;
    c8:=2*D*c7;
    r:=#deltaList;
    hdeltaList:=[AbsoluteLogarithmicHeight(delta) : delta in deltaList];
    hstarprod:=&*[(4*hs^2 + pi^2/D^2)^(1/2) : hs in hdeltaList];
    hstarprod:=hstarprod*( 4*c7^2 + pi^2/D^2 )^(1/2);
    c9:=6*30^(r+5)*(r+2)^(5.5)*D^(r+3)*Log(e*D)*hstarprod;
    c10:=c8 + c9*Log(e*(r+1));
    hdaggerprod:=&*[Max(2*hs,1/(16*e^2*D^2)) : hs in hdeltaList];
    hdaggerprod:=hdaggerprod*Max(2*c7,1/(16*e^2*D^2));
    relD:=D div d;
    c11:=0;
    for p in T do
	// Recall T is the set of rational primes below the primes in S.
	for uu in [1..relD] do
	    for vv in [1..Floor(relD/uu)] do
		// Looking at all pairs u, v such that uv <= D/d.
		c11:= Max(c11,uu^(r+1)*p^vv/(vv*Log(p)));
	    end for;
	end for;
    end for;
    c12:=(16*e*D)^(2*(r+1)+2)*(r+1)^(5/2)*Log(2*(r+1)*D)*Log(2*D);
    c12:=c12*c11*hdaggerprod ;
    embeds:=InfinitePlaces(K);
    assert d eq u+2*v;
    assert #embeds eq u+v;
    c13:=(u+v+2*#S)/d;
    c14:=2*htau + c13*Max(c8,c10);
    c15:=c13*Max(c9,c12);
    c16:=c14 + (3/2)*Log(2) + 2*htheta + htau;

    c17list:=[];
    for i in [1..#S + #embeds] do
	if i le #S then
	    shortS:=[S[j] : j in [1..#S] | j ne i];
	    shortembeds:=embeds;
	else
	    shortS:=S;
 	    shortembeds:=[embeds[j] : j in [1..#embeds] | j ne (i - #S)];
	end if;
	M:=[];
	for delta in deltaList do
	    M1:=[LocalDegree(j)*Log(AbsoluteValue(Evaluate(delta,j)))
		 : j in shortembeds];
	    M2:=[Log(Norm(P)^(-Valuation(delta*OK,P))) : P in shortS];
	    Append(~M, M1 cat M2);
	end for;
	M:=Transpose(Matrix(RealField(),M));
	assert NumberOfRows(M) eq r;
	assert NumberOfColumns(M) eq r;
	Minv:=M^(-1);
	Minv:=Rows(Minv);
	Minv:=[Eltseq(r) : r in Minv];
	Minv:=&cat(Minv);
	assert #Minv eq r^2;
	Minv:=[AbsoluteValue(m) : m in Minv];
	c17:=Max(Minv);
	Append(~c17list,c17);
    end for;
    c17:=Min(c17list);
    c18:=2*d*c17*c16;
    c19:=2*d*c17*c15;
    c20:=2*c18 + Max(4*e^2,2*c19*Log(c19));
    return c17,c20;
end function;

boundConstants:=function(K,theta,tau)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      the constants c22, c23 needed to reduce the bound for
      B = max(|b_1|,...,|b_r|), derived using the bounds of Matveev and Yu.

      Parameters
          K: FldNum
              The number field K = Q(theta) of degree at least 3.
          theta: RngOrdElt
              A generator of K.
          tau: FldNumElt
      Returns
          c22: FldReElt
          c23: FldReElt
   */
    d:=Degree(K);
    assert d ge 3;
    u,v:=Signature(K);
    pls:=InfinitePlaces(K);
    realPls:=[pl : pl in pls | LocalDegree(pl) eq 1];
    cmxPls:=[pl : pl in pls | LocalDegree(pl) eq 2];
    assert #realPls eq u;
    assert #cmxPls eq v;
    c22:=0;
    for pl in cmxPls do
	c22:=c22+2*Log(Max(1,AbsoluteValue(Evaluate(tau,pl)/
					   Imaginary(Evaluate(theta,pl)))));
    end for;
    realtheta:=[Evaluate(theta,pl) : pl in realPls];
    realtau:=[Evaluate(tau,pl) : pl in realPls];
    if u le 1 then
	c23:=1; // We do not use c23 when #realPls is 0.
    else
	c23:=Min([(AbsoluteValue(realtheta[i]-realtheta[j]))/
		  (AbsoluteValue(realtau[i])+AbsoluteValue(realtau[j]))
		  : i,j in [1..u] | i lt j]);
    end if;
    return c22,c23;
end function;

nonUnitExps:=function(K,S,deltaList,vecB)
    /*
      Compute the inverse of the matrix M_0 and initial vector u'' required to
      reduce the exponents of the non-unit delta_i.

      Parameters
          K: FldNum
              The number field K = Q(theta) of degree at least 3.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_r|.
      Returns
          absM0inv: AlgMatElt
              The inverse of the matrix M_0 = (Log|delta_j|_v), where each
	      entry is taken in absolute value.
          vecUpp: ModMatFldElt
              The most up-to-date vector u'' = (Log|eps|_v).
   */
    d:=Degree(K);
    OK:=MaximalOrder(K);
    r:=#deltaList;
    u,v:=Signature(K);
    s:=#S;
    assert d eq u+2*v;
    assert s eq (r-(u+v-1));
    assert &and[Abs(Norm(delta)) eq 1 : delta in deltaList[1..u+v-1]];
    assert &and[Abs(Norm(delta)) ne 1 : delta in deltaList[u+v..r]];
    if (s eq 0) then
	return [],[];
    else
	M0:=Matrix([[Log(Norm(P)^(-Valuation(deltaList[i]*OK,P)))
		     : i in [u+v..r]] : P in S]);
	assert NumberOfRows(M0) eq s;
	assert NumberOfColumns(M0) eq s;
	M0inv:=1/Determinant(M0)*Adjoint(M0);
	// Use the adjoint matrix of M0 to compute the matrix inverse
	// in order to avoid precision errors.
	absM0:=Matrix([[Abs(M0[i,j]) : j in [1..s]] : i in [1..s]]);
	absM0inv:= Matrix([[Abs(M0inv[i,j]) : j in [1..s]] : i in [1..s]]);
	vecUpp:=absM0*Matrix(RealField(),s,1,vecB[u+v..r]);
    end if;
    return absM0inv,vecUpp;
end function;

valuationBound:=function(tau,deltaList,fp,cB2)
    /*
      Determine an upper bound for the valuation of ord_{fp}(a_0 X - theta Y),
      as per Proposition 7.2.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          fp: RngOrdIdl
              A prime ideal of K, contained in S.
          cB2: FldReElt
              The smallest known upper bound for the L^2 norm of the vector
	      (b_1,...,b_r).
      Returns
          km1: RngIntElt
              An upper bound for ord_{fp}(a_0 X - theta Y).
   */
    K:=NumberField(Universe([tau] cat deltaList));
    theta:=K.1;
    OK:=MaximalOrder(K);
    d:=Degree(K);
    assert d ge 3;
    r:=#deltaList;
    p:=Characteristic(ResidueClassField(fp));
    assert Degree(fp) eq 1;
    assert Valuation(p*OK,fp) eq 1;
    fa:=p*OK/fp;
    for fact in Factorisation(fa) do
	assert &and[Valuation(delta*OK,fact[1]) eq 0 : delta in deltaList];
    end for;
    if cB2 lt 1 then
        return Valuation(tau*OK,fp);
    end if;
    k:=Floor(r*Log(cB2)/((d-2)*Log(p)));
    repeat
	k:=k+1;
	mpol:=MinimalPolynomial(theta);
	lcf:=LeadingCoefficient(mpol);
	assert lcf eq 1; // Surely theta is an algebraic integer!
	hprec:=2*Valuation(Discriminant(mpol),p)+k+1;
	rts:=Roots(mpol,pAdicRing(p,hprec));
	rts:=[Integers()!r[1] : r in rts];
	rts:=[r : r in rts | Valuation((r-theta)*OK,fp) ge k];
	assert #rts ge 1;
	theta0:=rts[1]; // Thus ord_{fp}(theta-theta_0) \ge k.
	taufacts:=Factorisation(tau*OK);
	taunumfacts:=[fact : fact in taufacts | fact[2] gt 0];
	// The factorization of the numerator ideal of tau*OK.
	if #taunumfacts eq 0 then
	    taunum:=1*OK;
	else
	    taunum:=&*[fact[1]^(fact[2]) : fact in taunumfacts];
	end if;
	// Now taunum is the numerator ideal of tau.
	if (taunum+fa^k) ne ((theta-theta0)*OK+fa^k) then
	    // Condition (i) of the proposition.
	    return k-1;
	end if;
	fb:=fa^k/(taunum+fa^k);
	tau0:=(theta0-theta)/tau;
	for fact in Factorisation(fb) do
	    assert Valuation(tau0*OK,fact[1]) eq 0; // Sanity check!
	end for;
	if (fb eq ideal<OK|1>) then
	    kprime:= 0;
	else
	    kprime:=Max([Ceiling(fact[2]/Valuation(p*OK,fact[1])) :
			 fact in Factorisation(fb)]);
	end if;
	pkprime:=p^kprime;
	assert Denominator(pkprime*OK/fb) eq 1;
	if kprime ge 1 then
	    assert Denominator((p^(kprime-1)*OK)/fb) ne 1;
	end if;
	// kprime really is the smallest positive integer such
	// that fb divides p^kprime.
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
	// of tau0, the generators of (\Z/p^\prime)^*, and the \delta_i.
	tau0im:=alphaList[1];
	alphaList:=alphaList[2..#alphaList];
	gensim:=alphaList[1..#gens];
	alphaList:=alphaList[(#gens+1)..#alphaList];
	assert #alphaList eq #deltaList;
	G,pi:=quo<D | sub<D | gensim>>;
	tau0im:=pi(tau0im);
	alphaList:=[pi(alpha) : alpha in alphaList];
	assert #alphaList eq r;
	Zr:=FreeAbelianGroup(r);
	phi:=hom<Zr->G | alphaList>;
	// This the map \phi in the paper.
	if tau0im in Image(phi) eq false then
	    // Condition (ii) of the proposition.
	    return k-1;
	end if;
	w:=tau0im@@phi;
	L:=Kernel(phi);
	ZZr:=StandardLattice(r);
	LL:=sub<ZZr | [ ZZr!(Eltseq((Zr!l))) : l in OrderedGenerators(L)]>;
	ww:=ZZr!Eltseq(Zr!w);
	// ClosestVectors is error-prone in Magma V2.26-10.
	// If there are no vectors in P, then there are no vectors v in w+L
	// such that |v+w| < cB2.
	P:=CloseVectorsProcess(LL,-ww,Integers()!Floor(cB2^2));
	// until DLwsq gt cB2^2; // We keep increasing k until condtion (iii)
    until (IsEmpty(P) eq true);
    // Condition (iii) of the proposition.
    // We keep increasing k until condition (i), (ii), or (iii) are satisfied.
    return k-1;
end function;

c21Func:=function(tau,deltaList,S,absM0inv,vecUpp,vecB)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      the constant c21 needed to reduce the bound for
      B = max(|b_1|,...,|b_r|), derived using the bounds of Matveev and Yu.
      This function also computes lower and upper bounds for the fp-valuation
      of eps = delta_1^b_1 ... delta_r^b_r, as well as updates the vector u''
      and the list of smallest known bounds for |b_1|,...,|b_r|.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
          absM0inv: AlgMatElt
              The inverse of the matrix M_0 = (Log|delta_j|_v), where each
	      entry is taken in absolute value.
          vecUpp: ModMatFldElt
              The most up-to-date vector u'' = (Log|eps|_v).
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_r|.
      Returns
          c21: FldReElt
	  expSbds: SeqEnum
	  vecUpp: ModMatFldElt
              The most up-to-date vector u'' = (Log|eps|_v).
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_r|.
   */
    K:=NumberField(Universe([tau] cat deltaList));
    theta:=K.1;
    OK:=MaximalOrder(K);
    d:=Degree(K);
    assert d ge 3;
    r:=#deltaList;
    s:=#S;
    u,v:=Signature(K);
    assert d eq u+2*v;
    assert s eq (r-(u+v-1));
    assert &and[Abs(Norm(delta)) eq 1 : delta in deltaList[1..u+v-1]];
    assert &and[Abs(Norm(delta)) ne 1 : delta in deltaList[u+v..r]];

    expSbds:=[];
    c21:=0;
    cB2:=Sqrt(&+[i^2 : i in vecB]); // The L^2 norm of the vector b.
    for i in [1..s] do
	fp:=S[i];
	km1:=valuationBound(tau,deltaList,fp,cB2); // This is k-1.
	kfpp:=Valuation(tau*OK,fp);
	kfppp:=km1-kfpp;
	Append(~expSbds,<-kfpp,kfppp>);
	vecUpp[i,1]:=Abs(Log(Norm(fp)))*Max([Abs(kfpp),Abs(kfppp)]);
	vecBnew:=Eltseq(absM0inv*vecUpp);
	for j in [1..s] do
	    if vecBnew[j] lt vecB[u+v-1+j] then
		vecB[r-s+j]:=Round(vecBnew[j]); // Update the vector b.
	    end if;
	end for;
	cB2:=Sqrt(&+[i^2 : i in vecB]); // Update cB2.
	// <-kfpp,kfppp> are lower and upper bounds for the fp-valuation
	// of eps = delta_1^b_1 ... delta_r^b_r.
	c21:=c21+Max(0,kfppp)*Log(Norm(fp));
    end for;
    return c21,expSbds,vecUpp,vecB;
end function;

totallyComplexRed:=function(tau,deltaList,S,consts : verb:=false)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      a reduced bound for B = max(|b_1|,...,|b_r|) in the totally complex case.
      This reduction is done succesively, and in fact yields a bound for each
      |b_1|,...,|b_r|.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
          consts: SeqEnum
              The list of positive real constants c17,c20,c22,c23,c26.
          verb: BoolElt
              A true/false value. If set to true, this function returns status
	      updates as it proceeds.
      Returns
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_r|.
          expSbds: SeqEnum
              The mode up-to-date lower and upper bounds for the fp-valuation
              of eps.
   */
    SetVerbose("User1",verb);
    K:=Universe([tau] cat deltaList);
    K:=NumberField(K);
    theta:=K.1;
    OK:=MaximalOrder(K);
    s:=#S;
    r:=#deltaList;
    d:=Degree(K);
    assert d ge 3;
    u,v:=Signature(K);
    pls:=InfinitePlaces(K);
    realPls:=[pl : pl in pls | LocalDegree(pl) eq 1];
    cmxPls:=[pl : pl in pls | LocalDegree(pl) eq 2];
    assert #realPls eq u;
    assert #cmxPls eq v;
    assert d eq u+2*v;
    assert u eq 0;

    c17,c20,c22,c23,c26:=Explode(consts);
    vprintf User1: "We're in the totally complex case.\n";
    vprintf User1: "The initial bound is %o.\n",c20;
    cBinf:=c20;
    vecB:=[cBinf : i in [1..r]];
    absM0inv,vecUpp:=nonUnitExps(K,S,deltaList,vecB);

    finished:=false;
    repeat
	c21,expSbds,vecUpp,vecB:=c21Func(tau,deltaList,S,absM0inv,vecUpp,vecB);
	if (s gt 0) then
	    vprintf User1: "The exponent bounds are %o.\n",expSbds;
	end if;
	cBinfNew:=2*c17*(c21+c22);
	if cBinfNew lt cBinf then
	    vprintf User1: "The reduction process gives a new bound of %o.\n",
			   Floor(cBinfNew);
	    cBinf:=Floor(cBinfNew);
	    for i in [1..r] do
		if cBinf lt vecB[i] then
		    vecB[i]:=cBinf;
		end if;
	    end for;
	else
	    vprintf User1: "The reduction process gives a worse bound of %o.\n",
			   Floor(cBinfNew);
	    finished:=true;
	end if;
    until finished;
    return vecB,expSbds;
end function;

boundConstantsDivc25:=function(K,theta,tau,c23,sigma)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      the modified list of constants c27, c28, and c29, all divided by c25,
      as well as the real embedding sigma_2, needed to reduce the bound for
      B = max(|b_1|,...,|b_r|).

      Parameters
          K: FldNum
              The number field K = Q(theta) of degree at least 3.
          theta: RngOrdElt
              A generator of K.
          tau: FldNumElt
          c23: FldReEld
          sigma: PlcNumElt
              A fixed real embedding of K.
      Returns
	  constsDivc25: Tup
              The list (c27Divc25,c28Divc25,c29Divc25) of constants which,
              when multiplied by c25, yield c27, c28, and the list of constants
              c29.
          sigma2: PlcNumElt
              The fixed real embedding sigma_2 != sigma of K.
   */
    d:=Degree(K);
    assert d ge 3;
    u,v:=Signature(K);
    pls:=InfinitePlaces(K);
    realPls:=[pl : pl in pls | LocalDegree(pl) eq 1];
    cmxPls:=[pl : pl in pls | LocalDegree(pl) eq 2];
    assert #realPls eq u;
    assert #cmxPls eq v;

    otherRealPls:=[pl : pl in realPls | pl ne sigma];
    otherPls:=otherRealPls cat cmxPls;
    assert #otherPls eq (u+v-1);
    sigma2:=otherPls[1];
    tau1:=AbsoluteValue(Evaluate(tau,sigma));
    if sigma2 in realPls then
	otherRealPls:=[pl : pl in otherRealPls | pl ne sigma2];
	assert #otherRealPls eq (u-2);
	tau2:=AbsoluteValue(Evaluate(tau,sigma2));
	tauconjAbs:=[AbsoluteValue(Evaluate(tau,pl)) : pl in otherRealPls];
	thetaconjIms:=[Imaginary(Evaluate(theta,pl)) : pl in cmxPls];
	c27divc25:=tau1/(Min(tauconjAbs cat [tau2])*c23);
	if IsEmpty(thetaconjIms) then
	    c28divc25:=1;
	else
	    c28divc25:=tau1/(Min(thetaconjIms));
	end if;
	c29divc25:=[0] cat [tau1/(tau2*c23)] cat
		   [tau1/tau*c23 : tau in tauconjAbs] cat
		   [tau1/thetaj : thetaj in thetaconjIms];
    else
	assert IsEmpty(otherRealPls);
	otherCmxPls:=[pl : pl in cmxPls | pl ne sigma2];
	assert #otherCmxPls eq (v-1);
	theta2:=Imaginary(Evaluate(theta,sigma2));
	thetaconjIms:=[Imaginary(Evaluate(theta,pl)) : pl in otherCmxPls];
	c27divc25:=1;
	c28divc25:=tau1/(Min(thetaconjIms cat [theta2]));
	c29divc25:=[0] cat [tau1/theta2] cat
		   [tau1/thetaj : thetaj in thetaconjIms];
    end if;
    assert #c29divc25 eq (u+v);
    return <c27divc25,c28divc25,c29divc25>,sigma2;
end function;

approxLattice:=function(tau,deltaList,S,sigma,sigma2,C)
    /*
      Given a real embedding sigma of K such that |sigma(eps)| is minimal, and
      a positive integer C, determine the matrix M and vector w used to reduce
      the bound for B = max(|b_1|,...,|b_r|). Here, the lattice L (as in
      Proposition 9.1) is spanned by the rows of M.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
          sigma: PlcNumElt
              A fixed real embedding of K.
          sigma2: PlcNumElt
              The fixed real embedding sigma_2 != sigma of K.
          C: RngIntElt
      Returns
          LatMat: AlgMatElt
              The matrix M whose rows define the lattice L.
          ww: SeqEnum
              The vector w.
   */
    K:=Universe([tau] cat deltaList);
    K:=NumberField(K);
    theta:=K.1;
    d:=Degree(K);
    u,v:=Signature(K);
    pls:=InfinitePlaces(K);
    realPls:=[pl : pl in pls | LocalDegree(pl) eq 1];
    cmxPls:=[pl : pl in pls | LocalDegree(pl) eq 2];
    assert sigma in realPls;
    assert LocalDegree(sigma) eq 1;
    assert #realPls eq u;
    assert #cmxPls eq v;
    s:=#S;
    r:=#deltaList;
    w:=u+v-2;
    n:=d+s-1;
    assert d eq u+2*v;
    assert r eq (u+v+s-1);
    assert n eq r+v;
    prec:=2*Ceiling(Log(C)/Log(10))+2000;
    theta1:=Evaluate(theta,sigma : Precision:=prec);
    pi:=Pi(RealField(prec));
    saj:=C*pi;
    _,xp:=MantissaExponent(saj);
    assert xp lt -100;
    // This a test to make sure that we have enough precision to compute
    // the nearest integer correctly.
    Cpi:=Round(saj);
    matblock1:=ZeroMatrix(Integers(),w,r+v);
    matblock2:=HorizontalJoin(IdentityMatrix(Integers(),s+1),
			      ZeroMatrix(Integers(),s+1,d-2));
    matblock3:=ZeroMatrix(Integers(),v,r+v);
    LatMat:=VerticalJoin(<matblock1,matblock2,matblock3>);
    assert NumberOfColumns(LatMat) eq NumberOfRows(LatMat);
    assert NumberOfColumns(LatMat) eq ((s+1)+(d-2));

    for i in [1..v] do
	LatMat[r+i,s+1+w+i]:=Cpi;
    end for;
    ww:=[0 : j in [1..(s+1)]];
    otherRealPls:=[pl : pl in realPls | pl ne sigma];
    otherPls:=otherRealPls cat cmxPls;
    assert #otherPls eq u+v-1;
    otherPls:=[pl : pl in otherPls | pl ne sigma2];
    assert #otherPls eq w;
    theta2:=Evaluate(theta,sigma2 : Precision:=prec);
    tau2:=Evaluate(tau,sigma2 : Precision:=prec);
    for j in [1..w] do
	sigma3:=otherPls[j];
	theta3:=Evaluate(theta,sigma3 : Precision:=prec);
	tau3:=Evaluate(tau,sigma3 : Precision:=prec);
	betaj:=Log(AbsoluteValue((theta2-theta1)*tau3))
	       - Log(AbsoluteValue((theta3-theta1)*tau2));
	saj:=C*betaj;
	_,xp:=MantissaExponent(saj);
	assert xp lt -100;
	saj:=Round(saj);
	ww:=ww cat [saj];
	for i in [1..r] do
	    delta:=deltaList[i];
	    ld2:=Log(AbsoluteValue(Evaluate(delta,sigma2 : Precision:=prec)));
	    ld3:=Log(AbsoluteValue(Evaluate(delta,sigma3 : Precision:=prec)));
	    saj:=C*(ld3-ld2);
	    _,xp:=MantissaExponent(saj);
	    assert xp lt -100;
	    saj:=Round(saj);
	    LatMat[i,s+1+j]:=saj;
	end for;
    end for;
    for j in [1..v] do
	CC:=ComplexField(prec);
	sigma2:=cmxPls[j];
	theta2:=CC!Evaluate(theta,sigma2 : Precision:=prec);
	tau2:=CC!Evaluate(tau,sigma2 : Precision:=prec);
	beta:=Imaginary(Log((theta1-theta2)/tau2));
	saj:=C*beta;
	_,xp:=MantissaExponent(saj);
	assert xp lt -100;
	saj:=Round(saj);
	ww:=ww cat [saj];
	for i in [1..r] do
	    delta:=deltaList[i];
	    alphai:=-Imaginary(Log(CC!Evaluate(delta,sigma2
					       : Precision:=prec)));
	    saj:=C*alphai;
	    _,xp:=MantissaExponent(saj);
	    assert xp lt -100;
	    saj:=Round(saj);
	    LatMat[i,s+1+w+j]:=saj;
	end for;
    end for;
    assert #ww eq n;
    return LatMat,ww;
end function;

distanceLBsq:=function(LL,ww)
   /*
      Determines the shortest distance from a vector in the lattice L to the
      vector w.

      Parameters
	  LL: Lat
	      The lattice L of Z^n given by the rows of the matrix M.
          ww: LatElt
              The vector w in Z^n.
      Returns
          DLwsq: RngIntElt
              The squared shortest distance from a vector in L to w, D(L,w)^2.
   */
    ww:=AmbientSpace(Parent(ww))!ww;
    // Changing the coefficient ring of ww to the rationals.
    n:=Rank(LL);
    assert n eq Degree(LL);
    LL:=LLL(LL);
    B:=Basis(LL);
    b1:=B[1];
    if ww in LL then
	return (2^(1-n))*Norm(b1);
    end if;
    B:=[Eltseq(b) : b in B];
    B:=[[Rationals()!a : a in b] : b in B];
    B:=Matrix(B);
    Binv:=B^-1;
    sigma:=Eltseq(ww*B^-1);
    sigma:=[AbsoluteValue(t-Round(t)) : t in sigma];
    sigma:=[t : t in sigma | t ne 0];
    sigmai0:=sigma[#sigma];
    return (2^(1-n))*sigmai0^2*Norm(b1);
end function;

distanceLBsq2:=function(LL,ww,cB5sq)
    /*
      Determines whether the shortest distance from a vector in the lattice L
      to the vector w is less than cB5, and outputs the shortest distance if
      true.

      Parameters
	  LL: Lat
	      The lattice L of Z^n given by the rows of the matrix M.
          ww: LatElt
              The vector w in Z^n.
          cB5sq: FldReElt
              The bound cB5^2.
      Returns
          shortvecs: SeqEnum
              A list of vectors v of the translated lattice w + L having squared
	      L^2 norm <= cB5^2.
          DLwsq: RngIntElt
              The squared shortest distance from a vector in L to w, D(L,w)^2.
   */
    ww:=AmbientSpace(Parent(ww))!ww;
    // Changing the coefficient ring of ww to the rationals.
    n:=Rank(LL);
    assert n eq Degree(LL);
    l:=Integers()!Floor(cB5sq);
    u:=5*l;
    if ww in LL then
	P:=ShortVectorsProcess(LL,l,u);
	while IsEmpty(P) do
	    u:=5*u;
	    P:=ShortVectorsProcess(LL,l,u);
	end while;
	// Find at least 1 vector v in L such that cB5^2 <= Norm(v+w) <= mult*cB5^2.
	P:=ShortVectorsProcess(LL,u);
	// Find all vectors v in L such that Norm(v+w) <= mult*cB5^2.
	vecs:=[];
	while (IsEmpty(P) eq false) do
	    Append(~vecs,NextVector(P));
	end while;
	shortvecs:=[xx : xx in vecs | Norm(xx) le l];
	DLwsq:=Min([Norm(xx) : xx in vecs | xx notin shortvecs]);
	shortvecs:=[ChangeRing(xx-ww,Rationals()) : xx in shortvecs];
	assert DLwsq gt cB5sq;
	return shortvecs,DLwsq;
    end if;
    P:=CloseVectorsProcess(LL,-ww,l,u);
    while IsEmpty(P) do
	u:=5*u;
	P:=CloseVectorsProcess(LL,-ww,l,u);
    end while;
    // Find at least 1 vector v in L such that cB5^2 <= Norm(v+w) <= mult*cB5^2.
    P:=CloseVectorsProcess(LL,-ww,u);
    // Find all vectors v in L such that Norm(v+w) <= mult*cB5^2.
    vecs:=[];
    while (IsEmpty(P) eq false) do
	Append(~vecs,NextVector(P));
    end while;
    shortvecs:=[vv : vv in vecs | Norm(vv+ww) le l];
    DLwsq:=Min([Norm(vv+ww) : vv in vecs | vv notin shortvecs]);
    shortvecs:=[ChangeRing(vv,Rationals()) : vv in shortvecs];
    assert DLwsq gt cB5sq;
    return shortvecs,DLwsq;
end function;

fixedRealEmbeddingRed:=function(tau,deltaList,S,consts,sigma : verb:=false)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      a reduced bound for B = max(|b_1|,...,|b_r|) under the assumption that
      sigma is a real embedding of K such that |sigma(eps)| is minimal. This
      reduction is done succesively, and in fact yields a bound for each
      |b_1|,...,|b_r|.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
          consts: SeqEnum
              The list of positive real constants c17,c20,c22,c23,c26.
          sigma: PlcNumElt
              A fixed real embedding of K.
          verb: BoolElt
              A true/false value. If set to true, this function returns status
	      updates as it proceeds.
      Returns
          vecs: SetEnum
              A list of vectors v of the translated lattice w + L having squared
	      L^2 norm <= cB5^2.
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_r|.
          expSbds: SeqEnum
              The mode up-to-date lower and upper bounds for the fp-valuation
              of eps.
   */
    SetVerbose("User1",verb);
    K:=Universe([tau] cat deltaList);
    K:=NumberField(K);
    theta:=K.1;
    OK:=MaximalOrder(K);
    d:=Degree(K);
    r:=#deltaList;
    u,v:=Signature(K);
    w:=u+v-2;
    s:=#S;
    n:=r+v;
    assert d eq u+2*v;
    assert s eq (r-(u+v-1));
    assert &and[Abs(Norm(delta)) eq 1 : delta in deltaList[1..u+v-1]];
    assert &and[Abs(Norm(delta)) ne 1 : delta in deltaList[u+v..r]];

    c17,c20,c22,c23,c26:=Explode(consts);
    constsDivc25,sigma2:=boundConstantsDivc25(K,theta,tau,c23,sigma);
    c27divc25,c28divc25,c29divc25:=Explode(constsDivc25);
    vprintf User1: "The initial bound is %o.\n",c20;
    cBinf:=c20;
    vecB:=[cBinf : i in [1..r]];
    absMinv,vecUpp:=nonUnitExps(K,S,deltaList,vecB);
    vecs:={};

    finished:=false;
    repeat
	c21,expSbds,vecUpp,vecB:=c21Func(tau,deltaList,S,absMinv,vecUpp,vecB);
	if (s gt 0) then
	    vprintf User1: "The exponent bounds are %o.\n",expSbds;
	end if;
	c24:=c21+c22+(u-1)*Log(Max(1,1/c23));
	c25:=Exp(c24);
	c27:=c27divc25*c25;
	c28:=c28divc25*c25;
	c29:=[c29divc25[j]*c25 : j in [1..(u+v)]];
	if u eq 1 then
	    c30:=Max([2*c17*c24,Log(2*c28)/c26]);
	elif v eq 0 then
	    c30:=Max([2*c17*c24,Log(2*c27)/c26]);
	else
	    c30:=Max([2*c17*c24,Log(2*c27)/c26,Log(2*c28)/c26]);
	end if;

	cBinf:=Max(vecB);
	cB1:=&+[Abs(i) : i in vecB];
	cB2:=Sqrt(&+[i^2 : i in vecB]);
	cA1:=(1+cB1)/2;
	cA2:=(1+cB1) + 1/(2*Pi(RealField()));
	cB3:=0;
	cB4:=0;
	for j in [1..w] do
	    cB3:=cB3 + (c29[2] + c29[j+2])^2;
	    cB4:=cB4 + cA1*(c29[2] + c29[j+2]);
	end for;
	for j in [1..v] do
	    cB3:=cB3 + (c29[u+j])^2;
	    cB4:=cB4 + cA2*(c29[u+j]);
	end for;
	cB5:=(cB2^2 - w*cBinf^2 + w*cA1^2 + v*cA2^2)^(1/2);
	C:=Ceiling(cB5^(n/(d-2)));
	// We follow the procedure in the section "Reduction of Bounds"
	// to repeatedly reduce the initial bound.
	iter:=1;
	repeat
	    ZZn:=StandardLattice(n);
	    M,ww:=approxLattice(tau,deltaList,S,sigma,sigma2,C);
	    while (Rank(M) ne n) do
		vprintf User1: "Increasing C.\n";
		C:=10*C;
		M,ww:=approxLattice(tau,deltaList,S,sigma,sigma2,C);
	    end while;
	    ww:=ZZn!ww;
	    LL:=Lattice(M);
	    if (ww in LL) then
		vv:=ChangeRing(-ww,Rationals());
		vecs:=vecs join
			   {(vv*ChangeRing(M,Rationals())^(-1))[1]};
	    end if;
	    DLwsq:=distanceLBsq(LL,-ww);
	    // This is D(L,w)^2 in the notation of Proposition 9.1.
	    tf:=DLwsq gt cB5^2;
	    if (tf eq false) then
		if (iter eq 20) then
		    shortvecs,DLwsq:=distanceLBsq2(LL,ww,cB5^2);
		    vecs:=vecs join {(vv*ChangeRing(M,Rationals())^(-1))[1]
				     : vv in shortvecs};
		    tf:=DLwsq gt cB5^2;
		else
		    vprintf User1: "Increasing C.\n";
		    C:=10*C;
		    iter:=iter+1;
		end if;
	    end if;
	until tf;
	denom:=(cB3*(DLwsq-cB5^2) + cB4^2)^(1/2) - cB4;
	cBinfNew:=(1/c26)*Log((2*C*cB3)/denom);
	cBinfNew:=Max(c30,cBinfNew);
	if cBinfNew lt 10^10 then
	    cBinfNew:=Floor(cBinfNew);
	end if;
	if cBinfNew lt cBinf then
	    vprintf User1: "The new bound is %o.\n",cBinfNew;
	    cBinf:=cBinfNew;
	    for i in [1..r] do
		if cBinf lt vecB[i] then
		    vecB[i]:= cBinf;
		end if;
	    end for;
	else
	    finished:=true;
	end if;
	if (cBinf eq 0) then
	    finished:=true;
	end if;
    until finished;
    return vecs,vecB,expSbds;
end function;

reducedBound:=function(tau,deltaList : verb:=false)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      the final bound for B = max(|b_1|,...,|b_r|) after successively reducing
      the initial, very large bound.

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          verb: BoolElt
              A true/false value. If set to true, this function returns status
	      updates as it proceeds.
      Returns
          vecs: SeqEnum
              A list of vectors v of the translated lattice w + L having squared
	      L^2 norm <= cB5^2.
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_r|.
          S: SeqEnum
              A list of prime ideals fp_1,...,fp_s of K.
          expSbds: SeqEnum
              The mode up-to-date lower and upper bounds for the fp-valuation
              of eps.
   */
    SetVerbose("User1",verb);
    K:=Universe([tau] cat deltaList);
    K:=NumberField(K);
    theta:=K.1;
    OK:=MaximalOrder(K);
    S:=&join[Support(delta*OK) : delta in deltaList];
    S:=SetToSequence(S);
    fn:=func<P1,P2 | Norm(P2)-Norm(P1) >;
    // Sort S so that the prime ideals with largest norm come first.
    Sort(~S,fn);
    s:=#S;
    assert &and[fact[1] in S : fact in Factorisation(tau*OK) | fact[2] lt 0];
    // This is a sanity check. According to the paper, the denominator ideal
    // of tau is supported on S.
    c17,c20:=initialBound(tau,deltaList,S);
    c22,c23:=boundConstants(K,theta,tau);
    c26:=1/(2*c17);
    consts:=[c17,c20,c22,c23,c26];
    r:=#deltaList;
    d:=Degree(K);
    assert d ge 3;
    u,v:=Signature(K);
    pls:=InfinitePlaces(K);
    realPls:=[pl : pl in pls | LocalDegree(pl) eq 1];
    cmxPls:=[pl : pl in pls | LocalDegree(pl) eq 2];
    assert #realPls eq u;
    assert #cmxPls eq v;
    assert d eq u+2*v;
    ZZn:=StandardLattice(r+v);
    vecs:={};

    if u eq 0 then
	vecB,expSbds:=totallyComplexRed(tau,deltaList,S,consts : verb:=verb);
	// The totally complex case.
	for i in [1..s] do
	    expSbds[i,1]:=expSbds[i,1]+Valuation(tau,S[i]);
	    expSbds[i,2]:=expSbds[i,2]+Valuation(tau,S[i]);
	end for;
	return [Integers()!b : b in vecB],S,expSbds;
    end if;
    vprintf User1: "We're carrying out the reduction process for each real ";
    vprintf User1: "embedding separately.\n";
    vecBList:=[];
    expSbdsList:=[];
    for i in [1..u] do
	vprintf User1: "++++++++++++++++++++++++++\n";
	if (i mod 10) eq 1 then
	    vprintf User1: "Dealing with the %o-st real embedding.\n",i;
	elif (i mod 10) eq 2 then
	    vprintf User1: "Dealing with the %o-nd real embedding.\n",i;
	elif (i mod 10) eq 3 then
	    vprintf User1: "Dealing with the %o-rd real embedding.\n",i;
	else
	    vprintf User1: "Dealing with the %o-th real embedding.\n",i;
	end if;
	sigma:=realPls[i];
	vecsig,vecBsig,expSbds:=
	    fixedRealEmbeddingRed(tau,deltaList,S,consts,sigma : verb:=verb);
	Append(~vecBList,vecBsig);
	Append(~expSbdsList,expSbds);
	vecs:=vecs join vecsig;
    end for;

    assert (#vecBList eq u) and (#expSbdsList eq u);
    expSbds:=expSbdsList[1];
    for expSbds2 in expSbdsList[2..u] do
	assert #expSbds2 eq s;
	for i in [1..s] do
	    expSbds[i,1]:=Min(expSbds[i,1], expSbds2[i,1]);
	    expSbds[i,2]:=Max(expSbds[i,2], expSbds2[i,2]);
	end for;
    end for;
    for i in [1..s] do
	expSbds[i,1]:=expSbds[i,1]+Valuation(tau,S[i]);
	expSbds[i,2]:=expSbds[i,2]+Valuation(tau,S[i]);
    end for;
    vecB:=vecBList[1];
    for vecB2 in vecBList[2..u] do
	assert #vecB eq r;
	for i in [1..r] do
	    vecB[i]:=Max(vecB[i],vecB2[i]);
	end for;
    end for;
    vecs:=[Eltseq(ZZn!vv) : vv in vecs];
    return vecs,[Integers()!b : b in vecB],S,expSbds;
end function;
