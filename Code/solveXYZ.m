/*
solveXYZ.m

These functions solve X+Y=Z for non-zero, positive coprime integers X, Y, and Z
whose prime divisors all lie in {p_1,...,p_s}.

Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    17 November 2022
*/

load "./parseIO.m";
SetAutoColumns(false);
SetColumns(235);

pAdicApprox:=function(x,p,mu)
    /*
      Given a p-adic number x, determine the positive rational integer x^{mu}
      such that x^{mu} < p^mu, with ord_p(x-x^{mu})>=mu.

      Parameters
          x: FldPadElt
	  p: RngIntElt
          mu: RngIntElt

      Returns
          xmu: RngIntEl
	      The rational integer congruent to x modulo p^{mu}.
   */
    xmu:=(IntegerRing()!x) mod p^mu;
    if xmu lt 0 then
        xmu:=xmu+p^mu;
    end if;
    return xmu;
end function;

initialBound:=function(primelist)
    /*
      For |XYZ| = p_1^{b_1} ... p_s^{b_s}, determine a bound for each entry in
      the vector vecB = (b_1,...,b_s) using the bounds of von Kanel and Matshke.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
      Returns
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_s|.
   */
    assert &and[IsPrime(p) : p in primelist];
    assert 2 in primelist;
    N:=(2^4)*(&*primelist);
    s0N:=N*(3/4);
    sInfN:=2;
    d:=N*(&*[1+1/p : p in primelist]);
    l:=(N/6)*(&*[p+1 : p in primelist]);
    lStar:=d/6;

    g:=1+d/12;
    // Upper bound for the genus of the corresponding modular curve.
    m:=(s0N+7)/12-sInfN/2;
    // Upper bound for the number of newforms of level dividing N.
    betaBar:=(m/2)*Log(m)+(5/8)*m*(18+Log(l));
    betaPrime:=m*(Log(m)/2+(3/4)*Log(l)+Log(8.5));
    betaBarStar:=(g/2)*Log(g*lStar)+(lStar/2)*Log(4+4*Log(lStar));
    alphaBar:=Min([betaBar,betaBarStar,betaPrime]);
    heightBound:=(6/5)*alphaBar+28;
    vecB:=[Ceiling(heightBound/Log(p)) : p in primelist];
    return vecB;
end function;

latticeType:=function(primelist,p)
    /*
      Given p in p_1,...,p_s, find the prime p0 in the same set such that
      m0 = ord_p(log_p(p0)) is minimal and, where possible, such that the
      sublattice Gamma_mu^{*} may be used in place of Gamma_mu.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
          p: RngIntElt
      Returns
          plist: SeqEnum
              A list of rational primes p_1,...,p_{s-2},p_0.
          m0: RngIntElt
              The value ord_p(log_p(p0)).
          isReduced: BoolElt
              A true/false value. For p >= 5, if p0 can be chosen such that
	      ord_p(log_p(p0)) is minimal and p0 is a primitive (p-1)th root of
	      unity, return true, otherwise return false.
   */
    s:=#primelist;
    assert &and[IsPrime(q) : q in primelist];
    assert 2 in primelist;
    pind:=Index(primelist,p);
    plist0:=Remove(primelist,pind);

    prec:=20;
    Qp:=pAdicField(p,prec);
    vals:=[Valuation(Log(Qp!pi)) : pi in plist0];
    while true do
	// Increase precision until ord_p(log_p(q)) is constant.
	prec:=prec+10;
	Qp:=pAdicField(p,prec);
	valsNew:=[Valuation(Log(Qp!pi)) : pi in plist0];
	if (valsNew eq vals) then
	    break;
	else
	    vals:=valsNew;
	end if;
    end while;
    p0s:=[plist0[i] : i in [1..s-1] | vals[i] eq Min(vals)];
    m0:=Min(vals);
    isReduced:=false;
    if (p ge 5) and (PrimitiveRoot(p) in p0s) then
	p0:=PrimitiveRoot(p);
	assert (p0 in p0s);
	isReduced:=true;
    else
	p0:=p0s[1];
    end if;
    plist:=Exclude(plist0,p0) cat [p0];
    return plist,m0,isReduced;
end function;

approxLattice:=function(primelist,p,newBound)
    /*
      Given a rational prime p, determine the matrix B whose columns define the
      lattice Gamma_mu (or Gamma_mu^{*}) containing the vectors such that
      ord_p(Z) >= newBound+1.
      If newBound is chosen such that ord_p(Z) >= 1/(p-1) does not hold, or such
      that m<=0, return false,[].

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
          p: RngIntElt
          newBound: RngIntElt
      Returns
          tf: BoolElt
              A true/false value. If newBound is chosen such that
	      ord_p(Z) >= 1/(p-1) does not hold or such that mu <= 0, return
	      false, otherwise return true.
          plist: SeqEnum
              A list of rational primes p_1,...,p_{s-2},p_0.
          matB: AlgMatElt
	      The matrix whose columns define the lattice Gamma_mu
	      (or Gamma_mu^{*)).
   */
    s:=#primelist;
    assert &and[IsPrime(q) : q in primelist];
    assert 2 in primelist;
    plist,m0,isReduced:=latticeType(primelist,p);
    p0:=plist[s-1];
    Exclude(~plist,p0);
    mu:=newBound+1-m0;
    if ((mu+m0) lt 1) or (mu le 0) then
	return false,[],[];
    end if;
    Qp:=pAdicField(p,mu+5);
    thetaList:=[-Log(Qp!pi)/Log(Qp!p0) : pi in plist];
    thetaMuList:=[pAdicApprox(theta,p,mu) : theta in thetaList];
    assert (#thetaMuList eq s-2);
    matB:=IdentityMatrix(Integers(),s-1);
    matB[s-1]:=Vector(Integers(),thetaMuList cat [p^mu]);

    if (p ge 5) and isReduced then
	// Generate the matrix associated to the lattice Gamma_mu^{*}.
	alphaList:=[];
        for pi in plist do
	    alpha:=1;
	    while (Modexp(pi,1,p) ne Modexp(p0,alpha,p)) do
		alpha:=alpha+1;
	    end while;
	    assert alpha in [1..p-1];
	    Append(~alphaList,alpha);
	end for;
	gamma0:=Integers()!((p-1)/2);
	gammaList:=[(alphaList[i] + thetaMuList[i]) mod gamma0 : i in [1..s-2]];
	for i in [1..s-2] do
	    AddColumn(~matB,-gammaList[i],s-1,i);
	end for;
        MultiplyColumn(~matB,gamma0,s-1);
    end if;
    return true,plist cat [p0],matB;
end function;

pReductionTest:=function(primelist,p,vecB,newBound)
    /*
      Determine whether the lattice Gamma_mu (or Gamma_mu^{*}) contains vectors
      of length less than cB2p^2, where cB2p is the smallest upper bound on the
      L^2-norm of |b_1|,...,|b_{s-2}|,|b_0|.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
          p: RngIntElt
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_s|.
          newBound: RngIntElt
      Returns
          tf: BoolElt
              A true/false value. If Gamma_mu (or Gamma_mu^{*}) does not contain
	      any vectors of length less than cB2p^2, return true, otherwise
	      return false.
   */
    s:=#primelist;
    assert &and[IsPrime(q) : q in primelist];
    assert 2 in primelist;
    tf,plist,matB:=approxLattice(primelist,p,newBound);
    if (tf eq false) then
	// Cannot construct the lattice, newBound is too small.
	return false;
    end if;
    pind:=Index(primelist,p);
    cB2psq:=Integers()!Floor(&+[vecB[i]^2 : i in [1..s] | i ne pind]);
    LL:=Lattice(Transpose(matB));
    P:=ShortVectorsProcess(LL,cB2psq);
    if IsEmpty(P) then
	return true;
    else
	return false;
    end if;
end function;

reducedBound:=function(primelist)
    /*
      For |XYZ| = p_1^{b_1} ... p_s^{b_s}, determine the final bound for each
      of |b_1|,...,|b_s| after successively reducing the initial, very
      large bounds.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
      Returns
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_s|.
   */
    s:=#primelist;
    assert &and[IsPrime(q) : q in primelist];
    assert 2 in primelist;
    assert primelist eq Sort(primelist);
    vecB:=initialBound(primelist);

    finished:=false;
    repeat
	newVecB:=vecB;
	for i in [s..1 by -1] do
	    p:=primelist[i];
	    if vecB[i] gt 10 then
		bound:=10/1.3;
		// Set bound so that the first iteration begins at bound:=10.
	    else
		bound:=1/1.3;
		// Set bound so that the first iteration begins at bound:=5.
	    end if;
	    repeat
		bound:=bound*1.3;
		newBound:=Floor(bound);
		succeeded:=pReductionTest(primelist,p,newVecB,newBound);
	    until (succeeded or (newBound ge vecB[i]));
	    newVecB[i]:=newBound;
	end for;
	if (newVecB eq vecB) then
	    finished:=true;
	end if;
	if &and[newVecB[i] le vecB[i] : i in [1..s]] then
	    vecB:=newVecB;
	else
	    finished:=true;
	end if;
    until finished;
    return vecB;
end function;

solutionVectors:=function(primelist,plist,vecs)
    /*
      Given a rational number x/y = p_0^{b_0} ... p_{s-2}^{b_{s-2}}, determine
      whether the prime divisors of x+y or x-y are a subset of p_1,...,p_s, in
      which case extract X, Y and Z such that X + Y = Z.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
          plist: SeqEnum
              A list of rational primes p_1,...,p_{s-2},p_0.
	  vecs: SeqEnum
              A list of vectors of the lattice Gamma_mu (or Gamma_mu^{*}) having
	      squared L^2 norm <= cB2p^2.
      Returns
          sols: SetEnum
              A list of solutions [X,Y,Z].
   */
    if #vecs eq 0 then
	return {};
    end if;
    sols:={};
    s:=#primelist;
    assert plist subset primelist;
    assert #plist ge s-1;
    for ww in vecs do
	xy:=&*[plist[i]^ww[1,i] : i in [1..#plist]];
	zplus:=Abs(Numerator(xy)+Denominator(xy));
	zpfacs:=&*[p^(Valuation(zplus,p)) : p in primelist];
	if (zplus div zpfacs) eq 1 then
	    X:=Max(Numerator(xy),Denominator(xy));
            Y:=Min(Numerator(xy),Denominator(xy));
	    sols:=sols join {[X,Y,X+Y]};
	end if;
	zminus:=Abs(Numerator(xy)-Denominator(xy));
	zmfacs:=&*[p^(Valuation(zminus,p)) : p in primelist];
	if (zminus ne 0) and (zminus div zmfacs) eq 1 then
	    if (Numerator(xy) lt Denominator(xy)) then
                X:=Max(Numerator(xy),Abs(Numerator(xy)-Denominator(xy)));
                Y:=Min(Numerator(xy),Abs(Numerator(xy)-Denominator(xy)));
                sols:=sols join {[X,Y,Denominator(xy)]};
	    else
                X:=Max(Denominator(xy),Abs(Numerator(xy)-Denominator(xy)));
                Y:=Min(Denominator(xy),Abs(Numerator(xy)-Denominator(xy)));
                sols:=sols join {[X,Y,Numerator(xy)]};
	    end if;
	end if;
    end for;
    return sols;
end function;

latticeVolume:=function(primelist,vecB)
    /*
      For each rational prime p_1,...,p_s, generate a smaller bound on the
      correspondng exponent |b_1|,...,|b_s| and determine the volume of the
      associated lattice Gamma_mu, returning the rational prime and bound of the
      smallest such lattice.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_s|.
      Returns
          tf: BoolElt
              A true/false value. If the lattice Gamma_mu (or Gamma_mu^{*})
	      cannot be generated for any prime in p_1,...,p_s, return false,
	      otherwise return true.
          p: RngIntElt
          newBound: RngIntElt
          cB2psq: RntIntElt
              The smallest upper bound on the L^2-norm of
	      |b_1|,...,|b_{s-2}|,|b_0|.
   */
    s:=#primelist;
    assert &and[IsPrime(q) : q in primelist];
    assert 2 in primelist;
    assert primelist eq Sort(primelist);

    volInfo:=[];
    for i in [1..s] do
	p:=primelist[i];
        bi:=vecB[i];
	cB2psq:=Integers()!Floor(&+[vecB[j]^2 : j in [1..s] | j ne i]);
        if bi ne 0 then
	    if bi gt 15 then
		newBound:=bi-Ceiling(bi/3);
	    else
		newBound:=bi-Ceiling(0.125*bi);
	    end if;
	    assert newBound lt bi;
            tf,plist,matB:=approxLattice(primelist,p,newBound);
	    if tf then
		vol:=Ceiling((cB2psq^((s-1)/2))/Determinant(matB));
		Append(~volInfo,[p,vol,newBound,cB2psq]);
	    end if;
	end if;
    end for;
    if IsEmpty(volInfo) then
	return false,0,0,0;
    end if;
    mvol,ind:=Min([v[2] : v in volInfo]);
    p,vol,newBound,cB2psq:=Explode(volInfo[ind]);
    return true,p,newBound,cB2psq;
end function;

latticeSieve:=function(primelist,vecB)
    /*
      Determine all vectors in the lattice Gamma_mu (or Gamma_mu^{*}) having
      squared L^2 norm less than cB2p^2 and output the corresponding solutions
      X, Y and Z where possible.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_s|.
      Returns
          vecB: SeqEnum
              The list of smallest known bounds for |b_1|,...,|b_s|.
          sols: SetEnum
              A list of solutions [X,Y,Z].
   */
    s:=#primelist;
    assert &and[IsPrime(q) : q in primelist];
    assert 2 in primelist;
    assert primelist eq Sort(primelist);

    tf1,p,newBound,cB2psq:=latticeVolume(primelist,vecB);
    if (tf1 eq false) then
	return vecB,{};
    end if;
    pind:=Index(primelist,p);
    tf2,plist,matB:=approxLattice(primelist,p,newBound);
    assert tf2;
    LL:=Lattice(Transpose(matB));
    shortvecs:=ShortVectors(LL,cB2psq);
    vecBp:=[];
    for pi in plist do
	Append(~vecBp,vecB[Index(primelist,pi)]);
    end for;
    vecs:=[ww : ww in shortvecs |
	   &and[Abs(ww[1,i]) le vecBp[i] : i in [1..s-1]]];
    // Keep only those vectors less than the smallest known bounds for
    // |b_1|,...,|b_s|
    sols:=solutionVectors(primelist,plist,vecs);
    vecB[pind]:=newBound;
    return vecB,sols;
 end function;

solveXYZ:=function(primelist)
    /*
      Solve X+Y=Z for non-zero, positive coprime integers X, Y, and Z whose prime
      divisors all lie in {p_1,...,p_s}.

      Parameters
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_s.
      Returns
          sols: SetEnum
              A list of solutions [X,Y,Z].
   */
    assert &and[IsPrime(p) : p in primelist];
    assert 2 in primelist;
    Sort(~primelist);
    assert #primelist ge 2;
    s:=#primelist;

    vecB:=reducedBound(primelist);
    sols:={};
    finished:=false;
    repeat
	newVecB,psols:=latticeSieve(primelist,vecB);
	sols:=sols join psols;
	if (newVecB eq vecB) then
	    finished:=true;
	end if;
	if &and[newVecB[i] le vecB[i] : i in [1..s]] then
	    vecB:=newVecB;
	else
	    finished:=true;
	end if;
    until finished;
    vecs:=[[c] : c in CartesianProduct([[-bi..bi] : bi in vecB])];
    sols:=sols join solutionVectors(primelist,primelist,vecs);
    return sols;
end function;


Nlist,primelist:=extractForm(set);
if #primelist ge 3 then
    time sols:=SUnitXYZ(primelist);
    out:="dWTest/" cat seqEnumToString(primelist) cat ".txt";
    fprintf out, "%o\n",sols;
end if;
exit;
