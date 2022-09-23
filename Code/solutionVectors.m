/*
solutionVectors.m

These functions sift through all potential solution vectors and output the
corresponding valid solutions as [X,Y,z_1,...,z_v].

Authors
    Adela Gherga <adelagherga@gmail.com>
    Samir Siksek <samir.siksek@gmail.com>
Created
    18 August 2022
*/

vectorTests:=function(alist,a,primelist,a0xy)
    /*
      Given an element of K, tau delta_1^{b_1} ... delta_r^{b_r}, determine
      whether it is of the form a_0 X - theta Y, in which case extract X, Y and
      determine whether it is a valid solution of
      a_0 X^d + ... + a_d Y^d = a * p_1^{z_1} ... p_v^{z_v}.

      Parameters
          alist: SeqEnum
              A list of integer coefficients a_0, a_1,...,a_d.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
          a0xy: SeqEnum
              The list of coefficients of tau delta_1^{b_1} ... delta_r^{b_r}.
      Returns
          sol: SetEnum
              A list of solutions [X,Y,z_1,...,z_v].
   */
    d:=#alist-1;
    ZUV<U,V>:=PolynomialRing(Integers(),2);
    F:=&+[alist[i+1]*U^(d-i)*V^i : i in [0..d]];
    a0:=alist[1];

    if (&and[a0xy[i] eq 0 : i in [3..d]] eq false) then
	// A valid solution must be of the form a0 X - theta Y.
	return {};
    end if;
    if (a0xy[1] notin Integers()) or (a0xy[2] notin Integers()) then
	return {};
    end if;
    a0xy:=[Integers()!a0xy[1],Integers()!a0xy[2]];
    if (IsDivisibleBy(a0xy[1],a0) eq false) then
	return {};
    end if;
    sol:=[a0xy[1] div a0,-a0xy[2]];
    if (GCD(sol[1],sol[2]) ne 1) or (GCD(a0,sol[2]) ne 1) then
	return {};
    end if;
    Fsol:=Evaluate(F,sol);
    if (IsDivisibleBy(Fsol,a) eq false) then
	return {};
    end if;
    Fsol:=Fsol div a;
    zlist:=[];
    for p in primelist do
	z:=Valuation(Fsol,p);
	Append(~zlist,z);
	Fsol:=Fsol div p^z;
    end for;
    if Fsol eq 1 then
	if IsEven(d) then
	    return {sol cat zlist,[-sol[1],-sol[2]] cat zlist};
	else
	    return {sol cat zlist};
	end if;
    elif Fsol eq -1 then
	if IsOdd(d) then
	    return {[-sol[1],-sol[2]] cat zlist};
	end if;
    end if;
end function;

solutionVectors:=function(alist,a,primelist,tau,deltaList,vecs)
    /*
      Given an element of K, tau delta_1^{b_1} ... delta_r^{b_r}, determine
      whether it is of the form a_0 X - theta Y, in which case extract X, Y and
      determine whether it is a valid solution of
      a_0 X^d + ... + a_d Y^d = a * p_1^{z_1} ... p_v^{z_v}.

      Parameters
          alist: SeqEnum
              A list of integer coefficients a_0, a_1,...,a_d.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
	  vecs: SeqEnum
              A list of vectors b of the translated lattice wc + Lc having
	      squared L^2 norm <= cBfsq.
      Returns
          sol: SetEnum
              A list of solutions [X,Y,z_1,...,z_v].
   */
    if #vecs eq 0 then
	return {};
    end if;
    eqSols:={};
    r:=#deltaList;
    v:=#primelist;
    d:=#alist-1;
    K:=NumberField(Universe([tau] cat deltaList));
    assert d eq Degree(K) and d ge 3;
    ZUV<U,V>:=PolynomialRing(Integers(),2);
    F:=&+[alist[i+1]*U^(d-i)*V^i : i in [0..d]];
    a0:=alist[1];
    for ww in vecs do
	a0xy:=tau*&*[deltaList[i]^ww[i] : i in [1..r]];
	a0xy:=Eltseq(K!a0xy);
	eqSols:=eqSols join vectorTests(alist,a,primelist,a0xy);
    end for;
    if (Abs(a0) eq 1) and (Abs(a) eq 1) then
	if (a0 eq a) then
	    assert (Evaluate(F,[1,0]) eq a0);
	    eqSols:=eqSols join {[1,0] cat [0 : i in [1..v]]};
	elif (a0*(-1)^d eq a) then
	    assert (Evaluate(F,[-1,0]) eq a);
	    eqSols:=eqSols join {[-1,0] cat [0 : i in [1..v]]};
	end if;
    end if;
    return eqSols;
end function;
