//SUnitXYZ2.m

/*
INPUT:
    S:= [p_1,...,p_s], p_i primes in S

OUTPUT:
    FinalSolns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = z_i^2, where
        x_i, y_i:= prod_{i:= 1 to s} p_i^{a_i} for rational integers x_i, y_i such that
            gcd(x_i, y_i) is squarefree
            x_i >= y_i and x_i >= 0
        z_i:= a rational integer, > 0

COMMENTS:
    Computes all S-unit equations x + y = z^2
    This algorithm uses the Fincke-Pohst algorithm and LLL

    Based on the algorithm found in Chapter 7 of the Reference

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

*/

ChangeDirectory("./Code/XYZ2");
load "./parseIO.m";
SetAutoColumns(false);
SetColumns(235);

Attach("./XYZ2Intrinsics.m");
load "./solveXYZ.m";
load "./FinckePohst.m";
load "./ExponentsXYZ2.m";
load "./DecompositionOfPrimesXYZ2.m";
load "./IdealExponentsXYZ2.m";
load "./IdealConjugatesXYZ2.m";
load "./AlphasXYZ2.m";
load "./FundamentalUnitXYZ2.m";
load "./IUFactorsXYZ2.m";
load "./SymmetricCaseXYZ2.m";
load "./SplitPrimePropertiesXYZ2.m";
load "./IIPrimeXYZ2.m";
load "./nExponentsXYZ2.m";
load "./C1pnXYZ2.m";
load "./LambdaBoundXYZ2.m";
load "./KappaBoundXYZ2.m";
load "./Kappa_BoundXYZ2.m";
load "./C12BoundXYZ2.m";
load "./MaximalC12BoundXYZ2.m";
load "./LambdaLogpXYZ2.m";
load "./LambdaInitialSortXYZ2.m";
load "./LambdaLatticeXYZ2.m";
load "./HenselLiftXYZ2.m";
load "./KappaInitialSortXYZ2.m";
load "./Kappa_InitialSortXYZ2.m";
load "./KappaLatticeXYZ2.m";
load "./Kappa_LatticeXYZ2.m";
load "./SmallestLatticeBoundXYZ2.m";
load "./FPRestrictionsKappaXYZ2.m";
load "./FPRestrictionsKappa_XYZ2.m";
load "./FPRestrictionsLambdaXYZ2.m";
load "./FPRestrictionsXYZ2.m";
load "./FPParametersXYZ2.m";
load "./FinalSearchXYZ2.m";

seqEnumToString:=function(X : quotes:=false)
    /*
      Converts a SeqEnum into a string without whitespace, enclosed by "[ ]" for
      .csv input

      Parameters
          X: SeqEnum
          quotes: BoolElt
              A true/false vale. If set to true, encloses the output in
	      quotations.
      Returns
          stX: MonStgElt
	      The set X as a string without whitespace.
   */
    strX:= "[";
    for i in [1..#X] do
	if X[i] in Integers() then
	    strX:=strX cat IntegerToString(Integers()!X[i]);
	elif X[i] in Rationals() then
	    strX:=strX cat IntegerToString(Numerator(X[i])) cat "/" cat
		  IntegerToString(Denominator(X[i]));
	end if;
	if (i ne #X) then
	    strX:=strX cat ",";
	end if;
    end for;
    strX:=strX cat "]";
    if quotes then
	strX:="\"" cat strX cat "\"";
    end if;
    return strX;
end function;

solTest:=function(u,y1,y2)
    X:=u^2;
    Y:=y1*y2;
    Z:=Integers()!((y1+y2)/2);
    if (X+Y eq Z^2) and (X ge Y) and (Z gt 0) then
	return {[X,Y,Z]};
    end if;
    return {};
end function;

convertToXYZ2:=function(xyzsols)
    /*
      <description>

      Parameters
          <param>: <param type>
	      <param description>
      Returns
          <param>: <param type>
	      <param description>
   */

    sols:={};
    for s in xyzsols do
	a,b,c:=Explode(s);
	if IsEven(a) then
	    sols:=sols join solTest(Integers()!(a/2),b,c);
	    sols:=sols join solTest(b,2*c,2*a);
	    sols:=sols join solTest(c,2*a,-2*b);
	elif IsEven(b) then
            sols:=sols join solTest(a,2*b,2*c);
	    sols:=sols join solTest(Integers()!(b/2),c,a);
	    sols:=sols join solTest(c,2*a,-2*b);
	else
	    assert IsEven(c);
            sols:=sols join solTest(a,2*b,2*c);
	    sols:=sols join solTest(b,2*c,2*a);
	    sols:=sols join solTest(Integers()!(c/2),a,-b);
	end if;
    end for;
    return sols;
end function;

function SUnitXYZ2(S)
    Sort(~S);
    assert 2 in S;
    s:=#S;
    assert s ge 2;
    Z:= IntegerRing();
    D0:= ExponentsXYZ2(S,0);    // generates all possible values for D:= ( (p_1)^(b_1) )* /cdots *( (p_n)^(b_n) ), where S:= [p_1, ..., p_n], b_i in {0,1}

    time xyzsols:=solveXYZ(S);       // computes all [x,y,z] where x + y = z
    sols:={[1,-1,0]} join convertToXYZ2(xyzsols);
    Exclude(~D0,1);

    for D in D0 do
        print "D:", D;
        b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
        if IsEmpty(SplitPrimes) then
            A0:= AlphasXYZ2([],b0,SplitPrimes,NonSplitPrimes,S,D);  // generates all alphas in the case I, I_ == []
            for A in A0 do
                sols:=sols join SymmetricCaseXYZ2(A,S,D);   // generates [x,y,z] where x + y = z in the symmetric case
            end for;
        else
            J:= IIPrimeXYZ2(SplitPrimes,D);         // generates all possible I, I_
            for II_ in J do
                A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);     // generates all alphas in the case I, I_ != []
                FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);    // computes maximal upper bound on U0, M0, M0_, absn
                U0, M0, M0_, absn:= SmallestLatticeBoundXYZ2(II_,FA,S);     // generates reduced bounds U0, M0, M0_, absn
                for A in A0 do
                    a:= A[3];
                    U01:= U0;
                    M01:= M0;
                    M0_1:= M0_;
                    absn1:= absn;
                    if (#M0 + #M0_) gt 2 then
                        U01, M01, M0_1, absn1, eqns:=
			    FPParametersXYZ2(FA,U01,M01,M0_1,absn1,b0,A,S);   // reduces U0, M0, M0_, absn via Fincke-Pohst, when more than 1 prime splits in K
                        sols:=sols join eqns;   // appends solutions that may have come from the Fincke-Pohst reduction
                    end if;
                    sols:=sols join FinalSearchXYZ2(U01,M01,M0_1,absn1,A,S,D);       // computes all [x,y,z] where x + y = z^2 below the reduced bounds U0, M0, M0_, absn
                end for;
            end for;
        end if;
    sols:=sols join {[D,-D,0]};    // appends the trivial solution [x,y,z]:= [D,-D,0]
    end for;

    FinalSolns:={};
    for s in sols do
        sqfree,sq:= Squarefree(GCD(Z!s[1],Z!s[2]));   // computes the squarefree integer sqfree as well as an integer sq, such that GCD(x,y) = (sqfree)*(sq^2)
	FinalSolns:=FinalSolns join {[s[1]/(sq^2), s[2]/(sq^2), s[3]/sq]};

//if sq ne 1 then         // if GCD(x,y) is not squarefree
//            FinalSolns:=FinalSolns join {[s[1]/(sq^2), s[2]/(sq^2), s[3]/sq]};
// appends only the reduced solution, [x,y,z] with GCD(x,y) squarefree, to FinalSolns
//        elif (s in FinalSolns) eq false then    // appends the solution, [x,y,z] to FinalSolns; this solution is reduced, ie. GCD(x,y) is squarefree
//            Append(~FinalSolns, s);
//        end if;
    end for;

    return FinalSolns;
end function;

//SUnitDXYZ2.m

/*
INPUT:
    S:= [p_1,...,p_s], p_i primes in S
    D:= the value ( (p_1)^(b_1) )* /cdots *( (p_{n_i})^(b_{n_i}) ) where b_i in {0,1} corresponding to S
    Ind:= the index of II_ where
        II_:= [I,I_], where, assuming the primes p_i of S which split in Q(sqrt(D)) are p_1,...,p_t, for some t in [0,s],
            I:= { i | 1 <= i <= t, a_i > b_i} and
            I_:= { i | 1 <= i <= t, a_i < b_i}
        Nb. if D == 1 or there are no split primes (ie. II_ == []), Ind is 0

OUTPUT:
    FinalSolns:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], all S-unit equations x_i + y_i = z_i^2, where
        x_i, y_i:= prod_{i:= 1 to s} p_i^{a_i} for rational integers x_i, y_i such that
            gcd(x_i, y_i) is squarefree
            x_i >= y_i and x_i >= 0
        z_i:= a rational integer, > 0

COMMENTS:
    Computes all S-unit equations x + y = z^2
    This algorithm uses the Fincke-Pohst algorithm and LLL

    Based on the algorithm found in Chapter 7 of the Reference

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5,7,11];
    > D:= 22;
    > Ind:= 3;
    > FinalSolns:= SUnitDXYZ2(S,D,Ind);
    > FinalSolns;
    [
        [ 126720000, 49, 11257 ],
        [ 117649, 12672, 361 ],
        [ 38808, 1, 197 ],
        [ 5632, -7, 75 ],
        [ 792, 49, 29 ],
        [ 352, -343, 3 ],
        [ 19800, 2401, 149 ],
        [ 88, -7, 9 ],
        [ 28512, 49, 169 ],
        [ 16038, -4802, 106 ],
        [ 198, -2, 14 ],
        [ 198, -98, 10 ],
        [ 66550, 14, 258 ],
        [ 22, 14, 6 ],
        [ 1078, 11, 33 ],
        [ 198, -77, 11 ],
        [ 550, 539, 33 ],
        [ 123750, 154, 352 ],
        [ 7546, 198, 88 ],
        [ 1782, 154, 44 ],
        [ 22, -22, 0 ]
    ]

*/

function SUnitDXYZ2(S,D,Ind)
    Sort(~S);
    assert 2 in S;
    s:=#S;
    assert s ge 2;

    Z:= IntegerRing();
    AllSolns:= [];      // stores all [x,y,z] where x + y = z^2

    if Ind eq 0 then    // corresponds to D == 1 or (I == [] and I_ == [])
        if D eq 1 then
            xyzsol:= SUnitXYZ(S);       // computes all [x,y,z] where x + y = z
            sol:= [];   // stores [x,y,z] where x + y = z^2 coming from SUnitXYZ.m
            for s in xyzsol do
                sol:= sol cat SUnitXYZtoSUnitXYZ2(s);   // converts xyz solutions to [x,y,z] where x + y = z^2
            end for;
            AllSolns:= AllSolns cat sol;        // appends solutions from SUnitXYZ.m
        else
            b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
            if (IsEmpty(SplitPrimes) eq false) then     // verification that indeed SplitPrimes == []; if false, there is an error in nSetsXYZ2.m in defining Ind
                print "Something is wrong in nSetsXYZ2: Ind == 0 but SplitPrimes != [].";
            end if;
            A0:= AlphasXYZ2([],b0,SplitPrimes,NonSplitPrimes,S,D);      // generates all alphas in the case I, I_ == []
            for A in A0 do
                AllSolns:= AllSolns cat SymmetricCaseXYZ2(A,S,D);       // generates [x,y,z] where x + y = z in the symmetric case
            end for;
        end if;
    else
        b0, SplitPrimes, NonSplitPrimes:= DecompositionOfPrimesXYZ2(S,D);
        J:= IIPrimeXYZ2(SplitPrimes,D);         // generates all possible I, I_
        II_:= J[Ind];   // selects II_ based on Ind
        A0:= AlphasXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);         // generates all alphas in the case I, I_ != []
        FA:= MaximalC12BoundXYZ2(II_,b0,SplitPrimes,NonSplitPrimes,S,D);        // computes maximal upper bound on U0, M0, M0_, absn
        U0, M0, M0_, absn:= SmallestLatticeBoundXYZ2(II_,FA,S);         // generates reduced bounds U0, M0, M0_, absn
        for A in A0 do
            a:= A[3];
            U01:= U0;
            M01:= M0;
            M0_1:= M0_;
            absn1:= absn;
            if (#M0 + #M0_) gt 2 then
                U01, M01, M0_1, absn1, eqns:= FPParametersXYZ2(FA,U01,M01,M0_1,absn1,b0,A,S);   // reduces U0, M0, M0_, absn via Fincke-Pohst, when more than 1 prime splits in K
                AllSolns:= AllSolns cat eqns;   // appends solutions that may have come from the Fincke-Pohst reduction
            end if;
            xyz2:= FinalSearchXYZ2(U01,M01,M0_1,absn1,A,S,D);   // computes all [x,y,z] where x + y = z^2 below the reduced bounds U0, M0, M0_, absn
            for s in xyz2 do
                if (s in AllSolns) eq false then
                    Append(~AllSolns, s);
                end if;
            end for;
        end for;
    end if;
    Append(~AllSolns, [D,-D,0]);        // appends the trivial solution [x,y,z]:= [D,-D,0]

    FinalSolns:= [];
    for s in AllSolns do
        sqfree,sq:= Squarefree(GCD(Z!s[1],Z!s[2]));     // computes the squarefree integer sqfree as well as an integer sq, such that GCD(x,y) = (sqfree)*(sq^2)
        if sq ne 1 then         // if GCD(x,y) is not squarefree
            Reduceds:= [s[1]/(sq^2), s[2]/(sq^2), s[3]/sq];
            if (Reduceds in FinalSolns) eq false then
                Append(~FinalSolns, Reduceds);  // appends only the reduced solution, [x,y,z] with GCD(x,y) squarefree, to FinalSolns
            end if;
        elif (s in FinalSolns) eq false then    // appends the solution, [x,y,z] to FinalSolns; this solution is reduced, ie. GCD(x,y) is squarefree
            Append(~FinalSolns, s);
        end if;
    end for;

    return FinalSolns;
end function;

extractForm:=function(set)
    /*
      Extracts Nlist,alist,a,primelist,[i,j] from the string set.

      Parameters
          set: MonStgElt
              A string in the format "Nlist,alist,a,primelist,[i,j]".
      Returns
          Nlist: SeqEnum
              A list of conductors.
          alist: SeqEnum
              A list of coefficients a_0, a_1,...,a_3.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
          ij: SeqEnum
              The index (i,j) of the corresponding S-unit equation.
   */
    bracketSplit:=Split(set,"[]");
    assert (#bracketSplit eq 3);
    Nlist:=[StringToInteger(N) : N in Split(bracketSplit[1],",")];
    primelist:=[StringToInteger(p) : p in Split(bracketSplit[3],",")];
    return Nlist,primelist;
end function;


Nlist,primelist:=extractForm(set);
time sols:=SUnitXYZ2(primelist);
out:="../../Data/Test/" cat seqEnumToString(primelist) cat ".txt";
fprintf out, "%o\n",sols;
exit;
