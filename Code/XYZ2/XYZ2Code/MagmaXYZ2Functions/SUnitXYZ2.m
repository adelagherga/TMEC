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

ChangeDirectory("./Code");
load "./parseIO.m";
load "./solveXYZ.m";
SetAutoColumns(false);
SetColumns(235);

ChangeDirectory("./XYZ2");
Attach("./XYZ2Code/XYZ2Intrinsics.m");
load "./FinckePohstCode/F3Approx.m";
load "./FinckePohstCode/MyCholesky.m";
load "./FinckePohstCode/ShortVectors.m";
load "./FinckePohstCode/FinckePohst.m";

load "./XYZ2Code/MagmaXYZ2Functions/ExponentsXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/DecompositionOfPrimesXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/IdealExponentsXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/IdealConjugatesXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/AlphasXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FundamentalUnitXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/IUFactorsXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/SymmetricCaseXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/SplitPrimePropertiesXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/IIPrimeXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/nExponentsXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/C1pnXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/LambdaBoundXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/KappaBoundXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/Kappa_BoundXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/C12BoundXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/MaximalC12BoundXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/LambdaLogpXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/LambdaInitialSortXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/LambdaLatticeXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/HenselLiftXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/KappaInitialSortXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/Kappa_InitialSortXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/KappaLatticeXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/Kappa_LatticeXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/SmallestLatticeBoundXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsKappaXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsKappa_XYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsLambdaXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FPRestrictionsXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FPParametersXYZ2.m";
load "./XYZ2Code/MagmaXYZ2Functions/FinalSearchXYZ2.m";

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
    D0:= ExponentsXYZ2(S,0);    // generates all possible values for D:= ( (p_1)^(b_1) )* /cdots *( (p_n)^(b_n) ), where S:= [p_1, ..., p_n], b_i in {0,1}

    time xyzsols:=solveXYZ(S);       // computes all [x,y,z] where x + y = z
    sols:=convertToXYZ2(xyzsols);
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
    end if;
    sols:=sols join {[D,-D,0]};    // appends the trivial solution [x,y,z]:= [D,-D,0]


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
//if #primelist ge 3 then
    time sols:=solveXYZ(primelist);
    out:="../Data/Test/" cat seqEnumToString(primelist) cat ".txt";
    fprintf out, "%o\n",sols;
//end if;
exit;
