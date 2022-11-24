/*
convertXYZ2ToEllipticCurves.m

Thie function determines all elliptic curves over Q of conductor N associated to
the solutions of X+Y=Z^2, where X, Y, and N are made up of rational primes
{p_1,...,p_s}. Here, all resulting elliptic curves have nontrivial rational
2-torsion points.

Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    23 November 2022
*/

convertXYZ2ToEllipticCurves:=function(N,sols)
    /*
      Determine elliptic curves over Q of conductor N associated to the
      solutions of X+Y=Z^2, where X, Y, and N are made up of rational primes
      {p_1,...,p_s}. Here, all elliptic curves have nontrivial rational 2-torsion
      points.

      Parameters
          N: RngIntElt
          sols: SetEnum
              A list of solutions [X,Y,Z] to X+Y=Z^2.
      Returns
          relECs: SetEnum
              A list of sets (N,aInvs,[X,Y,Z]), where aInvs are the a-invariants
	      defining an elliptic curve of minimal model and conductor N,
	      determined by the solution [X,Y,Z] of the equation X+Y=Z^2 where
	      X, Y, and N are composed of rational primes {p_1,...,p_s}.
   */
    assert &and[IsPrime(p) : p in primelist];
    assert 2 in primelist;
    Sort(~primelist);
    s:=#primelist;
    assert s ge 2;

    divs6N:=[D : D in Divisors(6*N) | IsSquarefree(D)];
    assert 1 in divs6N;
    ECs:=[];
    for sol in sols do
        X,Y,Z:=Explode(sol);
	a2:=Z;
	a4List:=[X/4,Y/4];
	for a4 in a4List do
            E:=[0,a2,0,a4,0];
            minE:=MinimalModel(EllipticCurve(E));
            condE:=Conductor(minE);
            Append(~ECs,<condE,aInvariants(minE),sol>);
	    for D in divs6N do
		minE2:=MinimalModel(QuadraticTwist(minE,D));
		condE2:=Conductor(minE2);
		Append(~ECs,<condE2,aInvariants(minE2),sol>);
		twistMinE2:=MinimalModel(QuadraticTwist(minE2,-1));
		twistCondE2:=Conductor(twistMinE2);
		Append(~ECs,<twistCondE2,aInvariants(twistMinE2),sol>);
	    end for;
	end for;
    end for;
    Sort(~ECs);
    relECs:={};
    for E in ECs do
	if E[1] eq N then
	    relECs:=relECs join {E};
	end if;
    end for;
    return relECs;
end function;
