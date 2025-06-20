/*
optimalForm.m

These functions generate cubic forms that are GL2(Z)-equivalent to
F(X,Y) = a_0 X^3 + ... + a_3 Y^3 and determine the form yielding the least
number of S-unit equations for F(X,Y) = a * p_1^{z_1} ... p_v^{z_v}.

Parameters
    set: MonStgElt
        A string in the format "Nlist,alist,a,primelist".
    dir: MonStgElt
        A string in the format "Data/[N1,N2,...]i/TM" which serves as the name of
	the directory to which all output files are printed.
Returns
    OutFile: MonStgElt
        A .csv file named Nlist,alist,a,primelist.csv where alist defines a
        GL2(Z)-optimal cubic form associated to a, primelist. Each row is
	written as "Nlist,alist,a,primelist,[i,j]", with i,j iterating through
	the indices of the associated S-unit equations to solve. This function
	returns "Nlist,alist,a,primelist" in the case of a Thue equation, and
	omits forms with no S-unit equations.
Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    27 September 2022
*/

SetOutputFile(dir cat "/" cat set cat "tmp.txt");
ChangeDirectory("./Code/TM");
load "./solveThueMahler.m";

findGL2Zactions:=function(a,c)
    /*
      Given a, c, determine integers b, d such that the 2x2 matrix with entries
      [a b c d] defines an element of GL2(Z).

      Parameters
          a: RngIntElt
          c: RngIntElt
      Returns
          actions: SetEnum
              A set of 0 or 1 lists of elements [a,b,c,d] such that the 2x2 matrix with entries
              [a b c d] defines an element of GL2(Z).
   */
    g,d,b:=XGCD(a,c);
    if g eq 1 then
       return {[a,-b,c,d]};
    else
       return {};
    end if;
end function;

equivForm:=function(alist)
    /*
      Generate at most 10 GL2(Z)-equivalent forms to a_0 X^3 + ... + a_3 Y^3,
      sorted by the number of prime divisors of a0.

      Parameters
          alist: SeqEnum
              A list of coefficients a_0, a_1,...,a_3.
      Returns
          GL2Zalists: SeqEnum
	      A list of elements (alist_i), where each alist_i defines a cubic
              form GL2(Z)-equivalent to alist.
   */
    QUV<U,V>:=PolynomialRing(Rationals(),2);
    Qx<x>:= PolynomialRing(Integers());
    assert &and[a_i in Integers() : a_i in alist];
    assert (#alist-1) eq 3;
    F:=&+[alist[i+1]*U^(3-i)*V^i : i in [0..3]];
    assert IsHomogeneous(F);
    GL2Zactions:={};
    GL2Zalists:=[];

    ThueF:=Thue(Evaluate(F,[x,1]));
    testset:=PrimesInInterval(1,200) cat [1,4,9,25];
    Sort(~testset);
    for i in testset do
	// Determine GL2(Z) action yielding a0 in testset using the Thue solver.
        sols := Solutions(ThueF,i);
        if not IsEmpty(sols) then
	    a:=sols[1][1];
	    c:=sols[1][2];
	    GL2Zactions:=GL2Zactions join findGL2Zactions(a,c);
	end if;
    end for;
    absMax:=25;
    for a,c in [-absMax..absMax] do
	// Determine GL2(Z) actions under which a0 is small, avoinding the
	// Thue solver.
       g,d,b:=XGCD(a,c);
       assert a*d+b*c eq g;
       if g eq 1 then
          GL2Zactions:=GL2Zactions join {[a,-b,c,d]};
       end if;
    end for;
    for action in GL2Zactions do
	a,b,c,d:=Explode(action);
	GL2ZF:=Evaluate(F,[a*U+b*V,c*U+d*V]);
	newalist:=[MonomialCoefficient(GL2ZF,U^(3-i)*V^i) : i in [0..3]];
	newalist:=[Integers()!a_i : a_i in newalist];
	if (newalist[1] lt 0) then
	    newalist:=[-a_i : a_i in newalist];
	end if;
	a0:=newalist[1];
	if (#Divisors(a0) le 3) and (a0 le 5000) then
	    if (newalist notin GL2Zalists) then
		Append(~GL2Zalists,newalist);
	    end if;
	end if;
    end for;
    a0Eq1:=[newalist : newalist in GL2Zalists | newalist[1] eq 1];
    a0IsPrime:=[newalist : newalist in GL2Zalists | IsPrime(newalist[1])];
    a0Other:=[newalist : newalist in GL2Zalists |
	      newalist notin a0Eq1 and newalist notin a0IsPrime ];
    GL2Zalists:=Sort(a0Eq1) cat Sort(a0IsPrime) cat Sort(a0Other);
    // Sort resulting GL2(Z)-equivalent forms by number of prime factors of a0.
    if alist in GL2Zalists then
	Exclude(~GL2Zalists,alist);
    end if;
    if #GL2Zalists lt 10 then
	// Return at most 11 candidate GL2(Z)-equivalent forms.
	return [alist] cat GL2Zalists;
    else
	return [alist] cat GL2Zalists[1..10];
    end if;
end function;

optimalForm:=function(alist,a,primelist)
    /*
      Generate and test cubic forms GL2(Z)-equivalent to
      F(X,Y) = a_0 X^3 + ... + a_3 Y^3 and determine the form yielding the
      least number of S-unit equations for F(X,Y) = a * p_1^{z_1} ... p_v^{z_v},
      along with the indices to iterate over all associated S-unit equations.

      Parameters
          alist: SeqEnum
              A list of coefficients a_0, a_1,...,a_3.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
      Returns
          alist_i: SeqEnum
	      A list of coefficients a_0, a_1,...,a_3 defining an optimal
	      GL2(Z)-equivalent form to alist. That is, alist_i generates the
	      least number of S-unit equations of the GL2(Z)-equivalent forms
	      tested.
          caseNo_i: SeqEnum
              A list of elements [i,j] corresponding to alist_i, a, primelist,
              where i denotes the index of the potential monic Thue--Mahler
	      equation associated to alist_i, and j denotes the number of S-unit
	      equations corresponding to that monic equation.
   */
  GL2Zalists:=[alist]; // equivForm(alist);
    caseNo:=[];
    for i in [1..#GL2Zalists] do
	alist:=GL2Zalists[i];
	assert &and[IsPrime(p) : p in primelist];
	assert &and[a_i in Integers() : a_i in alist];
	a0:=Integers()!alist[1];
	assert a0 ne 0;
	d:=#alist-1;
	assert d ge 3;
	QUV<U,V>:=PolynomialRing(Rationals(),2);
	Qx<x>:=PolynomialRing(Rationals());
	F:=&+[alist[j+1]*U^(d-j)*V^j : j in [0..d]];
	assert IsHomogeneous(F);
	f:=a0^(d-1)*Evaluate(F,[x/a0,1]);
	assert IsMonic(f);
	assert Degree(f) eq d;
	assert IsIrreducible(f);
	falist:=Reverse(Coefficients(f));
	assert &and[a_i in Integers() : a_i in falist];
	falist:=[Integers()!a_i : a_i in falist];
	newablist:=makeMonic(alist,a,primelist);
	count:=[];
	no:=0;
	for j in [1..#newablist] do
            new_a:=Integers()!newablist[j][1][1];
            blist:=newablist[j][2];
	    assert &and[Valuation(new_a,p) eq 0 : p in primelist];
	    tauDeltaList:=equationsInK(falist,new_a,primelist);
	    Append(~count,[j,#tauDeltaList]);
	    no:=no+#tauDeltaList;
	end for;
	caseNo[i]:=<no,count>;
    end for;
    min,ind:=Min([c[1] : c in caseNo]);
    return GL2Zalists[ind],caseNo[ind][2];
end function;

OutFile:="../../" cat dir cat "/" cat set cat ".csv";
Nlist,alist,a,primelist,_:=extractForm(set);
if IsEmpty(primelist) then
    fprintf OutFile, "%o\n",set;
else
    alist,caseNo:=optimalForm(alist,a,primelist);
    for i in [1..#caseNo] do
	assert i eq caseNo[i][1];
	if (caseNo[i][2] ne 0) then
	    for j in [1..caseNo[i][2]] do
		fprintf OutFile, "%o,%o,%o,%o,%o\n",seqEnumToString(Nlist),
			seqEnumToString(alist),IntegerToString(a),
			seqEnumToString(primelist),seqEnumToString([i,j]);
	    end for;
	end if;
    end for;
end if;
exit;
