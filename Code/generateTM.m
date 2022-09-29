/*
generateTM.m

<description>

Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    27 September 2022
*/

findGL2Zactions:=function(a,c)
    /*
      Given a, c, determines integers b, d such that 2x2 matrix with entries
      [a b c d] defines an element of GL2(Z).

      Parameters
          a: RngIntElt
          c: RngIntElt
      Returns
          actions: SetEnum
              A list of elements [a,b,c,d] such that the 2x2 matrix with entries
              [a b c d] defines an element of GL2(Z).
   */
    actions:={};
    if ((a eq 0) and (Abs(c) eq 1)) then
	b:=-1/c;
	d:=0;
	assert (a*d-b*c eq 1);
	actions:=actions join {[a,b,c,d]};
    elif ((c eq 0) and (Abs(a) eq 1)) then
	d:=1/a;
	b:=0;
	assert (a*d-b*c eq 1);
	actions:=actions join {[a,b,c,d]};
    elif ((a ne 0) and (c ne 0) and (GCD(a,c) eq 1)) then
	g,d,b:=XGCD(a,c);
	b:=-b;
	assert g eq 1;
	assert (a*d-b*c eq 1);
	actions:=actions join {[a,b,c,d]};
    end if;
    return actions;
end function;

equivForm:=function(alist)

    /*
    // generate possible GL2(Z) actions under which c0 is small, avoiding Thue solver

     Description: generate all possible GL2(Z) actions under which c0 lies in [1..20]
     Input: clist:= [c_0, \dots, c_n], the coefficients of F(X,Y)
     Output: GL2Zclists:= all possible coefficients of F(X,Y) under GL2(Z) action under which
                          c0 lies in the interval [1..20]
     Example:
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
	if IsEmpty(Solutions(ThueF,i)) eq false then
	    a:=Solutions(ThueF,i)[1][1];
	    c:=Solutions(ThueF,i)[1][2];
	    GL2Zactions:=GL2Zactions join findGL2Zactions(a,c);
	end if;
    end for;
    absMax:=25;
    for a,c in [-absMax..absMax] do
	// Determine GL2(Z) actions under which a0 is small, avoinding the
	// Thue solver.
	if GCD(a,c) eq 1 then
	    GL2Zactions:=GL2Zactions join findGL2Zactions(a,c);
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
    // Sort resulting GL2(Z)-equivalent forms // LEFT OFF HERE
    if alist in GL2Zalists then
	Exclude(~GL2Zalists,alist);
    end if;
    if #GL2Zalists lt 10 then
	return [alist] cat GL2Zalists;
    else
	return [alist] cat GL2Zalists[1..10];
    end if;
end function;

optimalForm:=function(alist,a,primelist)
    /*
      <description>

      Parameters
          <param>: <param type>
	      <param description>
      Returns
          <param>: <param type>
	      <param description>
   */
    GL2Zalists:=equivForm(alist);
    if #GL2Zalists eq 1 then
	return GL2Zalists[1];
    end if;
    caseNo:=[0 : i in [1..#GL2Zalists]];
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
	no:=0;
	for j in [1..#newablist] do
            new_a:=Integers()!newablist[j][1][1];
            blist:=newablist[j][2];
	    assert &and[Valuation(new_a,p) eq 0 : p in primelist];
	    no:=no+#equationsInK(falist,new_a,primelist);
	end for;
	caseNo[i]:=no;
    end for;
    min,ind:=Min(caseNo);
    return GL2Zalists[ind];
end function;

// this should be it's own function in the same file as seqenumtostring
// need also to sort the code folder and the directories (ie. where we should be, etc)

commaSplit:=Split(set,","); // Split bash input by ",".
bracketSplit:=Split(set,"[]"); // Split bash input by "[" and "]".
assert (commaSplit[2][1] eq "[") and (commaSplit[5][#commaSplit[5]] eq "]");
assert (#bracketSplit eq 4);
N:=StringToInteger(commaSplit[1]);
alist:=[StringToInteger(a_i) : a_i in Split(bracketSplit[2],",")];
a:=StringToInteger(commaSplit[6]);
primelist:=[StringToInteger(p) : p in Split(bracketSplit[4],",")];


alist:=optimalForm(alist,a,primelist);

assert &and[IsPrime(p) : p in primelist];
assert &and[a_i in Integers() : a_i in alist];
a0:=Integers()!alist[1];
assert a0 ne 0;
d:=#alist-1;
assert d eq 3;
QUV<U,V>:=PolynomialRing(Rationals(),2);
Qx<x>:=PolynomialRing(Rationals());
F:=&+[alist[i+1]*U^(d-i)*V^i : i in [0..d]];
assert IsHomogeneous(F);
f:=a0^(d-1)*Evaluate(F,[x/a0,1]);
assert IsMonic(f);
assert Degree(f) eq d;
assert IsIrreducible(f);
falist:=Reverse(Coefficients(f));
assert &and[a_i in Integers() : a_i in falist];
falist:=[Integers()!a_i : a_i in falist];
newablist:=makeMonic(alist,a,primelist);

for i in [1..#newablist] do
    new_a:=Integers()!newablist[i][1][1];
    assert &and[Valuation(new_a,p) eq 0 : p in primelist];
    tauDeltaList:=equationsInK(falist,new_a,primelist);
    for j in [1..#tauDeltaList] do
	print N,alist,a,primelist,i,j;
    end for;
end for;
print "-========================";
end for;
