/*
computeEllipticCurvesThue.m

This function solves the Thue equations defined in the string "set",
computes the corresponding elliptic curves, and outputs the logfile and elliptic
curve data in seperate files. Note that these output files may contain multiple
Thue form solutions, one for each integer rhs in "RHSlist".

Parameters
    set: MonStgElt
        A string in the format
	"N,\"form\",\"optimal form\",\"RHSlist\"".
Returns
    Outfile: MonStgElt
        A .csv file named "N,[optimal form]Out.csv" containing, for each
	elliptic curve E, the row
	"N,aInvariants(E),optimal form,rhs,sol". If there are no
	elliptic curves corresponding to this set, no such file is created.
    Logfile: MonStgElt
        A .txt file named "N,[optimal form]Log.txt" containing
	all output, timings, and solutions.
Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    24 August 2022
*/

ChangeDirectory("./Code");
load "./convertTMToEllipticCurves.m";

seqEnumToString:=function(X : quotes:=false)
    /*
      Convert a SeqEnum into a string without whitespace, enclosed by "[ ]" for
      .csv input

      Parameters
          X: SeqEnum
          quotes: BoolElt
              A true/false vale. If set to true, encloses the output in
	      quotations.
      Returns
          stX: MonStgElt
	      the set X as a string without whitespace.
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

CommaSplit:=Split(set,","); // Split bash input by ",".
RBracketSplit:=Split(set,"()"); // Split bash input by "(" and ")".
// Delimiter for form.
assert CommaSplit[2][2] eq "(" and CommaSplit[5][#CommaSplit[5]-1] eq ")";
// Delimiter for optimal form.
assert CommaSplit[6][2] eq "(" and CommaSplit[9][#CommaSplit[9]-1] eq ")";
assert (#RBracketSplit eq 7);
N:=StringToInteger(CommaSplit[1]);
alist:=[StringToInteger(i) : i in Split(RBracketSplit[4],",")];
rhsList:=[StringToInteger(i) : i in Split(RBracketSplit[6],",")];
OutFile:="../Data/ThueOutfiles/" cat
	 CommaSplit[1] cat "\,\[" cat RBracketSplit[4] cat "\]Out.csv";
LogFile:="../Data/ThueLogfiles/" cat
	 CommaSplit[1] cat "\,\[" cat RBracketSplit[4] cat "\]Log.txt";
SetLogFile(LogFile);
assert &and[a_i in Integers() : a_i in alist];
a0:=Integers()!alist[1];
assert a0 ne 0;
d:=#alist-1;
assert d ge 3;
ZUV<U,V>:=PolynomialRing(Integers(),2);
Zx<x>:=PolynomialRing(Integers());
F:=&+[alist[i+1]*U^(d-i)*V^i : i in [0..d]];
assert IsHomogeneous(F);
ThueF:=Thue(Evaluate(F,[x,1]));
for rhs in rhsList do
    printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    printf "N:=%o; alist:=%o; rhs:=%o; \n",N,alist,rhs;
    printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    time sols:=Solutions(ThueF,rhs);
    printf "sols:=%o\n",sols;
    ECs:=convertTMToEllipticCurves(N,alist,sols);
    printf "%o\n",ECs;
    printf "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    for E in ECs do
	assert E[1] eq N;
	fprintf OutFile, "%o, %o, %o, %o, %o\n",
		N,seqEnumToString(E[2]),seqEnumToString(alist),rhs,
		seqEnumToString(E[3]);
    end for;
end for;
UnsetLogFile();
