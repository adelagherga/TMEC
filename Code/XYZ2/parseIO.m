/*
parseIO.m

These functions parse terminal input/output data. Data intended for output is
converted into strings without whitespace for .csv files, and strings input from
the terminal are converted into appropriate magma structures.

Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    23 November 2022
*/

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

extractForm:=function(set)
    /*
      Extracts Nlist,primelist,[i,j] from the string set.

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
    assert (#bracketSplit in [3,5]);
    Nlist:=[StringToInteger(N) : N in Split(bracketSplit[1],",")];
    primelist:=[StringToInteger(p) : p in Split(bracketSplit[3],",")];
    if (#bracketSplit eq 5) then
	ij:=[StringToInteger(i) : i in Split(bracketSplit[5],",")];
    else
	ij:=[0,0];
    end if;
    return Nlist,primelist,ij;
end function;
