/*
parseString.m

These functions <description>

Authors
    Adela Gherga <adelagherga@gmail.com>
Created
    29 September 2022
*/

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
      Extracts Nlist,alist,a,primelist from the string set.

      Parameters
          set: MonStgElt
              A string in the format "Nlist,alist,a,primelist".
      Returns
          alist: SeqEnum
              A list of coefficients a_0, a_1,...,a_3.
          a: RngIntElt
          primelist: SeqEnum
              A list of rational primes p_1, p_2,...,p_v.
   */
    bracketSplit:=Split(set,"[]"); // Split bash input by "[" and "]".
    assert (#bracketSplit eq 4) or (#bracketSplit eq 5);
    alist:=[StringToInteger(a_i) : a_i in Split(bracketSplit[3],",")];
    a:=StringToInteger(Split(bracketSplit[4],",")[1]);
    if (#bracketSplit eq 4) then
	primelist:=[];
    else
	primelist:=[StringToInteger(p) : p in Split(bracketSplit[5],",")];
    end if;
    return alist,a,primelist;
end function;
