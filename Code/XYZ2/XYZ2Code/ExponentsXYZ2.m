// ExponentsXYZ2.m

/*
INPUT:
    type:= 0 if 
        A:= [a_1,...,a_n], ie. there is only one data type in A
        
    type:= 1 if 
        A:= [[a_1,b_1],...,[a_n,b_n]], ie. there are 2 data types in A
        
OUTPUT:
    if type:= 0:
        D:= [c_1,...,c_m], where
            m:= the number of all possible products of the a_i with exponent 0 or 1
            c_i: one of all possible products of the a_i with exponent 0 or 1 
        output displayed in order from least integer to greatest

    if type:= 1:
        D:= [[c_1,d_1],...,[c_m,d_m]], where
            m:= the number of all possible products of the a_i with exponent 0 or 1
            c_i:= one of all possible products of the a_i with exponent 0 or 1
            d_i:= one of all possible products of the b_i with exponent 0 or 1

EXAMPLE:
    > A:= [2,3,5,7];
    > type:= 0;
    > ExponentsXYZ2(A,type);
    [ 1, 2, 3, 5, 6, 7, 10, 14, 15, 21, 30, 35, 42, 70, 105, 210 ]

    > A:= [[1,2],[3,4],[5,6],[7,8]];
    > type:= 1;
    > ExponentsXYZ2(A,type);
    [ <1, 1>, <1, 2>, <3, 4>, <3, 8>, <5, 6>, <5, 12>, <15, 24>, <15, 48>, <7, 8>, 
    <7, 16>, <21, 32>, <21, 64>, <35, 48>, <35, 96>, <105, 192>, <105, 384> ]
        
*/


function ExponentsXYZ2(A,type)
    
    D:= [];
    s:= #A;     // computes the length of A

    if type eq 0 then
        for n in [0..(2^s)-1] do
            iv := Intseq(n,2); 
            while #iv lt s do  
                Append(~iv,0);     // computes the combinatorial class of integer vectors of 0's and 1's of length s
            end while;    
            Append(~D,&*{A[j]^iv[j] : j in [1..s]});        // computes the product of elements of A raised to the respective exponent 1,0 of iv
        end for;
        Sort(~D);           // sorts the elements if type = 0
    
    elif type eq 1 then
        for n in [0..(2^s)-1] do
            iv := Intseq(n,2); 
            while #iv lt s do  
                Append(~iv,0);     // computes the combinatorial class of integer vectors of 0's and 1's of length s
            end while;
            Append(~D,<&*{A[j][1]^iv[j] : j in [1..s]},&*{A[j][2]^iv[j] : j in [1..s]}>);        // computes the product of elements of A raised to the respective exponent 1,0 of iv
        end for;
    end if;
    
    return D;
end function;
