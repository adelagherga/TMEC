// ShortVectors.m

/*
INPUT:
    A:= a square, symmetric, positive-definite matrix, whose columns form the basis of the lattice L
    C:= a constant > 0
    Bound:= a constant > 0 determining the maximum number of short vectors to be computed

OUTPUT:
    if #sol > Bound:
        maxReached:= "true", indicates that the number of short vector solutions is larger than Bound; 
        sol:= the incomplete set of coefficient vectors x_i of the lattice L with norm <= C
              ie. Q(x_i) <= C, where Q is the binary quadratic form Q associated to the lattice L
              Equivalently, x_i*A*Transpose(x_i) = [y_i], where y_i <= C
                
    if #sol <= Bound:
        maxReached:= "false", indicates that the number of short vector solutions is smaller than Bound; 
        sol:= the set of all coefficient vectors x_i of the lattice L with norm <= C
              ie. Q(x_i) <= C, where Q is the binary quadratic form Q associated to the lattice L
              Equivalently, x_i*A*Transpose(x_i) = [y_i], where y_i <= C
    
COMMENTS:
    There are 2*(#sol) total solutions, including multiplication by a sign
    
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.
    
    R. Meissen. Lattice Isometries and Short Vector Enumeration. MsC thesis, University of Calgary, 2014.

EXAMPLE:
    > A:= Matrix(3,3,[2,1,1,1,2,1,1,1,2]);
    > C:= 5;
    > Bound:= 9;
    > sol, maxReached:= ShortVectors(A,C,Bound);
    > sol;    
    [
        [-1  1  1],

        [0 0 1],

        [-1  0  1],

        [ 1 -1  1],

        [ 0 -1  1],

        [-1 -1  1],

        [0 1 0],

        [-1  1  0],

        [1 0 0]
    ]
    > maxReached;
    false
    
    > A:= Matrix(3,3,[2,1,1,1,2,1,1,1,2]);
    > C:= 5;
    > Bound:= 8;
    > sol, maxReached:= ShortVectors(A,C,Bound);
    > sol;    
    [
        [-1  1  1],

        [0 0 1],

        [-1  0  1],

        [ 1 -1  1],

        [ 0 -1  1],

        [-1 -1  1],

        [0 1 0],

        [-1  1  0],

        [1 0 0]
    ]
    > maxReached;
    true

*/


function ShortVectors(A,C,Bound)
    bound:= C + 0.99;   // adds 0.99 to bound to avoid roundoff errors in Ceiling, Floor() which may result in missing vectors
    Q,R:= MyCholesky(A);        // applies Cholesky Decomposition to A, generating the quadratic form Q
    n:= NumberOfColumns(A);

    T:= ZeroMatrix(RealField(100), 1, n);
    U:= ZeroMatrix(RealField(100), 1, n);
    L:= ZeroMatrix(IntegerRing(), 1, n);
    x:= ZeroMatrix(IntegerRing(), 1, n);
    
    sol:= [];   // stores coefficient vectors
    i:= n;
    T[1,i]:= bound;
    stillEnumerating:= true;
    numSol:= 0;
    maxReached:= false;

    while (stillEnumerating eq true) do
        Z:= Sqrt(T[1,i]/Q[i,i]);
        L[1,i]:= Floor(Z-U[1,i]);
        x[1,i]:= Ceiling(-Z-U[1,i])-1;

        updatingEntries:= true;
        while updatingEntries eq true do
            x[1,i]:= x[1,i] + 1;

            if x[1,i] gt L[1,i] then
                i:= i + 1;      // updatingEntires still 'true': returns back up to x[i] = x[i]+1 until x[i] <= L[i]
            else
                if i eq 1 then
                    if x eq ZeroMatrix(IntegerRing(), 1, n) then       // if x_j = 0 for j in range(0,n+1)
                        stillEnumerating:= false;       // updatingEntries is 'false': exits subloop
                        updatingEntries:= false;        // stillEnumerating is 'false': terminates algorithm
                    else
                        Append(~sol,-x);
                        numSol:= numSol + 1;
                        if numSol gt Bound then
                            maxReached:= true;      // maxReached is 'true': maximum number of solutions (Bound) is reached before the algorithm has terminated
                            stillEnumerating:= false;       // updatingEntries is 'false': exits subloop
                            updatingEntries:= false;        // stillEnumerating is 'false': terminates algorithm
                        end if;
                    end if;    
                else
                    T[1,i-1]:= T[1,i] - Q[i,i]*(x[1,i]+U[1,i])^2;
                    i:= i - 1;
                    U[1,i]:= &+[Q[i,j]*x[1,j]: j in [i+1..n]];  // recomputes U_i
                    updatingEntries:= false;                    // updatingEntries = 'false': returns to StillEnumerating start
                end if;
            end if;
        end while;
    end while;
    
    return sol, maxReached;     // returns the coefficients of the short vector of L; whether Bound was reached before the algorithm terminated, maxReached
    
end function;    