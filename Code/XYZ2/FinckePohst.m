// F3Approx.m

/*
INPUT:
    B:= a 2 x 2 square matrix whose columns {xy[1],xy[2]} form the basis of a lattice Gamma_mu

OUTPUT:
    The smallest norm in the lattice
        ie. the norm of the point in the lattice Gamma_mu of minimal length

COMMENTS:
    p-adic Approximation Algorithm as in Figure 3 of Page 62 of the Reference

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > B:= Matrix(2,2,[1,2,4,1]);
    > F3Approx(B);
    2.23606797749978969640917366873

    > B:= Matrix(2,2,[101,2,43,11]);
    > F3Approx(B);
    11.1803398874989484820458683437

*/


function F3Approx(B)
    xy:= [(Transpose(B))[1], (Transpose(B))[2]]; ;       // stores first and second column of matrix B

    if Sqrt(Norm(xy[1])) gt Sqrt(Norm(xy[2])) then
        xy:= [xy[2],xy[1]];     // swaps xy[1],xy[2]
    end if;
    while Sqrt(Norm(xy[1])) lt Sqrt(Norm(xy[2])) do
        u:= Round(InnerProduct(xy[1],xy[2])/InnerProduct(xy[1],xy[1]));         // computes integer u such that |y-u*x| is minimal
        xy[2]:= xy[2] - u*xy[1];        // updates y
        xy:= [xy[2],xy[1]];     // swaps xy[1],xy[2]
    end while;

    xy:= [xy[2],xy[1]];         // swaps xy[1],xy[2] to have |xy[1]| < |xy[2]|
    return Sqrt(Norm(xy[1]));   // computes norm of point with minimal lenght in lattice
end function;

// MyCholesky.m

/*
INPUT:
    A:= a square, symmetric, positive-definite matrix

OUTPUT:
    Q:= matrix [q_{i,j}], where
        Q(x):= sum_{i=1}^{n} q_{i,i}(x_i + sum_{j = i+1}^m q_{i,j}*x_j)^2, for the quadratic form, Q(x)
    R:= upper triangular matrix such that Transpose(R)*R = A

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > A:= Matrix(3,3,[2,1,1,1,2,1,1,1,2]);
    > Q,R:= MyCholesky(A);
    > Q;
    [2.00000000000000000000000000000 0.500000000000000000000000000000
        0.500000000000000000000000000000]
    [0.000000000000000000000000000000 1.50000000000000000000000000000
        0.333333333333333333333333333333]
    [0.000000000000000000000000000000 0.000000000000000000000000000000
        1.33333333333333333333333333333]
    > R;
    [1.41421356237309504880168872421 0.707106781186547524400844362105
        0.707106781186547524400844362105]
    [-0.000000000000000000000000000000 1.22474487139158904909864203735
        0.408248290463863016366214012451]
    [0.000000000000000000000000000000 -0.000000000000000000000000000000
        1.15470053837925152901829756100]

*/


function MyCholesky(A)
    n:= NumberOfColumns(A);     // computes the size of A
    Q:= ZeroMatrix(RealField(),n,n);
    R:= Transpose(Cholesky(A));         // returns upper triangular matrix R such that A = Transpose(R)*R via built-in function "Cholesky"
    for i in [1..n] do          // contructs Q from R, only requiring upper triangular portion of R
        Q[i,i]:= (R[i,i])^2;
        for j in [i+1..n] do
            Q[i,j]:= R[i,j]/R[i,i];
        end for;
    end for;
    return Q,R;
end function;

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

// FinckePohst.m

/*
INPUT:
    B:= [b_1,...,b_n], a matrix over ZZ with columns b_i,
        {b_1,...,b_n} form a basis of a lattice, L
    C:= a constant > 0
    Bound:= a constant > 0 determining the maximum number of short vectors to be computed

OUTPUT:
    if #z > Bound:
        maxReached:= "true", indicates that the number of short vector solutions is larger than Bound;
        M:= associated lattice matrix such that y_i:= Transpose(M*Transpose(x_i)) is a vector if the lattice with norm <= C, where
            x_i:= a lattice coefficient vector
        z:= the incomplete set of lattice coefficient vectors x_i of the lattice L such that y_i:= Transpose(M*Transpose(x_i)) has norm <= C,
            ie. Q(y_i) <= C, where Q is the binary quadratic form Q associated to the lattice L
            Equivalently, (y_i)*A*Transpose(y_i) = [z_i], where z_i <= C

    if #z <= Bound:
        maxReached:= "false", indicates that the number of short vector solutions is smaller than Bound;
        M:= associated lattice matrix such that y_i:= Transpose(M*Transpose(x_i)) is a vector if the lattice with norm <= C, where
            x_i:= a lattice coefficient vector
        z:= the set of all lattice coefficient vectors x_i of the lattice L such that y_i:= Transpose(M*Transpose(x_i)) has norm <= C,
            ie. Q(y_i) <= C, where Q is the binary quadratic form Q associated to the lattice L
            Equivalently, (y_i)*A*Transpose(y_i) = [z_i], where z_i <= C

COMMENTS:
    There are 2*(#z) total such vectors, including multiplication by a sign

REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

    U. Fincke and M. Pohst. Improved Methods For Calculating Vectors of Short Length in a Lattice,
        Including a Complexity Analysis. Mathematics of Computation, vol. 44, pp. 463–471, April 1985.

EXAMPLE:
    > B:= Matrix(3,3,[2,1,1,1,2,1,1,1,2]);
    > C:= 5;
    > Bound:= 100;
    > M,z,maxReached:= FinckePohst(B,C,Bound);
    > M;
    [ 1  0  1]
    [ 0  1  1]
    [-1 -1  2]
    > z;
    [
        [0 1 0],

        [-1  1  0],

        [1 0 0]
    ]
    > maxReached;
    false

    > B:= Matrix(3,3,[2,1,1,1,2,1,1,1,2]);
    > C:= 5;
    > Bound:= 2;
    > M,z,maxReached:= FinckePohst(B,C,Bound);
    > M;
    [ 1  0  1]
    [ 0  1  1]
    [-1 -1  2]
    > z;
    [
        [0 1 0],

        [-1  1  0],

        [1 0 0]
    ]
    > maxReached;
    true

*/


function FinckePohst(B,C,Bound)
    A:= Transpose(B)*B;         // symmetric, positive-definite matrix B
    n:= NumberOfColumns(A);     // computes dimension of square input matrix A
    A1:= ZeroMatrix(IntegerRing(),n,n);        // creates matrix A1, intentionally not positive-definite
    A1[2,1]:= -1;      // changes matrix A1 so that it is not symmetric; nb. n >= 2 since S contains at least 3 primes
    prc0:= Max( [ Max([ A[i,j] : j in [1..n] ]) : i in [1..n] ] );      // computes the largest integer in the matrix A
    prc:= Max(100,#IntegerToString(prc0));      // sets precision for the Cholesky decomposition based on the number of digits in the largest integer in the matrix A

    while IsSymmetric(A1) eq false do
        R:= Transpose(Cholesky(A, RealField(prc))); // applies Cholesky Decomposition to the input matrix
        t,Ri:= IsInvertible(R);             // computes R^{-1} (easy to compute since R is upper triangular)
        Si,Ui:= LLL(Ri);            // computes row-reduced version S^{-1} of R^{-1}, along with U^{-1}, where S^{-1} = U^{-1}*R^{-1}

        t,U:=IsInvertible(Ui);      // computes inverse U of U^{-1}
        S:= R*U;

        norms:= [Sqrt(Norm(Si[i])) : i in [1..n]];        // computes norms of rows of S^{-1}
        order:= Reverse(Sort(norms));                     // sorts norms by decreasing size
        per:= [Index(norms, i) : i in order];       // computes the permutation on (1..n) for which the n rows of S^{-1} have decreasing norms
        for i in [1..n] do
            a:=  &+[1: x in per | x eq i];    // counts repetition of index i
            if a gt 1 then                    // corrects repetition in per, ie. in the event any two norms are equal
                for j in [1..a-1] do
                    Ind:= Index(per,i);       // computes index of i in per
                    per[Index(per[Ind + 1..n],i) + Ind]:= i + 1;    // replaces repeated index, if any
                end for;
            end if;
        end for;

        P:= Transpose(PermutationMatrix(IntegerRing(),per));        // computes the permutation matrix corresponding to the permutation, per
        S1:= S*P;           // rewrites S with columns rearranged by the permutation per, above
        A10:= Transpose(S1)*S1;     // symmetric, positive-definite matrix A1, with entries significantly smaller than those of A
        A1:= Matrix(IntegerRing(),n,n,[IntegerRing()! Round(A10[i][j]) : j in [1..n], i in [1..n]]);        // converts matrix parent (RealField()) of A10 to IntegerRing()
        prc:= prc + 5;  // increases precision on RealField until A1 is Symmetric
    end while;

    while IsPositiveDefinite(A1) eq false do
        R:= Transpose(Cholesky(A, RealField(prc))); // applies Cholesky Decomposition to the input matrix
        t,Ri:= IsInvertible(R);             // computes R^{-1} (easy to compute since R is upper triangular)
        Si,Ui:= LLL(Ri);            // computes row-reduced version S^{-1} of R^{-1}, along with U^{-1}, where S^{-1} = U^{-1}*R^{-1}

        t,U:=IsInvertible(Ui);      // computes inverse U of U^{-1}
        S:= R*U;

        norms:= [Sqrt(Norm(Si[i])) : i in [1..n]];        // computes norms of rows of S^{-1}
        order:= Reverse(Sort(norms));                     // sorts norms by decreasing size
        per:= [Index(norms, i) : i in order];       // computes the permutation on (1..n) for which the n rows of S^{-1} have decreasing norms
        for i in [1..n] do
            a:=  &+[1: x in per | x eq i];    // counts repetition of index i
            if a gt 1 then                    // corrects repetition in per, ie. in the event any two norms are equal
                for j in [1..a-1] do
                    Ind:= Index(per,i);       // computes index of i in per
                    per[Index(per[Ind + 1..n],i) + Ind]:= i + 1;    // replaces repeated index, if any
                end for;
            end if;
        end for;

        P:= Transpose(PermutationMatrix(IntegerRing(),per));        // computes the permutation matrix corresponding to the permutation, per
        S1:= S*P;           // rewrites S with columns rearranged by the permutation per, above
        A10:= Transpose(S1)*S1;     // symmetric, positive-definite matrix A1, with entries significantly smaller than those of A
        A1:= Matrix(IntegerRing(),n,n,[IntegerRing()! Round(A10[i][j]) : j in [1..n], i in [1..n]]);        // converts matrix parent (RealField()) of A10 to IntegerRing()
        prc:= prc + 5;  // increases precision on RealField until A1 is PositiveDefinite
    end while;

    assert IsSymmetric(A1);     // verification that indeed A1 is a symmetric matrix; if false, there is an error in FinckePohst.m
    assert IsPositiveDefinite(A1);      // verification that indeed A1 is a positive definite matrix; if false, there is an error in FinckePohst.m

    //maxReached:= false;     // TEMPORARY TEST!!!!!
    //L:= LatticeWithGram(A1);
    //z1:= ShortVectorsMatrix(L,C);
    //delete L;
    //z:= [Matrix(z1[i]) : i in [1..NumberOfRows(z1)]];       // END TEST

    z,maxReached:= ShortVectors(A1,C,Bound);     // computes coefficient matrix of lattice points x in L with norm(x) <= C. ie. x*A1*Transpose(x) = [y], where y <= C; maxReached
    M:= B*U*P;      // computes permutation matrix B*U*P such that y:= Transpose(B*U*P*Transpose(z[i])) is a short lattice vector, for z[i] a short coefficient vector

    return M, z, maxReached;    // returns associated lattice matrix, M; short coefficient vectors, z; whether Bound was reached before the algorithm terminated, maxReached

end function;
