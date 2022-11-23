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