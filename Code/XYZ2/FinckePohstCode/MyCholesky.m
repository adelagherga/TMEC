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