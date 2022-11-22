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