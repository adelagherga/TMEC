//SUnitXYZtoSUnitXYZ2.m

/*
INPUT:
    s:= [x,y,z], an S-unit equation x + y = z, where
        S:= [p_1,...,p_s], p_i primes in S
        x, y, z:= prod_{p_i in S} p_i^{a_i}, for pairwise relatively prime rational integers x, y, z

OUTPUT:
    Sol:= [[x_1,y_1,z_1],...,[x_n,y_n,z_n]], where
        [x_i, y_i, z_i]:= an S-unit equations x_i + y_i = z_i^2 pertaining to s:= [x,z,z] such that 
            x_i, y_i:= prod_{i:= 1 to s} p_i^{a_i} for rational integers x_i, y_i such that
                gcd(x_i, y_i) is squarefree
                x_i >= y_i and x_i >= 0
            z_i:= a rational integer, > 0 

COMMENTS:
    Corresponds to the Case D:= 1 in Chapter 7 of the Reference 
    
REFERENCE:
    B.M.M.De Weger. Algorithms For Diophantine Equations. PhD thesis, University of Leiden, 1988.

EXAMPLE:
    > S:= [2,3,5];
    > s:= [ 27, 5, 32 ]; 
    > SUnitXYZtoSUnitXYZ2(s);
    [
        [ 729, 640, 37 ],
        [ 256, -135, 11 ]
    ]
    
    > S:= [2,3,5,7];
    > s:= [420,21,441];
    > SUnitXYZtoSUnitXYZ2(s);
    [
        [ 44100, 9261, 231 ],
        [ 194481, -35280, 399 ]
    ]

*/


function SUnitXYZtoSUnitXYZ2(s)
    xyz:= [];   // stores all solutions corresponding to s
    Z:= IntegerRing();
    a:= Z!s[1];         // converts s[1] into an integer
    b:= Z!s[2];         // converts s[2] into an integer
    c:= Z!s[3];         // converts s[3] into an integer
    
    if a mod 2 eq 0 then        // if 2|a
        if ((a/2)^2 ge b*c) and (c-(a/2) gt 0) then     // tests if x >= y and z > 0
            Append(~xyz, [Z!(a/2)^2,Z!b*c,Z!c-(a/2)]);  // stores [x,y,z] where x + y = z^2 under the condition that x >= y and z > 0
        end if;
        if ((b)^2 ge 4*c*a) and (2*c-b gt 0) then
            Append(~xyz, [Z!(b)^2,Z!4*c*a,Z!2*c-b]);
        end if;
        if ((c)^2 ge -4*a*b) and (2*a-c gt 0) then
            Append(~xyz, [(c)^2,-4*a*b,2*a-c]);
        end if;
    elif b mod 2 eq 0 then      // if 2|b
        if ((a)^2 ge 4*b*c) and (2*c-a gt 0) then
            Append(~xyz, [(a)^2,4*b*c,2*c-a]);
        end if;
        if ((b/2)^2 ge c*a) and (c-(b/2) gt 0) then 
            Append(~xyz, [(b/2)^2,c*a,c-(b/2)]);
        end if;
        if ((c)^2 ge -4*a*b) and (2*a-c gt 0) then
            Append(~xyz,[(c)^2,-4*a*b,2*a-c]);
        end if;
    else                        // if 2|c
        if ((a)^2 ge 4*b*c) and (2*c-a gt 0) then
            Append(~xyz,[(a)^2,4*b*c,2*c-a]);
        end if;
        if ((b)^2 ge 4*c*a) and (2*c-b gt 0) then
            Append(~xyz,[(b)^2,4*c*a,2*c-b]);
        end if;
        if ((c/2)^2 ge -a*b) and (a-(c/2) gt 0) then
            Append(~xyz,[(c/2)^2,-a*b,a-(c/2)]);
        end if;
    end if;

    Sol:= [];   // stores the relevant solutions of xyz; ie. without redundancy and such that x >= y 
    for sol in xyz do
        if (sol in Sol eq false) and (sol[1] gt sol[2]) then
            Append(~Sol, sol);
        elif (sol in Sol eq false) and (sol[1] le sol[2]) then
            Append(~Sol, [sol[2],sol[1],sol[3]]);
        end if;
    end for;
    return Sol;
end function;
    