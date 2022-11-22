// NonEmptyMaxMinProductSum.m

/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]: 
        The maximal element of C, Max(C)
    if C = [ ]:
        Returns 0
    
COMMENTS:
    Built-in MAGMA function Max(C) returns error if C == [ ] 

EXAMPLE:
    > C:= [ ];
    > Max(C);

    >> Max(C);
          ^
    Runtime error in 'Max': Argument 1 is not non-empty
    
    > NonEmptyMax(C);
    0
    
*/


intrinsic NonEmptyMax(C::SeqEnum) -> .  // any type 
    {Returns Max(C) if C is nonempty, otherwise returns 0.}
    if IsEmpty(C) eq true then
        MaxC:= 0;       // consistent with SAGE output
    else
        MaxC:= Max(C);
    end if;
    return MaxC;        
end intrinsic;


/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]: 
        The minimal element of C, Min(C)
    if C = [ ]:
        Returns 0
    
COMMENTS:
    Built-in MAGMA function Min(C) returns error if C == [ ] 

EXAMPLE:
    > C:= [ ];
    > Min(C);

    >> Min(C);
          ^
    Runtime error in 'Min': Argument 1 is not non-empty
    
    > NonEmptyMin(C);
    0
    
*/


intrinsic NonEmptyMin(C::SeqEnum) -> .  // any type 
    {Returns Min(C) if C is nonempty, otherwise returns 0.}
    if IsEmpty(C) eq true then
        MinC:= 0;       // consistent with SAGE output
    else
        MinC:= Min(C);
    end if;
    return MinC;        
end intrinsic;


/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]: 
        The product of all element of C, &*(C)
    if C = [ ]:
        Returns 1
    
COMMENTS:
    Built-in MAGMA function &*(C) returns error if C == [ ] 

EXAMPLE:
    > C:= [ ];
    > &*(C);

    >> &*(C);
       ^
    Runtime error in '&*': Illegal null sequence
    
    > NonEmptyProduct(C);
    1
    
*/


intrinsic NonEmptyProduct(C::SeqEnum) -> .      // any type
    {Returns the product of elements of C, &*(C), if C is nonempty, otherwise returns 1.}
    if IsEmpty(C) eq true then
        ProductC:= 1;   // consistent with SAGE output
    else
        ProductC:= &*C;
    end if;
    return ProductC;        
end intrinsic;


/*
INPUT:
    C:= [c_1,...,c_n], an arbitrary sequence of elements c_i, potentially empty
        ie. C:= [ ]

OUTPUT:
    if C != [ ]: 
        The sum of all element of C, &*(C)
    if C = [ ]:
        Returns 0
    
COMMENTS:
    Built-in MAGMA function &+(C) returns error if C == [ ] 

EXAMPLE:
    > C:= [ ];
    > &+(C);

    >> &+(C);
       ^
    Runtime error in '&+': Illegal null sequence
    
    > NonEmptySum(C);
    0
    
*/


intrinsic NonEmptySum(C::SeqEnum) -> .
    {Returns the sum of elements of C, &+(C), if C is nonempty, otherwise returns 0.}
    if IsEmpty(C) eq true then
        SumC:= 0;       // consistent with SAGE output
    else
        SumC:= &+C;
    end if;
    return SumC;        
end intrinsic;