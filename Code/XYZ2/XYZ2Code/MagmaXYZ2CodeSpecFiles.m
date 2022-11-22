/* 
To Run XYZ2Code: 
    1. Import SUnitCode Directory
    2. Change current working directory to SUnitCode
    3. Open Magma
    4. To run SUnitXYZ.m:
        load "./XYZ2Code/MagmaXYZFunctions/SUnitXYZ.m";
    5. To run SUnitXYZ2.m:
        i. load "./XYZ2Code/MagmaXYZ2Functions/SUnitXYZ2.m";
            Nb. This code calls upon SUnitXYZ.m, hence loading this algorithm automatically loads SUnitXYZ.m as well
        ii. S:= [p_1,...,p_s];
            here, p_i is a rational prime
        iii. Sol:= SUnitXYZ2(S);

*/


{ XYZ2Code
    { MagmaXYZFunctions
        { SUnitXYZ.m
            { SmallestLatticeBoundXYZ.m
                { LatticeBoundXYZ.m }
                { F3Approx.m }  // contained in MagmaCode/FinckePohstCode
                { LatticeXYZ.m 
                    { ConvertpAdic.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { InitialSortXYZ.m
                        { pAdicLog.m } } } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
            { FPRestrictionsXYZ.m 
                { SFactors.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { InitialSortXYZ.m
                    { pAdicLog.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { LatticeXYZ.m 
                    { ConvertpAdic.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { InitialSortXYZ.m
                        { pAdicLog.m } } }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { FinckePohst.m         // contained in MagmaCode/FinckePohstCode
                    { ShortVectors.m }
                    { MyCholesky.m } } } } }    
    { MagmaXYZ2Functions
        { SUnitXYZ2.m
            { ExponentsXYZ2.m }
            { SUnitXYZ.m          // contained in MagmaCode/MagmaXYZCode
                { SmallestLatticeBoundXYZ.m
                    { LatticeBoundXYZ.m }
                    { F3Approx.m }      // contained in MagmaCode/FinckePohstCode
                    { LatticeXYZ.m 
                        { ConvertpAdic.m }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { InitialSortXYZ.m
                            { pAdicLog.m } } } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { FPRestrictionsXYZ.m 
                    { SFactors.m }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { InitialSortXYZ.m
                        { pAdicLog.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { LatticeXYZ.m 
                        { ConvertpAdic.m }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { InitialSortXYZ.m
                            { pAdicLog.m } } }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions 
                    { FinckePohst.m     // contained in MagmaCode/FinckePohstCode
                        { ShortVectors.m }
                        { MyCholesky.m } } } }
            { SUnitXYZtoSUnitXYZ2.m }
            { DecompositionOfPrimesXYZ2.m }
            { AlphasXYZ2.m
                { ExponentsXYZ2.m }
                { IdealExponentsXYZ2.m }
                { IdealConjugatesXYZ2.m } }
            { SymmetricCaseXYZ2.m
                { IUFactorsXYZ2.m 
                    { FundamentalUnitXYZ2.m }
                    { IdealExponentsXYZ2.m } }
                { FundamentalUnitXYZ2.m }
                { SymmetricCaseZerosXYZ2.m }
                { BinaryRecurrenceSequence.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { BinareRecurrenceSequencePeriod.m }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { SFactors.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
            { IIPrimeXYZ2.m 
                { SplitPrimePropertiesXYZ2.m } }
            { MaximalC12BoundXYZ2.m 
                { AlphasXYZ2.m
                    { ExponentsXYZ2.m }
                    { IdealExponentsXYZ2.m }
                    { IdealConjugatesXYZ2.m } }
                { C12BoundXYZ2.m 
                    { FundamentalUnitXYZ2.m }
                    { QsqrtDPrecision.m }       // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { IUFactorsXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { IdealExponentsXYZ2.m } }
                    { nExponentsXYZ2.m 
                        { SFactors.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { LambdaBoundXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { C1pnXYZ2.m } }
                    { Ordp.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { NonEmptyMaxMinProductSum.m }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { KappaBoundXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { C1pnXYZ2.m } }
                    { Kappa_BoundXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { C1pnXYZ2.m } } } }
            { SmallestLatticeBoundXYZ2.m 
                { IUFactorsXYZ2.m 
                    { FundamentalUnitXYZ2.m }
                    { IdealExponentsXYZ2.m } }
                { LambdaLatticeXYZ2.m 
                    { LambdaInitialSortXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { LambdaLogpXYZ2.m 
                            { Ordp.m } } }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { LambdaLogpXYZ2.m 
                        { Ordp.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { ConvertpAdic.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { F3Approx.m }  // contained in MagmaCode/FinckePohstCode
                { NonEmptyMaxMinProductSum.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { KappaLatticeXYZ2.m 
                    { HenselLiftXYZ2.m }
                    { KappaInitialSortXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { HenselLiftXYZ2.m }
                        { pAdicLog.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { pAdicLog.m }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { ConvertpAdic.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { Kappa_LatticeXYZ2.m 
                    { HenselLiftXYZ2.m }
                    { Kappa_InitialSortXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { HenselLiftXYZ2.m }
                        { pAdicLog.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { pAdicLog.m }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { ConvertpAdic.m } } }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
            { FPParametersXYZ2.m 
                { FPRestrictionsXYZ2.m 
                    { IUFactorsXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { IdealExponentsXYZ2.m } }
                    { nExponentsXYZ2.m 
                        { SFactors.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { LambdaInitialSortXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { LambdaLogpXYZ2.m 
                            { Ordp.m } } }      // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { LambdaLatticeXYZ2.m 
                        { LambdaInitialSortXYZ2.m 
                            { FundamentalUnitXYZ2.m }
                            { LambdaLogpXYZ2.m 
                                { Ordp.m } } }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions 
                        { LambdaLogpXYZ2.m 
                            { Ordp.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions  
                        { ConvertpAdic.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { KappaInitialSortXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { HenselLiftXYZ2.m }
                        { pAdicLog.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { KappaLatticeXYZ2.m 
                        { HenselLiftXYZ2.m }
                        { KappaInitialSortXYZ2.m 
                            { FundamentalUnitXYZ2.m }
                            { HenselLiftXYZ2.m }
                            { pAdicLog.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { pAdicLog.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { ConvertpAdic.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { Kappa_InitialSortXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { HenselLiftXYZ2.m }
                        { pAdicLog.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { Kappa_LatticeXYZ2.m 
                        { HenselLiftXYZ2.m }
                        { Kappa_InitialSortXYZ2.m 
                            { FundamentalUnitXYZ2.m }
                            { HenselLiftXYZ2.m }
                            { pAdicLog.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { pAdicLog.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { ConvertpAdic.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { FinckePohst.m     // contained in MagmaCode/FinckePohstCode
                        { ShortVectors.m }
                        { MyCholesky.m } }      
                    { FPRestrictionsKappaXYZ2.m
                        { FundamentalUnitXYZ2.m }
                        { NonEmptyMaxMinProductSum.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { SFactors.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { FPRestrictionsKappa_XYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { NonEmptyMaxMinProductSum.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { SFactors.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { FPRestrictionsLambdaXYZ2.m 
                        { FundamentalUnitXYZ2.m }
                        { NonEmptyMaxMinProductSum.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                        { SFactors.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                    { NonEmptyMaxMinProductSum.m } }    // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { NonEmptyMaxMinProductSum.m } }        // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
            { FinalSearchXYZ2.m 
                { FundamentalUnitXYZ2.m }
                { NonEmptyMaxMinProductSum.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { SFactors.m }  // contained in MagmaCode/MagmaXYZ2Code/MagmaXYZ2IntrinsicFunctions
                { IdealExponentsXYZ2.m } } } } }