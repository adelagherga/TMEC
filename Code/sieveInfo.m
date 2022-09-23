/*
sieveInfo.m

For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, these functions
determine information needed for the final sieve with prime ideals outside
S. In particular, these functions determine a list of rational primes q,
coprime to the supports of tau, delta_1,...,delta_r, as well as the
corresponding set Tq and map phiq: Z^r --> Bq (as in Proposition 10.3).

Authors
    Adela Gherga <adelagherga@gmail.com>
    Samir Siksek <samir.siksek@gmail.com>
Created
    11 August 2022
*/

smallSieveInfo:=function(smallInfs,a0,theta,qBound)
    /*
      Determine generic information needed for the final sieve with prime ideals
      outside S (as in Section 10.2) for all prime ideals above q less than
      qBound. This information is generic in that no equation to solve is given,
      and hence no set of prime ideals S is given.

      Parameters
          smallInfs: List
              The list of elements (q,Bq,imdelFunc,SqShifted) already computed
	      under this function; this is [* *] in the first iteration.
	  a0: RngIntElt
          theta: RngOrdElt
              A generator of K.
          qBound: RngIntElt
              A positive integer determining the highest rational prime to use.
      Returns
          smallInfs: List
	      A list of elements (q,Bq,Piq,SqShifted), where q is a
	      rational prime less than qBound, Bq is the associated group
	      Aq/mu(Fq^*), and the associated subset Sq is obtained by
	      subtracting the image of tau under Piq from SqShifted. Here,
	      Piq is the map Oq^* --> Bq.
   */
    if IsEmpty(smallInfs) then // Determine primes already used.
	qMax:=0;
    else
	qMax:=smallInfs[#smallInfs][1];
    end if;
    if (qBound lt qMax+1) or IsEmpty(PrimesInInterval(qMax+1,qBound)) then
	return smallInfs;
    end if;
    // Run SmallSieveInfo only for those primes remaining.
    qlist:=[q : q in PrimesInInterval(qMax+1,qBound) | Valuation(a0,q) eq 0];
    printf "Getting the small sieve information";
    K:=NumberField(Parent(theta));
    assert K.1 eq theta;
    OK:=MaximalOrder(K);

    for q in qlist do
	qfacts:=[fq[1] : fq in Factorisation(q*OK)];
	if &and[Norm(fq) lt 10^10 : fq in qfacts] then
	    AqOrd:=&*[Norm(fq)-1 : fq in qfacts];
	    assert IsDivisibleBy(AqOrd,q-1);
	    BqOrd:=AqOrd div (q-1);
	    OModq,modq:=quo<OK | q*OK>;
	    Aq,eta:=UnitGroup(OModq);
	    Zqs,fromZqs:=UnitGroup(Integers(q)); // This is isomorphic to Fq^*.
	    ZqsGens:=[Integers()!fromZqs(z) : z in Generators(Zqs)];
	    // These are the generators of mu(Fq^*).
	    Bq,toBq:=quo<Aq | sub<Aq | [(modq(OK!z))@@eta : z in ZqsGens]>>;
	    // Bq is the group Aq/mu(Fq^*).
	    if (IsDivisibleBy(Discriminant(OK),q) eq false) then
		assert #Bq eq BqOrd;
	    end if;
	    Rq:=[(a0*u-theta) : u in [0..(q-1)]] cat [a0];
	    Rqstar:=[rho : rho in Rq |
		     &and[Valuation(rho*OK,fq) eq 0 : fq in qfacts]];
	    // Rqstar is Rq intersected with Oq^*.
	    Piq:=func<del | toBq((residueClassImage(del,modq,q*OK))@@eta)>;
	    // Piq is the map Oq^* --> Bq.
	    SqShifted:=[Piq(rho) : rho in Rqstar];
	    // Subtract the image of tau from SqShifted to get Sq.
	    Append(~smallInfs,[* q,Bq,Piq,SqShifted *]);
	end if;
    end for;
    printf "...Done.\n";
    return smallInfs;
end function;

bigSieveInfo:=function(tau,deltaList,smallInfs)
    /*
      For a_0 X - theta Y = tau delta_1^{b_1} ... delta_r^{b_r}, determine
      information needed for the final sieve with prime ideals
      outside S (as in Section 10.2).

      Parameters
          tau: FldNumElt
          deltaList: SeqEnum
              The list delta_1,...,delta_r.
          smallInfs: List
	      The list of elements (q,Bq,Piq,SqShifted), where q is a
	      rational prime less than qBound, Bq is the associated group
	      Aq/mu(Fq^*), and the associated subset Sq is obtained by
	      subtracting the image of tau under Piq from SqShifted. Here,
	      Piq is the map Oq^* --> Bq.
      Returns
          Zr: GrpAb
              The abelian group isomorphic to Z^r.
	  bigInfs: List
	      For each entry (q,Bq,Piq,SqShifted) of smallInfs, this is the
	      pair (phiq,Tq), where phiq defines the map Z^r --> Bq, and Tq
	      is the set Sq intersected with piq(Z^r).
   */
    K:=NumberField(Universe([tau] cat deltaList));
    OK:=MaximalOrder(K);
    supp:=&join[Support(delta*OK) : delta in deltaList];
    supp:=supp join Support(tau*OK);
    r:=#deltaList;
    Zr:=FreeAbelianGroup(r);
    bigInfs:=[* *];
    for I in smallInfs do
	q:=I[1];
	Bq:=I[2];
	Piq:=I[3];
	SqShifted:=I[4];
	qfacts:=[fq[1] : fq in Factorisation(q*OK)];
	if #(Set(qfacts) meet supp) eq 0 then
	    imsDelta:=[Piq(delta) : delta in deltaList];
	    Cq:=sub<Bq | imsDelta>;
	    phiq:=hom<Zr->Cq | [Cq!delta : delta in imsDelta]>;
	    // This is the map Z^r --> Bq.
	    imtau:=Piq(tau);
	    Sq:=[rho-imtau : rho in SqShifted];
            Tq:=[];
            for j in [1..#Sq] do
                rho:=Sq[j];
                if rho in Cq then
                    Append(~Tq,Cq!rho);
                end if;
            end for;
	    Append(~bigInfs,[* phiq,Tq *]);
	end if;
    end for;
    return Zr,bigInfs;
end function;
