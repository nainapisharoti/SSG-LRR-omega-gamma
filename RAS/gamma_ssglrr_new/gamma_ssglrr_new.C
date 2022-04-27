/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "gamma_ssglrr_new.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr_new<BasicTurbulenceModel>::F1
(
    const volScalarField& CD
) const
{
    const volScalarField Ry(this->y_*sqrt(this->k_)/this->nu());
    const volScalarField F3(exp(-pow(Ry/120.0, 8)));

    return max(ssglrr_omega<BasicTurbulenceModel>::F1(CD), F3);
} 


template<class BasicTurbulenceModel>
tmp<volSymmTensorField> gamma_ssglrr_new<BasicTurbulenceModel>::P() const
{
 
    const volVectorField& U = this->U_;
    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const volScalarField k = this->k_;
    const volScalarField omega = this->omega_;

    volScalarField S(sqrt(2*magSqr(symm(gradU))));
    volScalarField Omega(sqrt(2*magSqr(skew(gradU))));
    volScalarField G(this->GName(), 0.5*mag(tr(ssglrr_omega<BasicTurbulenceModel>::P())));
   
    tgradU.clear();

    const volScalarField& nut = this-> nut();

    volScalarField Rev(sqr(this->y_)*S/nu);

    volScalarField F_lim = min(max(Rev/(Conset_*ReThetac_lim_)-1.0,scalar(0.0)),scalar(3.0));

    dimensionedScalar nu_dim = dimensionedScalar("0.0",sqr(dimLength)/dimTime,0.0);


    volScalarField P_lim
    (
        5.0*Ck_*max(gammaInt_-0.2,0.0)*(1.0-gammaInt_)*F_lim*max(3.0*Csep_*nu-nut,nu_dim)*G*(omega/k)
    );

    return  this->gammaInt_*ssglrr_omega<BasicTurbulenceModel>::P() + P_lim*I;
}

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> gamma_ssglrr_new<BasicTurbulenceModel>::D() const
{

    return max(this->gammaInt_,0.1)*ssglrr_omega<BasicTurbulenceModel>::D();
}

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> gamma_ssglrr_new<BasicTurbulenceModel>::Pi() const
{
    return  this->gammaInt_* ssglrr_omega<BasicTurbulenceModel>::Pi();
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> gamma_ssglrr_new<BasicTurbulenceModel>::Pi2() const
{
   const volScalarField::Internal gammaInt(this->gammaInt_); 
   return gammaInt*ssglrr_omega<BasicTurbulenceModel>::Pi2();
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
gamma_ssglrr_new<BasicTurbulenceModel>::ReThetac() const
{
    tmp<volScalarField::Internal> tReThetac
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::groupName("ReThetac", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            this->mesh_,
            dimless
        )
    );

    volScalarField::Internal& ReThetac = tReThetac.ref();

    const volVectorField& U = this->U_;
    const volVectorField& n_ = (wallDist::New(this->mesh_).n());

    const volScalarField::Internal& nu = this->nu();
    const volScalarField::Internal& k = this->k_;
    const volScalarField::Internal& omega = this->omega_();
    const volScalarField::Internal& y = this->y_();

    const volScalarField::Internal dVdy((fvc::grad(U & n_) & n_));
    
    const volScalarField::Internal TUL(min(100.0*(sqrt(2*k/3)/(y*omega)),scalar(100)));

    volScalarField::Internal lambdaTL1 = ((-0.00757*dVdy*(sqr(y)/nu) + 0.0128));
    volScalarField::Internal lambdaTL = (min(max(lambdaTL1,-1.0),1.0));

    volScalarField FPG(gammaInt_);

    forAll(gammaInt_,celli)
    { 

        if (lambdaTL[celli] >= 0)
        {
            FPG[celli] = min(1.0+CPG1_.value()*lambdaTL[celli],1.5);
        }
        else
        {
            FPG[celli] = min(1.0+CPG2_.value()*lambdaTL[celli]+CPG3_.value()*min(lambdaTL[celli]+0.0681,0.0),3.0);
        }

        FPG[celli]=(max(FPG[celli],0.0));

        ReThetac[celli] = CTU1_.value()+CTU2_.value()*exp(-CTU3_.value()*TUL[celli]*FPG[celli]);

    }
    return tReThetac;
} 

	
template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> gamma_ssglrr_new<BasicTurbulenceModel>::Fonset
(
    const volScalarField::Internal& Rev,
    const volScalarField::Internal& ReThetac,
    const volScalarField::Internal& RT
) const
{
    const volScalarField::Internal Fonset1(Rev/(Conset_*ReThetac));

    const volScalarField::Internal Fonset2
    (
        min(Fonset1, scalar(2.0))
    );

    const volScalarField::Internal Fonset3(max(1 - pow3(RT/3.5), scalar(0)));

    return tmp<volScalarField::Internal>
    (
        new volScalarField::Internal
        (
            IOobject::groupName("Fonset", this->U_.group()),
            max(Fonset2 - Fonset3, scalar(0))
        )
    );
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal> gamma_ssglrr_new<BasicTurbulenceModel>::Clam() const
{
   const volScalarField::Internal ReThetac(this->ReThetac());
  // const volScalarField::Internal& y = this->y_();
  // const volScalarField::Internal& k = this->k_;
  // const volScalarField::Internal& nu = this->nu();

    tmp<volScalarField::Internal> tClam
    (
        new volScalarField::Internal
        (
            IOobject
            (
                IOobject::groupName("Clam", this->U_.group()),
                this->runTime_.timeName(),
                this->mesh_
            ),
            this->mesh_,
            dimless
        )
    );
    volScalarField::Internal& Clam = tClam.ref();

    forAll(gammaInt_, celli)
    {
        const scalar gammaInt = gammaInt_[celli];
	//const scalar ReThetac = ReThetac[celli];
        const scalar c_omega(sqr(exp(-pow((420/(ReThetac[celli]+500)),4))));
        //const scalar Ry(y[celli]*sqrt(k[celli])/nu[celli]);
        //const scalar c_omega(sqr(exp(-pow((420/Ry),4))));

        if (c_omega>gammaInt)
        {
            Clam[celli] = 0;
        }
        else
        {
           Clam[celli] = (gammaInt-c_omega)/(1-c_omega);
          // Clam[celli] = gammaInt;
        }
    }

    return tClam;
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
gamma_ssglrr_new<BasicTurbulenceModel>::Pomega() const
{
   // const volScalarField k = this->k_;
   // const volScalarField omega = this->omega_;
   volScalarField G(this->GName(), 0.5*mag(tr(this->P())));

    const volScalarField CD
        (
            (fvc::grad(this->k_) & fvc::grad(this->omega_))
        );

    const volScalarField F1(this->F1(CD));
    const volScalarField::Internal alphaOmega(this->alphaOmega(F1));

    const volScalarField::Internal Clam(this->Clam());
    //const volScalarField::Internal gammaInt(this->gammaInt_);
    return Clam*alphaOmega*G*this->omega_()/(this->k_());
    //return alphaOmega*G*this->omega_()/(this->k_());
}

template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
gamma_ssglrr_new<BasicTurbulenceModel>::P2omega() const
{

   // const volScalarField k = this->k_;
   // const volScalarField omega = this->omega_;
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const volScalarField::Internal S(sqrt(2*magSqr(symm(gradU))));

    const volScalarField::Internal G(this->GName(), 0.5*mag(tr(this->P())));

    const volScalarField CD
    (
        (fvc::grad(this->k_) & fvc::grad(this->omega_))
    );

    const volScalarField F1(this->F1(CD));
    const volScalarField::Internal alphaOmega(this->alphaOmega(F1));

    const volScalarField::Internal Clam(this->Clam());
    //const volScalarField::Internal gammaInt(this->gammaInt_);
    return (1-Clam)*alphaOmega*sqr(S);
   // return alphaOmega*sqr(S);
}

/*template<class BasicTurbulenceModel>
tmp<volScalarField::Internal>
gamma_ssglrr_new<BasicTurbulenceModel>::Pomega() const
{
    const volScalarField k = this->k_;
    const volScalarField omega = this->omega_;
    const volScalarField::Internal G(this->GName(), 0.5*mag(tr(this->P())));

    const volScalarField CD
    (
        (fvc::grad(k) & fvc::grad(omega))
    );

    volScalarField F1(this->F1(CD));
    volScalarField::Internal alphaOmega(this->alphaOmega(F1));

    return alphaOmega*G*(this->omega_())/(this->k_());
}*/


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
gamma_ssglrr_new<BasicTurbulenceModel>::gamma_ssglrr_new
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    ssglrr_omega<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    ca2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ca2",
            this->coeffDict_,
            0.06
        )
    ),
    ce2_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "ce2",
            this->coeffDict_,
            50
        )
    ),
    Flength_
    (
        dimensionedScalar::lookupOrAddToDict
        (
            "Flength",
            this->coeffDict_,
            100
        )
    ),
    ReThetac_lim_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "ReThetac_lim",
            this->coeffDict_,
            1100
        )
    ),
    CTU1_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "CTU1",
            this->coeffDict_,
            100.0
        )
    ),
    CTU2_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "CTU2",
            this->coeffDict_,
            1000.0
        )
    ),
    CTU3_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "CTU3",
            this->coeffDict_,
            1.0
        )
    ),
    CPG1_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "CPG1",
            this->coeffDict_,
            14.68
        )
    ),
    CPG2_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "CPG2",
            this->coeffDict_,
            -7.34
        )
    ),
    CPG3_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "CTU3",
            this->coeffDict_,
            0.0
        )
    ),
    Conset_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "Conset",
            this->coeffDict_,
            2.2
        )
    ),
    Ck_
    (
        dimensionedScalar::lookupOrAddToDict
        (
           "Ck",
           this->coeffDict_,
           1.0
        )   
    ),
    Csep_
    (
        dimensionedScalar::lookupOrAddToDict
        (
           "Csep",
           this->coeffDict_,
           1.0
        )   
    ),


    deltaU_("deltaU", dimVelocity, SMALL),

    gammaInt_
    (
        IOobject
        (
            IOobject::groupName("gammaInt", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    
    Clam_
    (
	IOobject
	(
	    "Clam",
	    this->runTime_.timeName(),
	    this->mesh_,
	    IOobject::NO_READ,
	    IOobject::AUTO_WRITE
	),
	this->mesh_,
	dimensionedScalar("0.0",dimless,0.0)
    )
{
 
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool gamma_ssglrr_new<BasicTurbulenceModel>::read()
{
    if (gamma_ssglrr_new<BasicTurbulenceModel>::read())
    {
        ca2_.readIfPresent(this->coeffDict());
        ce2_.readIfPresent(this->coeffDict());
        Flength_.readIfPresent(this->coeffDict());
        ReThetac_lim_.readIfPresent(this->coeffDict());
        CTU1_.readIfPresent(this->coeffDict());
        CTU2_.readIfPresent(this->coeffDict());
        CTU3_.readIfPresent(this->coeffDict());
        CPG1_.readIfPresent(this->coeffDict());
        CPG2_.readIfPresent(this->coeffDict());
        CPG3_.readIfPresent(this->coeffDict());
        Conset_.readIfPresent(this->coeffDict());
        Ck_.readIfPresent(this->coeffDict());
        Csep_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class BasicTurbulenceModel>
void gamma_ssglrr_new<BasicTurbulenceModel>::correctGammaInt()
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const tmp<volScalarField> tnu = this->nu();
    const volScalarField::Internal& nu = tnu()();
    const volScalarField::Internal& y = this->y_();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Fields derived from the velocity gradient
    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal Omega(sqrt(2*magSqr(skew(tgradU()()))));
    const volScalarField::Internal S(sqrt(2*magSqr(symm(tgradU()()))));
    const volScalarField::Internal Us(max(mag(U()), deltaU_));
    const volScalarField::Internal dUsds((U() & (U() & tgradU()()))/sqr(Us));
    tgradU.clear();

    const volScalarField::Internal ReThetac(this->ReThetac());
    const volScalarField::Internal Rev(sqr(y)*S/nu);
    const volScalarField::Internal RT(k()/(nu*omega()));

    // volScalarField plim = this->gammaIntEff();

    // P_Lim_ = plim;

    {
        const volScalarField::Internal Pgamma
        (
            alpha()*rho()
            *Flength_*S*gammaInt_()*Fonset(Rev, ReThetac, RT)
        );

        const volScalarField::Internal Fturb(exp(-pow4(0.5*RT)));

        const volScalarField::Internal Egamma
        (
            alpha()*rho()*ca2_*Omega*Fturb*gammaInt_()
        );

        
        // Intermittency equation
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
          + fvm::div(alphaRhoPhi, gammaInt_)
          - fvm::laplacian(alpha*rho*DgammaInt(), gammaInt_)
        ==
            Pgamma - fvm::Sp(Pgamma, gammaInt_)
          + Egamma - fvm::Sp(ce2_*Egamma, gammaInt_)
          + fvOptions(alpha, rho, gammaInt_)
        );

        gammaIntEqn.ref().relax();
        fvOptions.constrain(gammaIntEqn.ref());
        solve(gammaIntEqn);
        fvOptions.correct(gammaInt_);
        bound(gammaInt_, 0);
    }
}


template<class BasicTurbulenceModel>
void gamma_ssglrr_new<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Correct gammaInt
    correctGammaInt();

    // Correct R and omega
    ssglrr_omega<BasicTurbulenceModel>::correct();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam


// ************************************************************************* //
