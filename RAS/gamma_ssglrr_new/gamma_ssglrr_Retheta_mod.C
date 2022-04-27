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

#include "gamma_ssglrr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::F1
(
    const volScalarField& CD
) const
{
    const volScalarField Ry(this->y_*sqrt(this->k_)/this->nu());
    const volScalarField F3(exp(-pow(Ry/120.0, 8)));

    return max(ssglrr_omega<BasicTurbulenceModel>::F1(CD), F3);
}  

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::gammaIntEff() const
{

	const volVectorField& U = this->U_;
    
    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

	tmp<volTensorField> tgradU(fvc::grad(U));
	const volTensorField& gradU = tgradU();
	volScalarField S(sqrt(2*magSqr(symm(gradU))));

    volScalarField ReThetac(max(this->ReThetac(),1e-6));

    volScalarField Rev(sqr(this->y_)*S/nu);
    volScalarField RT(this->k_/(nu*this->omega_));

	volScalarField Freattach(exp(-pow4(RT/20.0)));

    volScalarField F_thetat();

    volScalarField Re_omega((this->omega_*sqr(this->y_))/nu);
    volScalarField Fwake(exp(-sqr(Re_omega/1e5)));
    volScalarField Fturb(exp(-pow4(0.5*RT)));

    volScalarField F_bl(min(
                            max(Fwake,
                                1.0-sqr((this->gammaInt_-1/ce2_)/(1.0-1/ce2_))),1.0));

    volScalarField gammaSep(min(2*max(Rev/(3.235*(ReThetac)) - 1, scalar(0))*Freattach, scalar(2))*F_bl);

    return max(gammaSep,this->gammaInt_);

}

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::Pgamma() const
{

    const volVectorField& U = this->U_;
    const volScalarField& gammaInt = this->gammaInt_;
    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;

    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();
    volScalarField S(sqrt(2*magSqr(symm(gradU))));

    const volScalarField Rev(sqr(this->y_)*S/nu);
    const volScalarField RT(k/(nu*omega));

    volScalarField ReThetac(max(this->ReThetac(),1e-6));

    const volScalarField fonset(Fonset(Rev, ReThetac, RT));

    return Flength_*S*gammaInt*fonset;
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::Egamma() const
{

    const volVectorField& U = this->U_;

    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volScalarField& gammaInt = this->gammaInt_;

    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const volScalarField Omega(sqrt(2*magSqr(skew(gradU))));

    const volScalarField RT(k/(nu*omega));

    const volScalarField Fturb(exp(-pow4(0.5*RT)));

    return ca2_*Omega*Fturb*gammaInt;
}

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> gamma_ssglrr<BasicTurbulenceModel>::P() const
{

    volScalarField gammaInt_eff = gammaIntEff();

    const volVectorField& U = this->U_;
    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();
    volScalarField S(sqrt(2*magSqr(symm(gradU))));
    volScalarField Omega(sqrt(2*magSqr(skew(gradU))));

    tgradU.clear();

    const volScalarField& nut = this-> nut();

    volScalarField Rev(sqr(this->y_)*S/nu);

    volScalarField F_lim = min(max(Rev/(Conset_*ReThetac_lim_)-1.0,scalar(0.0)),scalar(3.0));

    dimensionedScalar nu_dim = dimensionedScalar("0.0",sqr(dimLength)/dimTime,0.0);

    volScalarField CD 
    (
        (fvc::grad(this->k_) & fvc::grad(this->omega_))
    );

    volScalarField gammaSep = gammaIntEff();

    volScalarField P_lim
    (
        5.0*Ck_*max(this->gammaInt_-0.2,0.0)*(gammaSep-this->gammaInt_)*F_lim*max(3.0*Csep_*nu-nut,nu_dim)*Omega*S //  *(1.0-this->gammaInt_)*F_lim*max(3.0*Csep_*nu-nut,nu_dim)*Omega*S
    );

    return  this->gammaInt_*ssglrr_omega<BasicTurbulenceModel>::P() + P_lim*I;
}


template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::P2omega() const
{

    volScalarField ReThetac(max(this->ReThetac(),1.0));

    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();
    volScalarField S(sqrt(2*magSqr(symm(gradU))));
    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    volScalarField Rev(sqr(this->y_)*S/nu);

    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volScalarField RT(k/(nu*omega));

   const volScalarField fonset(Fonset(Rev, ReThetac, RT));

    return pow((1-fonset),m_blend_) //(1.0 - pow(mag(1.0 - this->gammaInt_),n_blend_)) //min(pow((1-F1(CD)),m_blend_)+(1.0 - pow(mag(1.0 - this->gammaInt_),n_blend_)),1.0) 

   				* ssglrr_omega<BasicTurbulenceModel>::P2omega();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::Pomega() const
{
    volScalarField ReThetac(max(this->ReThetac(),1.0));

    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();
    volScalarField S(sqrt(2*magSqr(symm(gradU))));
    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    volScalarField Rev(sqr(this->y_)*S/nu);

    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volScalarField RT(k/(nu*omega));

   const volScalarField fonset(Fonset(Rev, ReThetac, RT));

   return  (1 - pow((1-fonset),m_blend_)) //pow(mag(1.0 - this->gammaInt_),n_blend_) //max((1-pow((1.0-F1(CD)),m_blend_)) + pow(mag(1.0 - this->gammaInt_),n_blend_) - 1.0 ,0.0)

  			 * ssglrr_omega<BasicTurbulenceModel>::Pomega();
}


template<class BasicTurbulenceModel>
tmp<volSymmTensorField> gamma_ssglrr<BasicTurbulenceModel>::D() const
{

    return min(max(this->gammaInt_,0.1),1.0)*ssglrr_omega<BasicTurbulenceModel>::D();
}

template<class BasicTurbulenceModel>
tmp<volSymmTensorField> gamma_ssglrr<BasicTurbulenceModel>::Pi() const
{
    return  this->gammaInt_ * ssglrr_omega<BasicTurbulenceModel>::Pi();
}

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::Pi2() const
{
    return  this->gammaInt_* ssglrr_omega<BasicTurbulenceModel>::Pi2();
}

template<class BasicTurbulenceModel>
tmp<volScalarField>
gamma_ssglrr<BasicTurbulenceModel>::ReThetac() const
{
    // tmp<volScalarField> tReThetac
    // (
    //     new volScalarField
    //     (
    //         IOobject
    //         (
    //             IOobject::groupName("ReThetac", this->U_.group()),
    //             this->runTime_.timeName(),
    //             this->mesh_,
    //             IOobject::NO_READ,
    //             IOobject::AUTO_WRITE
    //         ),
    //         this->mesh_,
    //         dimless //dimensionedScalar("1e6",dimless,1e6)
    //     )
    // );

    // volScalarField ReThetac = tReThetac.ref();

    const volVectorField& U = this->U_;
    const volVectorField& n_ = (wallDist::New(this->mesh_).n());

    tmp<volScalarField> tnu = this->nu();
    volScalarField nu = tnu;

    const volScalarField& k = this->k_;
    const volScalarField& omega = this->omega_;
    const volScalarField& y = this->y_;

    const volScalarField dVdy((fvc::grad(U & n_) & n_));
    
    const volScalarField TUL(min(100.0*(sqrt(2*k/3)/(y*omega)),scalar(100)));

    volScalarField lambdaTL1 = ((-0.00757*dVdy*(sqr(y)/nu) + 0.0128));
    volScalarField lambdaTL = (min(max(lambdaTL1,-1.0),1.0));

    volScalarField FPG(lambdaTL);

    volScalarField rethetac(lambdaTL);

    forAll(FPG,celli)
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

        rethetac[celli] = CTU1_.value()+CTU2_.value()*exp(-CTU3_.value()*TUL[celli]*FPG[celli]);

    }

    // forAll(ReThetac.boundaryField(), patchi)
    // {
    //     fvPatchScalarField& pT = ReThetac.boundaryFieldRef()[patchi];

    //     forAll(pT, facei)
    //     {
    //         pT[facei] = max(pT[facei],1e6);
    //     }
    // }

    return ReThetac;
} 

template<class BasicTurbulenceModel>
tmp<volScalarField> gamma_ssglrr<BasicTurbulenceModel>::Fonset
(
    const volScalarField& Rev,
    const volScalarField& ReThetac,
    const volScalarField& RT
) const
{

    const volScalarField Fonset1(Rev/(Conset_*ReThetac));

    const volScalarField Fonset2
    (
        min(Fonset1, scalar(2.0))
    );

    const volScalarField Fonset3(max(1 - pow3(RT/3.5), scalar(0)));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject::groupName("Fonset", this->U_.group()),
            max(Fonset2 - Fonset3, scalar(0))
        )
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
gamma_ssglrr<BasicTurbulenceModel>::gamma_ssglrr
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
    n_blend_
        (
        dimensionedScalar::lookupOrAddToDict
        (
            "n_blend",
            this->coeffDict_,
            0.75
        )
    ),
    m_blend_
    (
        dimensionedScalar::lookupOrAddToDict
        (
           "m_blend",
           this->coeffDict_,
           0.75
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
    gammaIntEff_
    (
        IOobject
        (
            "gammaIntEff_",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0",dimless,0.0)
    ),
    P_Lim_
    (
        IOobject
        (
            "P_Lim",
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
bool gamma_ssglrr<BasicTurbulenceModel>::read()
{
    if (gamma_ssglrr<BasicTurbulenceModel>::read())
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
        n_blend_.readIfPresent(this->coeffDict());
        m_blend_.readIfPresent(this->coeffDict());
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
void gamma_ssglrr<BasicTurbulenceModel>::correctGammaInt()
{
    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    // const volVectorField& U = this->U_;
    // const volScalarField& k = this->k_;
    // const volScalarField& omega = this->omega_;
    // const tmp<volScalarField> tnu = this->nu();
    // const volScalarField::Internal& nu = tnu()();
    // const volScalarField::Internal& y = this->y_();
    fv::options& fvOptions(fv::options::New(this->mesh_));

    // Fields derived from the velocity gradient
    //tmp<volTensorField> tgradU = fvc::grad(U);
    //const volScalarField::Internal Omega(sqrt(2*magSqr(skew(tgradU()()))));
    //const volScalarField::Internal S(sqrt(2*magSqr(symm(tgradU()()))));
    //const volScalarField::Internal Us(max(mag(U()), deltaU_));
    //const volScalarField::Internal dUsds((U() & (U() & tgradU()()))/sqr(Us));
    //tgradU.clear();

    //const volScalarField::Internal ReThetac(this->ReThetac());
    //const volScalarField::Internal Rev(sqr(y)*S/nu);
    //const volScalarField::Internal RT(k()/(nu*omega()));

    const volScalarField ReThetac(this->ReThetac());

    P_Lim_ = ReThetac;

    {        
        // Intermittency equation
        tmp<fvScalarMatrix> gammaIntEqn
        (
            fvm::ddt(alpha, rho, gammaInt_)
          + fvm::div(alphaRhoPhi, gammaInt_)
          - fvm::laplacian(alpha*rho*DgammaInt(), gammaInt_)
        ==
            alpha*rho*this->Pgamma() - fvm::Sp(alpha*rho*this->Pgamma(), gammaInt_)
          + alpha*rho*this->Egamma() - fvm::Sp(ce2_*alpha*rho*this->Egamma(), gammaInt_)
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
void gamma_ssglrr<BasicTurbulenceModel>::correct()
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
