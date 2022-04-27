/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015-2016 OpenFOAM Foundation
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

#include "ssglrr_omega_base.H"
#include "fvOptions.H"
#include "wallFvPatch.H"
#include "bound.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::ssglrr_omega::F1
(
    const volScalarField& CD
) const
{
    tmp<volScalarField> CDPlus = (sigma_d_epsilon_/omega_)*max
    (
        CD,
        dimensionedScalar("1.0e-10", dimless/pow3(dimTime), 1.0e-10)
    );

    tmp<volScalarField> arg1 = min
    (
        min
        (
            max
            (
                (scalar(1)/Cmu_)*sqrt(k_)/(omega_*y_),
                (scalar(500)*this->nu())/(sqr(y_)*omega_)
            ),
            (4*sigma_epsilon_)*k_/(CDPlus*sqr(y_))
        ),
        scalar(10)
    );

    return tanh(pow4(arg1));
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volSymmTensorField>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::P() const
{
    const volVectorField& U = this->U_;
    const volSymmTensorField& R = this->R_;
    
     tmp<volTensorField> tgradU(fvc::grad(U));
     const volTensorField& gradU = tgradU();

    return -twoSymm(R & gradU);
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volSymmTensorField>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::D() const
{
    const volSymmTensorField& R = this->R_;
    const volScalarField epsilon = this->epsilon_;
    const volScalarField omega = this->omega_;

    return ((2.0/3.0)*epsilon*I)+((omega*R)-(omega*R));
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volSymmTensorField>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::Pi() const
{
    const volScalarField epsilon = this->epsilon_;
    const volScalarField omega = this->omega_;
    const volScalarField k = this->k_;
    const volVectorField& U = this->U_;
    const volSymmTensorField& R = this->R_;

    const volScalarField CD 
    (
        (fvc::grad(k) & fvc::grad(omega))
    );

    const volScalarField F1(this->F1(CD));

    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const volScalarField G(this->GName(), 0.5*mag(tr(this->P())));
    const volSymmTensorField b((R/this->k_)-(2.0/3.0)*I);
    const volSymmTensorField S(symm(gradU));
    const volTensorField Omega(skew(gradU));

    const volScalarField C1(this->C1(F1));
    const volScalarField C1s(this->C1s(F1));
    const volScalarField C2(this->C2(F1));
    const volScalarField C3(this->C3(F1));
    const volScalarField C3s(this->C3s(F1));
    const volScalarField C4(this->C4(F1));
    const volScalarField C5(this->C5(F1));

    return 
            (((2.0/3.0)*I)*(C1*epsilon+C1s*G) 
          + (C2*epsilon)*dev(innerSqr(b))
          + k*(
                   (C3 - C3s*mag(b))*dev(S)
                +  C4*dev(twoSymm(b&S))
                +  C5*twoSymm(b&Omega)
            )
            );
}


template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::Pi2() const
{
   // const volScalarField epsilon = this->epsilon_;
   // const volScalarField k = this->k_;
   // const volScalarField omega = this->omega_;
    const volScalarField G(this->GName(), 0.5*mag(tr(this->P())));

    const volScalarField CD 
    (
        (fvc::grad(k_) & fvc::grad(omega_))
    );

    const volScalarField F1(this->F1(CD));
    const volScalarField C1(this->C1(F1));
    const volScalarField C1s(this->C1s(F1));

    const volScalarField::Internal C(C1);
    const volScalarField::Internal Cs(C1s);

    return (C*epsilon_() + Cs*G)/k_();
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::Pomega() const
{

    //const volScalarField k = this->k_;
    //const volScalarField omega = this->omega_;
    const volScalarField G(this->GName(), 0.5*mag(tr(this->P())));   

    const volScalarField CD 
    (
        (fvc::grad(k_) & fvc::grad(omega_))
    );

     volScalarField F1(this->F1(CD));
     volScalarField::Internal alphaOmega(this->alphaOmega(F1));
   // volScalarField::Internal F(F1);

    return (1)*alphaOmega*G*omega_()/(k_());
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::P2omega() const
{

    //const volScalarField k = this->k_;
    //const volScalarField omega = this->omega_;
    const volVectorField& U = this->U_;
    tmp<volTensorField> tgradU(fvc::grad(U));
    const volTensorField& gradU = tgradU();

    const volScalarField::Internal S(sqrt(2*magSqr(symm(gradU))));

    const volScalarField G(this->GName(), 0.5*mag(tr(this->P())));   

    const volScalarField CD 
    (
        (fvc::grad(k_) & fvc::grad(omega_))
    );

    volScalarField F1(this->F1(CD));
    volScalarField::Internal alphaOmega(this->alphaOmega(F1));
    //volScalarField::Internal F(F1);

    return (1-1)*alphaOmega*sqr(S);
}

/*template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField::Internal>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::Pomega() const
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


template<class TurbulenceModel, class BasicTurbulenceModel>
void ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::correctNut()
{
	this->nut_ = this->k_/this->omega_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::ssglrr_omega
(
    const word& type,
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName
)
:
    ReynoldsStress<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    alpha_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha_omega",
            this->coeffDict_,
            0.5556
        )
    ),
    beta_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta_omega",
            this->coeffDict_,
            0.075
        )
    ),
    sigma_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigma_omega",
            this->coeffDict_,
            0.5
        )
    ),
    sigma_d_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigma_d_omega",
            this->coeffDict_,
            0.0
        )
    ),
    C1_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1_omega",
            this->coeffDict_,
            1.8
        )
    ),
    C1s_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1s_omega",
            this->coeffDict_,
            0.0
        )
    ),
    C2_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2_omega",
            this->coeffDict_,
            0.0
        )
    ),
    C3_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3_omega",
            this->coeffDict_,
            0.8
        )
    ),
    C3s_omega_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3s_omega",
            this->coeffDict_,
            0.0
        )
    ),
    C2_LRR_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2_LRR",
            this->coeffDict_,
            0.52
        )
    ),
    alpha_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alpha_epsilon",
            this->coeffDict_,
            0.44
        )
    ),
    beta_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta_epsilon",
            this->coeffDict_,
            0.0828
        )
    ),
    sigma_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigma_epsilon",
            this->coeffDict_,
            0.856
        )
    ),
    sigma_d_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigma_d_epsilon",
            this->coeffDict_,
            1.712
        )
    ),
    C1_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1_epsilon",
            this->coeffDict_,
            1.7
        )
    ),
    C1s_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C1s_epsilon",
            this->coeffDict_,
            0.9
        )
    ),
    C2_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2_epsilon",
            this->coeffDict_,
            1.05
        )
    ),
    C3_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3_epsilon",
            this->coeffDict_,
            0.8
        )
    ),
    C3s_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C3s_epsilon",
            this->coeffDict_,
            0.65
        )
    ),
    C4_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C4_epsilon",
            this->coeffDict_,
            0.625
        )
    ),
    C5_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C5_epsilon",
            this->coeffDict_,
            0.2
        )
    ),
    D_epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "D_epsilon",
            this->coeffDict_,
            0.22
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        0.5*tr(this->R_)
    ),
    omega_
    (
        IOobject
        (
            "omega",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            "epsilon",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        Cmu_*this->k_*this->omega_
    ),
    F1_
    (
        IOobject
        (
            "F1",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0",dimless,0.0)
    ),
    D_
    (
        IOobject
        (
            "D",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0",dimless,0.0)
    ),
    Dreff_
    (
        IOobject
        (
            "Dreff",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("0.0",dimless,0.0)
    ),

    y_(wallDist::New(this->mesh_).y())
{

        this->boundNormalStress(this->R_);
        bound(omega_, this->omegaMin_);
        k_ = 0.5*tr(this->R_);
        epsilon_ = Cmu_*this->k_*this->omega_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class TurbulenceModel, class BasicTurbulenceModel>
bool ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::read()
{
    if (ReynoldsStress<RASModel<BasicTurbulenceModel>>::read())
    {
        Cmu_.readIfPresent(this->coeffDict());

        alpha_omega_.readIfPresent(this->coeffDict());
        beta_omega_.readIfPresent(this->coeffDict());
        sigma_omega_.readIfPresent(this->coeffDict());
        sigma_d_omega_.readIfPresent(this->coeffDict());
        C1_omega_.readIfPresent(this->coeffDict());
        C1s_omega_.readIfPresent(this->coeffDict());
        C2_omega_.readIfPresent(this->coeffDict());
        C3_omega_.readIfPresent(this->coeffDict());
        C3s_omega_.readIfPresent(this->coeffDict());
        C2_LRR_.readIfPresent(this->coeffDict());

        alpha_epsilon_.readIfPresent(this->coeffDict());
        beta_epsilon_.readIfPresent(this->coeffDict());
        sigma_epsilon_.readIfPresent(this->coeffDict());
        sigma_d_epsilon_.readIfPresent(this->coeffDict());
        C1_epsilon_.readIfPresent(this->coeffDict());
        C1s_epsilon_.readIfPresent(this->coeffDict());
        C2_epsilon_.readIfPresent(this->coeffDict());
        C3_epsilon_.readIfPresent(this->coeffDict());
        C3s_epsilon_.readIfPresent(this->coeffDict());
        C4_epsilon_.readIfPresent(this->coeffDict());
        C5_epsilon_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volSymmTensorField>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::DREff(const volScalarField& F1) const
{
    volScalarField D(this->D(F1));

    return tmp<volSymmTensorField>
    (
        new volSymmTensorField
        (
            "DREff",
            (D*(this->k_/this->epsilon_))*this->R_ + I*this->nu()
        )
    );
}

template<class TurbulenceModel, class BasicTurbulenceModel>
tmp<volScalarField>
ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::DomegaEff(const volScalarField& F1)
const
{
    volScalarField sigmaOmega(this->sigmaOmega(F1));

    return tmp<volScalarField>
    (
        new volScalarField
        (
            "DomegaEff",
            sigmaOmega*(this->k_/this->omega_) + this->nu()
        )
    );
}

template<class TurbulenceModel, class BasicTurbulenceModel>
void ssglrr_omega<TurbulenceModel, BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    volSymmTensorField& R = this->R_;
    fv::options& fvOptions(fv::options::New(this->mesh_));

    ReynoldsStress<RASModel<BasicTurbulenceModel>>::correct();

     tmp<volTensorField> tgradU(fvc::grad(U));
     const volTensorField& gradU = tgradU();

    volSymmTensorField P(-twoSymm(R & gradU));
    volScalarField G(this->GName(), 0.5*mag(tr(this->P())));

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField CD 
    (
        (fvc::grad(k_) & fvc::grad(omega_))
    );

    volScalarField CDomega
    (
        max(CD, dimensionedScalar("1.0e-10", dimless/pow3(dimTime), 1.0e-10))/omega_
    );

    volScalarField F1(this->F1(CD));

    F1_ = F1;


    dimensionedScalar P_dim = dimensionedScalar("1.0",pow3(dimTime)/sqr(dimLength),1.0);

    volScalarField Dest(tr(this->D())*P_dim);

    D_ = Dest;

    dimensionedScalar D_dim = dimensionedScalar("1.0",dimTime/sqr(dimLength),1.0);

    volScalarField Diff(tr(this->DREff(F1))*D_dim);

    Dreff_ = Diff;

    {
	    volScalarField alphaOmega(this->alphaOmega(F1));
	    volScalarField betaOmega(this->betaOmega(F1));
	    volScalarField sigmaD(this->sigmaD(F1));

	    // Turbulence specific dissipation rate equation - blended
	    tmp<fvScalarMatrix> omegaEqn
	    (
	        fvm::ddt(alpha, rho, omega_)
	      + fvm::div(alphaRhoPhi, omega_)
	      - fvm::laplacian(alpha*rho*DomegaEff(F1), omega_)
	     ==
	        alpha()*rho()*(Pomega() + P2omega())
	      - fvm::Sp(betaOmega*alpha*rho*omega_, omega_)
	      + fvm::SuSp(alpha*rho*sigmaD*CDomega()/omega_(),omega_)
	      + fvOptions(alpha, rho, omega_)
	    );

	    omegaEqn.ref().relax();
	    fvOptions.constrain(omegaEqn.ref());
	    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
	    solve(omegaEqn);
	    fvOptions.correct(omega_);
	    bound(omega_, this->omegaMin_);
    }

    // Correct the trace of the tensorial production to be consistent
    // with the near-wall generation from the wall-functions
    const fvPatchList& patches = this->mesh_.boundary();

    forAll(patches, patchi)
    {
        const fvPatch& curPatch = patches[patchi];

        if (isA<wallFvPatch>(curPatch))
        {
            forAll(curPatch, facei)
            {
                label celli = curPatch.faceCells()[facei];
                P[celli] *= min
                (
                    G[celli]/(0.5*mag(tr(P[celli])) + SMALL),
                    1.0
                );
            }
        }
    }

	{	    
	    // Reynolds stress equation - Blended
	    tmp<fvSymmTensorMatrix> REqn
	    (
	        fvm::ddt(alpha, rho, R)
	      + fvm::div(alphaRhoPhi, R)
	      - fvm::laplacian(alpha*rho*DREff(F1), R)
          + fvm::Sp(this->Pi2(), R)
	     ==
	        alpha()*rho()*this->P()
          - alpha()*rho()*this->D()
          + alpha()*rho()*this->Pi()
	      + fvOptions(alpha, rho, R)
	    );
	    

	    REqn.ref().relax();
	    fvOptions.constrain(REqn.ref());
	    solve(REqn);
	    fvOptions.correct(R);

	    this->boundNormalStress(R);
	}

    k_ = 0.5*tr(R);

    epsilon_ = Cmu_*this->k_*this->omega_;

    correctNut();

    // Correct wall shear-stresses when applying wall-functions
    this->correctWallShearStress(R);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
