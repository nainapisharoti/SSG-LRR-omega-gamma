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
Class
    Foam::RASModels::ssglrr_omega
Group
    grpRASTurbulence
Description
    
    Blended SSG/LRR Reynolds stress turbulence model from Eisfeld

    Speziale, Sarkar and Gatski Reynolds-stress turbulence model for
    incompressible and compressible flows.
    Reference:
    \verbatim
        Speziale, C. G., Sarkar, S., & Gatski, T. B. (1991).
        Modelling the pressurestrain correlation of turbulence:
        an invariant dynamical systems approach.
        Journal of Fluid Mechanics, 227, 245-272.
    \endverbatim
    Including the generalized gradient diffusion model of
    Daly and Harlow:
    \verbatim
        Daly, B. J., & Harlow, F. H. (1970).
        Transport equations in turbulence.
        Physics of Fluids (1958-1988), 13(11), 2634-2649.
    \endverbatim
    The default model coefficients are:
    \verbatim
        ssglrr_omegaCoeffs
        {
            Cmu             0.09;
            C1              3.4;
            C1s             1.8;
            C2              4.2;
            C3              0.8;
            C3s             1.3;
            C4              1.25;
            C5              0.4;
            Ceps1           1.44;
            Ceps2           1.92;
            Cs              0.25;
            Ceps            0.15;
        }
    \endverbatim
SourceFiles
    ssglrr_omega_base.C
\*---------------------------------------------------------------------------*/

#ifndef ssglrr_omega_base_H
#define ssglrr_omega_base_H

//#include "RASModel.H"
//#include "ReynoldsStress.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ssglrr_omega Declaration
\*---------------------------------------------------------------------------*/

template<class TurbulenceModel, class BasicTurbulenceModel>
class ssglrr_omega
:
    public TurbulenceModel//ReynoldsStress<RASModel<BasicTurbulenceModel>>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        ssglrr_omega(const ssglrr_omega&);
        void operator=(const ssglrr_omega&);


protected:

    // Protected data

        // Model coefficients

            dimensionedScalar Cmu_;

            dimensionedScalar alpha_omega_;
            dimensionedScalar beta_omega_;
            dimensionedScalar sigma_omega_;
            dimensionedScalar sigma_d_omega_;
            dimensionedScalar C1_omega_;
            dimensionedScalar C1s_omega_;
            dimensionedScalar C2_omega_;
            dimensionedScalar C3_omega_;
            dimensionedScalar C3s_omega_;
            dimensionedScalar C2_LRR_;

            dimensionedScalar alpha_epsilon_;
            dimensionedScalar beta_epsilon_;
            dimensionedScalar sigma_epsilon_;
            dimensionedScalar sigma_d_epsilon_;
            dimensionedScalar C1_epsilon_;
            dimensionedScalar C1s_epsilon_;
            dimensionedScalar C2_epsilon_;
            dimensionedScalar C3_epsilon_;
            dimensionedScalar C3s_epsilon_;
            dimensionedScalar C4_epsilon_;
            dimensionedScalar C5_epsilon_;
            dimensionedScalar D_epsilon_;        

        // Fields

            volScalarField k_;
            volScalarField omega_;
            volScalarField epsilon_;

            volScalarField F1_;
            volScalarField D_;
            volScalarField Dreff_;

        //- Wall distance
            const volScalarField& y_;

    // Protected Member Functions

        virtual tmp<volScalarField> F1(const volScalarField& CDkOmega) const;

        //- Return R production rate
        virtual tmp<volSymmTensorField> P() const;

        //- Return R destruction rate
        virtual tmp<volSymmTensorField> D() const;

        //- Return Pressure strain correlation
        virtual tmp<volSymmTensorField> Pi() const;

        //- Return implicit part of pressure strain correlation;
        virtual tmp<volScalarField::Internal> Pi2() const;

        //- Return modified diffusion production term
        virtual tmp<volScalarField::Internal> Pomega() const;

        //- Return modified diffusion production term
        virtual tmp<volScalarField::Internal> P2omega() const;

        tmp<volScalarField> blend
        (
            const volScalarField& F1,
            const dimensionedScalar& psi1,
            const dimensionedScalar& psi2
        ) const
        {
             return F1*(psi1 - psi2) + psi2;
        }

        tmp<volScalarField> alphaOmega(const volScalarField& F1) const
        {
            return blend(F1, alpha_omega_, alpha_epsilon_);
        }
        tmp<volScalarField> betaOmega(const volScalarField& F1) const
        {
            return blend(F1, beta_omega_, beta_epsilon_);
        }
        tmp<volScalarField> sigmaOmega(const volScalarField& F1) const
        {
            return blend(F1, sigma_omega_, sigma_epsilon_);
        }
        tmp<volScalarField> sigmaD(const volScalarField& F1) const
        {
            return blend(F1, sigma_d_omega_, sigma_d_epsilon_);
        }
        tmp<volScalarField> C1(const volScalarField& F1) const
        {
            return blend(F1, C1_omega_, C1_epsilon_);
        }
        tmp<volScalarField> C1s(const volScalarField& F1) const
        {
            return blend(F1, C1s_omega_, C1s_epsilon_);
        }
        tmp<volScalarField> C2(const volScalarField& F1) const
        {
            return blend(F1, C2_omega_, C2_epsilon_);
        }
        tmp<volScalarField> C3(const volScalarField& F1) const
        {
            return blend(F1, C3_omega_, C3_epsilon_);
        }
        tmp<volScalarField> C3s(const volScalarField& F1) const
        {
            return blend(F1, C3s_omega_, C3s_epsilon_);
        }
        tmp<volScalarField> C4(const volScalarField& F1) const
        {
            dimensionedScalar C4_omega_ = 0.5*(18.0*C2_LRR_+12.0)/11.0;

            return blend(F1, C4_omega_, C4_epsilon_);
        }
        tmp<volScalarField> C5(const volScalarField& F1) const
        {
            dimensionedScalar C5_omega_ = 0.5*(-14.0*C2_LRR_+20.0)/11.0;
 
            return blend(F1, C5_omega_, C5_epsilon_);
        }
        tmp<volScalarField> D(const volScalarField& F1) const
        {
            dimensionedScalar D_omega_ = 0.75*Cmu_;

            return blend(F1, D_omega_, D_epsilon_);
        }

        //- Update the eddy-viscosity
        virtual void correctNut();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;

    // Constructors

        //- Construct from components
        ssglrr_omega
        (
            const word& type,
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName
        );


    //- Destructor
    virtual ~ssglrr_omega()
    {}


    // Member Functions

        //- Read model coefficients if they have changed
        virtual bool read();

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const
        {
            return k_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const
        {
            return epsilon_;
        }

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> omega() const
        {
            return omega_;
        }

        //- Return the effective diffusivity for R
        tmp<volSymmTensorField> DREff(const volScalarField& F1) const;

        //- Return the effective diffusivity for omega
        tmp<volScalarField> DomegaEff(const volScalarField& F1) const;

        //- Solve the turbulence equations and correct eddy-Viscosity and
        //  related properties
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "ssglrr_omega_base.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
