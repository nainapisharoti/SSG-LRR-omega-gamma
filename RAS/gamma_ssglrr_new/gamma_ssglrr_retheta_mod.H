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

Class
    Foam::RASModels::gamma_ssglrr

Group
    grpLESTurbulence

Description
    Langtry-Menter 4-equation transitional SST model
    based on the k-omega-SST RAS model.

    References:
    \verbatim
        Langtry, R. B., & Menter, F. R. (2009).
        Correlation-based transition modeling for unstructured parallelized
        computational fluid dynamics codes.
        AIAA journal, 47(12), 2894-2906.

        Menter, F. R., Langtry, R., & Volker, S. (2006).
        Transition modelling for general purpose CFD codes.
        Flow, turbulence and combustion, 77(1-4), 277-303.

        Langtry, R. B. (2006).
        A correlation-based transition model using local variables for
        unstructured parallelized CFD codes.
        Phd. Thesis, Universität Stuttgart.
    \endverbatim

    The model coefficients are
    \verbatim
        gamma_ssglrrCoeffs
        {
            // Default SST coefficients
            alphaK1     0.85;
            alphaK2     1;
            alphaOmega1 0.5;
            alphaOmega2 0.856;
            beta1       0.075;
            beta2       0.0828;
            betaStar    0.09;
            gamma1      5/9;
            gamma2      0.44;
            a1          0.31;
            b1          1;
            c1          10;
            F3          no;

            // Default gammma_kOmegaSST coefficients
            ca2         0.06;
            ce2         50;
            Flength     100;
            cThetat     0.03;
            sigmaThetat 2;

            lambdaErr   1e-6;
            maxLambdaIter 10;
        }
    \endverbatim

SourceFiles
    gamma_ssglrr.C

\*---------------------------------------------------------------------------*/

#ifndef gamma_ssglrr_H
#define gamma_ssglrr_H

#include "ssglrr_omega.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                     Class gamma-ssglrr Declaration
\*---------------------------------------------------------------------------*/

template<class BasicTurbulenceModel>
class gamma_ssglrr
:
    public ssglrr_omega<BasicTurbulenceModel>
{
    // Private Member Functions

        // Disallow default bitwise copy construct and assignment
        gamma_ssglrr(const gamma_ssglrr&);
        void operator=(const gamma_ssglrr&);


protected:

    // Protected data

        // Model constants
            dimensionedScalar ca2_;
            dimensionedScalar ce2_;
            dimensionedScalar Flength_;
            dimensionedScalar ReThetac_lim_;
            
            dimensionedScalar CTU1_;
            dimensionedScalar CTU2_;
            dimensionedScalar CTU3_;

            dimensionedScalar CPG1_;
            dimensionedScalar CPG2_;
            dimensionedScalar CPG3_;
            dimensionedScalar Conset_;
            dimensionedScalar n_blend_;
            dimensionedScalar m_blend_;

            dimensionedScalar Ck_;
            dimensionedScalar Csep_;

            //- Stabilization for division by the magnitude of the velocity
            const dimensionedScalar deltaU_;


        // Fields

            //- Intermittency
            volScalarField gammaInt_;

            volScalarField gammaIntEff_;

            volScalarField P_Lim_;



    // Protected Member Functions

        //- Modified form of the k-omega SST F1 function
        virtual tmp<volScalarField> F1(const volScalarField& CD) const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volScalarField> gammaIntEff() const;

        //- gamma equation production term
        virtual tmp<volScalarField> Pgamma() const;

        //- gamma equation destruction term
        virtual tmp<volScalarField> Egamma() const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volSymmTensorField> P() const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volScalarField> Pomega() const;

        //- Modified form of the k-omega SST k production rate
        virtual tmp<volScalarField> P2omega() const;

        //- Modified form of the k-omega SST k destruction rate
        virtual tmp<volSymmTensorField> D() const;

        //- Modified form of the ssglrr_omega pressure strain correlation
        virtual tmp<volSymmTensorField> Pi() const;

        //- Modified form of the ssglrr_omega pressure strain correlation
        virtual tmp<volScalarField> Pi2() const;

        //- Empirical correlation for critical Reynolds number where the
        //  intermittency first starts to increase in the boundary layer
        virtual tmp<volScalarField> ReThetac() const;

        //- Transition onset location control function
        tmp<volScalarField> Fonset
        (
            const volScalarField& Rev,
            const volScalarField& ReThetac,
            const volScalarField& RT
        ) const;

        //- Solve the turbulence equations and correct the turbulence viscosity
        void correctGammaInt();


public:

    typedef typename BasicTurbulenceModel::alphaField alphaField;
    typedef typename BasicTurbulenceModel::rhoField rhoField;
    typedef typename BasicTurbulenceModel::transportModel transportModel;


    //- Runtime type information
    TypeName("gamma_ssglrr");


    // Constructors

        //- Construct from components
        gamma_ssglrr
        (
            const alphaField& alpha,
            const rhoField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~gamma_ssglrr()
    {}


    // Member Functions

        //- Re-read model coefficients if they have changed
        virtual bool read();

        //- Access function to intermittency
        const volScalarField& gammaInt() const
        {
            return gammaInt_;
        }

        //- Return the effective diffusivity for intermittency
        tmp<volScalarField> DgammaInt() const
        {
            return tmp<volScalarField>
            (
                new volScalarField
                (
                    "DgammaInt",
                    this->nut_ + this->nu()
                )
            );
        }

        //- Solve the turbulence equations and correct the turbulence viscosity
        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "gamma_ssglrr.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
