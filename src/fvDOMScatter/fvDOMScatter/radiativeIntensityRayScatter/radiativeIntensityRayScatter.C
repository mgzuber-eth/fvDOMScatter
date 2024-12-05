/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "radiativeIntensityRayScatter.H"
#include "fvm.H"
#include "fvDOMScatter.H"
#include "constants.H"

using namespace Foam::constant;

const Foam::word Foam::radiationModels::radiativeIntensityRayScatter::intensityPrefix
(
    "ILambda"
);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiationModels::radiativeIntensityRayScatter::radiativeIntensityRayScatter
(
    const fvDOMScatter& dom,
    const fvMesh& mesh,
    const scalar phi,
    const scalar theta,
    const scalar deltaPhi,
    const scalar deltaTheta,
    const label nLambda,
    const absorptionEmissionModel& absorptionEmission,
    const blackBodyEmissionScatter& blackBody,
    const label rayId
)
:
    dom_(dom),
    mesh_(mesh),
    absorptionEmission_(absorptionEmission),
    blackBody_(blackBody),
    I_
    (
        IOobject
        (
            "I" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qr_
    (
        IOobject
        (
            "qr" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qin_
    (
        IOobject
        (
            "qin" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    qem_
    (
        IOobject
        (
            "qem" + name(rayId),
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimMass/pow3(dimTime), 0)
    ),
    d_(Zero),
    dAve_(Zero),
    theta_(theta),
    phi_(phi),
    omega_(0.0),
    nLambda_(nLambda),
    ILambda_(nLambda),
    myRayId_(rayId)
{
    scalar sinTheta = Foam::sin(theta);
    scalar cosTheta = Foam::cos(theta);
    scalar sinPhi = Foam::sin(phi);
    scalar cosPhi = Foam::cos(phi);

    omega_ = 2.0*sinTheta*Foam::sin(deltaTheta/2.0)*deltaPhi;
    d_ = vector(sinTheta*sinPhi, sinTheta*cosPhi, cosTheta);
    dAve_ = vector
    (
        sinPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        cosPhi
       *Foam::sin(0.5*deltaPhi)
       *(deltaTheta - Foam::cos(2.0*theta)
       *Foam::sin(deltaTheta)),
        0.5*deltaPhi*Foam::sin(2.0*theta)*Foam::sin(deltaTheta)
    );

    // Transform directions so that they fall inside the bounds of reduced
    // dimension cases
    if (mesh_.nSolutionD() == 2)
    {
        vector meshDir(vector::zero);
        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (mesh_.geometricD()[cmpt] == -1)
            {
                meshDir[cmpt] = 1;
            }
        }
        const vector normal(vector(0, 0, 1));

        const tensor coordRot = rotationTensor(normal, meshDir);

        dAve_ = coordRot & dAve_;
        d_ = coordRot & d_;
    }
    else if (mesh_.nSolutionD() == 1)
    {
        vector meshDir(vector::zero);
        for (direction cmpt=0; cmpt<vector::nComponents; cmpt++)
        {
            if (mesh_.geometricD()[cmpt] == 1)
            {
                meshDir[cmpt] = 1;
            }
        }
        const vector normal(vector(1, 0, 0));

        dAve_ = (dAve_ & normal)*meshDir;
        d_ = (d_ & normal)*meshDir;
    }


    autoPtr<volScalarField> IDefaultPtr;

    forAll(ILambda_, lambdaI)
    {
        IOobject IHeader
        (
            intensityPrefix + "_" + name(rayId) + "_" + name(lambdaI),
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        );

        // Check if field exists and can be read
        if (IHeader.typeHeaderOk<volScalarField>(true))
        {
            ILambda_.set
            (
                lambdaI,
                new volScalarField(IHeader, mesh_)
            );
        }
        else
        {
            // Demand driven load the IDefault field
            if (!IDefaultPtr.valid())
            {
                IDefaultPtr.reset
                (
                    new volScalarField
                    (
                        IOobject
                        (
                            "IDefault",
                            mesh_.time().timeName(),
                            mesh_,
                            IOobject::MUST_READ,
                            IOobject::NO_WRITE
                        ),
                        mesh_
                    )
                );
            }

            // Reset the MUST_READ flag
            IOobject noReadHeader(IHeader);
            noReadHeader.readOpt() = IOobject::NO_READ;

            ILambda_.set
            (
                lambdaI,
                new volScalarField(noReadHeader, IDefaultPtr())
            );
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiationModels::radiativeIntensityRayScatter::~radiativeIntensityRayScatter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::scalar Foam::radiationModels::radiativeIntensityRayScatter::correct
(
	const scalar nRay_, //fvDOMScatter
	PtrList<radiativeIntensityRayScatter>& allRays, //fvDOMScatter
	const scalar currentRay //fvDOMScatter
)
{
    // Reset boundary heat flux to zero
    qr_.boundaryFieldRef() = 0.0;

    scalar maxResidual = -great;

    const surfaceScalarField Ji(dAve_ & mesh_.Sf());

	//fvDOMScatter Step 1: attentuation by scattering
	//make variable for scattering coefficient
	const dimensionedScalar sigmaScatteringCoefficient
	(
		"sigmaScatteringCoefficient",
		dimensionSet(0, -1, 0, 0, 0, 0, 0),
		scalar(0.0)//12.5
	);
	// end Step 1
	
	//fvDOMScatter Step 2: initialization
	//~ const dimensionedScalar augScatter
	//~ (
		//~ "augScatter",
		//~ dimensionSet(1,0,-3,0,0,0,0),
		//~ scalar(0.0)
	//~ );
	//dimensionedScalar augScatter;
	//end initialization
	
    forAll(ILambda_, lambdaI)
    {
        const volScalarField& k = dom_.aLambda(lambdaI);

		//fvDOMScatter Step 2: augmentation by scattering
		//loops stolen from fvDOM-Ru -- start
		// Compute total incident radiation within frequency band
		tmp<DimensionedField<scalar, volMesh>> augScatterLoop
		(
			allRays[0].ILambda(lambdaI)()*allRays[0].omega()
		);

		for (label rayI=1; rayI < nRay_; rayI++)
		{
			augScatterLoop.ref() += allRays[rayI].ILambda(lambdaI)()*allRays[rayI].omega();
		}
		
		//& is dot product
		// is magnitude of vector
		//muS = u&v / (magu *magv) //cosine of angle between two vectors
		//phaseFunction = 0.63*mus*mus - 1.43*mu + 0.79;
		
		volScalarField scatteringAugmentation ("scatteringAugmentation", 0.*blackBody_.bLambda(0));
		//scalar cosAngle;
		scalar muS;
		scalar phaseFunc;
		//scalar radians;
		for (label rayI=0; rayI < nRay_; rayI++)
		{
			//cosAngle = allRays[currentRay].dAve_ & allRays[rayI].dAve_;
			muS = (allRays[currentRay].dAve_ & allRays[rayI].dAve_)/(mag(allRays[currentRay].dAve_)*mag(allRays[rayI].dAve_));
			phaseFunc = 0.63*sqr(muS) + 1.43*muS + 0.79;
			scatteringAugmentation.ref() += allRays[rayI].ILambda(lambdaI)()*allRays[rayI].omega()*phaseFunc;
			//Info << "phaseFunc: " << phaseFunc << nl;
			//Info << "omega: " << allRays[rayI].omega() << nl;
			//radians = acos(muS);
		}
		
		//Info << "scatteringAugmentation: " << scatteringAugmentation << nl;
		//Info << "augScatterLoop: " << augScatterLoop << nl;
		//Info << "blackBody: " << blackBody_.bLambda(lambdaI) << nl;
		//Info << "Ilambdy: " << allRays[1].ILambda(0)() << nl;
		//dimensionedScalar augScatter
		//(
		//	sum(augScatterLoop)
		//);
		//const volDimensionedField& aug = augScatter();		
		//augScatter = sum(augScatterLoop);
		//loops stolen from fvDOM-Ru -- end
		
        fvScalarMatrix IiEq
        (
            fvm::div(Ji, ILambda_[lambdaI], "div(Ji,Ii_h)")
          + fvm::Sp(k*omega_, ILambda_[lambdaI])
		  + fvm::Sp(sigmaScatteringCoefficient*omega_, ILambda_[lambdaI]) //fvDOMScatter
        ==
            1.0/constant::mathematical::pi*omega_
           *(
                // Remove aDisp from k
                (k - absorptionEmission_.aDisp(lambdaI))
               *blackBody_.bLambda(lambdaI)

              + absorptionEmission_.E(lambdaI)/4
            )\
		  + sigmaScatteringCoefficient*omega_/(4.0*constant::mathematical::pi)*scatteringAugmentation //fvDOMScatter
        );

        IiEq.relax();

        const solverPerformance ILambdaSol = solve(IiEq, "Ii");

        const scalar initialRes =
            ILambdaSol.initialResidual()*omega_/dom_.omegaMax();

        maxResidual = max(initialRes, maxResidual);
    }

    return maxResidual;
}


void Foam::radiationModels::radiativeIntensityRayScatter::addIntensity()
{
    I_ = dimensionedScalar(dimMass/pow3(dimTime), 0);

    forAll(ILambda_, lambdaI)
    {
        I_ += ILambda_[lambdaI];
    }
}


// ************************************************************************* //
