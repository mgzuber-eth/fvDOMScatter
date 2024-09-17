/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2008-2011 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
                            | Copyright (C) 2011-2017 OpenFOAM Foundation
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

inline const Foam::radiation::radiativeMGZIntensityRay&
Foam::radiation::fvMGZDOM::IRay(const label rayI) const
{
    return  IRay_[rayI];
}


inline const Foam::volScalarField&
Foam::radiation::fvMGZDOM::IRayLambda
(
    const label rayI,
    const label lambdaI
) const
{
    return IRay_[rayI].ILambda(lambdaI);
}


inline Foam::label Foam::radiation::fvMGZDOM::nTheta() const
{
    return nTheta_;
}


inline Foam::label Foam::radiation::fvMGZDOM::nPhi() const
{
    return nPhi_;
}


inline Foam::label Foam::radiation::fvMGZDOM::nRay() const
{
    return nRay_;
}


inline Foam::label Foam::radiation::fvMGZDOM::nLambda() const
{
    return nLambda_;
}


inline Foam::label Foam::radiation::fvMGZDOM::nBands() const
{
    return nLambda_;
}

inline const Foam::volScalarField& Foam::radiation::fvMGZDOM::a() const
{
    return a_;
}


inline const Foam::volScalarField& Foam::radiation::fvMGZDOM::aLambda
(
    const label lambdaI
) const
{
    return aLambda_[lambdaI];
}


inline const Foam::volScalarField& Foam::radiation::fvMGZDOM::G() const
{
    return G_;
}


inline const Foam::volScalarField& Foam::radiation::fvMGZDOM::qr() const
{
    return qr_;
}

inline const Foam::volScalarField& Foam::radiation::fvMGZDOM::qin() const
{
    return qin_;
}


inline const Foam::volScalarField& Foam::radiation::fvMGZDOM::qem() const
{
    return qem_;
}


inline const Foam::radiation::blackMGZBodyEmission&
Foam::radiation::fvMGZDOM::blackBody() const
{
    return blackBody_;
}


inline Foam::scalar Foam::radiation::fvMGZDOM::omegaMax() const
{
    return omegaMax_;
}


inline Foam::vector Foam::radiation::fvMGZDOM::meshOrientation() const
{
    return meshOrientation_;
}


inline bool Foam::radiation::fvMGZDOM::useSolarLoad() const
{
    return useSolarLoad_;
}

// ************************************************************************* //
