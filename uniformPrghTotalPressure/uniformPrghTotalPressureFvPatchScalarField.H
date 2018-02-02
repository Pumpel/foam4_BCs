/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::uniformPrghTotalPressureFvPatchScalarField

Group
    grpInletBoundaryConditions grpOutletBoundaryConditions

Description
    This boundary condition provides a time-varying form of the uniform total
    pressure boundary condition Foam::totalPressureFvPatchField.

    For incompressible and subsonic flow, the patch pressure is set while
    accounting for the hydrostatic pressure.

        \f[
            p_rgh = p - \rho g.(h - hRef)
        \f]

        \f[
            p = p0 - 0.5 \rho |U|^2
        \f]

    where
    \vartable
        p_rgh   | Pseudo hydrostatic pressure [Pa]
        p       | Static pressure [Pa]
        p0      | Total pressure [Pa]
        h       | Height in the opposite direction to gravity
        hRef    | Reference height in the opposite direction to gravity
        \rho    | Density
        g       | Acceleration due to gravity [m/s^2]
    \endtable

Usage
    \table
        Property     | Description                | Required    | Default value
        U            | Velocity field name        | no          | U
        phi          | Flux field name            | no          | phi
        rho          | Density field name         | no          | rho
        psi          | Compressibility field name | no          | none
        gamma        | (Cp/Cv)                    | no          | 1
        p0           | Total pressure as a function of time | yes |
        value	$internalField;
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            uniformTotalPressure;
        p0              uniform 1e5;
    }
    \endverbatim

    The \c p0 entry is specified as a Function1 type, able to describe
    time varying functions.

See also
    Foam::Function1Types
    Foam::uniformFixedValueFvPatchField
    Foam::totalPressureFvPatchField

SourceFiles
    uniformTotalPressureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef uniformTotalPressureFvPatchScalarField_H
#define uniformTotalPressureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class uniformTotalPressureFvPatchField Declaration
\*---------------------------------------------------------------------------*/

class uniformPrghTotalPressureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

        //- Name of the velocity field
        word UName_;

        //- Name of the flux transporting the field
        word phiName_;

        //- Name of the density field used to normalise the mass flux
        //  if neccessary
        word rhoName_;

        //- Name of the compressibility field used to calculate the wave speed
        word psiName_;

        //- Heat capacity ratio
        scalar gamma_;

        //- Table of time vs total pressure, including the bounding treatment
        autoPtr<Function1<scalar>> p0_;


public:

    //- Runtime type information
    TypeName("uniformPrghTotalPressure");


    // Constructors

        //- Construct from patch and internal field
    	uniformPrghTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
    	uniformPrghTotalPressureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
    	uniformPrghTotalPressureFvPatchScalarField
        (
            const uniformPrghTotalPressureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
    	uniformPrghTotalPressureFvPatchScalarField
        (
            const uniformPrghTotalPressureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformPrghTotalPressureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        uniformPrghTotalPressureFvPatchScalarField
        (
            const uniformPrghTotalPressureFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new uniformPrghTotalPressureFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Access

            //- Return the name of the velocity field
            const word& UName() const
            {
                return UName_;
            }

            //- Return reference to the name of the velocity field
            //  to allow adjustment
            word& UName()
            {
                return UName_;
            }

            //- Return the heat capacity ratio
            scalar gamma() const
            {
                return gamma_;
            }

            //- Return reference to the heat capacity ratio to allow adjustment
            scalar& gamma()
            {
                return gamma_;
            }


        // Evaluation functions

            //- Inherit updateCoeffs from fixedValueFvPatchScalarField
            using fixedValueFvPatchScalarField::updateCoeffs;

            //- Update the coefficients associated with the patch field
            //  using the given patch velocity field
            virtual void updateCoeffs(const vectorField& Up);

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
