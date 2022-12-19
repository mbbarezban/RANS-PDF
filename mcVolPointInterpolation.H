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
    Foam::mcVolPointInterpolation

Description
    Interpolate from cell centres to points (vertices) using inverse distance
    weighting

SourceFiles
    mcVolPointInterpolation.C
    mcVolPointInterpolate.C

\*---------------------------------------------------------------------------*/

#ifndef mcVolPointInterpolation_H
#define mcVolPointInterpolation_H

#include "scalarList.H"
#include "volFields.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class pointMesh;

/*---------------------------------------------------------------------------*\
                       Class mcVolPointInterpolation Declaration
\*---------------------------------------------------------------------------*/

class mcVolPointInterpolation
{
    // Private data

    // Boundary handling

    const fvMesh& mesh_;

    //- IDs for inlet/outlet/wall patches
    labelList patchIDs_;

    //- points on inlet/outlet/wall boundaries
    labelList InlOutWallBndPoints_;

    //- face labels of all points on inlet/outlet/wall patches
    labelListList InlOutWallBndPointFaces_;
    scalarListList InlOutWallBndPointWeights_;

    //- point labels of wedge/cyclic patches (except points on axis)
    labelList WedCycPointList1_;
    labelList WedCycPointList2_;

    //- axis points for both wedge/cyclic patches
    labelList axisPoints_;
    labelListList axisPointCells_;
    scalarListList axisPointWeights_;

    //- all inlet/outlet/wall/axis points
    labelList allBoundaryPoints_;

    // Private Member Functions

    //- Construct addressing over all boundary faces
    void calcBoundaryAddressing();

    //- Make weights for points on uncoupled patches
    void makeBoundaryWeights(scalarField& sumWeights);

    //- Construct all point weighting factors
    void makeWeights();

    //- Get boundary field in same order as boundary faces. Field is
    //  zero on all coupled and empty patches
    template<class Type>
    tmp<Field<Type>> flatBoundaryField
    (
            const GeometricField<Type, fvPatchField, volMesh>& vf
    ) const
    {
        const fvMesh& mesh = vf.mesh();
        const fvBoundaryMesh& bm = mesh.boundary();

        tmp<Field<Type>> tboundaryVals
        (
            new Field<Type>(mesh.nFaces(), Zero)
        );
        Field<Type>& boundaryVals = tboundaryVals.ref();

        forAll(patchIDs_, patchI)
        {
            label patchi = patchIDs_[patchI];
            const Field<Type>& vfB = vf.boundaryField()[patchi];
            const fvPatch& fvp = bm[patchi];

            label facei = fvp.start();
            forAll(fvp, fI)
            {
                boundaryVals[facei] = vfB[fI];
                facei++;
            }
        }

        return tboundaryVals;
    }

    //- Disallow default bitwise copy construct
    mcVolPointInterpolation(const mcVolPointInterpolation&);

    //- Disallow default bitwise assignment
    void operator=(const mcVolPointInterpolation&);


public:

    // Declare name of the class and its debug switch
    ClassName("mcVolPointInterpolation");


    // Constructors

    //- Constructor given fvMesh and pointMesh.
    explicit mcVolPointInterpolation(const fvMesh&);


    //- Destructor
    ~mcVolPointInterpolation();


    // Member functions

    // Interpolation functions

    //- Interpolate boundary field without applying constraints/boundary
    //  conditions
    void interpolateBoundaryScalarField
    (
            const GeometricField<scalar, fvPatchField, volMesh>& vf,
            GeometricField<scalar, pointPatchField, pointMesh>& pf
    ) const;

    void interpolateBoundaryTensorField
    (
            const GeometricField<symmTensor, fvPatchField, volMesh>& vf,
            GeometricField<symmTensor, pointPatchField, pointMesh>& pf
    ) const;

    void correctGradUpointField
    (
            GeometricField<tensor, pointPatchField, pointMesh>& pGradU
    ) const;

    void correctDivRpointField
    (
            GeometricField<vector, pointPatchField, pointMesh>& pDivR
    ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#endif
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
