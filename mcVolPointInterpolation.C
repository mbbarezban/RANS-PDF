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

\*---------------------------------------------------------------------------*/

#include "mcVolPointInterpolation.H"
#include "fvMesh.H"
#include "volFields.H"
#include "pointFields.H"
#include "demandDrivenData.H"
#include "surfaceFields.H"
#include "wedgePolyPatch.H"
#include "wallPolyPatch.H"
#include "cyclicPolyPatch.H"
#include "syncTools.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

    defineTypeNameAndDebug(mcVolPointInterpolation, 0);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mcVolPointInterpolation::calcBoundaryAddressing()
{
    const polyBoundaryMesh& pbm = mesh_.boundaryMesh();
    const pointField& points = mesh_.points();
    const faceUList& faces = mesh_.faces();
    const cellList& cells = mesh_.cells();

    boolList meshPointCheck(mesh_.nPoints(), false);
    labelListList meshPointFaces(mesh_.nPoints());

    //- assign all point labels on inlet/outlet/wall patches to InlOutWallBndPoints_
    // loop over all patch to finding inlet/outlet/wall patches
    forAll(pbm, patchI)
    {
        const polyPatch& pp = pbm[patchI];

        //- inlet/outlet patch or wall patch
        if( pbm.types()[patchI] == "patch" || isA<wallPolyPatch>(pp) )
        {
            patchIDs_.append(patchI);

            //- the first face of this patch
            label facei = pp.start();

            //- loop over all faces of this patch
            forAll(pp, i)
            {
                //- face object
                const face& thisFace = faces[facei];

                //- loop over all vertices of this face
                forAll(thisFace, verI)
                {
                    label pointi = thisFace[verI];
                    meshPointFaces[pointi].append(facei);

                    if(!meshPointCheck[pointi])
                    {
                        InlOutWallBndPoints_.append(pointi);
                        meshPointCheck[pointi] = true;
                    }
                }
                //- the next face of this patch
                facei++;
            }
        }
    }

    //- assign face labels of all points on inlet/outlet/wall patches to InlOutWallBndPointFaces_
    InlOutWallBndPointFaces_.setSize(InlOutWallBndPoints_.size());
    forAll(InlOutWallBndPointFaces_, pI)
    {
        labelList& pointFaces = InlOutWallBndPointFaces_[pI];
        label pointi = InlOutWallBndPoints_[pI];
        pointFaces.setSize(meshPointFaces[pointi].size());
        pointFaces = meshPointFaces[pointi];
    }

    //- for cases with wedge/empty patch
    if(mesh_.nGeometricD() <= 2)
    {
        labelHashSet wedgesID = pbm.findPatchIDs<wedgePolyPatch>();
        label wedge1 = wedgesID.toc()[0];
        label wedge2 = wedgesID.toc()[1];
        const polyPatch& pp1 = pbm[wedge1];
        const polyPatch& pp2 = pbm[wedge2];

        // const wedgePolyPatch& wpp =
        //        static_cast<const wedgePolyPatch&>(pp1);

        // vector axis = wpp.axis();

        //- fill the axisPoints_ and WedCycPointList1_
        //- loop over al points of this wedge patch
        forAll(pp1.meshPoints(), pointI)
        {
            label pointi = pp1.meshPoints()[pointI];

            //- if y component of this point is 0 then it is on the axis
            //-(for cases with axis along to the x direction)
            if((points[pointi][1]) == 0.)
            {
                axisPoints_.append(pointi);
            }
            else
                WedCycPointList1_.append(pointi);
        }

        //- fill the WedCycPointList2_
        WedCycPointList2_.setSize(WedCycPointList1_.size());

        forAll(WedCycPointList1_, pI)
        {
            label point1 = WedCycPointList1_[pI];
            const point& thisPoint1 = points[point1];

            forAll(pp2.meshPoints(), pointI)
            {
                label point2 = pp2.meshPoints()[pointI];
                const point& thisPoint2 = points[point2];

                //- if x and y component of these two points are equal then they are coresponding!
                //-(for cases with axis along to the x direction)
                if((thisPoint1[0] == thisPoint2[0]) && (thisPoint1[1] == thisPoint2[1]))
                    WedCycPointList2_[pI] = point2;
            }
        }

        //- fill the axisPointCells_
        axisPointCells_.setSize(axisPoints_.size());

        forAll(axisPoints_, pI)
        {
            label pointi = axisPoints_[pI];
            labelList& axisPointCells = axisPointCells_[pI];

            forAll(cells, cellI)
            {
                const cell& thisCell = cells[cellI];
                const labelList& cellVertices = thisCell.labels(faces);
                forAll(cellVertices, verI)
                {
                    if(pointi == cellVertices[verI])
                    {
                        axisPointCells.append(cellI);
                        break;
                    }
                }
            }
        }

    }
    //- for cases with cyclic patch
    else
    {
        labelHashSet cyclicID = pbm.findPatchIDs<cyclicPolyPatch>();
        label cyclic1 = cyclicID.toc()[0];
        label cyclic2 = cyclicID.toc()[1];
        const polyPatch& pp1 = pbm[cyclic1];
        const polyPatch& pp2 = pbm[cyclic2];

        // const cyclicPolyPatch& wpp =
        //        static_cast<const cyclicPolyPatch&>(pp1);

        // vector axis = wpp.rotationAxis();

        //- fill the axisPoints_ and WedCycPointList1_
        //- loop over al points of this cyclic patch
        forAll(pp1.meshPoints(), pointI)
        {
            label pointi = pp1.meshPoints()[pointI];

            //- if y component of this point is 0 then it is on the axis
            //-(for cases with axis along to the x direction)
            if((points[pointi][1]) == 0.)
            {
                axisPoints_.append(pointi);
            }
            else
                WedCycPointList1_.append(pointi);
        }

        //- fill the WedCycPointList2_
        WedCycPointList2_.setSize(WedCycPointList1_.size());
        forAll(WedCycPointList1_, pI)
        {
            label point1 = WedCycPointList1_[pI];
            const point& thisPoint1 = points[point1];

            forAll(pp2.meshPoints(), pointI)
            {
                label point2 = pp2.meshPoints()[pointI];
                const point& thisPoint2 = points[point2];

                //- if x and y component of these two points are equal then they are coresponding!
                //-(for cases with axis along to the x direction)
                if((thisPoint1[0] == thisPoint2[0]) && (thisPoint1[1] == thisPoint2[1]))
                    WedCycPointList2_[pI] = point2;
            }
        }

        //- fill the axisPointCells_
        axisPointCells_.setSize(axisPoints_.size());
        forAll(axisPoints_, pI)
        {
            label pointi = axisPoints_[pI];
            labelList& axisPointCells = axisPointCells_[pI];

            forAll(cells, cellI)
            {
                const cell& thisCell = cells[cellI];
                const labelList& cellVertices = thisCell.labels(faces);
                forAll(cellVertices, verI)
                {
                    if(pointi == cellVertices[verI])
                    {
                        axisPointCells.append(cellI);
                        break;
                    }
                }
            }
        }

    }

    //- fill the allBoundaryPoints_ (inlet + outlet + wall + axis points)
    allBoundaryPoints_.setSize(InlOutWallBndPoints_.size() + axisPoints_.size());
    forAll(InlOutWallBndPoints_, pI)
    {
        allBoundaryPoints_[pI] = InlOutWallBndPoints_[pI];
    }
    forAll(axisPoints_, pI)
    {
        allBoundaryPoints_[InlOutWallBndPoints_.size() + pI] = axisPoints_[pI];
    }
}

void Foam::mcVolPointInterpolation::makeBoundaryWeights
(
        scalarField& sumWeights
)
{
    const pointField& points = mesh_.points();
    const pointField& faceCentres = mesh_.faceCentres();
    const pointField& cellCentres = mesh_.cellCentres();

    axisPointWeights_.clear();
    axisPointWeights_.setSize(axisPoints_.size());

    //- make axis point weights
    forAll(axisPoints_, pI)
    {
        const label& pointi = axisPoints_[pI];
        const labelList& pCells = axisPointCells_[pI];

        scalarList& apw = axisPointWeights_[pI];
        apw.setSize(pCells.size());

        forAll(pCells, cI)
        {
            label celli = pCells[cI];

            apw[cI] = 1.0/mag(points[pointi] - cellCentres[celli]);
            sumWeights[pointi] += apw[cI];
        }
    }

    //- make inlet/outlet/wall point weights
    InlOutWallBndPointWeights_.clear();
    InlOutWallBndPointWeights_.setSize(InlOutWallBndPoints_.size());

    forAll(InlOutWallBndPoints_, pI)
    {
        const label& pointi = InlOutWallBndPoints_[pI];
        const labelList& pFaces = InlOutWallBndPointFaces_[pI];

        scalarList& bpw = InlOutWallBndPointWeights_[pI];
        bpw.setSize(pFaces.size());

        sumWeights[pointi] = 0.0;

        forAll(pFaces, fI)
        {
            label facei = pFaces[fI];

            bpw[fI] = 1.0/mag(points[pointi] - faceCentres[facei]);
            sumWeights[pointi] += bpw[fI];
        }
    }
}

void Foam::mcVolPointInterpolation::makeWeights()
{
    // Update addressing over all boundary faces
    calcBoundaryAddressing();

    // Running sum of weights
    pointScalarField sumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh_.polyMesh::instance(),
            mesh_
        ),
        pointMesh::New(mesh_),
        dimensionedScalar("zero", dimless, 0.)
    );

    // Create boundary weights; override sumWeights
    makeBoundaryWeights(sumWeights);

    //- Normalise inlet/outlet/wall boundary weights
    forAll(InlOutWallBndPoints_, pI)
    {
            const label& pointi = InlOutWallBndPoints_[pI];

        // Normalise boundary weights
            scalarList& bpw = InlOutWallBndPointWeights_[pI];
            forAll(bpw, i)
            {
                bpw[i] /= sumWeights[pointi];
            }
    }

    //- Normalise axis boundary weights
    forAll(axisPoints_, pI)
    {
            const label& pointi = axisPoints_[pI];

        // Normalise boundary weights
            scalarList& apw = axisPointWeights_[pI];
            forAll(apw, i)
            {
                apw[i] /= sumWeights[pointi];
            }
    }
}

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

Foam::mcVolPointInterpolation::mcVolPointInterpolation(const fvMesh& vm)
:
    mesh_(vm),
    patchIDs_(0),
    InlOutWallBndPoints_(0),
    InlOutWallBndPointFaces_(0),
    InlOutWallBndPointWeights_(0),
    WedCycPointList1_(0),
    WedCycPointList2_(0),
    axisPoints_(0),
    axisPointCells_(0),
    axisPointWeights_(0),
    allBoundaryPoints_(0)
{
    makeWeights();

    /*if(Pstream::myProcNo() == 0)
    {
        Pout<< "\npatchIDs: " << endl;
        Pout<< patchIDs_ << endl;
        Pout<< "\nboundaryPoints: " << endl;
        Pout<< InlOutWallBndPoints_ << endl;
        Pout<< "\nboundaryPointFaces: " << endl;
        Pout<< InlOutWallBndPointFaces_ << endl;
        Pout<< "\nboundaryPointWeights: " << endl;
        Pout<< InlOutWallBndPointWeights_ << endl;
        Pout<< "\nwpl1: " << endl;
        Pout<< WedCycPointList1_ << endl;
        Pout<< "\nwpl2: " << endl;
        Pout<< WedCycPointList2_ << endl;
        Pout<< "\naxisPoints: " << endl;
        Pout<< axisPoints_ << endl;
        Pout<< "\naxisPointCells: " << endl;
        Pout<< axisPointCells_ << endl;
        Pout<< "\naxisPointWeights: " << endl;
        Pout<< axisPointWeights_ << endl;
    }*/
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::mcVolPointInterpolation::~mcVolPointInterpolation()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::mcVolPointInterpolation::interpolateBoundaryScalarField
(
    const GeometricField<scalar, fvPatchField, volMesh>& vf,
    GeometricField<scalar, pointPatchField, pointMesh>& pf
) const
{
    Field<scalar>& pfi = pf.primitiveFieldRef();

    // Get face data in flat list
    tmp<Field<scalar>> tboundaryVals(flatBoundaryField(vf));
    const Field<scalar>& boundaryVals = tboundaryVals();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(InlOutWallBndPoints_, pI)
    {
        label pointi = InlOutWallBndPoints_[pI];

        const labelList& pFaces = InlOutWallBndPointFaces_[pI];
        const scalarList& pWeights = InlOutWallBndPointWeights_[pI];

        scalar& val = pfi[pointi];

        val = Zero;
        forAll(pFaces, j)
        {
            val += pWeights[j] * boundaryVals[pFaces[j]];
        }
    }

    if(WedCycPointList1_.size() > 0)
    {
        forAll(WedCycPointList1_, pI)
        {
            pfi[WedCycPointList1_[pI]] = 0.5 * ( pfi[WedCycPointList1_[pI]] + pfi[WedCycPointList2_[pI]] );
            pfi[WedCycPointList2_[pI]] = pfi[WedCycPointList1_[pI]];
        }
    }
}

void Foam::mcVolPointInterpolation::interpolateBoundaryTensorField
(
    const GeometricField<symmTensor, fvPatchField, volMesh>& vf,
    GeometricField<symmTensor, pointPatchField, pointMesh>& pf
) const
{
    Field<symmTensor>& pfi = pf.primitiveFieldRef();

    // Get face data in flat list
    tmp<Field<symmTensor>> tboundaryVals(flatBoundaryField(vf));
    const Field<symmTensor>& boundaryVals = tboundaryVals();

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    forAll(axisPoints_, pI)
    {
        label pointi = axisPoints_[pI];

        symmTensor& val = pfi[pointi];

        val.xy() = 0.;
        val.xz() = 0.;
        val.yz() = 0.;
    }

    forAll(InlOutWallBndPoints_, pI)
    {
        label pointi = InlOutWallBndPoints_[pI];

        const labelList& pFaces = InlOutWallBndPointFaces_[pI];
        const scalarList& pWeights = InlOutWallBndPointWeights_[pI];

        symmTensor& val = pfi[pointi];

        val = Zero;
        forAll(pFaces, j)
        {
            val += pWeights[j] * boundaryVals[pFaces[j]];
        }
    }

    if(WedCycPointList1_.size() > 0)
    {
        forAll(WedCycPointList1_, pI)
        {
            pfi[WedCycPointList1_[pI]] = 0.5 * ( pfi[WedCycPointList1_[pI]] + pfi[WedCycPointList2_[pI]] );
            pfi[WedCycPointList2_[pI]] = pfi[WedCycPointList1_[pI]];
        }
    }
}


void Foam::mcVolPointInterpolation::correctGradUpointField
(
        GeometricField<tensor, pointPatchField, pointMesh>& pGradU
) const
{
    Field<tensor>& pfi = pGradU.primitiveFieldRef();

    forAll(axisPoints_, pI)
    {
        label pointi = axisPoints_[pI];

        tensor& val = pfi[pointi];

        val.xy() = 0.;
        val.xz() = 0.;
        val.yx() = 0.;
        val.zx() = 0.;
    }
}

void Foam::mcVolPointInterpolation::correctDivRpointField
(
        GeometricField<vector, pointPatchField, pointMesh>& pDivR
) const
{
    Field<vector>& pfi = pDivR.primitiveFieldRef();

    forAll(axisPoints_, pI)
    {
        label pointi = axisPoints_[pI];

        vector& val = pfi[pointi];

        val[1] = 0.;
        val[2] = 0.;
    }
}
// ************************************************************************* //
