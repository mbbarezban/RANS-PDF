/*---------------------------------------------------------------------------*\
                pdfFoam: General Purpose PDF Solution Algorithm
                   for Reactive Flow Simulations in OpenFOAM

 Copyright (C) 2012 M.B.Barezban, Heng Xiao, Patrick Jenny,
                    Institute of Fluid Dynamics, ETH Zurich
-------------------------------------------------------------------------------
License
    This file is part of pdfFoam.

    pdfFoam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 2 of the License, or
    (at your option) version 3 of the same License.

    pdfFoam is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with pdfFoam.  If not, see <http://www.gnu.org/licenses/>.

Application
    reactPdfSimpleFoam

Description
    a consistent hybrid RANS/flactuating velocity-turbulent frecuency-composition
    joint PDF method for simulation of turbulent reacting flow.

Author
    M.B.Barezban

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiReactionThermo.H"
#include "mcCompTurbulenceModel.H"
#include "CombustionModel.H"
#include "radiationModel.H"
#include "reactPdfCloud.H"
#include "simpleControl.H"
#include "pressureControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"

    Info<< "\nStarting time loop\n" << endl;

    // starting of moving time average for FV fields including U.
    label j0MTA = pdfCloud.solution().j0MTA();
    scalar NTA  = 5;
    scalar alpha = 1. / NTA;

    while (simple.loop(runTime))
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // ************************** FV inner loop (7 loops) **********************************
        // calculate relaxation factor for U
        label j = pdfCloud.solution().iter();

        label numFVLoop = pdfCloud.solution().numFVLoop();  // -number of FV loop
        for(int i = 1; i <= numFVLoop; i++)
        {
            Info<< "simple iteration: " << i << endl;

            #include "UEqn.H" // update U from phi, rho, p and R/divR field
            #include "pEqn.H" // update p, phi, rho = p / RT from U field

            // time averaging for U
            if( j > j0MTA)
                U.relax(alpha);

        } // end FV loop

        // ************************** MC inner loop (3 Loops) ***********************************

        //- calculate R, h, Y, CMomega, T, psi, mu, alpha
        pdfCloud.solve(); // update T, Y, rho = p / RT from U and P field

        // **************************************************************************************

        turbulence->correct();  // update k, divR = div(rho , R) and nut = Cmu*k/CMomega

        Info<< "min/max(T) = "
            << min(T).value() << ", " << max(T).value() << endl;

        //- update heat release rate Qdot = RR * Hc
        //        Qdot = reaction->Qdot();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;

    } // end runTime.run()

    Info<< "End\n" << endl;

    return 0;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
