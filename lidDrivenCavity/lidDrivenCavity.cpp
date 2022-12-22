#include "fileWriter.h"
#include "fvOperations.h"
#include "extrafunc.h"
#include "FiniteVolume.h"
#include "Equation.h"
#include "MatrixClass.h"
#include "SipSolSolver.h"
#include "Simple.h"
#include <memory>

using namespace std;

int main()
{

    int N_ = 5;
    int M_ = 5;
    double length1 = 1.0;
    Simple* useSimple = new Simple(N_, M_, length1);
    useSimple->assembleSolveMomentum();
    delete useSimple;
    useSimple = nullptr;

    //////setting valori iniziali
    //int N1 = 47;
    //double length1 = 1.0;
    //Grid myGrid_(N1, N1, length1);

    //Solution sol_;
    //sol_.dt = 0.001; sol_.R = 1.0; sol_.visc = 1e-03; sol_.density = 1.0;
    //sol_.nsteps = 1; sol_.maxit = 2;

    //int NI = myGrid_.pNI();
    //int NJ = myGrid_.pNJ();

    //Fields FieldOper;
    //Fields::vectorField U(NI, Fields::vec1dField(NJ));
    //Fields::vectorField UO(NI, Fields::vec1dField(NJ));
    //Fields::vectorField UOO(NI, Fields::vec1dField(NJ));
    //Fields::vectorField V(NI, Fields::vec1dField(NJ));
    //Fields::vectorField VO(NI, Fields::vec1dField(NJ));
    //Fields::vectorField VOO(NI, Fields::vec1dField(NJ));
    //Fields::vectorField P(NI, Fields::vec1dField(NJ));
    //Fields::vectorField PP(NI, Fields::vec1dField(NJ));
    //Fields::vectorField DPX(NI, Fields::vec1dField(NJ));
    //Fields::vectorField DPY(NI, Fields::vec1dField(NJ));
    //Fields::vectorField F1(NI, Fields::vec1dField(NJ));
    //Fields::vectorField F2(NI, Fields::vec1dField(NJ));

    //FieldOper.getGridInfoPassed(U, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(UO, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(UOO, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(V, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(VO, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(P, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(PP, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(DPX, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(DPY, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(F1, myGrid_, sol_);
    //FieldOper.getGridInfoPassed(F2, myGrid_, sol_);

    //FiniteMatrix::finiteMat AE(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat AW(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat AS(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat AN(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat AP(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat SU(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat SV(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat APU(NI, vector<FiniteMatrix>(NJ));
    //FiniteMatrix::finiteMat APV(NI, vector<FiniteMatrix>(NJ));

    ////variabili per la definizione del problema specifico
    //string EastB = "East"; string Northb = "North";
    //string WestB = "West"; string Southb = "South";
    //string timeMethod = "EULER2"; // either EULER2 or EULER3
    //double ULID = 1.0; // LID VELOCITY 

    //// Qui si inizia a risolvere i problemi
    //int Liter = 5;
    //SipSolSolver SIP(NI, NJ, Liter);

    //// initial boundary on Lid
    //FieldOper.dirichletBoundary(U, Northb, ULID);

    ////elementi che vengono modificati dentro le funzioni che le utilizzano
    //int pressureGEast = 1;
    //int pressureGNorth = 2;
    //int solveIter = 1;
    //double blendConvec = 1.0;
    //Equation* UEqn = nullptr;
    //Equation* VEqn = nullptr;
    //Equation* PEqn = nullptr;
    //for (int iter = 0; iter < 100; iter++)
    //{
    //    cout << " ITER -> " << iter + 1 << endl;

    //    // 1. ASSEMBLE AND SOLVE MOMENTUM EQUATION 

    //    FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, Northb);
    //    FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, Southb);
    //    FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, WestB);
    //    FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, EastB);

    //    blendConvec = 1.0;
    //    FiniteMatrix::finiteMat fmU = fvm::convectionTerm(U, F1, F2, blendConvec) + fvm::diffusiveTerm(U) + fvm::pressureGrad(P, pressureGEast);
    //    UEqn = new Equation(fmU);
    //    UEqn->noWallShearXBoundaryConditions(U);
    //    UEqn->assembleEquation();
    //    UEqn->relax(U);
    //    FiniteMatrix::finiteMat SU(UEqn->sourceRelaxed);

    //    Fields::vectorField Utemp = UEqn->solve(U, SU, sol_, solveIter);
    //    FieldOper.copyInternalField(U, Utemp);
    //    APU = UEqn->rAP;

    //    //FieldOper.print2dfield(U);
    //    //string timeMethod_ = "EULER3";
    //    FiniteMatrix::finiteMat fmV = fvm::convectionTerm(V, F1, F2, blendConvec) + fvm::diffusiveTerm(V) + fvm::pressureGrad(P, pressureGNorth);
    //    VEqn = new Equation(fmV);
    //    VEqn->EqnName = "V-Eqn";
    //    VEqn->noWallShearYBoundaryConditions(V);
    //    VEqn->assembleEquation();
    //    VEqn->relax(V);
    //    FiniteMatrix::finiteMat SV(VEqn->sourceRelaxed);
    //    solveIter = 1;
    //    Fields::vectorField Vtemp = VEqn->solve(V, SV, sol_, solveIter);
    //    FieldOper.copyInternalField(V, Vtemp);
    //    APV = VEqn->rAP;

    //    FieldOper.computeCellCenterPressureGrad(P, DPX, DPY);

    //    FiniteMatrix FiniteClass_;
    //    Fields::vectorField PgradPtoE(FieldOper.interpolatedFieldEast(DPX, myGrid_));    // DPX okay DPY okay    // PgradPtoE okay
    //    Fields::vectorField PgradPtoN(FieldOper.interpolatedFieldNorth(DPY, myGrid_));   //  okay
    //    Fields::vectorField UVelPtoE(FieldOper.interpolatedFieldEast(U, myGrid_));       // this seems okay
    //    Fields::vectorField VVelPtoN(FieldOper.interpolatedFieldNorth(V, myGrid_));              //okay
    //    FiniteMatrix::finiteMat APUPtoE(FiniteClass_.interpolatedFieldEast(APU, myGrid_));       //okay
    //    FiniteMatrix::finiteMat APVPtoN(FiniteClass_.interpolatedFieldNorth(APV, myGrid_));      // okay

    //    // Cell Face Gradients
    //    Fields::vectorField DPXE(FieldOper.cellFaceGradientEast(P, myGrid_));    //okay
    //    Fields::vectorField DPYN(FieldOper.cellFaceGradientNorth(P, myGrid_));   //okay now
    //    //FieldOper.print2dfield(test2);
    //    //FiniteClass_.print2dfield(APVPtoN);

    //    // corrected Cell Face Velocity
    //    Fields::vectorField UE(FiniteClass_.correctedFaceVelocityEast(UVelPtoE, DPXE, PgradPtoE, APUPtoE, myGrid_));    //okay
    //    Fields::vectorField VN(FiniteClass_.correctedFaceVelocityNorth(VVelPtoN, DPYN, PgradPtoN, APVPtoN, myGrid_));   //okay

    //    // 2. COMPUTE MASS FLUXES USING RHIE-CHOW INTERPOLATION
    //    FieldOper.computeEastMassFluxes(F1, UE);     //okay
    //    FieldOper.computeNorthMassFluxes(F2, VN);    //okay

    //    // 3. ASSEMBLE AND SOLVE PRESSURE CORRECTION EQUATION
    //    FiniteMatrix::finiteMat fmP = fvm::HTerm(APUPtoE, APVPtoN, U, V) + fvm::divPhi(F1, F2);
    //    PEqn = new Equation(fmP);
    //    PEqn->assembleEquation();
    //    forAllInternal(PP)
    //    {
    //        PP[i][j].value = 0.0;
    //    }
    //    FiniteMatrix::finiteMat SP(PEqn->sourceFinal);
    //    PEqn->EqnName = "PP-Eqn";
    //    PEqn->URF = 0.2;
    //    solveIter = 6;
    //    FiniteMatrix::finiteMat PAE(PEqn->AE);
    //    FiniteMatrix::finiteMat PAN(PEqn->AN);
    //    PP = PEqn->solve(PP, SP, sol_, solveIter);

    //    // Up to here now okay
    //    FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, Northb);
    //    FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, Southb);
    //    FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, WestB);
    //    FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, EastB);
    //    cout << " Pressure at Ref location 2,2 " << PP[4][4].value << endl;

    //    // 4. CORRECT FLUXES, U, V, P AND UPDATE VARIABLES FOR NEXT ITERATION
    //    FiniteClass_.correctEastMassFluxes(F1, PP, PAE);  //F1 okay
    //    FiniteClass_.correctNorthMassFluxes(F2, PP, PAN); //F1 okay
    //    // Finally correct U,V,P at Cell Center
    //    forAllInternal(U)
    //    {
    //        double DX = (myGrid_.X[i] - myGrid_.X[i - 1]);
    //        double DY = (myGrid_.Y[j] - myGrid_.Y[j - 1]);
    //        double RP = 1.0;

    //        double PPE = (PP[i + 1][j].value * myGrid_.FX[i]) + (PP[i][j].value * (1.0 - myGrid_.FX[i]));
    //        double PPW = (PP[i][j].value * myGrid_.FX[i - 1]) + (PP[i - 1][j].value * (1.0 - myGrid_.FX[i - 1]));
    //        double PPN = (PP[i][j + 1].value * myGrid_.FY[j]) + (PP[i][j].value * (1.0 - myGrid_.FY[j]));
    //        double PPS = (PP[i][j].value * myGrid_.FY[j - 1]) + (PP[i][j - 1].value * (1.0 - myGrid_.FY[j - 1]));

    //        U[i][j].value = U[i][j].value - (PPE - PPW) * DY * RP * APU[i][j].value;
    //        V[i][j].value = V[i][j].value - (PPN - PPS) * DX * RP * APV[i][j].value;
    //        P[i][j].value = P[i][j].value + sol_.URFPressure * (PP[i][j].value - PP[4][4].value);
    //    }


    //     fileWriter fileWriterr;
    //     string results="OOPLid";
    //     int time1=iter;

    //    if( time1%200 == 0)
    //    {
    //        fileWriterr.writeUVP(results,time1,myGrid_,U,V,P);
    //    }

    //    if (UEqn)
    //    {
    //        delete UEqn;
    //        UEqn = nullptr;
    //    }
    //    if (VEqn)
    //    {
    //        delete VEqn;
    //        VEqn = nullptr;
    //    }
    //    if (PEqn)
    //    {
    //        delete PEqn;
    //        PEqn = nullptr;
    //    }

    //}   // end outer loop
    //
    //if (UEqn)
    //{
    //    delete UEqn;
    //    UEqn = nullptr;
    //}
    //if (VEqn)
    //{
    //    delete VEqn;
    //    VEqn = nullptr;
    //}
    //if (PEqn)
    //{
    //    delete PEqn;
    //    PEqn = nullptr;
    //}


    return 0;
}
