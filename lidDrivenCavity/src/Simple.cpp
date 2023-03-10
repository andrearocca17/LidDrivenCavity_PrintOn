#include "Simple.h"

Simple::Simple(int& N_, int& M_, double& l_):myGrid_(N_, M_, l_)
{
    NI = myGrid_.pNI();
    NJ = myGrid_.pNJ();
    U   = Fields::vectorField(NI, Fields::vec1dField(NJ));
    UO  = Fields::vectorField(NI, Fields::vec1dField(NJ));
    UOO = Fields::vectorField(NI, Fields::vec1dField(NJ));
    V   = Fields::vectorField(NI, Fields::vec1dField(NJ));
    VO  = Fields::vectorField(NI, Fields::vec1dField(NJ));
    VOO = Fields::vectorField(NI, Fields::vec1dField(NJ));
    P   = Fields::vectorField(NI, Fields::vec1dField(NJ));
    PP  = Fields::vectorField(NI, Fields::vec1dField(NJ));
    DPX = Fields::vectorField(NI, Fields::vec1dField(NJ));
    DPY = Fields::vectorField(NI, Fields::vec1dField(NJ));
    F1  = Fields::vectorField(NI, Fields::vec1dField(NJ));
    F2  = Fields::vectorField(NI, Fields::vec1dField(NJ));
    AE = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    AW = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    AS = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    AN = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    AP = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    SU = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    SV = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    APU = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    APV = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    APUPtoE = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    APVPtoN = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    PAE = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
    PAN = FiniteMatrix::finiteMat(NI, vector<FiniteMatrix>(NJ));
 }

Simple::~Simple()
{
    
}

void Simple::passGrid() {
    sol_.dt = 0.001; sol_.R = 1.0; sol_.visc = 1e-03; sol_.density = 1.0;
    sol_.nsteps = 1; sol_.maxit = 2;
   
    FieldOper.getGridInfoPassed(U, myGrid_, sol_);
    FieldOper.getGridInfoPassed(UO, myGrid_, sol_);
    FieldOper.getGridInfoPassed(UOO, myGrid_, sol_);
    FieldOper.getGridInfoPassed(V, myGrid_, sol_);
    FieldOper.getGridInfoPassed(VO, myGrid_, sol_);
    FieldOper.getGridInfoPassed(P, myGrid_, sol_);
    FieldOper.getGridInfoPassed(PP, myGrid_, sol_);
    FieldOper.getGridInfoPassed(DPX, myGrid_, sol_);
    FieldOper.getGridInfoPassed(DPY, myGrid_, sol_);
    FieldOper.getGridInfoPassed(F1, myGrid_, sol_);
    FieldOper.getGridInfoPassed(F2, myGrid_, sol_);

}

void Simple::resolveU() {
    FiniteMatrix::finiteMat fmU = fvm::convectionTerm(U, F1, F2, blendConvec) + fvm::diffusiveTerm(U) + fvm::pressureGrad(P, pressureGEast);
    UEqn = new Equation(fmU);
    UEqn->noWallShearXBoundaryConditions(U);
    UEqn->assembleEquation();
    UEqn->relax(U);
    FiniteMatrix::finiteMat SU(UEqn->sourceRelaxed);
    solveIter = 1;
    Fields::vectorField Utemp = UEqn->solve(U, SU, sol_, solveIter);
    FieldOper.copyInternalField(U, Utemp);
    APU = UEqn->rAP;

}
void Simple::resolveV() {
    FiniteMatrix::finiteMat fmV = fvm::convectionTerm(V, F1, F2, blendConvec) + fvm::diffusiveTerm(V) + fvm::pressureGrad(P, pressureGNorth);
    VEqn = new Equation(fmV);
    VEqn->EqnName = "V-Eqn";
    VEqn->noWallShearYBoundaryConditions(V);
    VEqn->assembleEquation();
    VEqn->relax(V);
    FiniteMatrix::finiteMat SV(VEqn->sourceRelaxed);
    solveIter = 1;
    Fields::vectorField Vtemp = VEqn->solve(V, SV, sol_, solveIter);
    FieldOper.copyInternalField(V, Vtemp);
    APV = VEqn->rAP;
}
void Simple::resolveP() {
    FiniteMatrix::finiteMat fmP = fvm::HTerm(APUPtoE, APVPtoN, U, V) + fvm::divPhi(F1, F2);
    PEqn = new Equation(fmP);
    PEqn->assembleEquation();
    forAllInternal(PP)
    {
        PP[i][j].value = 0.0;
    }
    FiniteMatrix::finiteMat SP(PEqn->sourceFinal);
    PEqn->EqnName = "PP-Eqn";
    PEqn->URF = 0.2;
    solveIter = 6;
    PAE = PEqn->AE;
    PAN = PEqn->AN;
    FiniteMatrix::finiteMat PAW = PEqn->AW;
    FiniteMatrix::finiteMat PAS = PEqn->AS;
    FiniteMatrix::finiteMat PAP = PEqn->AP;


    //FM.print2dfield(PAP);
    PP = PEqn->solve(PP, SP, sol_, solveIter);
    
}

void Simple::assembleSolveMomentum()
{
    passGrid();
    FieldOper.dirichletBoundary(U, Northb, ULID);

    for (int iter = 1; iter < 3; iter++)
    {
        cout << " ITER -> " << iter  << endl;

        // 1. ASSEMBLE AND SOLVE MOMENTUM EQUATION 
        FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, Northb);
        FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, Southb);
        FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, WestB);
        FieldOper.extaraPolateZeroGrad(P, myGrid_.FX, myGrid_.FY, EastB);
        cout << endl;
        cout << "------------------- P -------------------- " << endl;
        FieldOper.print2dfield(P);
    
        //FieldOper.print2dfield(P.FX);
        resolveU();
        resolveV();
        cout << endl;
        cout << "------------------- U-------------------- "<< endl;
        FieldOper.print2dfield(U);
        cout << endl;
        cout << "------------------- V-------------------- " << endl;
        FieldOper.print2dfield(V);
        cout << endl;

        FieldOper.computeCellCenterPressureGrad(P, DPX, DPY);

        FiniteMatrix FiniteClass_;
        Fields::vectorField PgradPtoE(FieldOper.interpolatedFieldEast(DPX, myGrid_));    // DPX okay DPY okay    // PgradPtoE okay
        Fields::vectorField PgradPtoN(FieldOper.interpolatedFieldNorth(DPY, myGrid_));   //  okay
        Fields::vectorField UVelPtoE(FieldOper.interpolatedFieldEast(U, myGrid_));       // this seems okay
        Fields::vectorField VVelPtoN(FieldOper.interpolatedFieldNorth(V, myGrid_));              //okay
        APUPtoE=FiniteClass_.interpolatedFieldEast(APU, myGrid_);       //okay
        APVPtoN=FiniteClass_.interpolatedFieldNorth(APV, myGrid_);      // okay

        // Cell Face Gradients
        Fields::vectorField DPXE(FieldOper.cellFaceGradientEast(P, myGrid_));    //okay
        Fields::vectorField DPYN(FieldOper.cellFaceGradientNorth(P, myGrid_));   //okay now
   
        // corrected Cell Face Velocity
        Fields::vectorField UE(FiniteClass_.correctedFaceVelocityEast(UVelPtoE, DPXE, PgradPtoE, APUPtoE, myGrid_));    //okay
        Fields::vectorField VN(FiniteClass_.correctedFaceVelocityNorth(VVelPtoN, DPYN, PgradPtoN, APVPtoN, myGrid_));   //okay

        // 2. COMPUTE MASS FLUXES USING RHIE-CHOW INTERPOLATION
        FieldOper.computeEastMassFluxes(F1, UE);     //okay
        FieldOper.computeNorthMassFluxes(F2, VN);    //okay
        cout << "-------------------F1-------------------- " << endl;
        FieldOper.print2dfield(F1);
        cout << "-------------------F2-------------------- " << endl;
        FieldOper.print2dfield(F2);
        cout << "----------------------------------------- " << endl;

        // 3. ASSEMBLE AND SOLVE PRESSURE CORRECTION EQUATION
        resolveP();
        cout << "-------------------PP-------------------- " << endl;
        FieldOper.print2dfield(PP);
        cout << endl;
        FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, Northb);
        FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, Southb);
        FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, WestB);
        FieldOper.extaraPolateZeroGrad(PP, myGrid_.FX, myGrid_.FY, EastB);
        cout << " Pressure at Ref location 2,2 " << PP[4][4].value << endl;
        cout << "-------------------PP extrapolate-------------------- " << endl;
        FieldOper.print2dfield(PP);
        cout << endl;
        //FieldOper.print2dfield(PP);
        // 4. CORRECT FLUXES, U, V, P AND UPDATE VARIABLES FOR NEXT ITERATION
        FiniteClass_.correctEastMassFluxes(F1, PP, PAE);  //F1 okay
        FiniteClass_.correctNorthMassFluxes(F2, PP, PAN); //F1 okay
        cout << "-------------------F1 NEW-------------------- " << endl;
        FieldOper.print2dfield(F1);
        cout << "-------------------F2 NEW-------------------- " << endl;
        FieldOper.print2dfield(F2);
        cout << "----------------------------------------- " << endl;
        cout << "-------------------P BEFORE CORRECTION-------------------- " << endl;
        FieldOper.print2dfield(P);
        // Finally correct U,V,P at Cell Center
        forAllInternal(U)
        {
            double DX = (myGrid_.X[i] - myGrid_.X[i - 1]);
            double DY = (myGrid_.Y[j] - myGrid_.Y[j - 1]);
            double RP = 1.0;

            double PPE = (PP[i + 1][j].value * myGrid_.FX[i]) + (PP[i][j].value * (1.0 - myGrid_.FX[i]));
            double PPW = (PP[i][j].value * myGrid_.FX[i - 1]) + (PP[i - 1][j].value * (1.0 - myGrid_.FX[i - 1]));
            double PPN = (PP[i][j + 1].value * myGrid_.FY[j]) + (PP[i][j].value * (1.0 - myGrid_.FY[j]));
            double PPS = (PP[i][j].value * myGrid_.FY[j - 1]) + (PP[i][j - 1].value * (1.0 - myGrid_.FY[j - 1]));

            U[i][j].value = U[i][j].value - (PPE - PPW) * DY * RP * APU[i][j].value;
            V[i][j].value = V[i][j].value - (PPN - PPS) * DX * RP * APV[i][j].value;
            P[i][j].value = P[i][j].value + sol_.URFPressure * (PP[i][j].value);// -PP[4][4].value);
        }
        cout << "-------------------U NEW-------------------- " << endl;
        FieldOper.print2dfield(U);
        cout << "-------------------V NEW-------------------- " << endl;
        FieldOper.print2dfield(V);
        cout << "-------------------P NEW-------------------- " << endl;
        FieldOper.print2dfield(P);

     
        int time1 = iter;

     /*   if (time1 % 200 == 0)
        {
            fileWriterr.writeUVP(results, time1, myGrid_, U, V, P);
        }*/
   /*     if (fileWriterr)
        {
            delete fileWriterr;
            fileWriterr = nullptr;
        }*/


       if (UEqn)
        {
            delete UEqn;
            UEqn = nullptr;
        }
        if (VEqn)
        {
            delete VEqn;
            VEqn = nullptr;
        }
        if (PEqn)
        {
            delete PEqn;
            PEqn = nullptr;
        }

    }   // end outer loop

}