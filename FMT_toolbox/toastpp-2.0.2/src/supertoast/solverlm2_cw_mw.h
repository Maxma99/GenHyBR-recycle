// ==========================================================================
// SolverLM: Levenberg-Marquardt

#ifndef __SOLVERLM2_CW_MW_H
#define __SOLVERLM2_CW_MW_H

#include "solver_cw.h"
#include "of.h"

class SolverLM2_CW_MW: public Solver_CW {
public:
    SolverLM2_CW_MW (ParamParser *_pp = NULL);
    ~SolverLM2_CW_MW ();
    SOLVER Type() { return SOLVER_LM; }
    void Solve (RFwdSolverMW &FWS, const Raster &raster,
       const Scaler *pscaler, const ObjectiveFunction &OF, const RVector &data,
       const RVector &sd, Solution &bsol, MWsolution &msol,
       const RCompRowMatrix &qvec, const RCompRowMatrix &mvec);
    void ReadParams (ParamParser &pp);
    void WriteParams (ParamParser &pp);

    enum LM_PRECON {       // Hessian inversion method
	LM_PRECON_NONE,    //   no preconditioner
	LM_PRECON_HDIAG,   //   diagonal of Hessian
	LM_PRECON_CH,      //   Cholesky factorisation using explicit Hessian
	LM_PRECON_ICH,     //   Incomplete Cholesky with sparse Hessian
	LM_PRECON_PCG,     //   PCG solver / implicit Hessian
	LM_PRECON_BICGSTAB,//   BiCGSTAB solver / implicit Hessian
	LM_PRECON_GMRES,   //   GMRES solver / implicit Hessian
	LM_PRECON_GMRES_JACOBIANFREE,// implicit Jacobian (using d+a fields)
	LM_PRECON_GMRES_DIRECT       // implicit Jacobian (direct method)
    } precon;
       
    enum LM_HESS_SCALING { // Hessian scaling method
	LM_HSCALE_NONE,       // no scaling
	LM_HSCALE_IMPLICIT,   // implicit scaling (requires adjoint fields)
	LM_HSCALE_EXPLICIT,   // explicit scaling
	LM_HSCALE_DISTANCE    // boundary distance scaling
    } hscale;

private:
    Regularisation *reg;   // regularisation instance
    PRIOR_OLD prior;           // regularisation method
    double lambda0;        // initial value of control parameter
    double lambda_scale;   // scaling factor for lambda adjustment
    double alpha_min;      // minimum step length for line search
    double alpha0;         // initial (DGN) or fixed (LM) step length
    double gn_tol;         // convergence criterion for (outer) GN solver
    double gmres_tol;      // convergence criterion for (inner) GMRES solver
    int nrmax;             // max GN iterations in outer loop
    int itmax;             // max LM iterations in inner loop
    bool do_linesearch;    // perform linesearch for parameter alpha (including
                           // prior)
};

#endif // !__SOLVERLM2_CW_MW_H
