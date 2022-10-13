
from libc.stdint cimport int32_t

cdef extern from "sundials/sundials_types.h":
  ctypedef double realtype;
  ctypedef int32_t sunindextype;

# N_Vector
cdef extern from "sundials/sundials_nvector.h":
  struct _generic_N_Vector:
    pass
  ctypedef _generic_N_Vector *N_Vector

  realtype *N_VGetArrayPointer(N_Vector v)
  void N_VDestroy(N_Vector v)

cdef extern from "nvector/nvector_serial.h":
  N_Vector N_VNew_Serial(sunindextype vec_length);
  N_Vector N_VMake_Serial(sunindextype vec_length, realtype *v_data)

# SUNMatrix
cdef extern from "sundials/sundials_matrix.h":
  struct _generic_SUNMatrix:
    pass
  ctypedef _generic_SUNMatrix *SUNMatrix

  void SUNMatDestroy(SUNMatrix A)

cdef extern from "sunmatrix/sunmatrix_dense.h":
  SUNMatrix SUNDenseMatrix(sunindextype M, sunindextype N)

# SUNLinearSolver
cdef extern from "sundials/sundials_linearsolver.h":
  struct _generic_SUNLinearSolver:
    pass
  ctypedef _generic_SUNLinearSolver *SUNLinearSolver

  int SUNLinSolFree(SUNLinearSolver S)

cdef extern from "sunlinsol/sunlinsol_dense.h":
  SUNLinearSolver SUNLinSol_Dense(N_Vector y, SUNMatrix A)
  
# CVode
cdef extern from "cvode/cvode_ls.h":
  int CVodeSetLinearSolver(void *cvode_mem, SUNLinearSolver LS, SUNMatrix A)

cdef extern from "cvode/cvode.h":
  ctypedef int (*CVRhsFn)(realtype t, N_Vector y, N_Vector ydot, void *user_data);
  ctypedef int (*CVRootFn)(realtype t, N_Vector y, realtype *gout, void *user_data);
  ctypedef void (*CVErrHandlerFn)(int error_code, const char *module, const char *function, char *msg, void *user_data);
  void *CVodeCreate(int lmm);
  int CVodeReInit(void *cvode_mem, realtype t0, N_Vector y0)
  int CVodeInit(void *cvode_mem, CVRhsFn f, realtype t0, N_Vector y0)
  int CVodeSStolerances(void *cvode_mem, realtype reltol, realtype abstol)
  int CVodeSetUserData(void *cvode_mem, void *user_data)
  int CVodeSetMaxNumSteps(void *cvode_mem, long int mxsteps)
  int CVodeSetMaxErrTestFails(void *cvode_mem, int maxnef)
  int CVodeSetMaxConvFails(void *cvode_mem, int maxncf)
  int CVode(void *cvode_mem, realtype tout, N_Vector yout, realtype *tret, int itask);
  int CVodeRootInit(void *cvode_mem, int nrtfn, CVRootFn g)
  int CVodeSetErrHandlerFn(void *cvode_mem, CVErrHandlerFn ehfun, void *eh_data)
  int CVodeGetRootInfo(void *cvode_mem, int *rootsfound)
  void CVodeFree(void **cvode_mem)
