from numpy cimport ndarray
cimport cvode_pxd as pxd

cdef class CVodeException(Exception):
  pass

cdef struct CVodeData:
  int neq; # neq
  int ng;
  void *rhsfcn_py; # rhs function in python
  void *rootfcn_py;
  void *args; # tuple of args to be passed to rhsfcn_py

  void *cvode_mem
  pxd.N_Vector y
  pxd.SUNMatrix A
  pxd.SUNLinearSolver LS

cdef void create_CVodeData(CVodeData *d):
  d.neq = 0
  d.ng = 0
  d.rhsfcn_py = NULL
  d.rootfcn_py = NULL
  d.args = NULL
  d.cvode_mem = NULL
  cdef void *y = <void *>d.y
  y = NULL
  cdef void *A = <void *>d.A
  A = NULL
  cdef void *LS = <void *>d.LS
  LS = NULL

cdef void destroy_CVodeData(CVodeData *d):
  if <void*> d.y:
    pxd.N_VDestroy(d.y);
  if <void*> d.A:
    pxd.SUNMatDestroy(d.A)
  if <void*> d.LS:
    pxd.SUNLinSolFree(d.LS);
  if d.cvode_mem:
    pxd.CVodeFree(&d.cvode_mem);  

cdef int CVrhs(double t, pxd.N_Vector y, pxd.N_Vector ydot, void *user_data):
  cdef CVodeData *d = <CVodeData *> user_data

  # Get pointer to nvector
  cdef double *y_1 = pxd.N_VGetArrayPointer(y)
  cdef double *ydot_1 = pxd.N_VGetArrayPointer(ydot)

  # Cython view of pointer
  cdef double[:] y_2 = <double[:d.neq]> y_1
  cdef double[:] ydot_2 = <double[:d.neq]> ydot_1

  # Call RHS runction
  return (<object>d.rhsfcn_py)(t, y_2, ydot_2, *<tuple>d.args)

cdef int CVroot(double t, pxd.N_Vector y, double *gout, void *user_data):
  cdef CVodeData *d = <CVodeData *> user_data

  # Get pointer to nvector
  cdef double *y_1 = pxd.N_VGetArrayPointer(y)

  # Cython view of pointer
  cdef double[:] y_2 = <double[:d.neq]> y_1
  cdef double[:] gout_2 = <double[:d.ng]> gout

  return (<object>d.rootfcn_py)(t, y_2, gout_2, *<tuple>d.args)

cdef void CVerrfcn(int error_code, const char *module, const char *function, char *msg, void *user_data):
  # this function prevents error messages from being printed
  pass
  
cdef class CVode():
  cdef CVodeData d;
  
  def __init__(self, object rhsfcn, double t0, ndarray[double,ndim=1] y0, double rtol = 1.0e-3, double atol = 1.0e-6, 
               tuple args = (), int ng = 0, object g = None, object verbose = False):
    cdef CVodeData *d = &self.d
    create_CVodeData(d)

    if not callable(rhsfcn):
      raise ValueError("rhsfcn must be callable")

    d.neq = y0.shape[0]
    d.ng = ng
    d.rhsfcn_py = <void *>rhsfcn
    d.args = <void *>args

    d.y = pxd.N_VNew_Serial(d.neq)
    if (<void *>d.y == NULL):
      raise CVodeException("N_VNew_Serial failed")
    cdef double *y_p = pxd.N_VGetArrayPointer(d.y);
    cdef int i;
    for i in range(d.neq):
      y_p[i] = y0[i]

    # 2 == BDF method
    d.cvode_mem = pxd.CVodeCreate(2);
    if d.cvode_mem == NULL:
      raise CVodeException("CVodeCreate failed")

    cdef int ret;
    ret = pxd.CVodeInit(d.cvode_mem, CVrhs, t0, d.y);
    if ret < 0:
      raise CVodeException("CVodeInit failed")

    ret = pxd.CVodeSStolerances(d.cvode_mem, rtol, atol);
    if ret < 0:
      raise CVodeException("CVodeSStolerances failed")
    
    ret = pxd.CVodeSetUserData(d.cvode_mem, d);
    if ret < 0:
      raise CVodeException("CVodeSetUserData failed")

    d.A = pxd.SUNDenseMatrix(d.neq, d.neq);
    if (<void *>d.A == NULL):
      raise CVodeException("SUNDenseMatrix failed")

    d.LS = pxd.SUNLinSol_Dense(d.y, d.A);
    if (<void *>d.LS == NULL):
      raise CVodeException("SUNLinSol_Dense failed")

    ret = pxd.CVodeSetLinearSolver(d.cvode_mem, d.LS, d.A)
    if ret < 0:
      raise CVodeException("CVodeSetLinearSolver failed")

    cdef long int mxsteps = 1000
    ret = pxd.CVodeSetMaxNumSteps(d.cvode_mem, mxsteps)
    if ret < 0:
      raise CVodeException("CVodeSetMaxNumSteps failed")

    # default is 7
    ret = pxd.CVodeSetMaxErrTestFails(d.cvode_mem, 15)
    if ret < 0:
      raise CVodeException("CVodeSetMaxErrTestFails failed")

    if ng > 0:
      if not callable(g):
        raise ValueError("g must be callable")
      d.rootfcn_py = <void *>g
      ret = pxd.CVodeRootInit(d.cvode_mem, ng, CVroot)
      if ret < 0:
        raise CVodeException("CVodeSetMaxNumSteps failed")

    if not verbose:
      ret = pxd.CVodeSetErrHandlerFn(d.cvode_mem, CVerrfcn, NULL)
      if ret < 0:
          raise CVodeException("CVodeSetErrHandlerFn failed")
  
  def __dealloc__(self):
    destroy_CVodeData(&self.d)

  def reinit(self, double t0, ndarray[double,ndim=1] y0):
    cdef CVodeData *d = &self.d

    if y0.shape[0] != d.neq:
      raise ValueError("y0 has the wrong shape")

    cdef double *y_p = pxd.N_VGetArrayPointer(d.y);
    cdef int i;
    for i in range(d.neq):
      y_p[i] = y0[i]
    cdef int ret = pxd.CVodeReInit(self.d.cvode_mem, t0, d.y)
    if ret < 0:
      raise CVodeException("CVodeReInit failed")

  def integrate(self, ndarray[double,ndim=0] t, ndarray[double,ndim=1] y):
    cdef CVodeData *d = &self.d

    if y.shape[0] != d.neq:
      raise ValueError("y has the wrong shape")

    cdef double tret;
    cdef int ret = pxd.CVode(d.cvode_mem, t.flat[0], d.y, &tret, 1);
    cdef double *y_p = pxd.N_VGetArrayPointer(d.y);
    cdef int i;
    for i in range(d.neq):
      y[i] = y_p[i]
      
    t.flat[0] = tret

    return ret

  def rootinfo(self, ndarray[int,ndim=1] roots):

    if self.d.rootfcn_py == NULL:
      raise CVodeException("No root function has been sent")

    if roots.shape[0] != self.d.ng:
      raise ValueError("roots has the wrong shape")

    cdef int ret;
    ret = pxd.CVodeGetRootInfo(self.d.cvode_mem, <int *>roots.data)
    if ret < 0:
      raise CVodeException("CVodeGetRootInfo failed")
