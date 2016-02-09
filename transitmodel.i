%module transimodel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}
%apply(double *INPLACE_ARRAY1, int DIM1){double *star_params, int np};
%apply(double *INPLACE_ARRAY1, int DIM1){double *planet_params, int np};
%apply(double *INPLACE_ARRAY1, int DIM1){double *phase, int np, *double *deficit_flux, int nf};

%{
#include "transitmodel.h"
%}

%include "transitmodel.h"
