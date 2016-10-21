%module simplemodel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}
%apply(double *INPLACE_ARRAY1, int DIM1){(double *star_params, int np)};
%apply(double *INPLACE_ARRAY1, int DIM1){(double *planet_params, int np)};
%apply(double *INPLACE_ARRAY1, int DIM1){(double *phase, int np), (double *flux_out, int npo)};

%{
#include "simplemodel.h"
#include "GDmodel.h"
#include "Laramodel.h"
#include "Zeipelmodel.h"
%}

%include "simplemodel.h"
%include "GDmodel.h"
%include "Laramodel.h"
%include "Zeipelmodel.h"
