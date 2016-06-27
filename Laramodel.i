%module LaraModel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "GDmodel.h"
#include "LaraModel.h"
%}

%include "GDmodel.h"
%include "LaraModel.h"
