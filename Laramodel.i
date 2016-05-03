%module LaraModel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "model.h"
#include "LaraModel.h"
%}

%include "model.h"
%include "LaraModel.h"
