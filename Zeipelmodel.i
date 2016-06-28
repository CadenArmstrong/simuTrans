%module ZeipelModel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "GDmodel.h"
#include "Zeipelmodel.h"
%}

%include "GDmodel.h"
%include "Zeipelmodel.h"
