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
#include "ZeipelModel.h"
%}

%include "GDmodel.h"
%include "ZeipelModel.h"
