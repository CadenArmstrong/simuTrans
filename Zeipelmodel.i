%module ZeipelModel
%{
#define SWIG_FILE_WITH_INIT
%}

%include "numpy.i"
%init %{
import_array();
%}


%{
#include "ZeipelModel.h"
%}

%include "ZeipelModel.h"
