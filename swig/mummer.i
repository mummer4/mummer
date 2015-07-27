%module(docstring="Mummer binding") mummer
%naturalvar; // Use const reference instead of pointers


%{
#ifdef SWIGPERL
#undef seed
#endif

#include <mummer/nucmer.hpp>
%}


%include "std_string.i"
%include "exception.i"
%include "std_except.i"
%include "typemaps.i"
%include "std_vector.i"

%include "nucmer.i"
