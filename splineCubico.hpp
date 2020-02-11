#ifndef splineCubico
#define splineCubico

#include "salvador_matriz.hpp" //replace w available matrix library

#ifdef salvador_matrices
using namespace salvador;
#endif

#ifdef _WIN32
#include <windows.h>
#endif

int splineCub(Matriz &tabla,const int & talla);

#endif
