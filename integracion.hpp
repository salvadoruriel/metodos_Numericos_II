#ifndef integracion
#define integracion

#include "salvador_matriz.hpp" //replace w available matrix library

#ifdef salvador_matrices
using namespace salvador;
#endif

#ifdef _WIN32
#include <windows.h>
#endif
int menintegracion(Matriz &tabla);

#endif
