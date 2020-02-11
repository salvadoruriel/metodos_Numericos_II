#ifndef interNewton
#define interNewton

#include "salvador_matriz.hpp" //replace w available matrix library

#ifdef salvador_matrices
using namespace salvador;
#endif

#ifdef _WIN32
#include <windows.h>
#endif
int factorial(int n);
int menu_interNewt(Matriz &puntos);

#endif
