/** \brief menu de materia de metodos numericos II
 *
 * Remember to g++ -std=c++11 menu_[actual].cpp salvador_matriz.cpp salvador_matriz.hpp
 *     		-o menu_metodos
 *  Creado Por Salvador Uriel Aguirre Andrade
 */
#include <cstdio>
#include <cstdlib>
#include <ctime> /*Utilizado para pruebas con random*/
/*#define _USE_MATH_DEFINES*/ /**< Se coloca por si acaso, pero no se recomienda usar */
#include <cmath>

#ifdef _WIN32
#include <windows.h>
#endif

#include <iostream>
#include <iomanip>
/*#include <limits>*//**< Not used, alternative solution to ignore incorrect input */
#include <new>

#include "salvador_matriz.hpp"

static const double C_Euler = std::exp(1.0); /*safer than using a macro M_E*/
char c; /*para control de ingreso incorrecto*/

int main(){
char opt;
salvador::Matriz Jacob(3,3),vecF(3,1),vecxk(3,1),vex(3,1);
vex.change_val(0,0,-4);
vex.change_val(1,0,-0.5);
vex.change_val(2,0,15.5);
double vx,vy,vz;
int iter{3};

    for(int i{0}; i<= iter; i++){
        vx=vex.get_val(0,0);
        vy=vex.get_val(1,0);
        vz=vex.get_val(2,0);
        Jacob.change_val(0,0, 2*vx -1 );
        Jacob.change_val(1,0, 5);
        Jacob.change_val(2,0, -2*vx );
        Jacob.change_val(0,1, 4*vy + vz);
        Jacob.change_val(1,1, -6 );
        Jacob.change_val(2,1, -2*vy );
        Jacob.change_val(0,2, vy );
        Jacob.change_val(1,2, 1);
        Jacob.change_val(2,2, 1 );

        vecF.change_val(0,0,vx*vx -vx +2*vy*vy +vy*vz -10  );
        vecF.change_val(1,0,5*vx -6*vy +vz );
        vecF.change_val(2,0,vz -vx*vx -vy*vy);

        vecxk= vex -Jacob.inversa()*vecF;
        std::cout << " J: " << Jacob << "J^-1" << Jacob.inversa() << "F: " << vecF
                    << "xk+1: " <<vecxk <<"Norma espectral: "
                    << std::fabs(vecxk.get_norma_ESPECTRAL() -vex.get_norma_ESPECTRAL());
        vex=vecxk;
    }

return 0;}
