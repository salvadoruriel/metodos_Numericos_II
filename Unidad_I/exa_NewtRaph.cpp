/** \brief
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
salvador::Matriz Jacob(3,3),vecF(3,1),vecxk(vecF),vex(vecF);
vex.change_val(0,0,3);
vex.change_val(1,0,2);
vex.change_val(2,0,2);
double vx,vy,vz;
int iter{6};

    for(int i{0}; i<= iter; i++){
        vx=vex.get_val(0,0);
        vy=vex.get_val(1,0);
        vz=vex.get_val(2,0);
        Jacob.change_val(0,0, 2*vx -4 );
        Jacob.change_val(1,0, 2*vx -1);
        Jacob.change_val(2,0, 6*vx -12 );
        Jacob.change_val(0,1, 2*vy);
        Jacob.change_val(1,1, -12 );
        Jacob.change_val(2,1, 2*vy );
        Jacob.change_val(0,2, 0 );
        Jacob.change_val(1,2, 0);
        Jacob.change_val(2,2, -6*vz );

        if(i != 0)
        std::cout <<"Norma espectral: "
                    << vecF.get_norma_ESPECTRAL() << std::endl;

        vecF.change_val(0,0,vx*vx -4*vx +vy*vy );
        vecF.change_val(1,0,vx*vx -vx -12*vy +1 );
        vecF.change_val(2,0,3*vx*vx -12*vx +vy*vy -3*vz*vz +8);

        vecxk= vex -Jacob.inversa()*vecF;
        std::cout << " J: " << Jacob << "J^-1" << Jacob.inversa() << "F: " << vecF
                    << "xk+"<<i<<": " <<vecxk ;
        
        vex=vecxk;
    };
vex.change_val(0,0,1);
vex.change_val(1,0,1);
vex.change_val(2,0,2);

std::cout << std::endl<< "--------- Con vec inicial 1,1,2" 
        <<std::endl;

    for(int i{0}; i<= iter; i++){
        vx=vex.get_val(0,0);
        vy=vex.get_val(1,0);
        vz=vex.get_val(2,0);
        Jacob.change_val(0,0, 2*vx -4 );
        Jacob.change_val(1,0, 2*vx -1);
        Jacob.change_val(2,0, 6*vx -12 );
        Jacob.change_val(0,1, 2*vy);
        Jacob.change_val(1,1, -12 );
        Jacob.change_val(2,1, 2*vy );
        Jacob.change_val(0,2, 0 );
        Jacob.change_val(1,2, 0);
        Jacob.change_val(2,2, -6*vz );

        if(i != 0)
        std::cout <<"Norma espectral: "
                    << vecF.get_norma_ESPECTRAL() << std::endl;

        vecF.change_val(0,0,vx*vx -4*vx +vy*vy );
        vecF.change_val(1,0,vx*vx -vx -12*vy +1 );
        vecF.change_val(2,0,3*vx*vx -12*vx +vy*vy -3*vz*vz +8);

        vecxk= vex -Jacob.inversa()*vecF;
        std::cout << " J: " << Jacob << "J^-1" << Jacob.inversa() << "F: " << vecF
                    << "xk+"<<i<<": " <<vecxk ;
        
        vex=vecxk;
    }

return 0;}
