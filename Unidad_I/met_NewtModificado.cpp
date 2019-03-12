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
salvador::Matriz vecx(2,1),tvec(2,1);
double sx{2.5},sy{3.5},vx{sx},vy{sy},temx;
double f1,f1d,f2,f2d;
int prec{6},espacio{prec+2+1+1};
int iter{7};

std::cout << "Iteraciones simultaneas: "<<std::endl;
    for(int i{0}; i<= iter; i++){
        vecx.change_val(0,0,vx);
        vecx.change_val(1,0,vy);
        /*
        f1= vx*vx -10*vx +vy*vy +6;
        f1d=2*vx -10;
        /*/
        f1= vx*vy*vy +vx -10*vy +4;
        f1d= vy*vy+1;
        //*/
        temx=vx;

        /*
        f2= vx*vy*vy +vx -10*vy +4;
        f2d= 2*vx*vy -10;
        /*/
        f2= vx*vx -10*vx +vy*vy +6;
        f2d=2*vy;
        //*/

        std::cout << std::setw(2)<<i << ":"<< std::setprecision(prec)
                    << std::fixed <<std::setw(espacio)
                    <<temx << " |" << std::setw(espacio) << f1 
                    <<" |" << std::setw(espacio) << f1d ;
        if(i!=0){
            std::cout << " |" <<std::setprecision(prec) << std::setw(espacio)
                    <<std::fabs( (vecx -tvec).get_norma_ESPECTRAL() );
        }
        std::cout << std::endl;
        std::cout << "  "<< std::setprecision(prec)<< std::fixed <<std::setw(espacio)
                    <<vy << " |" << std::setw(espacio) << f2 
                    <<" |" << std::setw(espacio) << f2d << std::endl;

        tvec.change_val(0,0,temx);
        tvec.change_val(1,0,vy);
        vx=vx -f1/f1d;
        vy=vy -f2/f2d;
    }
vx=sx;
vy=sy;
std::cout << "Iteraciones sucesivas: "<<std::endl;
    for(int i{0}; i<= iter; i++){
        vecx.change_val(0,0,vx);
        vecx.change_val(1,0,vy);
        /*
        f1= vx*vx -10*vx +vy*vy +6;
        f1d=2*vx -10;
        /*/
        f1= vx*vy*vy +vx -10*vy +4;
        f1d= vy*vy+1;
        //*/
        temx=vx;
        vx=vx -f1/f1d;

        /*
        f2= vx*vy*vy +vx -10*vy +4;
        f2d= 2*vx*vy -10;
        /*/
        f2= vx*vx -10*vx +vy*vy +6;
        f2d=2*vy;
        //*/

        std::cout << std::setw(2)<<i << ":"<< std::setprecision(prec)
                    << std::fixed <<std::setw(espacio)
                    <<temx << " |" << std::setw(espacio) << f1 
                    <<" |" << std::setw(espacio) << f1d ;
        if(i!=0){
            std::cout << " |" <<std::setprecision(prec) << std::setw(espacio)
                    <<std::fabs( (vecx -tvec).get_norma_ESPECTRAL() );
        }
        std::cout << std::endl;
        std::cout << "  "<< std::setprecision(prec)<< std::fixed <<std::setw(espacio)
                    <<vy << " |" << std::setw(espacio) << f2 
                    <<" |" << std::setw(espacio) << f2d << std::endl;

        tvec.change_val(0,0,temx);
        tvec.change_val(1,0,vy);
        vy=vy -f2/f2d;
    }

return 0;}
