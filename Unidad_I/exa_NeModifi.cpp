/** \brief menu de materia de metodos numericos II
 *
 * Remember to g++ -std=c++11 menu_[actual].cpp salvador_matriz.cpp salvador_matriz.hpp
 *          -o menu_metodos
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
salvador::Matriz vecx(3,1),tvec(3,1);
double sx{3.5},sy{1},sz{1},vx{sx},vy{sy},vz{sz},temx,temy;
double f1,f1d,f2,f2d,f3,f3d;
int prec{6},espacio{prec+2+1+1};
int iter{6};

std::cout << "Iteraciones simultaneas: "<<std::endl;
    for(int i{0}; i<= iter; i++){
        vecx.change_val(0,0,vx);
        vecx.change_val(1,0,vy);
        vecx.change_val(2,0,vz);
        //*
        f1= vx*vx -4*vx +vy*vy;
        f1d= 2*vx -4;
        /*/
        f1= vx;
        //*/
        temx=vx;

        //*
        f2= vx*vx -vx -12*vy +1;
        f2d= -12;
        /*/
        f2= vx6;
        f2d=2;
        //*/
        temy= vy;

        //*
        f3= vx +vy +vz -6;
        f3d= 1;
        /*/
        f2= vx6;
        f2d=2;
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
        std::cout << "   "<< std::setprecision(prec)<< std::fixed <<std::setw(espacio)
                    <<vy << " |" << std::setw(espacio) << f2 
                    <<" |" << std::setw(espacio) << f2d << std::endl;
        std::cout << "   "<< std::setprecision(prec)<< std::fixed <<std::setw(espacio)
                    <<vz << " |" << std::setw(espacio) << f3 
                    <<" |" << std::setw(espacio) << f3d << std::endl;

        tvec.change_val(0,0,temx);
        tvec.change_val(1,0,temy);
        tvec.change_val(2,0,vz);
        vx=vx -f1/f1d;
        vy=vy -f2/f2d;
        vz=vz -f3/f3d;
    }
vx=sx;
vy=sy;
vz=sz;
std::cout << "Iteraciones sucesivas: "<<std::endl;
    for(int i{0}; i<= iter; i++){
        vecx.change_val(0,0,vx);
        vecx.change_val(1,0,vy);
        vecx.change_val(2,0,vz);
        //*
        f1= vx*vx -4*vx +vy*vy;
        f1d= 2*vx -4;
        /*/
        f1= vx;
        //*/
        temx=vx;
        vx=vx- f1/f1d;

        //*
        f2= vx*vx -vx -12*vy +1;
        f2d= -12;
        /*/
        f2= vx6;
        f2d=2;
        //*/
        temy= vy;
        vy=vy -f2/f2d;

        //*
        f3= vx +vy +vz -6;
        f3d= 1;
        /*/
        f2= vx6;
        f2d=2;
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
        std::cout << "   "<< std::setprecision(prec)<< std::fixed <<std::setw(espacio)
                    <<temy << " |" << std::setw(espacio) << f2 
                    <<" |" << std::setw(espacio) << f2d << std::endl;
        std::cout << "   "<< std::setprecision(prec)<< std::fixed <<std::setw(espacio)
                    <<vz << " |" << std::setw(espacio) << f3 
                    <<" |" << std::setw(espacio) << f3d << std::endl;

        tvec.change_val(0,0,temx);
        tvec.change_val(1,0,temy);
        tvec.change_val(2,0,vz);
        vz=vz -f3/f3d;
    }

return 0;}
