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
salvador::Matriz Jacob(3,3),vecF(3,1),Ak(Jacob),Akpas(Jacob);
salvador::Matriz DeltaF(vecF),DeltaX(vecF),vecxk(vecF),vex(vecF),pasF(vecF);
vex.change_val(0,0,3); /*<** 1.vector inicial(cercano a la raiz) */
vex.change_val(1,0,1);
vex.change_val(2,0,1);
double vx,vy,vz; 
int iter{9}; //con x en 3 son 8, en 1 son 14

vx=vex.get_val(0,0);
vy=vex.get_val(1,0);
vz=vex.get_val(2,0);
/*<** 2.Jacob y F(x) */
Jacob.change_val(0,0, 2*vx -4); //f1/x
Jacob.change_val(1,0, 2*vx -1); //f2/x
Jacob.change_val(2,0, 6*vx -12); //f3/x
Jacob.change_val(0,1, 2*vy); //f1/y
Jacob.change_val(1,1, -12);
Jacob.change_val(2,1, 2*vy);
Jacob.change_val(0,2, 0); //f1/z
Jacob.change_val(1,2, 0);
Jacob.change_val(2,2, -6*vz);
pasF.change_val(0,0, vx*vx -4*vx +vy*vy);
pasF.change_val(1,0, vx*vx -vx -12*vy +1);
pasF.change_val(2,0, 3*vx*vx -12*vx +vy*vy -3*vz*vz +8);
/*<**3. 1era iteracion por Newton*/
vecxk= vex -Jacob.inversa()*pasF;

std::cout << "Vector X inicial: "<<vex << "Vector F: "<< pasF 
            << "Jacobiana" << Jacob << "Vector Xk:" << vecxk;
Akpas = Jacob.inversa();
    for(int i{2}; i<= iter; i++){
        /*<** 4.Obtener x0,x1,f0,f1,Dx,Df */
        vx=vecxk.get_val(0,0);
        vy=vecxk.get_val(1,0);
        vz=vecxk.get_val(2,0);

        vecF.change_val(0,0, vx*vx -4*vx +vy*vy);
        vecF.change_val(1,0, vx*vx -vx -12*vy +1);
        vecF.change_val(2,0, 3*vx*vx -12*vx +vy*vy -3*vz*vz +8);

        DeltaF= vecF -pasF;
        DeltaX= vecxk -vex;
        std::cout << "vecF" << vecF
                 << "DeltaX" << DeltaX << "DeltaF" <<DeltaF ;

        /*<** 5.Calcular A^-1, nota que aqui Ak sera la inversa */
        Ak = Akpas +  (1/(DeltaX.transpuesta()*Akpas*DeltaF).get_det() )
                        *(DeltaX - Akpas*DeltaF)*DeltaX.transpuesta()*Akpas;
        //*
        std::cout << "Matriz a sumar: "<<(DeltaX - Akpas*DeltaF)*DeltaX.transpuesta()*Akpas
                    << "Denominador: "<<(DeltaX.transpuesta()*Akpas*DeltaF).get_det()
                    << std::endl << "Ak: " << Ak;
        //*/

        /*<**7. acomodar los vals. */
        pasF= vecF;
        vex=vecxk;

        /*<** 6.Obtener x2=x1 -(A^-1)F1 */
        vecxk = vecxk -Ak*vecF;
        std::cout << "Vector aproximacion x" << i << " :"<<vecxk
                    <<  "Con aproximacion en la norma espectral: "
                    << (DeltaF).get_norma_ESPECTRAL()
                    << std::endl;

        /*<**8. repetir hasta alcanzar la tolerancia */
    }


return 0;}
