#ifndef newtonRaphson
#define newtonRaphson

#include <locale>
#include "salvador_matriz.hpp" //replace w available matrix library

#ifdef salvador_matrices
using namespace salvador;
#endif

#ifdef _WIN32
#include <windows.h>
#endif

/*!< Inciso 1*/
double f1_1(double x, double y);//x^2 + xy – 10
double Dxf1_1(double x, double y);
double Dyf1_1(double x, double y);
double f1_2(double x, double y);//y + 3xy^2 – 50
double Dxf1_2(double x, double y);
double Dyf1_2(double x, double y);
/*!< Inciso 2 */
double f2_1(double x, double y);//f1(x,y)=  x^2 + y^2 – 9
double Dxf2_1(double x, double y);
double Dyf2_1(double x, double y);
double f2_2(double x, double y);//f2(x,y)=  -e^x - 2y - 3
double Dxf2_2(double x, double y);
double Dyf2_2(double x, double y);
/*!< Inciso 3 */
double f3_1(double x, double y, double z);//f1=  2x^2 – 4x + y^2 + 3z^2 + 6z + 2
double Dxf3_1(double x, double y, double z);
double Dyf3_1(double x, double y, double z);
double Dzf3_1(double x, double y, double z);
double f3_2(double x, double y, double z);// f2= x^2 + y^2 – 2y + 2z^2 – 5
double Dxf3_2(double x, double y, double z);
double Dyf3_2(double x, double y, double z);
double Dzf3_2(double x, double y, double z);
double f3_3(double x, double y, double z);// f3= 3x^2 – 12x + y^2 - 3z^2 + 8
double Dxf3_3(double x, double y, double z);
double Dyf3_3(double x, double y, double z);
double Dzf3_3(double x, double y, double z);
/*!< Inciso 4 */
double f4_1(double x, double y, double z);//f1= x^2 – 4x + y^2
double Dxf4_1(double x, double y, double z);
double Dyf4_1(double x, double y, double z);
double Dzf4_1(double x, double y, double z);
double f4_2(double x, double y, double z); // f2= x^2 – x – 12y + 1
double Dxf4_2(double x, double y, double z);
double Dyf4_2(double x, double y, double z);
double Dzf4_2(double x, double y, double z);
double f4_3(double x, double y, double z);// f3= 3x^2 – 12x + y^2 - 3z^2 + 8
double Dxf4_3(double x, double y, double z);
double Dyf4_3(double x, double y, double z);
double Dzf4_3(double x, double y, double z);

void metodoNewtonRaphson1(Matriz &Sol, const double &TOL, const int &N);
void metodoNewtonRaphson2(Matriz &Sol, const double &TOL, const int &N);
void metodoNewtonRaphson3(Matriz &Sol, const double &TOL, const int &N);
void metodoNewtonRaphson4(Matriz &Sol, const double &TOL, const int &N);
int menu_NewtRaph();

#endif
