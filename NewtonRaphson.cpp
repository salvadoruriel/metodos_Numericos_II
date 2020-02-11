#include "NewtonRaphson.hpp"
using namespace std;

/*!< Inciso 1*/
double f1_1(double x, double y){
	return x*x +x*y -10; } //x^2 + xy – 10
double Dxf1_1(double x, double y){
	return 2*x +y; }
double Dyf1_1(double x, double y){
	return x;
}
double f1_2(double x, double y){
	return y + 3*x*y*y -50; } //y + 3xy^2 – 50
double Dxf1_2(double x, double y){
	return 3*y*y; }
double Dyf1_2(double x, double y){
	return 1 +6*x*y;
}
/*!< Inciso 2 */
double f2_1(double x, double y){
	return ((x*x)+(y*y)-9); } //f1(x,y)=  x^2 + y^2 – 9
double Dxf2_1(double x, double y){
	return (2*x); }
double Dyf2_1(double x, double y){
	return (2*y);
}
double f2_2(double x, double y){
	return (-exp(x)-(2*y)-3); }//f2(x,y)=  -e^x - 2y - 3
double Dxf2_2(double x, double y){
	return (-exp(x)); }
double Dyf2_2(double x, double y){
	return (-2);
}
/*!< Inciso 3 */
double f3_1(double x, double y, double z){
	return 2*x*x -4*x +y*y +3*z*z +6*z +2; }//f1=  2x^2 – 4x + y^2 + 3z^2 + 6z + 2
double Dxf3_1(double x, double y, double z){
	return 4*x -4; }
double Dyf3_1(double x, double y, double z){
	return 2*y;}
double Dzf3_1(double x, double y, double z){
	return 6*z +6;
}
double f3_2(double x, double y, double z){
	return x*x + y*y -2*y +2*z*z -5; } // f2= x^2 + y^2 – 2y + 2z^2 – 5
double Dxf3_2(double x, double y, double z){
	return 2*x ; }
double Dyf3_2(double x, double y, double z){
	return 2*y -2;}
double Dzf3_2(double x, double y, double z){
	return 4*z;
}
double f3_3(double x, double y, double z){
	return 3*x*x -12*x +y*y -3*z*z +8; } // f3= 3x^2 – 12x + y^2 - 3z^2 + 8
double Dxf3_3(double x, double y, double z){
	return 6*x -12; }
double Dyf3_3(double x, double y, double z){
	return 2*y;}
double Dzf3_3(double x, double y, double z){
	return -6*z;
}
/*!< Inciso 4 */
double f4_1(double x, double y, double z){
	return x*x -4*x +y*y; }//f1= x^2 – 4x + y^2
double Dxf4_1(double x, double y, double z){
	return 2*x -4; }
double Dyf4_1(double x, double y, double z){
	return 2*y;}
double Dzf4_1(double x, double y, double z){
	return 0;
}
double f4_2(double x, double y, double z){
	return x*x -x -12*y +1; } // f2= x^2 – x – 12y + 1
double Dxf4_2(double x, double y, double z){
	return 2*x -1; }
double Dyf4_2(double x, double y, double z){
	return -12;}
double Dzf4_2(double x, double y, double z){
	return 0;
}
double f4_3(double x, double y, double z){
	return 3*x*x -12*x +y*y -3*z*z +8; } // f3= 3x^2 – 12x + y^2 - 3z^2 + 8
double Dxf4_3(double x, double y, double z){
	return 6*x -12; }
double Dyf4_3(double x, double y, double z){
	return 2*y;}
double Dzf4_3(double x, double y, double z){
	return -6*z;
}

void metodoNewtonRaphson1(Matriz &Sol, const double &TOL, const int &N){
	Matriz F(2,1), Jacobiana(2,2), InversaJacobiana(2,2);
	int i=0;
	double vx = Sol.get_val(0,0), vy = Sol.get_val(1,0);

	while(i<N){

		cout << "\n\n X^" << i << Sol;

		F.change_val(0,0,f1_1(vx, vy ) );
		F.change_val(1,0,f1_2(vx, vy ) );

		if( F.get_norma_ESPECTRAL() < TOL || nearly_Equal(F.get_norma_ESPECTRAL(), TOL) ){
            cout << endl << "Se ha alcanzado la tolerancia "<<TOL << " en F(X^"
                << i <<"): " <<F.get_norma_ESPECTRAL();
			wcout << L"\n La solución del sistema de ecuaciones es: " ;
			cout << Sol;
			return;
		}
		Jacobiana.change_val(0,0,Dxf1_1(vx, vy ) );
		Jacobiana.change_val(1,0,Dxf1_2(vx, vy ) );
		Jacobiana.change_val(0,1,Dyf1_1(vx, vy ) );
		Jacobiana.change_val(1,1,Dyf1_2(vx, vy ) );

		InversaJacobiana=Jacobiana.inversa();

		cout << "\n\n F(X^" << i <<")" << F;
		cout << "\n\n J(X^" << i <<")" << Jacobiana;
		cout << "\n\n J(X^" << i <<")^-1" << InversaJacobiana;

		Sol=(Sol-(InversaJacobiana*F));
		vx = Sol.get_val(0,0);
        vy = Sol.get_val(1,0);

		i++;
	}

	wcout << L"\n\n Se llegó al limite de " << N << " iteraciones." << endl;
	cout << "\n\n La aproximacion al sistema de ecuaciones obtenida es: " << Sol;

return;}

void metodoNewtonRaphson2(Matriz &Sol, const double &TOL, const int &N){
	Matriz F(2,1), Jacobiana(2,2), InversaJacobiana(2,2);
	int i=0;
	double vx = Sol.get_val(0,0), vy = Sol.get_val(1,0);

	while(i<N){

		cout << "\n\n X^" << i <<Sol;

		F.change_val(0,0,f2_1(vx, vy ) );
		F.change_val(1,0,f2_2(vx, vy ) );

		if( F.get_norma_ESPECTRAL() < TOL || nearly_Equal(F.get_norma_ESPECTRAL(), TOL) ){
            cout << endl << "Se ha alcanzado la tolerancia "<<TOL << " en F(X^"
                << i <<"): " <<F.get_norma_ESPECTRAL();
			wcout << L"\n La solución del sistema de ecuaciones es: " ;
			cout << Sol;
			return;
		}
		Jacobiana.change_val(0,0,Dxf2_1(vx, vy ) );
		Jacobiana.change_val(1,0,Dxf2_2(vx, vy ) );
		Jacobiana.change_val(0,1,Dyf2_1(vx, vy ) );
		Jacobiana.change_val(1,1,Dyf2_2(vx, vy ) );

		InversaJacobiana=Jacobiana.inversa();

		cout << "\n\n F(X^" << i <<")" << F;
		cout << "\n\n J(X^" << i <<")" << Jacobiana;
		cout << "\n\n J(X^" << i <<")^-1" << InversaJacobiana;

		Sol=(Sol-(InversaJacobiana*F));
		vx = Sol.get_val(0,0);
        vy = Sol.get_val(1,0);

		i++;
	}

	wcout << L"\n\n Se llegó al limite de " << N << " iteraciones." << endl;
	cout << "\n\n La aproximacion al sistema de ecuaciones obtenida es: " << Sol;

return;}

void metodoNewtonRaphson3(Matriz &Sol, const double &TOL, const int &N){
	Matriz F(3,1), Jacobiana(3,3), InversaJacobiana(3,3);
	int i=0;
	double vx = Sol.get_val(0,0), vy = Sol.get_val(1,0), vz = Sol.get_val(2,0);

	while(i<N){

		cout << "\n\n X^" << i << Sol;

		F.change_val(0,0,f3_1(vx, vy, vz ) );
		F.change_val(1,0,f3_2(vx, vy, vz ) );
		F.change_val(2,0,f3_3(vx, vy, vz ) );

		if( F.get_norma_ESPECTRAL() < TOL || nearly_Equal(F.get_norma_ESPECTRAL(), TOL) ){
            cout << endl << "Se ha alcanzado la tolerancia "<<TOL << " en F(X^"
                << i <<"): " <<F.get_norma_ESPECTRAL();
			wcout << L"\n La solución del sistema de ecuaciones es: " ;
			cout << Sol;
			return;
		}
		Jacobiana.change_val(0,0,Dxf3_1(vx, vy, vz ) );
		Jacobiana.change_val(1,0,Dxf3_2(vx, vy, vz ) );
		Jacobiana.change_val(2,0,Dxf3_3(vx, vy, vz ) );
		Jacobiana.change_val(0,1,Dyf3_1(vx, vy, vz ) );
		Jacobiana.change_val(1,1,Dyf3_2(vx, vy, vz ) );
		Jacobiana.change_val(2,1,Dyf3_3(vx, vy, vz ) );
		Jacobiana.change_val(0,2,Dzf3_1(vx, vy, vz ) );
		Jacobiana.change_val(1,2,Dzf3_2(vx, vy, vz ) );
		Jacobiana.change_val(2,2,Dzf3_3(vx, vy, vz ) );

		InversaJacobiana=Jacobiana.inversa();

		cout << "\n\n F(X^" << i <<")" << F;
		cout << "\n\n J(X^" << i <<")" << Jacobiana;
		cout << "\n\n J(X^" << i <<")^-1" << InversaJacobiana;

		Sol=(Sol-(InversaJacobiana*F));
		vx = Sol.get_val(0,0);
        vy = Sol.get_val(1,0);
        vz = Sol.get_val(2,0);

		i++;
	}

	wcout << L"\n\n Se llegó al limite de " << N << " iteraciones." << endl;
	cout << "\n\n La aproximacion al sistema de ecuaciones obtenida es: " << Sol;

return;}

void metodoNewtonRaphson4(Matriz &Sol, const double &TOL, const int &N){
	Matriz F(3,1), Jacobiana(3,3), InversaJacobiana(3,3);
	int i=0;
	double vx = Sol.get_val(0,0), vy = Sol.get_val(1,0), vz = Sol.get_val(2,0);

	while(i<N){

		cout << "\n\n X^" << i << Sol;

		F.change_val(0,0,f4_1(vx, vy, vz ) );
		F.change_val(1,0,f4_2(vx, vy, vz ) );
		F.change_val(2,0,f4_3(vx, vy, vz ) );

		if( F.get_norma_ESPECTRAL() < TOL || nearly_Equal(F.get_norma_ESPECTRAL(), TOL) ){
            cout << endl << "Se ha alcanzado la tolerancia "<<TOL << " en F(X^"
                << i <<"): " <<F.get_norma_ESPECTRAL();
			wcout << L"\n La solución del sistema de ecuaciones es: " ;
			cout << Sol;
			return;
		}
		Jacobiana.change_val(0,0,Dxf4_1(vx, vy, vz ) );
		Jacobiana.change_val(1,0,Dxf4_2(vx, vy, vz ) );
		Jacobiana.change_val(2,0,Dxf4_3(vx, vy, vz ) );
		Jacobiana.change_val(0,1,Dyf4_1(vx, vy, vz ) );
		Jacobiana.change_val(1,1,Dyf4_2(vx, vy, vz ) );
		Jacobiana.change_val(2,1,Dyf4_3(vx, vy, vz ) );
		Jacobiana.change_val(0,2,Dzf4_1(vx, vy, vz ) );
		Jacobiana.change_val(1,2,Dzf4_2(vx, vy, vz ) );
		Jacobiana.change_val(2,2,Dzf4_3(vx, vy, vz ) );

		InversaJacobiana=Jacobiana.inversa();

		cout << "\n\n F(X^" << i <<")" << F;
		cout << "\n\n J(X^" << i <<")" << Jacobiana;
		cout << "\n\n J(X^" << i <<")^-1" << InversaJacobiana;

		Sol=(Sol-(InversaJacobiana*F));
		vx = Sol.get_val(0,0);
        vy = Sol.get_val(1,0);
        vz = Sol.get_val(2,0);

		i++;
	}

	wcout << L"\n\n Se llegó al limite de " << N << " iteraciones." << endl;
	cout << "\n\n La aproximacion al sistema de ecuaciones obtenida es: " << Sol;

return;}

int menu_NewtRaph(){
    char option,c;
    double TOL; //Tolerancia
    int N; //Iteraciones
    //setlocale(LC_ALL,"C");
    do{
        cout << endl<< "--Metodo de Newton Raphson--"<< endl;
        cout<< endl<<"Seleccione sistema a resolver: "<< endl
                 << "1 - f1(x,y)=  x^2 + xy - 10" << endl
                 << "    f2(x,y)=  y + 3xy^2 - 50" << endl << endl
                 << "2 - f1(x,y)=  x^2 + y^2 - 9" << endl
                 << "    f2(x,y)=  -e^x - 2y - 3 " << endl << endl
                 << "3 - f1(x,y,z)= 2x^2 - 4x + y^2 + 3z^2 + 6z + 2" << endl
                 << "    f2(x,y,z)= x^2 + y^2 - 2y + 2z^2 - 5" << endl
                 << "    f3(x,y,z)= 3x^2 - 12x + y^2 - 3z^2 + 8" << endl << endl
                 << "4 - f1(x,y,z)= x^2 - 4x + y^2" << endl
                 << "    f2(x,y,z)= x^2 - x - 12y + 1" << endl
                 << "    f3(x,y,z)= 3x^2 - 12x + y^2 - 3z^2 + 8" << endl << endl
                 << "m - Regresar al menu principal" << endl
                 << "--> ";
        while(scanf(" %c",&option) != 1)
                    while((c = getchar()) != '\n' && c != EOF);
        while((c = getchar()) != '\n' && c != EOF);

        switch(option){
            case '1':{
            	Matriz Sol(2,1);

                cout << endl << "Ingrese el vector inicial: ";
                Sol.IngresaMatriz();
                cout << "Ingrese la tolerancia (|F(X)|): ";
                cin >> TOL;
                wcout << L"Ingrese el numero máximo de iteraciones: " ;
                cin >> N;
                while((c = getchar()) != '\n' && c != EOF);
                metodoNewtonRaphson1(Sol,TOL,N);

                cout <<endl <<"Presione ENTER para seleccionar otro sistema... ";
                cin.get();
                    #ifdef _WIN32
                    system("cls");
                    #endif
                }
                break;

            case '2':{
                Matriz Sol(2,1);
                cout << endl << "Ingrese el vector inicial: ";
                Sol.IngresaMatriz();
                cout << "Ingrese la tolerancia (|F(X)|): ";
                cin >> TOL;
                wcout << L"Ingrese el numero máximo de iteraciones: " ;
                cin >> N;
                while((c = getchar()) != '\n' && c != EOF);
                metodoNewtonRaphson2(Sol,TOL,N);

                cout <<endl <<"Presione ENTER para seleccionar otro sistema... ";
                cin.get();
                    #ifdef _WIN32
                    system("cls");
                    #endif
                    }
                break;

            case '3':{
            	Matriz Sol(3,1);
                cout << endl <<"Ingrese el vector inicial: ";
                Sol.IngresaMatriz();
                cout << "Ingrese la tolerancia (|F(X)|): ";
                cin >> TOL;
                wcout << L"Ingrese el numero máximo de iteraciones: " ;
                cin >> N;
                while((c = getchar()) != '\n' && c != EOF);

                metodoNewtonRaphson3(Sol,TOL,N);

                cout <<endl <<"Presione ENTER para seleccionar otro sistema... ";
                cin.get();
                    #ifdef _WIN32
                    system("cls");
                    #endif
                }

                break;

            case '4':{
            	Matriz Sol(3,1);
                cout << endl << "Ingrese el vector inicial: " << endl;
                Sol.IngresaMatriz();
                cout << "Ingrese la tolerancia (|F(X)|): ";
                cin >> TOL;
                wcout << L"Ingrese el numero máximo de iteraciones: " ;
                cin >> N;
                while((c = getchar()) != '\n' && c != EOF);
                metodoNewtonRaphson4(Sol,TOL,N);

                cout <<endl <<"Presione ENTER para seleccionar otro sistema... ";
                cin.get();
                    #ifdef _WIN32
                    system("cls");
                    #endif
                }
                break;

            case 'M':
            case 'm':
                break;

            default:
                cout <<endl <<"Opcion incorrecta" << endl;
                cout <<endl <<"Presione ENTER para continuar... ";
                cin.get();
                system("cls"); /*WARNING: Comando especifico de WINDOWS*/
                break;
            }
    }while(option != 'm' && option != 'M');

return 0;}
