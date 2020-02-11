#include <iostream>
#include <locale>
#include "NewtonRaphson.hpp"
#include "interpolNewton.hpp"
#include "splineCubico.hpp"
#include "diferenciacion.hpp"
#include "integracion.hpp"

using namespace std;

int main(){
    char opt,c;
    setlocale(LC_ALL,"");
    do{
        cout<< endl <<"~~Programa para Metodos numericos II~~"<< endl;
        cout<< endl<<"Seleccione tema: "<< endl << endl
                 << "1 - Solucion de ecuaciones no lineales" << endl
                 << "2 - Interpolacion" <<endl
                 << "3 - Diferenciacion e Integracion" <<endl
                 << "q - salir" << endl
                 << "--> ";
        while(scanf(" %c",&opt) != 1)
					while((c = getchar()) != '\n' && c != EOF);
		while((c = getchar()) != '\n' && c != EOF);

	    switch(opt){
            case '1':
                menu_NewtRaph();
                //menu_SolEcNLineales(); solo tenemos 1 metodo ahorita
                break;

            case '2':{
                char interopt,intermenu,again;
                int pun;
                double temx,temy;
                Matriz *puntos;
                do{
                    do{
                        cout<< endl <<"Ingrese cuantos puntos estan en la tabla:"<< endl;
                        cin >> pun;
                        try{ /*Remember new throws exceptions, in cpp yes*/
                            puntos = new Matriz(pun,2);
                        }catch(const std::bad_alloc& e){
                            std::cerr <<std::endl <<"Error al intentar crear una nueva matriz, posible memoria insuficiente."
                                <<std::endl << "Error: " << e.what() << std::endl;
                            return -1;
                        }
                        cout<< endl <<"Ingrese los puntos en x,y :"<< endl;
                        for(int i=0; i<pun ; i++ ){
                            cout << "x" <<i <<": ";
                            cin >> temx;
                            (*puntos).change_val(i,0,temx);
                            cout << "y" <<i <<": ";
                            cin >> temy;
                            (*puntos).change_val(i,1,temy);
                        }
                        do{
                            cout<< *puntos << endl <<"Son correctos los datos?"<< endl
                                << "S - si" << endl
                                << "N - no" <<endl
                                << "--> ";
                            while(scanf(" %c",&interopt) != 1)
                                        while((c = getchar()) != '\n' && c != EOF);
                            while((c = getchar()) != '\n' && c != EOF);

                        }while(interopt!='S' && interopt!='s' &&
                               interopt!='N' && interopt!='n');
                        if(interopt == 'N' || interopt=='n') delete(puntos);
                    }while(interopt!='S' && interopt!='s');
                    do{
                        cout<< endl <<"Seleccione el metodo a usar"<< endl
                            << "1 - Interpolacion por Newton progresivo y regresivo" << endl
                            << "2 - Ajuste de curvas por Spline cubico" <<endl
                            << "--> ";
                        while(scanf(" %c",&intermenu) != 1)
                                    while((c = getchar()) != '\n' && c != EOF);
                        while((c = getchar()) != '\n' && c != EOF);
                        switch(intermenu){
                            case '1':
                                menu_interNewt(*puntos);
                                break;
                            case '2':
                                splineCub(*puntos, pun);
                                break;
                            default:
                                cout <<endl <<"Opcion incorrecta" << endl;
                                break;
                        }
                    }while(intermenu != '1' && intermenu != '2');

                    do{
                        cout<< endl <<"Desea realizar otra interpolacion con otra tabla?"<< endl
                            << "S - si" << endl
                            << "N - no" <<endl
                            << "--> ";
                        while(scanf(" %c",&again) != 1)
                                    while((c = getchar()) != '\n' && c != EOF);
                        while((c = getchar()) != '\n' && c != EOF);

                    }while(again!='S' && again!='s' &&
                           again!='N' && again!='n');
                }while(again!='N' && again!='n');

            }
                break;

            case '3':{
                char diferopt,difermenu,again;
                int pun;
                double temx,paso,temy;
                Matriz *puntos;
                do{
                    do{
                        cout<< endl <<"Ingrese cuantos puntos estan en la tabla:"<< endl;
                        cin >> pun;
                        try{ /*Remember new throws exceptions, in cpp yes*/
                            puntos = new Matriz(pun,2);
                        }catch(const std::bad_alloc& e){
                            std::cerr <<std::endl <<"Error al intentar crear una nueva matriz, posible memoria insuficiente."
                                <<std::endl << "Error: " << e.what() << std::endl;
                            return -1;
                        }
                        cout<< endl <<"De valor inicial x0: "<< endl;
                        cin >> temx;
                        (*puntos).change_val(0,0,temx);
                        cout<< endl <<"Ahora el tamanio de paso: "<< endl;
                        cin >> paso;

                        cout<< endl <<"Ingrese los puntos en f(x):"<< endl;
                        cout << "\t i" << "\t x " <<"\t f(x) "<<endl;
                        for(int i=0; i<pun ; i++ ){
                            if(i!= 0) (*puntos).change_val(i,0, (*puntos).get_val(i-1,0) + paso );
                            cout << "\t " <<i << "\t "<<(*puntos).get_val(i,0) <<"\t ";
                            cin >> temy;
                            (*puntos).change_val(i,1,temy);
                        }
                        do{
                            cout<< *puntos << endl <<"Son correctos los datos?"<< endl
                                << "S - si" << endl
                                << "N - no" <<endl
                                << "--> ";
                            while(scanf(" %c",&difermenu) != 1)
                                        while((c = getchar()) != '\n' && c != EOF);
                            while((c = getchar()) != '\n' && c != EOF);

                        }while(difermenu!='S' && difermenu!='s' &&
                               difermenu!='N' && difermenu!='n');
                        if(difermenu == 'N' || difermenu=='n') delete(puntos);
                    }while(difermenu!='S' && difermenu!='s');
                    do{
                        cout<< endl <<"Seleccione el metodo a usar"<< endl
                            << "1 - Diferenciacion centrada (primera y segunda derivada)" << endl
                            << "2 - Integracion" <<endl
                            << "--> ";
                        while(scanf(" %c",&difermenu) != 1)
                                    while((c = getchar()) != '\n' && c != EOF);
                        while((c = getchar()) != '\n' && c != EOF);
                        switch(difermenu){
                            case '1':
                                pun = mendiferenciacion(*puntos);
                                do{
                                    cout<< endl <<"Desea integrar con la misma tabla?"<< endl
                                        << "S - si" << endl
                                        << "N - no" <<endl
                                        << "--> ";
                                    while(scanf(" %c",&diferopt) != 1)
                                                while((c = getchar()) != '\n' && c != EOF);
                                    while((c = getchar()) != '\n' && c != EOF);

                                }while(diferopt!='S' && diferopt!='s' &&
                                       diferopt!='N' && diferopt!='n');
                                if( (diferopt == 's' || diferopt== 'S') && pun == 0) menintegracion(*puntos);
                                break;
                            case '2':
                                pun = menintegracion(*puntos);
                                do{
                                    cout<< endl <<"Desea diferenciar con la misma tabla?"<< endl
                                        << "S - si" << endl
                                        << "N - no" <<endl
                                        << "--> ";
                                    while(scanf(" %c",&diferopt) != 1)
                                                while((c = getchar()) != '\n' && c != EOF);
                                    while((c = getchar()) != '\n' && c != EOF);

                                }while(diferopt!='S' && diferopt!='s' &&
                                       diferopt!='N' && diferopt!='n');
                                if( (diferopt == 's' || diferopt== 'S') && pun == 0) mendiferenciacion(*puntos);
                                break;
                            default:
                                cout <<endl <<"Opcion incorrecta" << endl;
                                break;
                        }
                    }while(difermenu != '1' && difermenu != '2');

                    do{
                        cout<< endl <<"Desea operar con otra tabla?"<< endl
                            << "S - si" << endl
                            << "N - no" <<endl
                            << "--> ";
                        while(scanf(" %c",&again) != 1)
                                    while((c = getchar()) != '\n' && c != EOF);
                        while((c = getchar()) != '\n' && c != EOF);

                    }while(again!='S' && again!='s' &&
                           again!='N' && again!='n');
                }while(again!='N' && again!='n');

            }
                //menu_DifereInteg();
                break;

            case 'Q':
            case 'q':
                cout << endl << "Fin de programa" << endl
                    << "Integrantes: " << endl
                    << "Aguirre Andrade Salvador Uriel" << endl
                    << "Mejia Maldonado Jose Fernando" << endl
                    << "Cabrejos Reyes Eliseo Aldair" << endl
                    << "Benitez Vazquez Tania Denisse" << endl;
                break;

            default:
                cout <<endl <<"Opcion incorrecta" << endl;
                break;
	    }
        cout <<endl <<"Presione ENTER para continuar... " ;
        cin.get();
		#ifdef _WIN32
                system("cls"); /*WARNING: Comando especifico de WINDOWS*/
		#endif
	}while(opt != 'q' && opt != 'Q');

return 0;}
