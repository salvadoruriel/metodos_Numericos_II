#include "interpolNewton.hpp"
using namespace std;

int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

int menu_interNewt(Matriz &puntos){
    char redooption,c,recursive;
    //setlocale(LC_ALL,"C");

    if(puntos.get_renglones()<2 ){
        cout << endl<<"No hay suficientes puntos para el metodo.";
        return -1;
    }
    int grado;
    double point, dist,distpas=puntos.get_val(1,0) -puntos.get_val(0,0);

    cout << endl<< "--Interpolacion por Newton progresivo y regresivo--"<< endl;
    for(int i=0; i<puntos.get_renglones()-1; i++){
        if(puntos.get_val(i,0) >= puntos.get_val(i+1,0)){
            cout << "Los puntos de la tabla no estan ordenados"
                << " no se puede usar este metodo.";
            return -1;
        }
        dist = puntos.get_val(i+1,0) -puntos.get_val(i,0);
        if( !salvador::nearly_Equal(dist ,distpas) ){
            cout << i<< " "<<dist <<" = " << distpas <<"Los puntos de la tabla no estan igualmente espaciados,"
                << " no se puede usar este metodo.";
            return -1;
        }
        distpas=puntos.get_val(i+1,0) -puntos.get_val(i,0);
    }
    do{
        do{
            cout<< endl<<"Ingrese punto a interpolar (debe estar dentro del intervalo de la tabla): ";
            cin >> point;
            if(point >puntos.get_val(0,0) && point <puntos.get_val(puntos.get_renglones()-1,0)) break;
        }while(1);
        do{
            cout<< endl<<"Ingrese el grado del polinomio: ";
            cin >> grado;
            if(grado <0){
                cout<< endl<<"Ingrese un grado positivo";
                continue;
            }
            if(grado <=puntos.get_renglones()-1 ) break;
            cout<< endl<<"No se puede hacer un polinomio de ese grado con los puntos dados";
        }while(1);
        /*Para operar encontramos que parte de la tabla tomar*/
        int pos=1,arriba,abajo; /*arriba y abajo se refieren a los puntos disponibles
                    en esas direcciones al intentar ubicar el punto a evaluar en la tabla*/
        for(int i=1; i<puntos.get_renglones()-1; i++){
            if(puntos.get_val(i,0) < point ) pos++;
            else break;
        }
        arriba= pos;
        abajo = puntos.get_renglones() -pos;

        /*Paso 3 presentar resultado*/
        if(abajo >= arriba && abajo >=grado ){ //se hace newton progresivo
            Matriz difDiv(1+grado,2+grado);

            for(int i=0; i<difDiv.get_renglones(); i++){
                difDiv.change_val(i,0,puntos.get_val(pos-1 +i,0) );
                difDiv.change_val(i,1,puntos.get_val(pos-1 +i,1) );
            }
            for(int j=2; j<difDiv.get_columnas(); j++){
                for(int i=0; i<difDiv.get_renglones()-j+1; i++){
                    difDiv.change_val(i,j,difDiv.get_val(i+1,j-1)
                                      - difDiv.get_val(i,j-1) );
                }
            }
            cout << difDiv;
            double sol=difDiv.get_val(0,1),valS =(point - difDiv.get_val(0,0))/dist
            ,mult=valS;
            for(int j=2; j<difDiv.get_columnas(); j++){
                sol += mult*(difDiv.get_val(0,j)/factorial(j-1) ) ;
                cout << endl <<sol <<" , " << mult*(difDiv.get_val(0,j)/factorial(j-1) ) ;
                mult*=valS - (j-1);
            }
            cout << endl << "El valor encontrado con el polinomio de grado "
                << grado << " es: "<< sol;

        }else if( arriba > abajo && arriba>=grado){ //newton regresivo
            Matriz difDiv(1+grado,2+grado);

            for(int i=difDiv.get_renglones()-1; i>=0; i--){
                difDiv.change_val(i,0,
                                puntos.get_val(pos+1
                                               +(i-difDiv.get_renglones()),0) );
                difDiv.change_val(i,1,
                                puntos.get_val(pos+1
                                               +(i-difDiv.get_renglones()),1) );
            }
            for(int j=2; j<difDiv.get_columnas(); j++){
                for(int i=0; i<difDiv.get_renglones()-j+1; i++){
                    difDiv.change_val(i,j,difDiv.get_val(i+1,j-1)
                                      - difDiv.get_val(i,j-1) );
                }
            }
            cout << difDiv;
            double sol=difDiv.get_val(difDiv.get_renglones(),1),
                    valS =(point - difDiv.get_val(difDiv.get_renglones(),0))/dist
            ,mult=valS;
            for(int j=2; j<difDiv.get_columnas(); j++){
                sol += mult*(difDiv.get_val(difDiv.get_renglones()-j,j)/factorial(j-1) ) ;
                cout << endl <<sol <<" , " << mult*(difDiv.get_val(0,j)/factorial(j-1) ) ;
                mult*=valS + (j-1);
            }
            cout << endl << "El valor encontrado con el polinomio de grado "
                << grado << " es: "<< sol;

        }else {
            cout << endl<<"No hay suficientes puntos para obtener el polinomio"
                << " de interpolacion con el grado deseado, en el punto buscado.";
        }

        do{
            cout<< endl <<"Desea interpolar otro punto con la misma tabla?"<< endl
                << "S - si" << endl
                << "N - no" <<endl
                << "--> ";
            while(scanf(" %c",&redooption) != 1)
                        while((c = getchar()) != '\n' && c != EOF);
            while((c = getchar()) != '\n' && c != EOF);

        }while(redooption!='S' && redooption!='s' &&
               redooption!='N' && redooption!='n');
    }while(redooption!='N' && redooption!='n');

    do{
        cout<< endl <<" Quiere cambiar la tabla o regresar al submenú de interpolacion?"<< endl
            << "C - cambiar" << endl
            << "R - regresar" <<endl
            << "--> ";
        while(scanf(" %c",&recursive) != 1)
                    while((c = getchar()) != '\n' && c != EOF);
        while((c = getchar()) != '\n' && c != EOF);
        switch(recursive){
            case 'C':
            case 'c':{
                char interopt,intermenu,again;
                int pun;
                double temx,temy;
                Matriz *recurpuntos;
                do{
                    cout<< endl <<"Ingrese cuantos puntos estan en la tabla:"<< endl;
                    cin >> pun;
                    try{ /*Remember new throws exceptions, in cpp yes*/
                        recurpuntos = new Matriz(pun,2);
                    }catch(const std::bad_alloc& e){
                        std::cerr <<std::endl <<"Error al intentar crear una nueva matriz, posible memoria insuficiente."
                            <<std::endl << "Error: " << e.what() << std::endl;
                        return -1;
                    }
                    cout<< endl <<"Ingrese los puntos en x,y :"<< endl;
                    for(int i=0; i<pun ; i++ ){
                        cout << "x" <<i <<": ";
                        cin >> temx;
                        (*recurpuntos).change_val(i,0,temx);
                        cout << "y" <<i <<": ";
                        cin >> temy;
                        (*recurpuntos).change_val(i,1,temy);
                    }
                    do{
                        cout<< *recurpuntos << endl <<"Son correctos los datos?"<< endl
                            << "S - si" << endl
                            << "N - no" <<endl
                            << "--> ";
                        while(scanf(" %c",&interopt) != 1)
                                    while((c = getchar()) != '\n' && c != EOF);
                        while((c = getchar()) != '\n' && c != EOF);

                    }while(interopt!='S' && interopt!='s' &&
                           interopt!='N' && interopt!='n');
                    if(interopt == 'N' || interopt=='n') delete(recurpuntos);
                }while(interopt!='S' && interopt!='s');

                menu_interNewt(*recurpuntos);
            }
                break;
            case 'R':
            case 'r':
                cout <<endl <<"Presione ENTER para regresar... ";
                cin.get();
                    #ifdef _WIN32
                    system("cls");
                    #endif
                return 0;
                break;

            default:
                break;
        }

    }while(recursive!='C' && recursive!='c' &&
           recursive!='R' && recursive!='r');

return 0;}
