#include "splineCubico.hpp"
using namespace std;

int splineCub(Matriz &tabla,const int & talla){
    int i;
     if(tabla.get_renglones() <= 2 ){
        cout << endl<<"No hay suficientes puntos para el metodo.";
        return -1;
    }
    Matriz h(talla-1,1), diferenciasD(talla-1,1);
    //cout << "hi";
    Matriz soluciones(talla-2,1), splineN(talla,1);
    //cout << "hi2";
    Matriz vIndependiente(talla-2,1), sistemaEc(talla-2,talla-2),
        coeficientes(talla-1,3);
    //cout << "hi3";

     cout << endl<< "--Ajuste de curvas por spline cubico--"<< endl;


    for(i=0; i<talla-1;i++){
        h.change_val(i,0, tabla.get_val(i+1,0)-tabla.get_val(i,0));//Diferencias de los valores de x
        diferenciasD.change_val(i,0, (tabla.get_val(i+1,1)-tabla.get_val(i,1))/h.get_val(i,0) ); //Diferencias divididas
    }

    //cout << "    4hi";
    for(i=0; i<talla-2; i++){
        //Vector independiente del sistema de ecuaciones
        vIndependiente.change_val(i,0, 6*(diferenciasD.get_val(i+1,0)-diferenciasD.get_val(i,0)) );
    }
    //cout << "   5hi";
    //cout << h;

    //Matriz tridiagonal
    for(i=0;i<talla-3;i++){
        sistemaEc.change_val(i,i, 2*(h.get_val(i,0)+h.get_val(i+1,0)) );
        sistemaEc.change_val(i+1,i, h.get_val(i+1,0) );
        sistemaEc.change_val(i,i+1, h.get_val(i+1,0) );
    }
    //cout << endl<<"   llenandohi";
    sistemaEc.change_val(i,i, 2*(h.get_val(i,0)+h.get_val(i+1,0)) );
    //cout << endl<<"   sistemashi";
    //cout << sistemaEc;

    //Solucion del sistema de ecuaciones sin los ceros en S0 y Sn
    cout << " sistEc: "<< sistemaEc.get_renglones() << " " <<sistemaEc.get_columnas() <<endl
        << " vIndep: "<< vIndependiente.get_renglones() << " " <<vIndependiente.get_columnas() <<endl;
    cout << sistemaEc <<sistemaEc.menores();
    soluciones=sistemaEc.inversa()*vIndependiente;
    //cout << endl<<"   solucionesshi";

    //Solucion del sistema del sistema de ecuaciones con ceros en S0 y Sn
    for(i=1;i<talla-1;i++){
        splineN.change_val(i,0, soluciones.get_val(i-1,0) );
    }

    //Calculo de los coeficientes del spline cubico
    for(i=0;i<talla-1;i++){
        coeficientes.change_val(i,0, (splineN.get_val(i+1,0)-splineN.get_val(i,0))/(6*h.get_val(i,0)) );
        coeficientes.change_val(i,1, splineN.get_val(i,0)/2);
        coeficientes.change_val(i,2, diferenciasD.get_val(i,0)-((splineN.get_val(i+1,0)
                                                          +2*(splineN.get_val(i,0)))/6)*h.get_val(i,0) );
    }

    //Impresion del spline
    cout << endl << "---Los ajustadores de curvas (spline cubicos) son:" << endl;
    for(i=0;i<talla-1;i++)
    {
        cout<<"\n g_"<<i<<"(x)= "<<coeficientes.get_val(i,0)<<"(x-"<<tabla.get_val(i,0)<<")^3 + "
        <<coeficientes.get_val(i,1)<<"(x-"<<tabla.get_val(i,0)<<")^2 + "
        <<coeficientes.get_val(i,2)<<"(x-"<<tabla.get_val(i,0)<<") " << ( tabla.get_val(i,1)>0 ? "+" : " ")
        <<tabla.get_val(i,1);

        if( tabla.get_val(i, 0) < tabla.get_val(i+1, 0) ){
            cout<<"\n  \t"<<tabla.get_val(i,0)<<" <= x <= "<<tabla.get_val(i+1,0)<<endl;
        }
        else{
            cout<<"\n  \t"<<tabla.get_val(i+1,0)<<" <= x <= "<<tabla.get_val(i,0)<<endl;
        }
    }

    return 0;
}
