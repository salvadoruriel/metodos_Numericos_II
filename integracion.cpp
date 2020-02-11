#include "integracion.hpp"
using namespace std;

int menintegracion(Matriz &tabla){

	int pun,simpoctavos;
	double h;
	pun=tabla.get_renglones();
	if(pun < 3){
        cout << "No hay suficientes puntos para el metodo";
        return -1;
	}
	h =tabla.get_val(1,0) -tabla.get_val(0,0);

	simpoctavos = ( (pun-1)%2 == 0 ? 0 : 1 );
	if(simpoctavos == 1) pun-=3;
	double ressimp13= (tabla.get_renglones() == 4? 0 : tabla.get_val(0,1) + tabla.get_val(pun-1 ,1) );
    //cout << endl << ressimp13 << endl;
	for(int i=1; i<pun-1; i++){ /*Sumando por simpsion 1/3 */
        ressimp13 += (tabla.get_val(i,1) *(i%2 == 0 ? 2 : 4 ) );
        /*cout << tabla.get_val(i,1) <<" *" << (i%2 == 0 ? 2 : 4 )<< " = "
            << (tabla.get_val(i,1) *(i%2 == 0 ? 2 : 4 ) ) <<endl;
        */
	}
    //cout << endl << ressimp13 << endl;
	ressimp13 *= h/3;
    //cout << endl << ressimp13 << endl;

	double ressimp38 = (simpoctavos == 0 ? 0 : tabla.get_val(pun-1 ,1)
                     + 3*tabla.get_val(pun ,1) + 3*tabla.get_val(pun+1 ,1)
                     + tabla.get_val(pun+2 ,1)  );
	ressimp38 *= 3*h/8;


	cout <<endl <<"Valor de la integral por"
        <<( simpoctavos == 0 ? " regla de simpson 1/3"  : (
            tabla.get_renglones()==4 ? " regla de simpson 3/8"
            : " regla de simpson 1/3 y simpson 3/8" ) )
        << " desde x0= " << tabla.get_val(0,0)
        << ", hasta xn= "<< tabla.get_val(tabla.get_renglones()-1 ,0)
        << "\n Es:  "<< ressimp13 + ressimp38 <<endl;

	return 0;
}
