#include "diferenciacion.hpp"
using namespace std;

int mendiferenciacion(Matriz &tabla){

	int dim;
	double h;
	dim=tabla.get_renglones();
	if(dim < 3){
        cout << "No hay suficientes puntos para el metodo";
        return -1;
	}
	h=tabla.get_val(1,0)-tabla.get_val(0,0);
	Matriz diferencias(dim-2,4);

	for(int i=1; i<dim+1; i++){
		//Tabla que guarda x,f(x),f'(x),f''(x)
		diferencias.change_val(i-1,0,tabla.get_val(i,0));
		diferencias.change_val(i-1,1,tabla.get_val(i,1));
		diferencias.change_val(i-1,2,(1/(2*h))*(tabla.get_val(i+1,1)-tabla.get_val(i-1,1)));
		diferencias.change_val(i-1,3,(1/(h*h))*(tabla.get_val(i-1,1)-(2*tabla.get_val(i,1))+tabla.get_val(i+1,1)));
	}

	cout<<"Tabla de diferencias de primera y segunda derivada"<<endl;
	cout<<"Las columnas se imprimen en el siguiente orden: x,f(x),f'(x),f''(x)\n"<<endl
        << diferencias;
	return 0;
}
