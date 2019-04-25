#include <cstdlib> /*para calloc y free*/
#include <cmath> /*operaciones*/

#include <iostream>
#include <algorithm> /*for the max*/
#include <iomanip> /*Para dar formato a std con internal,left,right*/

/** \brief operaciones basicas con matrices dinamicas
 *
 * Se recomienda leer los comentarios al lado de cada funcion
 *  ya que algunas requieren detalles, de todas formas
 *  los campos deberian ser lo suficientemente especificos
 *  Autor Salvador Uriel Aguirre Andrade
 * WARNING: Usa algunos detalles de c++11.
 *
 */
#include "salvador_matriz.hpp"

namespace salvador{

    Matriz::Matriz(int oren=2,int ocol=2){
        ren = oren;
        col = ocol;
        det = 0; /*Ya que iniciamos la matriz en ceros*/

        if(ren > 0) mat = (double**) calloc(ren,sizeof(double));
        for(int i = 0; i < ren; i++){
            mat[i] = (double*) calloc(col,sizeof(double));
        }
    }

    /*! \brief creates a smaller matrix without the row and column selected
     *
     * \param A, la matriz a copiar
     * \param minren el renglon a excluir
     * \param mincol la columna a excluir
     * \note use -1 to avoid excluding any row or col
     *
     * \WARNING avoid vectores, they will make empty memory
     */
    Matriz::Matriz(const Matriz &A,int minren, int mincol){
        if(minren >= 0) ren = (A.get_renglones() -1 > 0 ? A.get_renglones()-1 : 0);
        else ren = A.get_renglones();

        if(mincol >= 0)col = (A.get_columnas() -1 > 0 ? A.get_columnas()-1 : 0);
        else col = A.get_columnas();

        if(ren > 0) mat = (double**) calloc(ren,sizeof(double));
        for(int i = 0; i < ren; i++){
            mat[i] = (double*) calloc(col,sizeof(double));
        }
        int avoidren{0}, avoidcol{0};
        for(int i=0; i<ren; i++){
            if(i == minren) avoidren = 1;
            avoidcol = 0;
            for(int j=0; j<col; j++){
                if(j == mincol) avoidcol = 1;
                mat[i][j] = A.mat[i+avoidren][j+avoidcol];
            }
        }
        det = (*this).get_det();
    }
    Matriz::Matriz(const Matriz &A){
        ren = A.ren;
        col = A.col;
        det = A.det; /*WARNING: we should still ask for det with get_det
                    this exception is cuz get_det does modify the matrix,
                    but we dont want that here, thats why we called it const*/
        if(ren > 0) mat = (double**) calloc(ren,sizeof(double));
        for(int i = 0; i < ren; i++){
            mat[i] = (double*) calloc(col,sizeof(double));
        }
        for(int i=0; i<ren; i++){
            for(int j=0; j<col; j++){
                mat[i][j] = A.mat[i][j];
            }
        }
    }

    Matriz::~Matriz(){
        int i;
        for(i = 0; i < ren; i++){
            free(mat[i]);
        }
        free(mat);
    }

    std::ostream& operator<<(std::ostream &out, const Matriz &A){
        int i,j,ren = A.get_renglones(), col = A.get_columnas()
                ,precision{ 7 },espacio{precision+1+2+1} ;
                /**< precision: n de digitos despues del punto */
                /**< setw(r) donde r= n +1(por el .) +2(por unidad a decena)
                                        +1(por el signo)*/
        double temp;
        out <<std::endl<<"r\\c";/*Ya que desconocemos el tamanio de la matriz*/
        for(i = 0; i < col ; i++){
            out <<std::setw(espacio) << i;
            if(i != col-1) out << " |";
        }
        out << std::endl;
        for(i = 0; i < ren ; i++){
                out  << std::setw(2) << i<<":";
            for (j = 0; j < col ; j++){
                temp = ( std::fabs(A.get_val(i,j))
                        < (1/std::pow(10, precision)) ? 0 : A.get_val(i,j) );
                out  <<std::setprecision(precision)
                    << std::fixed << std::setw(espacio)<< temp  ;
                if(j != col-1) out << " |";
            }
            if(j == col) out << std::endl;
        }
        out << std::endl;
    return out;}

    /**************************************************
                    Operaciones logicas
    ****************************************************/
    bool Matriz::es_cuadrada() const { return ren == col;}
    bool Matriz::es_simetrica() const{
        if( ren == 0 || col == 0 || ren != col) return false;

        for(int i =0; i < ren ; i++){
            for(int j =0; j < i; j++){
                if(mat[i][j] != mat[j][i] ) return false;

            }
        }
    return true;}

    bool Matriz::es_EDD() const{
        if( ren == 0 || col == 0) return false;

        double mayor,resultado;

        for(int i =0; i < ren ; i++){
            mayor = std::fabs(mat[i][i]);
            resultado=0;
            for(int j =0; j < col; j++){
                if(j == i) continue;
                resultado += std::fabs(mat[i][j]);
            }
            if(mayor <= resultado) return false;
        }
    return true;}

    bool Matriz::es_DD() const{
        if( ren == 0 || col == 0) return false;

        double mayor;

        for(int i =0; i < ren ; i++){
            mayor = std::fabs(mat[i][i]);
            for(int j =0; j < col; j++){ /*even if we go to i,i the evaluation should fail*/
                if(mayor < std::fabs(mat[i][j]) ) return false;
            }
        }
    return true;}

    /**************************************************
            Operaciones sobre la matriz o sus valores
    ****************************************************/
    int Matriz::replace_Mat(double data[],size_t datasize){ /*Ya que la matriz esta iniciada en ceros*/
        int temp_ren= 0,temp_col= 0, available =(int) (datasize/sizeof(double));
        int current=0,i,j,state=available - ren*col;

        if(state < 0){/*control en caso de tener menos datos que los necesarios*/
            while(available - 1 >=0){
                temp_col++;
                if(temp_col == col){
                    temp_ren++;
                    temp_col = 0;
                }
                available--;
            }
        }else{
            temp_ren = ren;
            temp_col = col;
        }

        for(i = 0; i < temp_ren; i++){
            for (j = 0; j < temp_col; j++){
                        /*
                        std::cout << "mat["<<i << "][" << j<<"] = data[" << current <<"] === " ;
                        printf("%f = %f",(mat[i][j]), (data[current]) );
                        std::cout  <<std::endl;
                        //*/           /*******DEBUGGING*/
                mat[i][j] = data[current];
                        /*
                        std::cout << "mat["<<i << "][" << j<<"] = data[" << current <<"] === " ;
                        printf("%f = %f",(mat[i][j]), (data[current]) );
                        std::cout  <<std::endl;
                        //*/           /*******DEBUGGING*/
                current++;
            }
        }
    return state;}
	Matriz Matriz::replace_Mat(double data){
		 for(int i = 0; i < ren; i++){
	            for (int j = 0; j < col; j++){
		                mat[i][j] = data;
	                        /*
	                        std::cout << "mat["<<i << "][" << j<<"] = data === " ;
	                        printf("%f = %f",(mat[i][j]), data );
	                        std::cout  <<std::endl;
	                        //*/           /*******DEBUGGING*/
	            }
       		 }

	return *this;}
	Matriz Matriz::replace_Mat(std::string type){
		if(type.compare("I") == 0 ){ /*IDENTITY matrix*/
			for(int i = 0; i < ren; i++){
		            for (int j = 0; j < col; j++){
		                        /*
		                        std::cout << "mat["<<i << "][" << j<<"] = i==j?1:0 === " ;
		                        printf("%f = %f",(mat[i][j]), i==j?1:0 );
		                        std::cout  <<std::endl;
		                        //*/           /*******DEBUGGING*/
			               mat[i][j] = ( i==j ? 1 : 0 );
		                        /*
		                        std::cout << "mat["<<i << "][" << j<<"] = i==j?1:0 === " ;
		                        printf("%f = %f",(mat[i][j]), i==j?1:0 );
		                        std::cout  <<std::endl;
		                        //*/           /*******DEBUGGING*/
		            }
	       		 }
		return *this;}
		Matriz::replace_Mat(std::nan(""));
	return *this;} /*function doesnt exist*/

    Matriz Matriz::mpow(Matriz a,int x){
        if(x <= 1) return a;
        return a * Matriz::mpow(a,x-1);
    }

    int Matriz::change_val(int ren_buscado,int col_buscada,double nuevo_val){/*contando desde 0*/
        if(ren_buscado >= ren || col_buscada >= col) return -1; /*error: value out of bounds*/

        mat[ren_buscado][col_buscada] = nuevo_val;
    return 0;}/*success!*/

    void Matriz::determinante_Mat(){ /*ya que el determinante es un valor nato, sera void la funcion*/
        /*Usando metodo de gauss, para triangular*/
        if( !((*this).es_cuadrada()) ) det = std::nan("");
        if( ren == 0 || col == 0) det = std::nan("");

        double resultado,mayor = mat[0][0], aux_res_operaciones=1, aux_temp
                        ,pivote=0,temp_first=0;
        unsigned int i ,j,ren_mayor=0,actual=0;
        Matriz temp(*this);
                    /*
                    std::cout << std::endl <<"A" << std::endl << (*this)
                                << std::endl <<"temp" << std::endl << temp;
                    //*/ /*******DEBUGGING*/
        do{
                //*
            aux_res_operaciones *= pow(-1,temp.hacer_DD());
                /*/
                        std::cout << std::endl <<"Despues de hacer diagonalmente dominante, A:" << std::endl
                                    << (*this) << std::endl <<"temp" << std::endl << temp;//this should go w/1st lol
            mayor = std::fabs(temp.mat[actual][actual]);
            ren_mayor = actual;
            for(i=actual; i< ren;i++){
                aux_temp = std::fabs(temp.mat[i][actual])
                if(aux_temp > mayor){
                    mayor = aux_temp;
                    ren_mayor = i;
                }
            }
            if(ren_mayor != actual){
                        std::cout <<std::endl<<"{DEBUGGING}intercambio: "<< actual <<","<<ren_mayor<<std::endl;
                temp.intercambia_rens(actual,ren_mayor);
                aux_res_operaciones *= -1;
                        std::cout << std::endl <<"{DEBUGGING}Despues de intercambiar rens, A:" << std::endl
                                    << (*this) << std::endl <<"temp" << std::endl << temp;
            }
                //*/ /**< ATTENTION: la primer funcion es superior y menos verbose que la segunda*/
            temp_first = temp.mat[actual][actual];
            if(temp_first == 0){
                actual++;
                continue;
            }
            for(i=actual+1;i<ren;i++){
                pivote = temp.mat[i][actual];
                for(j=actual;j<col;j++){
                    temp.mat[i][j] -= pivote*temp.mat[actual][j]/temp_first;
                }
            }
                        /*
                        std::cout << std::endl <<"Despues de operar, A:" << std::endl << (*this)
                                << std::endl <<"temp" << std::endl << temp;
                        //*/ /*******DEBUGGING*/
            actual++;
        }while(actual < ren );
        resultado = temp.mat[0][0];
        for(i=1;i<ren;i++){
            resultado *= temp.mat[i][i];
        }
                    /*
                    printf("\nresultado: %f\n",resultado);
                    std::cout << std::endl <<"Final A" << std::endl<< (*this)
                                << std::endl <<"Final temp" << std::endl << temp;
                    //*/ /*******DEBUGGING*/

    det = resultado*aux_res_operaciones;}

    void Matriz::intercambia_rens(unsigned int original,unsigned int cambio){
        if(original == cambio) return;
        double a_orig[col], b_camb[col];
        Matriz::da_ren(a_orig,original);
        Matriz::da_ren(b_camb,cambio);
        int i;
        /*Ya estamos dentro de la funcion de la clase Matriz, no hace
        falta llamar a la func de cambiar un cierto valor*/
        for(i= 0; i< col; i++){
            mat[original][i] = b_camb[i];
            mat[cambio][i] = a_orig[i];
        }
    }

    int Matriz::hacer_DD(){ /**<  Regresa la cantidad de renglones movidos*/
        int movimientos=0,actual=0,ren_mayor; /*actual es el renglon actual en el que operamos*/
        double mayor,aux;
        do{
            mayor = mat[actual][actual];
            ren_mayor = actual;
            for(int i=actual; i< ren;i++){
                aux = std::fabs(mat[i][actual]);
                if(aux > mayor){
                    mayor = aux;
                    ren_mayor = i;
                }
            }
            if(ren_mayor != actual){ /*Cambiamos si no es el mas grande*/
                        /*
                        std::cout <<std::endl<<"intercambio: "<< actual <<","<<ren_mayor<<std::endl;
                        //*/ /*******DEBUGGING*/
                (*this).intercambia_rens(actual,ren_mayor);
                movimientos++;
            }
            actual++;
        }while(actual < ren);
    return movimientos;}

    int Matriz::hacer_DD_con_VECTOR(Matriz &vec){ /**<  Regresa la cantidad de renglones movidos*/
        int movimientos=0,actual=0,ren_mayor; /*actual es el renglon actual en el que operamos*/
        if (ren != vec.get_renglones()) return -1; /*Error, vector must match*/
        double mayor,aux;
        do{
            mayor = mat[actual][actual];
            ren_mayor = actual;
            for(int i=actual; i< ren;i++){
                aux = std::fabs(mat[i][actual]);
                if(aux > mayor){
                    mayor = aux;
                    ren_mayor = i;
                }
            }
            if(ren_mayor != actual){ /*Cambiamos si no es el mas grande*/
                        /*
                        std::cout <<std::endl<<"intercambio: "<< actual <<","<<ren_mayor<<std::endl;
                        //*/ /*******DEBUGGING*/
                (*this).intercambia_rens(actual,ren_mayor);
                vec.intercambia_rens(actual,ren_mayor);
                movimientos++;
            }
            actual++;
        }while(actual < ren);
    return movimientos;}

    /**************************************************
                RETORNO de datos de la matriz
    ****************************************************/
    /*Imprime la matriz hasta la posicion ren,col, empezando desde 0*/
    int Matriz::print_Mat(int busca_ren,int busca_col) const{
        if(busca_ren >= ren) busca_ren = ren -1; /*Para evitar desbordes */
        if(busca_col >= col) busca_col = col -1; /* al buscar en memoria que no exista*/
        int i,j;
                        /*
                        std::cout << "\n busca_col = " << busca_col << ", busca_ren= " << busca_ren << std::endl;
                        //*/ /*******DEBUGGING*/
        std::cout << "   ";
        for(i = 0; i <= ((busca_ren>0) ? col-1 : busca_col ) ; i++){
            std::cout <<std::setw(8) << std::left << i;
            if(i != col-1) std::cout << std::setw(3)<< " | ";
        }
        std::cout << std::endl;
        for(i = 0; i <= busca_ren ; i++){
                std::cout  << std::left << i<<": ";
                if (i==busca_ren && busca_col<0) return 1;
            for (j = 0; j < col ; j++){
                std::cout  <<  std::setw(8) << std::left << std::showpos << mat[i][j]  << std::noshowpos ;
                if(j != col-1) std::cout << std::setw(3)<< " | ";
                        /*
                        std::cout <<"(" <<std::boolalpha  << "(busca_ren == ren -1)= " <<(busca_ren == ren -1)
                                    <<",(busca_col == col -1)= " << (busca_col == col -1)
                            <<",(busca_ren == ren -1) && (busca_col == col -1)= " << ((busca_ren == ren -1) && (busca_col == col -1))
                                    << std::noboolalpha<<")" << std::endl;
                            //*/ /*******DEBUGGING*/
                if ( !((busca_ren == ren -1) && (busca_col == col -1)) ){ /*If we don't search the full matrix*/
                    if (i==busca_ren && (j==busca_col || busca_col<0) ) return 1;} /*found spot*/
            }
            if(j == col) std::cout << std::endl;
        }
    return 0;}/*Printed whole matrix*/
    //*
    int Matriz::print_full_Mat() const{ return Matriz::print_Mat(ren,col);}
    //*/ /*DEPRECATED: use operator << and send it to stream*/
    double * Matriz::toArray(double data[]) const{ /*data is required to be of size ren*col atleast*/
        for(int i = 0; i < ren ; i++){
            for (int j = 0; j < col ; j++){
                data[i*col + j] = mat[i][j] ;
            }
        }
    return data;}

    double * Matriz::da_ren(double data[],unsigned int ren_buscado) const{ /*Warning, data is required to be of size col*/
        if(ren_buscado > ren){ /*Non existent row*/
            data[0] = std::nan("");
            return data;
        }
        if(ren_buscado == ren) ren_buscado--;
        int i;
        for(i = 0; i < col; i++){
            data[i] = mat[ren_buscado][i];
        }

    return data;}

    double Matriz::get_det(){/*WARNING: Siempre se debe buscar el determinante con esta funcion*/
        Matriz::determinante_Mat(); /*Para que se calcule en el momento requerido y tener mas velocidad*/
    return det;}
    double Matriz::get_Solo_det() const{ return det;} /*Solo si ya se uso get_det deberiamos usar esta func.*/

    double Matriz::get_val(int ren_buscado,int col_buscada) const{
        if(ren_buscado == ren) ren_buscado = ren -1; /*Para evitar desbordes */
        if(col_buscada == col) col_buscada = col -1;
        if(ren_buscado > ren) return std::nan("");
        if(col_buscada > col) return std::nan("");
    return mat[ren_buscado][col_buscada];}

    double Matriz::get_norma_ESPECTRAL() const{
        if( ren == 0 || col == 0) return std::nan("");

        double mayor = std::fabs(mat[0][0]), aux_temp;

        for(int i =0; i < ren ; i++){
            for(int j =0; j < col; j++){ /*even if we go to i,i the evaluation should fail*/
                aux_temp = std::fabs(mat[i][j]);
                if(mayor < aux_temp ) mayor = aux_temp;
            }
        }
    return mayor;}
    double Matriz::get_Biggest() const{
            if( ren == 0 || col == 0) return std::nan("");

            double mayor = mat[0][0];

            for(int i =0; i < ren ; i++){
                for(int j =0; j < col; j++){
                    if(std::fabs(mayor) < std::fabs(mat[i][j]) ) mayor = mat[i][j];
                }
            }
    return mayor;}

    double Matriz::get_Traza() const{
        if( !((*this).es_cuadrada()) ) return std::nan("");
        double traz{0};
        for(int i=0; i < ren; i++){
            traz += mat[i][i];
        }

    return traz;}

    int Matriz::get_renglones() const{ return ren;}
    int Matriz::get_columnas() const{ return col;}

	void Matriz::LU_Doolittle(Matriz& L,Matriz& U) const{
		if( Matriz::LU_DoolittlePROC(L,U) == -1){
			double oneCallOnly = std::nan("");
			L.replace_Mat(oneCallOnly );
			U.replace_Mat(oneCallOnly );
		}
	}
	int Matriz::LU_DoolittlePROC(Matriz& L,Matriz& U) const{
		if( L.get_renglones() != U.get_renglones()  /*Para factorizacion LU necesitamos matrices iguales*/
		|| U.get_renglones() != Matriz::get_renglones()
		|| L.get_columnas() != U.get_columnas()
		|| U.get_columnas() != Matriz::get_columnas()
        || !Matriz::es_cuadrada() ) return -1; /**TODO: quitar requisito de ser cuadrada*/
		L.replace_Mat("I"); /*En Doolittle la mat L tiene 1's en la diagonal*/
		U.replace_Mat(0.0); /*we make sure the other one is empty*/
		int actual{0};
		double sum;

		for(int end{0}; end < ren; end++){
            for(int i{actual}; i < col; i++){ /**primero operamos por renglon*/
                sum = 0;
                for(int j{0}; j < actual ; j++){
                        sum+= L.get_val(actual, j)*U.get_val(j,i);
                }
                /*
                std::cout <<"sum["<<i<<"]= "<<sum<<std::endl;
                /*/ /*!< Debugging */
                U.change_val(actual,i, (Matriz::get_val(actual,i) - sum) );
            }
        /*
        std::cout <<"U = " << U <<std::endl;
        /*/ /*!< Debugging */
            for(int j{actual+1}; j < ren ; j++ ){ /**despues por columna*/
                sum = 0;
                for(int i{0}; i < actual; i++){
                        sum+= L.get_val(j, i)*U.get_val(i,actual);
                }
                /*
                std::cout <<"sum["<<j<<"]= "<<sum<<std::endl;
                /*/ /*!< Debugging */
                L.change_val(j,actual, (Matriz::get_val(j,actual) - sum) /U.get_val(actual, actual) );
            }
            actual++;
        /*
        std::cout << "L = " << L <<std::endl;
        /*/ /*!< Debugging */

		}
	return 0;}

    Matriz Matriz::LU_Cholesky() const{
        if( !(*this).es_simetrica() ) return Matriz(0,0);
        Matriz temp(ren,col);
		int actual{0};
		double sum;

		for(int end{0}; end < ren; end++){
            for(int i{actual}; i < col; i++){
                sum = 0;
                for(int j{0}; j < actual ; j++){
                        sum+= temp.get_val(actual, j)*temp.get_val(i,j);
                }
                if( i == actual )
                temp.change_val(i,actual, sqrt( (*this).get_val(actual,i) - sum) );/*!<
                    La matriz a regresar va ser Lower, su transpuesta sera la Upper*/
                else
                temp.change_val(i,actual, ((*this).get_val(actual,i) - sum)
                                                    /temp.get_val(actual,actual) );
                /*
                std::cout<<"==sum["<<i<<"]= "<<sum <<std::endl
                        <<"A("<<actual<<","<<i<<"): "<< (*this).get_val(actual,i)<<std::endl
                        <<"A:" << temp <<std::endl;
                /*/ /*!< Debugging */
            }
            actual++;
		}

    return temp;}

    void Matriz::LU_Crout(Matriz &L,Matriz &U) const{
		if( Matriz::LU_CroutPROC(L,U) == -1){
			double oneCallOnly = std::nan("");
			L.replace_Mat(oneCallOnly );
			U.replace_Mat(oneCallOnly );
		}

    }
    int Matriz::LU_CroutPROC(Matriz &L, Matriz &U) const{
		if( L.get_renglones() != U.get_renglones()  /*necesitamos matrices iguales*/
		|| U.get_renglones() != Matriz::get_renglones()
		|| L.get_columnas() != U.get_columnas()
		|| U.get_columnas() != Matriz::get_columnas()
        || !Matriz::es_cuadrada() ) return -1;

		L.replace_Mat(0.0); /*we make sure the other one is empty*/
		U.replace_Mat("I"); /*En Doolittle la mat L tiene 1's en la diagonal*/
		int actual{0};
		double sum;

		for(int end{0}; end < ren; end++){
            for(int j{actual}; j < ren ; j++ ){ /**primero operamos por  columna*/
                sum = 0;
                for(int i{0}; i < actual; i++){
                        sum+= L.get_val(j, i)*U.get_val(i,actual);
                }
                /*
                std::cout <<"sum["<<j<<"]= "<<sum<<std::endl;
                /*/ /*!< Debugging */
                L.change_val(j,actual, (Matriz::get_val(j,actual) - sum) );
            }
        /*
        std::cout <<"L = " << L <<std::endl;
        /*/ /*!< Debugging */
            for(int i{actual+1}; i < col; i++){ /**despues por renglon*/
                sum = 0;
                for(int j{0}; j < actual ; j++){
                        sum+= L.get_val(actual, j)*U.get_val(j,i);
                }
                /*
                std::cout <<"sum["<<i<<"]= "<<sum<<std::endl;
                /*/ /*!< Debugging */
                U.change_val(actual,i, (Matriz::get_val(actual,i) - sum) /L.get_val(actual, actual) );
            }
            actual++;
        /*
        std::cout << "U = " << U <<std::endl;
        /*/ /*!< Debugging */

		}


    return 0;}


	Matriz Matriz::transpuesta() const {
	    if((*this).es_cuadrada() ){
            Matriz trans(*this);
            for(int i=0; i< ren; i++){
                for(int j=i+1; j< col; j++){
                    trans.change_val(i,j, (*this).get_val(j,i)  );
                    trans.change_val(j,i, (*this).get_val(i,j)  );
                }
            }
            return trans;
	    }
	    int new_ren{ (*this).get_columnas() }, new_col{ (*this).get_renglones() };
	    Matriz trans( new_ren, new_col);
            for(int i=0; i< new_ren; i++){
                for(int j=0; j< new_col; j++){
                    trans.change_val(i,j, (*this).get_val(j,i)  );
                    trans.change_val(j,i, (*this).get_val(i,j)  );
                }
            }

	return trans;}
	Matriz Matriz::menores() const {
	    Matriz mens(this->ren,this->col);
		for(int i=0; i< ren; i++){
			for(int j=0; j< col; j++){
                /*
                std::cout << std::endl << Matriz(*this,i,j) << Matriz(*this,i,j).get_det() << std::endl;
                /*/ /*!< DEBUGGING */
                mens.change_val(i,j, Matriz(*this,i,j).get_det() );
			}
		}

	return mens;}
	Matriz Matriz::cofactores() const {
	    Matriz cofacs( (*this).menores() );
	    int aux=1;
		for(int i=0; i< ren; i++){
			for(int j=0; j< col; j++){
                aux = ( (i+j)%2 == 0 ? 1: -1 );
                cofacs.change_val(i,j, aux*cofacs.get_val(i,j) );
			}
		}

	return cofacs;}
	Matriz Matriz::adjunta() const{
	    return (  ( (*this).cofactores() ).transpuesta()   ); }
    Matriz Matriz::inversa() const{
        /*
        std::cout << (*this) << std::endl
                <<(*this).menores() << std::endl
                <<(*this).cofactores() << std::endl
                << ((*this).cofactores()).transpuesta() << std::endl
                << (*this).adjunta() << std::endl
                << (1/Matriz(*this).get_det())*((*this).cofactores()).transpuesta() << std::endl;
        /*/ /*!< DEBUGGING */
        return (  (1/Matriz(*this).get_det()) * (*this).adjunta()  );
    }
    Matriz Matriz::transformacionHouseholder() const{
        Matriz transfHol(*this);
        double tH_G{0},tH_a21,tH_r ;
        int stop = transfHol.get_renglones() ;
        Matriz vec_w( stop,1 );

        for(int col_actual=0; col_actual < stop-2 ; col_actual++ ){
            tH_a21 = transfHol.get_val(col_actual+1,col_actual) ;

            tH_G = 0;
            for(int i=col_actual+1; i< stop; i++ ){
                tH_G += transfHol.get_val(i,col_actual) * transfHol.get_val(i,col_actual);
            }
            tH_G = sqrt(tH_G);
            tH_G *= (tH_a21>=0 ? 1:-1);

            tH_r = sqrt( (1.0/2)*tH_G*tH_G + (1.0/2)*tH_a21*tH_G );

            vec_w.replace_Mat(0);
            vec_w.change_val(col_actual+1,0, (tH_a21 + tH_G)/(2*tH_r) );
            for(int i=col_actual+2; i< stop; i++ ){
                vec_w.change_val(i,0, (transfHol.get_val(i,col_actual) )/(2*tH_r) );
            }
            Matriz th_P(  (Matriz(stop,stop).replace_Mat("I")) - (2.0*( vec_w*(vec_w.transpuesta()) ) ) );
            //*
            std::cout << "Hi" << std::endl << transfHol << std::endl
                        << std::endl << vec_w << std::endl << vec_w.transpuesta()
                        << std::endl << vec_w*vec_w.transpuesta()
                        << std::endl << th_P << std::endl
                        << std::endl << (th_P*transfHol*th_P)
                        << std::endl << (*this).get_Traza()
                        << std::endl << "G: " << tH_G << "   r: " << tH_r<< std::endl ;
            //*/ /*!< DEBUGGING */
            /*
            std::cout << "Hi" << std::endl << transfHol << std::endl
                      //  << std::endl << Matriz(stop,stop).replace_Mat("I")
                        << std::endl << vec_w << std::endl << vec_w.transpuesta()
                        << std::endl << 2.0*( vec_w*(vec_w.transpuesta()))
                        << std::endl << th_P << std::endl << ((th_P*transfHol)*th_P)
                        << std::endl << (th_P*transfHol*th_P)
                        << std::endl << (*this).get_Traza() << std::endl << (th_P*transfHol*th_P).get_Traza() ;
            //*/ /*!< DEBUGGING */
            transfHol = th_P*transfHol*th_P;
        }
    return transfHol ;}

    Matriz Matriz::valspropios() const{
        return  Matriz::valspropios(25 + ((*this).get_renglones()*25), 0.000000001, "" );
    }/**This direct use is arbitrarily aimed to optimal precision*/
    Matriz Matriz::valspropios(int iteraciones,double tolerancia,std::string arg) const{
        //*
        Matriz vals_QR((*this).transformacionHouseholder());
        /*/
        Matriz vals_QR((*this));
        //*/ /*!< 2nd one is without transformation to see the process*/
            //std::cout << "Hi" << std::endl << vals_QR; /*!< DEBUGGING */
        Matriz vals_Q( vals_QR.get_renglones(),vals_QR.get_columnas() );
        double costheta,sentheta, mayor_temp, aprox_toler;
        double traza_original{vals_QR.get_Traza() };
        double temp_ajj,temp_aij, raiz_aij2_ajj2; /*to avoid multiple calls*/
        int big_row,big_col, actual_iter{1};

        while(true){
            //std::cout << "Hi" << std::endl << vals_QR;
            mayor_temp = 0;
            for(int i = 0; i < ren; i++){
                for (int j = 0; j < i; j++){ /*checando debajo de la diagonal*/
                    if( std::fabs( vals_QR.get_val(i,j)) > std::fabs(mayor_temp)){
                        std::cout << std::endl << vals_QR.get_val(i,j) << " >"
                                << mayor_temp;
                        mayor_temp = vals_QR.get_val(i,j);
                        big_row = i;
                        big_col = j;
                    }
                }
            }
            temp_ajj = vals_QR.get_val(big_col, big_col);
            temp_aij = vals_QR.get_val(big_row, big_col);
            raiz_aij2_ajj2 = sqrt( temp_aij*temp_aij + temp_ajj*temp_ajj  );
            costheta = ( temp_ajj/raiz_aij2_ajj2 );
            sentheta = ( temp_aij/raiz_aij2_ajj2 );

            vals_Q.replace_Mat("I");
            vals_Q.change_val(big_row,big_row, costheta);
            vals_Q.change_val(big_col,big_col, costheta);

            vals_Q.change_val(big_row,big_col, sentheta);
            vals_Q.change_val(big_col,big_row, -1*sentheta);/*creamos la matriz Q*/
            //*
            std::cout << std::endl<<"Hi" << std::endl << big_row << " " <<big_col << std::endl
                        << std::endl << vals_Q << std::endl << vals_Q.transpuesta()
                        << std::endl << (vals_Q.transpuesta())*vals_QR
                        << ((vals_Q.transpuesta())*vals_QR*vals_Q).get_Traza()  ;
            //*/ /*!< DEBUGGING */

            vals_QR = (vals_Q.transpuesta())*vals_QR*vals_Q; /*matriz aproximacion*/

            aprox_toler = std::fabs( vals_QR.get_val(big_row,big_col ));

            if(arg.compare("v") == 0)
            std::cout << std::endl << "+++++++++++Matriz aproximacion en iteracion "
                        << actual_iter <<": " << std::endl << vals_QR
                        << "Con aproximacion en valor: "
                        << std::setprecision(15) << aprox_toler <<"  ++++"<<std::endl;
            if( ! nearly_Equal(vals_QR.get_Traza() ,traza_original ) ){
                if(arg.compare("v") == 0)
                std::cerr << std::endl << "El error de redondeo ha crecido mucho: " <<std::endl
                            << std::setprecision(15) << vals_QR.get_Traza() <<" != "
                            << traza_original <<std::endl ;
                break;
            }
            if( aprox_toler < tolerancia ){
                if(arg.compare("v") == 0)
                std::cout << std::endl << "Se alcanzo la tolerancia: " <<std::endl
                            << std::setprecision(15) << aprox_toler <<" < "<< tolerancia
                            <<std::endl ;
                break;
            }
            if( actual_iter >= iteraciones ){ /*should stop when its equal*/
                if(arg.compare("v") == 0)
                std::cout << std::endl << "Se alcanzo la  iteracion buscada: " << iteraciones<<std::endl ;
                break;
            }
            actual_iter++;
        }
        Matriz vals_propios(vals_QR.get_renglones(),1);
        for(int i = 0; i < ren; i++){
            vals_propios.change_val(i,0,vals_QR.get_val(i,i) );
        }
    return vals_propios;}


    /**************************************************
                        Sobrecarga de operadores
    ****************************************************/
                /**< operadores de asignacion */
        /* TODO: TEST 2nd case and replace in the other operators.
    Matriz &Matriz::operator=(const Matriz &B){
        if((*this).ren != B.ren ) return *this;
        if((*this).col != B.col) return *this;
        /*/
    Matriz Matriz::operator=(const Matriz &B){
        if((*this).ren != B.ren || (*this).col != B.col){
            (*this).replace_Mat(std::nan(""));
        return *this;}
        //*/ /*2nd case is optimal error return value, 1st case is what allows us to use &operator*/
        /*Los demas operadores de asignacion se definen como el segundo caso
            para tener un error rastreable en caso de suceder*/

        for(int i=0; i<(*this).ren; i++){
            for(int j=0; j<(*this).col; j++){
                (*this).mat[i][j] = B.mat[i][j];
            }
        }
    return *this;}

    Matriz Matriz::operator*=(const Matriz &B){ /*NECESITAN SER CUADRADAS y mismo tamanio*/
        if((*this).ren != B.ren || (*this).col != B.col){
            (*this).replace_Mat(std::nan(""));
        return *this;} /*Error: not same size*/

        Matriz aux(*this);
        (*this) = ((*this)* B);

    return *this;}

    Matriz Matriz::operator+=(const Matriz &B){
        if((*this).ren != B.ren) return Matriz(0,0); /*Error: not same size*/
        if((*this).col != B.col) return Matriz(0,0); /*Error: not same size*/
        (*this) = (*this) + B;

    return *this;}

    Matriz Matriz::operator-=(const Matriz &B){
        if((*this).ren != B.ren) return Matriz(0,0); /*Error: not same size*/
        if((*this).col != B.col) return Matriz(0,0); /*Error: not same size*/
        (*this) = (*this) - B;

    return *this;}

                /**< operadores aritmeticos */
    Matriz operator*(const double &c,const Matriz &A) {/*TODO: Template*/
        Matriz temp(A);
        for(int i=0; i< A.get_renglones() ; i++){
            for(int j=0; j<A.get_columnas() ; j++){
                temp.change_val(i,j, temp.get_val(i,j) * c);
            }
        }
    return temp;}

    Matriz operator*(const Matriz &A,const Matriz &B) {
        int A_ren = A.get_renglones(), A_col = A.get_columnas()
            ,B_ren = B.get_renglones(), B_col = B.get_columnas();
        if(A_col != B_ren) return Matriz(0,0); /*Error: not same size*/
        Matriz temp( A_ren,B_col );
        /*
        std::cout << std::endl << A <<" * " << B;
        //*/ /**< DEBUGGING */
        double resultado;
        for(int r=0; r< A_ren; r++){
            for(int c=0; c< B_col; c++){
                resultado =0;
                for(int i=0; i< A_col; i++){
                    resultado += A.get_val(r,i) * B.get_val(i,c);
                    /*
                        printf("\nA(%d,%d)= %f , ",r,i,A.get_val(r,i));
                        printf("B(%d,%d) = %f ",i,c,B.get_val(i,c));
                    std::cout << "  , A*B resultado: " <<resultado;
                    //*/ /**< DEBUGGING */
                }
                temp.change_val(r,c, resultado);
                    /*
                    std::cout << std::endl<<" ,r: " << r << " ,c: "<<c
                                <<" ,temp: "<< temp;
                    //*/ /**< DEBUGGING */
            }
        }

    return temp;}

    Matriz operator+(const Matriz &A,const Matriz &B) {
        if(A.get_renglones() != B.get_renglones()) return Matriz(0,0); /*Error: not same size*/
        if(A.get_columnas() != B.get_columnas()) return Matriz(0,0); /*Error: not same size*/
        Matriz temp(A.get_renglones(),A.get_columnas());
        for(int i=0; i<temp.get_renglones(); i++){
            for(int j=0; j<temp.get_columnas(); j++){
                temp.change_val(i,j, A.get_val(i,j) + B.get_val(i,j) );
            }
        }

    return temp;}

    Matriz operator-(const Matriz &A,const Matriz &B) {
        if( A.get_renglones() != B.get_renglones() ) return Matriz(0,0); /*Error: not same size*/
        if( A.get_columnas() != B.get_columnas() ) return Matriz(0,0); /*Error: not same size*/
        Matriz temp(A.get_renglones(),A.get_columnas());
        for(int i=0; i<temp.get_renglones(); i++){
            for(int j=0; j<temp.get_columnas(); j++){
                temp.change_val(i,j, A.get_val(i,j) - B.get_val(i,j) );
            }
        }

    return temp;}

                /**< operadores relacionales  */
    bool operator==(const Matriz &A,const Matriz &B) {
        if(A.get_columnas() != B.get_columnas()) return false; /*Error: not same size*/
        if(A.get_renglones() != B.get_renglones()) return false; /*Error: not same size*/

        bool ans = true;
        for(int i=0; i<A.get_renglones(); i++){
            for(int j=0; j<A.get_renglones(); j++){
             ans = (ans && nearly_Equal(A.get_val(i,j),B.get_val(i,j))  );
            }
        }

    return ans;}

    bool operator!=(Matriz &A,Matriz &B) {
    return !(A == B);}
}
