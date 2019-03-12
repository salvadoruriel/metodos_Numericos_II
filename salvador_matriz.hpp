#ifndef salvador_matrices
#define salvador_matrices
#include <cmath>
#include <string>
#include <iostream>
/*WARNING: Usa algunos detalles de c++11*/
namespace salvador{

    static constexpr double MAT_TOL_EPSILON{0.0000001}; /* Tolerancia para doubles,
                                        0.0000001 como valor esperado en uso academico comun*/
    static const double C_Euler = std::exp(1.0); /*safer than using a macro M_E*/

    static bool nearly_Equal(const double &left, const double &right){
        return (std::fabs(left - right)
                <= MAT_TOL_EPSILON * std::max(1.0,std::max(std::fabs(left), std::fabs(right)) )  );
    } /*a relative and absolute error aproach, for big and small numbers*/
    static bool nearly_Equal(const double &left, const double &right, const double &tolerancia){
        return (std::fabs(left - right)
                <= tolerancia * std::max(1.0,std::max(std::fabs(left), std::fabs(right)) )  );
    } /*a relative and absolute error aproach, for big and small numbers*/

	class Matriz{
		int ren,col; /*renglones,columnas, contando desde 1*/
		double** mat;
		double det; /*determinante*/

		public:/*sobrecarga de operadores al final*/
			Matriz(int oren,int ocol);
			Matriz(int oren,int ocol,double oval);
			Matriz(const Matriz &A,int minren, int mincol); /*Matriz reducida en minren y mincol*/
			Matriz(const Matriz &A);
			~Matriz();
            friend std::ostream& operator<<(std::ostream &,const Matriz &);
			/**Operaciones logicas*/
			bool es_cuadrada() const;
			bool es_simetrica() const;
			bool es_EDD() const; /*EDD: Estrictamente dominante diagonalmente */
			bool es_DD() const; /*Menciona si es dominante diagonalmente*/
            /*SUGGESTION: una sola funcion que retorna si es DD,EDD o ninguno,
                                                    como 1,2 o 0 respectivamente*/


            /**************************************************
                    Operaciones sobre la matriz o sus valores
            ****************************************************/
            /** \brief reemplaza una matriz con un vector de valores
             *
             * \param data[] double, vector con valores
             * \param datasize size_t, tamanio del vector
             * \return int,  retorna cantidad de elementos faltantes(negativo),
             *			    o sobrantes(positivo), si fueron exactos, el return es 0;
             */
			int replace_Mat(double data[],size_t datasize); /*to change it with a vector of values*/
			Matriz replace_Mat(double data);/*replace all slots with data, sirven para clases anonimas*/
			Matriz replace_Mat(std::string type);/*matrices especiales*/
			//*
            Matriz IngresaMatriz();
            //*/ /* \WARNING :non-standard, works only on CONSOLE. */
			Matriz mpow(Matriz a,int x); /*TODO: use double exp.*/
			int change_val(int ren_buscado,int col_buscada,double nuevo_val);
			void determinante_Mat();
			void intercambia_rens(unsigned int original,unsigned int cambio);
			int hacer_DD(); /*Intenta hacerla dominante diagonalmente, */
			int hacer_DD_con_VECTOR(Matriz &vec); /*Mueve tambien el vector independiente*/

            /**************************************************
                        RETORNO o impresion de datos de la matriz
            ****************************************************/
			int print_Mat(int busca_ren,int busca_col) const;
			//*
			int print_full_Mat() const;
			//*/ /*DEPRECATED: use operator << and send it to stream*/
			double * toArray(double data[]) const;
            /** \brief manda el renglon de numeros de la matriz al vector data
             *
             * \param data[] double: vector donde se alojaraon los valores del renglon
             * \param ren_buscado unsigned int
             * \return double* un puntero al vector
             * \TODO regresar un puntero a un vector de datos creado localmente en la func
             */
			double * da_ren(double data[], unsigned int ren_buscado) const;
			double get_det(); /*Unica excepcion de const, queremos que se actualice el det*/
			double get_Solo_det() const; /*Usarlo solo si ya se calculo por get_det*/
			double get_val(int ren_buscado,int col_buscada) const;
			double get_norma_ESPECTRAL() const; /*regresa el VALOR de mayor magnitud*/
			double get_Biggest() const; /*regresa el numero con SIGNO de mayor magnitud*/
			double get_Traza() const;

			int get_renglones() const;
			int get_columnas() const;

			void LU_Doolittle(Matriz &,Matriz &) const;
			int LU_DoolittlePROC(Matriz &, Matriz &) const;
			Matriz LU_Cholesky() const;
			void LU_Crout(Matriz &,Matriz &) const;
			int LU_CroutPROC(Matriz &, Matriz &) const;

			Matriz transpuesta() const;
			Matriz menores() const;
			Matriz cofactores() const;
			Matriz adjunta() const;
			Matriz inversa() const;
			Matriz transformacionHouseholder() const;
            /*! \brief da los valores propios de la matriz por metodo QR
             * \return el vector de los valores propios de la matriz
             */
			Matriz valspropios() const;
			Matriz valspropios(int iteraciones,double tolerancia,std::string arg) const;

            /**************************************************
                            Sobrecarga de operadores
            ****************************************************/
                /**< operadores de asignacion */
          /*Matriz &operator=(const Matriz &B);*//**< Leer comentario en el cpp */
			Matriz operator=(const Matriz &B);
			Matriz operator*=(const Matriz &B);
			Matriz operator+=(const Matriz &B);
			Matriz operator-=(const Matriz &B);
                /**< operadores aritmeticos */
			friend Matriz operator*(const double &c,const Matriz &A) ;/*TODO: template*/
			friend Matriz operator*(const Matriz &A,const Matriz &B) ;
			friend Matriz operator+(const Matriz &A,const Matriz &B) ;
			friend Matriz operator-(const Matriz &A,const Matriz &B) ;
                /**< operadores relacionales  */
			friend bool operator==(const Matriz &A,const Matriz &B) ;
			friend bool operator!=(Matriz &A,Matriz &B) ;
	};
}
#endif
