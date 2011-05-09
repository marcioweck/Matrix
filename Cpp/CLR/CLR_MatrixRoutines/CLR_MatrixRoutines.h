// CLR_MatrixRoutines.h

#pragma once

#include <iostream>

using namespace System;
using namespace std;

namespace CLR_MatrixRoutines {

	
	
	public ref class MatrixCommon
	{
	public:

		//static void Tranpstart(int rows, int cols);
		/*{
			double* d = new double[rows*cols];
			d[0] = 1;
			d[3] = 2;
			d[6] = 3;
			d[1] = 4;
			d[4] = 5;
			d[7] = 6;
			d[2] = 7;
			d[5] = 8;
			d[8] = 9;

			Class1::Transpose2(d, rows, cols);
		}*/
		static void TransposeNxM(double* mat, int rows, int cols);
		static void TransposeNxM(double* mat, double* matout, int rows, int cols);
		static void TransposeNxN(double* mat, int N);
		static void Plus(double* mat1, double* mat2, int length);
		static void Plus(double* mat, double num, int length);

		static void Minus(double* mat1, double* mat2, int length);
		static void Minus(double* mat, double num, int length);

		static void Times(double* mat1, double* mat2, int length);
		static void Times(double* mat, double num, int length);

		static void MTimes(double* mat1, int rows1, int cols1, double* mat2, int rows2, int cols2, double* outmat);
		//static void MTimes2(double* mat1, int rows1, int cols1, double* mat2, int rows2, int cols2, double* outmat);
		static void Interp2Linear(double* mat, int rows, int cols, double* matout);
		static void LUPDecomposition(double *mat, int rows, int cols, double* out_l, double* out_u, double* out_p);
		//static void Solve(double* l, int lrows, int lcols, double*u, double* p);
		//static void FSolve(double* l, int lrows, int lcols, double* p, int prows, int pcols, double* z);
		static void FSolve(double* lt, int ltrows, int ltcols, double* p, int prows, int pcols, double* z);
		static void BSolve(double* ut, int utrows, int utcols, double* z, int zrows, int zcols, double* x);

		static void Round(double* mat, int length, int decimals);

		static void Diag(int rows, int cols, double* outmat, double* diagvec);
		
		static void Determinant3x3(double* mat, double* det);
		static void Determinant2x2(double* mat, double* det);

		static void Householder(double* mat, int size, double* matout);
		template <typename T> static int sgn(T val);

		static void alert2();

		template <class T>
		static void alert(T in);
		
		//static void alert(double item);
		/*static void alert()
		{
			cout << "alert" << endl;
		}*/
		/*{
			double* ret = new double[rows*cols];

			for (int r = 0; r < rows; r++)
					for (int c = 0; c < cols; c++)
					{
						ret[c*rows+r] = mat[r*cols+c];
						//ret[c*rows
						//ret[r, c] = matrix[r, c];
					}
		 
			/*delete [] mat;
			mat = ret;*/
		/*	delete [] mat;
			//Buffer::BlockCopy(ret, 0, mat, 0, sizeof(double)*rows*cols);
			memcpy(mat, ret, sizeof(double)*rows*cols);
			/*System::BlockCopy(ret, 0, mat, 0, sizeof(double)*rows*cols);
			delete [] ret;*/
		//}
	};

	//generic <class ItemType>
	//template<typename T>
	/*public ref class Mat
	{
	public:
		template <typename T>
		static void alert(T item)
		//static void alert()
		{
			cout << "alert" << endl;
		}
		
		static void alert2()
		{
			cout << "alert2" << endl;
			alert(1.0);
		}
	};
	#pragma make_public(Mat::alert2)*/

	
	
}
