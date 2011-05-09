// This is the main DLL file.

#include "stdafx.h"

#include "CLR_MatrixRoutines.h"




namespace CLR_MatrixRoutines {
/*void  MatrixCommon::Tranpstart(int rows, int cols)
{
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

	//Transpose2(d, rows, cols);
	
}*/
void MatrixCommon::alert2()
{
	cout << "alert2" << endl;
}

template <class T>
void MatrixCommon::alert(T in)
{
	cout << sizeof(in) << endl;
}

void MatrixCommon::TransposeNxN(double* mat, int N)
{
	double tmp;
	
	for (int col=0; col< N; col++)
	{
		for (int row=col; row< N; row++)
		{
			tmp = mat[row+col*N];
			mat[row+col*N] = mat[col+row*N];
			mat[col+row*N] = tmp;
		}
	}
	
}

void MatrixCommon::TransposeNxM(double* mat, double* matout, int rows, int cols)
{
	int k;
	//double* matout = new double[rows*cols];
	for (int r = 0; r < rows; r++)
	{
		k = r*cols;
		for (int c = 0; c < cols; c++)
		{
			//ret[r*cols+c] = mat[c*rows+r];
			matout[k+c] = mat[c*rows+r];
		}		
	}
	
	//memcpy(mat, ret, sizeof(double)*rows*cols);
			
	//delete [] ret;	
}
void MatrixCommon::TransposeNxM(double* mat, int rows, int cols)
{
	double* matout = new double[rows*cols];
	TransposeNxM(mat, matout, rows, cols);
	memcpy(mat, matout, sizeof(double)*rows*cols);
	delete [] matout;
}

void MatrixCommon::Plus(double* mat1, double* mat2, int length)
{
	double* p1 = mat1;
	double* p2 = mat2;
	for(int i=0; i<length; i++, p1++, p2++)
		*p1 = (*p1)+(*p2);
	/*double* p1 = mat1;
	double* p2 = mat2;*/
	/*for (int i=0; i<length; i++)
	{
		//p1[l] = p1[l] + p2[l];
		mat1[i] = mat1[i]-mat2[i];
		/**p1 = *p1+*p2;
		p1++;
		p2++;*/
	//}
	//alert3();
	/*delete p1;
	delete p2;*/
}
void MatrixCommon::Plus(double* mat, double num, int length)
{
	double* p1 = mat;
	for(int i=0; i<length; i++, p1++)
		*p1 = (*p1)+num;
}
void MatrixCommon::Minus(double* mat1, double* mat2, int length)
{
	double* p1 = mat1;
	double* p2 = mat2;
	for(int i=0; i<length; i++, p1++, p2++)
		*p1 = (*p1)-(*p2);
}
void MatrixCommon::Minus(double* mat, double num, int length)
{
	double* p1 = mat;
	for(int i=0; i<length; i++, p1++)
		*p1 = (*p1)-num;
}
void MatrixCommon::Times(double* mat1, double* mat2, int length)
{
	double* p1 = mat1;
	double* p2 = mat2;
	for(int i=0; i<length; i++, p1++, p2++)
		*p1 = (*p1)*(*p2);
}

void MatrixCommon::Times(double* mat, double num, int length)
{
	double* p1 = mat;
	for(int i=0; i<length; i++, p1++)
		*p1 = (*p1)*num;
}

/*void MatrixCommon::MTimes2(double* mat1, int rows1, int cols1, double* mat2, int rows2, int cols2, double* outmat)
{
	double d;
	for (int c=0; c<cols2; c++)
		for (int r = 0; r<rows1; r++)
		{
			d=0;
			for (int t=0;t< cols1; t++)
			{
				d = d + mat1[rows1*t+r]*mat2[rows2*c+t];
			}
			outmat[c*rows1+r] = d;
		}
}*/

//void MatrixCommon::MTimes(double* mat1, int rows1, int cols1, double* mat2, int rows2, int cols2, double* outmat)
void MatrixCommon::MTimes(double* mat1t, int rows1, int cols1, double* mat2, int rows2, int cols2, double* outmat)
{
	/*
	Matrix multiplication
	A*B=C
	A = transpose(mat1);  to be able to access columns sequencially
	B = mat2;
	C = outmat;

	          -- cols1 --                -- cols2 --			               -- cols2 --
	   |   [ a11  a12  a13 ]      |   [ b11  b12  b13 ]       |   [ a11*b11+a21*b21+a31*b31      ...     ]
	rows1  [ a21  a22  a23 ]   rows2  [ b21  b22  b23 ]  = rows1  [ a12*b11+a22*b21+a32*b31      etc     ]
	   |   [ a31  a32  a33 ]      |   [ b31  b32  b33 ]       |   [ a13*b11+a23*b21+a33*b31      ...     ]

	A is accessed sequencially for each column of B and written sequencially to columns of C
	*/


	//double* mat1t = new double[rows1*cols1];
	//double* outmat = new double[rows1, cols2];
	//cout << "ver1" << endl;
	//TransposeNxM(mat1, mat1t,  rows1,  cols1);
	
	double* p1 = mat1t;
	double* p2 = mat2;
	double* o = outmat;
	//double* tmp;
	//double tmp2;
	//int i;
	int t, c, r;
	//for (int c=0; c<cols2; c++)
	for (c=0; c<cols2; c++)
	{		
		p1 = mat1t;
		//tmp = mat2+rows2*c;
		for ( r=0; r<rows1; r++)
		{			
			//p2 = mat2+rows2*c;
			//p2 = tmp;
			p2 = mat2+rows2*c;
			//tmp2 = 0;
			*o = 0;
			for ( t=0; t< rows2; t++, p1++, p2++)
			//i = p1+rows2;
			//while (p1<i)
			//for (p2=tmp; p2<tmp+rows2; p2++)
			{
				//*o = *o + (*p1++)*(*p2++);
				//*o += (*p1++)*(*p2++);
				*o += (*p1)*(*p2);
				//*o += (*p1++)*6.7;
				//tmp2 += (*p1++)*(*p2++);
				//*o += (*p1++)+(*p2);
				//*o += (*p1++)+p2[t];
				//p1++;
				//p2++;
			}
			//*o = tmp2;
			o++;			
		}
	}
	/*delete p1;
	delete p2;
	delete o;
	delete tmp;*/

	/*for (int c=0; c<cols1; c++)
	{
		for (int r=0; r<rows1; r++)
		{
			//p1 = mat1t+sizeof(double)*cols1*r;

		}
	}*/
}

void MatrixCommon::Interp2Linear(double* mat, int rows, int cols, double* matout)
{
	double *mp;
	double *op;
	for (int c=0; c<cols; c++)
	{
		mp = mat+c*rows;
		op = matout+(c*2)*(rows*2-1);
		for (int r=0; r<rows-1; r++)
		{
			*op = *mp;
			op++;
			*op = (*(mp+1)-*mp)/2.0+*mp;

			mp++; op++;
		}
		*op = *mp;

	}

	double* opl, *opr;
	for (int c=0; c<cols-1; c++)
	{
		opl = matout+(c*2)*(rows*2-1);
		op = opl+rows*2-1;
		opr = op + rows*2-1;

		for (int r=0; r<(rows*2-1); r++, op++, opl++, opr++)
		{
			*op = (*opr-*opl)/2+*opl;
		}
	}
}

void MatrixCommon::Determinant3x3(double* mat, double* det)
{
	*det = mat[0]*(mat[4]*mat[8]-mat[5]*mat[7])-mat[3]*(mat[1]*mat[8]-mat[2]*mat[7])+mat[6]*(mat[1]*mat[5]-mat[2]*mat[4]);
}

void MatrixCommon::Determinant2x2(double* mat, double* det)
{
	*det = mat[0]*mat[3]-mat[1]*mat[2];
}

void MatrixCommon::Householder(double* mat, int n, double* matout)
{
	for (int k=0; k< n-2;k++)
	{
		double alpha = 0;
		double* ap = mat +k*n+k+1; //A(k+1,k)
		for (int i=k+1; i<n; i++)
		{
			alpha = alpha + (*ap)*(*ap);
			ap++;
		}
		alpha = -sgn<double>(mat[k*n+k])*sqrt(alpha);

		double r;
		r = sqrt(0.5*alpha*(alpha-mat[k*n+k+1]));

		double* v = new double[n];
		double* vp = v;

		ap = mat +k*n; //A(1,k)		
		for (int i=0;i<n; i++, vp++, ap++)
		{
			if (i<k+1)
				*vp=0.0;
			else
				if (i==k+1)
					*vp=(*ap-alpha)/2/r;
				else
					*vp = (*ap)/2/r;
			
		}
		
		double* vmat = new double[n*n];
		MTimes(v,  n,  1, v,  1,  n, vmat);
		vp = vmat;
		for (int c=0; c<n; c++)
			for(int r=0; r<n; r++)
			{
				if (c!=r)
					*vp = -2*(*vp);
				else
					*vp = 1-2*(*vp);
				
				vp++;
			}
	    
		double* tmp = new double[n*n];
		// vmat is symmetric so vmat=vmat'
		MTimes(vmat, n, n, mat, n, n, tmp);
		TransposeNxN(tmp, n);
		MTimes(tmp,n,n, vmat, n, n, mat);

		//memcpy(mat, tmp, sizeof(double)*n*n);


		//alert3();
	}
}


/* LUDecomposition
	"mat" is shuffled
*/
void MatrixCommon::LUPDecomposition(double *mat, int rows, int cols, double* out_l, double* out_u, double* out_p)
{
	double* pl;
	double* pc;
	double* pm;

	for (int c=0;c<rows; c++)
	{
		out_p[c*rows+c] = 1;
	}

	for (int r=0; r< min(rows, cols); r++)
	{
		int index=r;
		pc = mat + r*rows+ r;
		double max = fabs(*pc);
		for (int t=r+1; t< rows; t++)
		{
			pc++;
			if (fabs(*pc)>max)
			{
				index = t;
				max = fabs(*pc);
			}
		}
		if (index != r)
		{
			double tmp;
			for (int c=0; c<cols; c++)
			{
				tmp = mat[c*rows+r];
				mat[c*rows+r] = mat[c*rows + index];
				mat[c*rows + index] = tmp;
			}
			for (int c=0; c< min(rows, cols); c++)
			{
				tmp = out_l[c*rows+r];
				out_l[c*rows+r] = out_l[c*rows + index];
				out_l[c*rows+index] = tmp;
			}
			for (int c=0; c<rows; c++)
			{
				tmp = out_p[c*rows+r];
				out_p[c*rows+r] = out_p[c*rows + index];
				out_p[c*rows+index]=tmp;
			}

		}

			

		for (int c=r; c<cols; c++)
		{			
			//for (int t=r+1; t<rows;t++)
			//pl = out_l+r*rows+r+1;
			pl = out_l+r*rows+r;
			*pl++ = 1;

			pm = mat+r+c*rows;

			pc = mat+c*rows+r+1;

			for (int t=r+1; t<rows;t++)
			{

				//pl = out_l+r*rows+t;
				//pm = mat+r+c*rows;
				//pc = mat+c*rows+t;
				if (c==r)
				{
					*pl++ = (*pc)/(*pm);
					*pc++ = 0.0;
				}
				else
				{
					*pc = *pc-(*pl++)*(*pm);
					pc++;
				}
			}
		}
	}

	/*double* pu = out_u;
	pc = mat;*/
	for (int c=0; c< cols;c++)
	{
		memcpy(out_u+min(rows, cols)*c, mat+rows*c, sizeof(double)*min(rows, cols));
	}
}

// Backward solve UX = Z
void MatrixCommon::BSolve(double* ut, int utrows, int utcols, double* z, int zrows, int zcols, double* x)
{
	int urows = utcols;
	int ucols = utrows;
	
	double* xm;
	xm = x+ucols*zcols-1;
	double* zp;
	zp = z+zcols*zrows-1;

	for (int xc =zcols-1; xc>=0; xc--)
	{
		for (int xr=ucols-1; xr>=0; xr--)
		{
			double* up;
			up = ut+ucols*xr+ucols-1;//xr;
			*xm = *zp;
			double* xp;
			xp = x+ucols*xc+ucols-1;
			//for (int ur=0; ur<xr; ur++)
			for (int ur=xr; ur<ucols-1; ur++)
			{
				*xm = *xm - (*up)*(*xp);
				xp--; up--;
			}
			*xm = (*xp)/(*up);
			xm--;
			zp--;
		}		
	}
}


/*	Forward solve LZ = P
	lt = Transpose(L)
	returns Z
	zrows = lcols = ltrows
	zcols = pcols
	*/
void MatrixCommon::FSolve(double* lt, int ltrows, int ltcols, double* p, int prows, int pcols, double* z)
{
	//double* lt = NULL;
	int lrows = ltcols;
	int lcols = ltrows;
	
	
	double *zm = z;
	double *pp = p;

	for (int zc =0; zc<pcols; zc++)
	{
		for (int zr=0; zr<lcols; zr++)
		{
			double* lp;
			lp = lt+lcols*zr;
			*zm = *pp;
			double* zp;
			zp = z+lcols*zc;
			for (int lr=0; lr<zr; lr++)
			{
				*zm = *zm - (*lp)*(*zp);
				zp++; lp++;
			}
			//*zp = (*zp)/(*lp);/// <-- Taka út??
			zm++;
			pp++;
		}		
	}
}

/*void MatrixCommon::FSolve(double* l, int lrows, int lcols, double* p, int prows, int pcols, double* z)
{
	double* lt = NULL;
	if (lrows==lcols)
	{
		TransposeNxN(l, lrows);
		lt = l;
	}
	else
		TransposeNxM(l, lt, lrows, lcols);
	
	// Forward solve LZ = P
	double *zm = z;
	double *pp = p;

	for (int zc =0; zc<pcols; zc++)
	{
		for (int zr=0; zr<lcols; zr++)
		{
			double* lp;
			lp = lt+lcols*zr;
			*zm = *pp;
			double* zp;
			zp = z+lcols*zc;
			for (int lr=0; lr<zr; lr++)
			{
				*zm = *zm - (*lp)*(*zp);
				zp++; lp++;
			}
			*zp = (*zp)/(*lp);
			zm++;
			pp++;
		}		
	}
}*/
/*void MatrixCommon::BSolve(double* l, int lrows, int lcols, double*u, double* p)
{

}*/
/*void MatrixCommon::LUPDecomposition2(double *mat, int rows, int cols, double* out_l, double* out_u, double* out_p)
{
	double* pmat = mat;
	double* pmrow = mat;
	double* pcrow = mat+1; 
	double* pl = out_l + 1;
	for (int r=0; r< min(rows, cols); r++)
	{
		//pl = out_l+1+rows*r+r;
		//pmat = pmat + c;
		//pmrow = mat + (double*)r*(double*)rows+(double*)r;
		//pl = out_l+(double*)r*rows+(double*)r+(double*)1;
		//pcrow = mat + (double*)r*rows+(double*)r+(double*)1;
		for (int t=r+1; t<rows;t++)
		{
			*pl = (*pcrow)/(*pmrow);
			*pcrow = 0;
			pl++;
			pcrow++;
		}

		pcrow = pcrow + 2 + r;
		
		for (int c=r+1; c<cols; c++)
		{
			pl = out_l + rows*r+r;
			pcrow = mat + c*rows + r;
			pmrow = mat + c*rows + r;
			for (int t=r+1; t< rows; t++)
			{
				pl++;
				pcrow++;
				*pcrow = *pcrow-(*pl)*(*pmrow);
				
			}
			
			
		}
	}
}*/

void MatrixCommon::Round(double* mat, int length, int decimals)
{
	for (int i=0; i<length; i++)
		mat[i] = Math::Round(mat[i], decimals);
}

void MatrixCommon::Diag(int rows, int cols, double* outmat, double* diagvec)
{
}


template <typename T>  int MatrixCommon::sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}

static void alert3()
{
	cout << "alert3" << endl;
}

}