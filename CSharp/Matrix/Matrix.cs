using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using CLR_MatrixRoutines;
using System.Globalization;

namespace EmsMatrix
{
    public partial class Matrix<T>
    {
        public Double[] Data;
        public Int32 Rows { get; private set; }
        public Int32 Cols { get; private set; }


        public Matrix(Int32 Rows, Int32 Cols)
        {
            if (typeof(T) != typeof(Double) && typeof(T) != typeof(Single))
                throw new Exception("Only supports doubles and floats");

            Data = new Double[Cols * Rows];
            this.Cols = Cols;
            this.Rows = Rows;
        }

        public Double this[Int32 Row, Int32 Col]
        {
            get
            {
                return Data[Row + Col * Rows];
            }
            set
            {
                Data[Row + Col * Rows] = value;
            }

        }
        public Matrix<T> Transpose()
        {
           
            Matrix<T> ret = this.Clone();

            if (typeof(T) == typeof(Double))
            {
                if (Rows == Cols)
                {
                    unsafe
                    {
                        fixed (double* pointer = ret.Data)
                        {
                            MatrixCommon.TransposeNxN(pointer, Rows);
                        }
                    }
                    Int32 tmp = Rows;
                    Rows = Cols;
                    Cols = tmp;
                }
                else
                {
                    unsafe
                    {
                        fixed (double* pointer = ret.Data)
                        {

                            MatrixCommon.TransposeNxM(pointer, Rows, Cols);

                        }
                    }
                    Int32 tmp = ret.Rows;
                    ret.Rows = ret.Cols;
                    ret.Cols = tmp;
                }
                
            }
            else
                throw new NotImplementedException();

            return ret;
            
        }
        public Matrix<T> Minus(Double d)
        {
            Matrix<T> ret = this.Clone();
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    //double* mat1; double* mat2;
                    fixed (double* mat1 = ret.Data)
                    {
                        MatrixCommon.Minus(mat1, d, Rows * Cols);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public Matrix<T> Minus(Matrix<T> mat)
        {

            Matrix<T> ret = this.Clone();
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    fixed (double* mat1 = ret.Data, mat2 = mat.Data)
                    {
                        MatrixCommon.Minus(mat1, mat2, Rows * Cols);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public Matrix<T> Plus(Double d)
        {
            Matrix<T> ret = this.Clone();
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    //double* mat1; double* mat2;
                    fixed (double* mat1 = ret.Data)
                    {

                        MatrixCommon.Plus(mat1, d, Rows * Cols);
                        //MatrixCommon.MinusC(mat1, 10.0, Rows * Cols);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public Matrix<T> Plus(Matrix<T> mat)
        {
            Matrix<T> ret = this.Clone();
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    //double* mat1; double* mat2;
                    fixed (double* mat1 = ret.Data, mat2 = mat.Data)
                    {

                        MatrixCommon.Plus(mat1, mat2, Rows * Cols);
                        //MatrixCommon.MinusC(mat1, 10.0, Rows * Cols);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }


        /// <summary>
        /// Matrix multiply.
        /// </summary>
        /// <param name="mat"></param>
        /// <returns></returns>
        public Matrix<T> MTimes(Matrix<T> mat)
        {
            if (this.Cols != mat.Rows)
                throw new Exception("MTimes: Inner matrix dimensions must agree.");

            Matrix<T> ret = new Matrix<T>(this.Rows, mat.Cols);
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    fixed (double* mat1 = this.Clone().Transpose().Data, mat2 = mat.Data, outmat = ret.Data)
                    {
                        MatrixCommon.MTimes(mat1, Rows, Cols, mat2, mat.Rows, mat.Cols, outmat);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public Matrix<T> Times(Double num)
        {
            Matrix<T> ret = this.Clone();
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    fixed (double* mat = ret.Data)
                    {
                        MatrixCommon.Times(mat, num, this.Rows * this.Cols);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public Matrix<T> Interp2Linear()
        {
            Matrix<T> ret = new Matrix<T>(this.Rows * 2 - 1, this.Cols * 2 - 1);
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    fixed (double* mat = this.Data, pret = ret.Data)
                    {
                        MatrixCommon.Interp2Linear(mat, this.Rows, this.Cols, pret);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public override string ToString()
        {
            IFormatProvider ifp = CultureInfo.InvariantCulture;
            String tmp;
            Int32 longest = 0;
            for (int r = 0; r < this.Rows; r++)
            {
                for (int c = 0; c < this.Cols; c++)
                {
                    tmp = this[r, c].ToString(ifp);
                    if (tmp.Length > longest)
                        longest = tmp.Length;
                }
            }

            StringBuilder sb = new StringBuilder();
            String[] lines = new String[this.Cols];
            for (int r = 0; r < this.Rows; r++)
            {

                for (int c = 0; c < this.Cols; c++)
                {
                    lines[c] = this[r, c].ToString(ifp).PadRight(longest, ' ');
                }
                sb.Append("[").Append(String.Join(" ", lines)).AppendLine("]");
            }
            //sb.Append(longest.ToString());
            return sb.ToString();

        }


        public static Matrix<T> Householder(Matrix<T> mat)
        {
            Matrix<T> ret = mat.Clone();
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    fixed (double* pm = mat.Data, pr = ret.Data)
                    {
                        MatrixCommon.Householder(pm, mat.Rows, pr);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }
        public Double Determinant()
        {
            Double ret;
            if (typeof(T) == typeof(Double))
            {
                unsafe
                {
                    fixed (double* pm = this.Data)
                    {
                        MatrixCommon.Determinant3x3(pm, &ret);
                    }
                }
            }
            else
                throw new NotImplementedException();
            return ret;
        }

        public Matrix<T> Inverse()
        {
            Matrix<T> l, u, p, x;
            Matrix<T>.LUDecomposition(this, out l, out u, out p);
            Matrix<T>.SolveLUP(l, u, p, out x);
            return x;
        }

        static public Matrix<T> Rand(Int32 Rows, Int32 Cols, Int32 Seed)
        {
            Random rand = new Random(Seed);
            Matrix<T> ret = new Matrix<T>(Rows, Cols);
            for (int c = 0; c < Cols; c++)
                for (int r = 0; r < Rows; r++)
                {
                    ret[r, c] = rand.NextDouble();
                }
            return ret;
        }
        static public Matrix<T> Rand(Int32 Rows, Int32 Cols)
        {
            return Matrix<T>.Rand(Rows, Cols, (int)DateTime.Now.Ticks);
        }

        static public void SolveLUP(Matrix<T> l, Matrix<T> u, Matrix<T> p, out Matrix<T> x)
        {
            Matrix<T> lt = l.Transpose();

            Matrix<T> z = new Matrix<T>(l.Rows, p.Cols);
            unsafe
            {
                fixed (double* plt = lt.Data, pu = u.Data, pp = p.Data, pz = z.Data)
                {
                    MatrixCommon.FSolve(plt, l.Rows, l.Cols, pp, p.Rows, p.Cols, pz);
                    //MatrixCommon.Solve(pl, l.Rows, l.Cols, pu, pp);
                }
            }


            Matrix<T> ut = u.Transpose();
            x = new Matrix<T>(u.Cols, z.Cols);
            unsafe
            {
                fixed (double* utp = ut.Data, zp = z.Data, px = x.Data)
                {
                    MatrixCommon.BSolve(utp, ut.Rows, ut.Cols, zp, z.Rows, z.Cols, px);
                }
            }

        }



        static public void LUDecomposition(Matrix<T> mat, out Matrix<T> l, out Matrix<T> u, out Matrix<T> p)
        {
            l = new Matrix<T>(mat.Rows, Math.Min(mat.Rows, mat.Cols));
            u = new Matrix<T>(Math.Min(mat.Rows, mat.Cols), mat.Cols);
            p = new Matrix<T>(mat.Rows, mat.Rows);
            Matrix<T> mat_clone = mat.Clone();
            unsafe
            {
                fixed (double* m = mat_clone.Data, ml = l.Data, mu = u.Data, mp = p.Data)
                {
                    MatrixCommon.LUPDecomposition(m, mat.Rows, mat.Cols, ml, mu, mp);
                }
            }
        }
        public Matrix<T> Clone()
        {
            Matrix<T> ret = new Matrix<T>(Rows, Cols);
            ret.Data = (Double[])Data.Clone();
            return ret;
        }
    }
}
