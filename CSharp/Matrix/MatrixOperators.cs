using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace EmsMatrix
{
    public partial class Matrix<T>
    {
        public static Matrix<T> operator +(Matrix<T> m1, Matrix<T> m2)
        {
            return m1.Plus(m2);
        }
        public static Matrix<T> operator +(Matrix<T> m1, Double d)
        {
            return m1.Plus(d);
        }
        public static Matrix<T> operator +(Double d, Matrix<T> m1)
        {            
            return m1.Plus(d);
        }

        public static Matrix<T> operator -(Matrix<T> m1, Matrix<T> m2)
        {
            return m1.Minus(m2);
        }
        public static Matrix<T> operator -(Matrix<T> m1, Double d)
        {
            return m1.Minus(d);
        }
        public static Matrix<T> operator -(Double d, Matrix<T> m1)
        {
            return m1.Minus(d);
        }

        public static Matrix<T> operator *(Matrix<T> m1, Matrix<T> m2)
        {
            return m1.MTimes(m2);
        }
        public static Matrix<T> operator *(Matrix<T> m, Double d)
        {
            return m.Times(d);
        }
        public static Matrix<T> operator *(Double d, Matrix<T> m)
        {
            return m.Times(d);
        }

        public static Matrix<T> operator /(Matrix<T> m1, Matrix<T> m2)
        {
            throw new NotImplementedException();
        }
        public static Matrix<T> operator /(Matrix<T> m, Double d)
        {
            return m.Times(1/d);
        }
        public static Matrix<T> operator /(Double d, Matrix<T> m)
        {
            return m.Times(1 / d);
        }
    }
}
