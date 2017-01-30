namespace Matrix.NET
{
    /*
        Matrix class in C#
        Written by Ivan Kuckir (ivan.kuckir@gmail.com, http://blog.ivank.net)
        Faculty of Mathematics and Physics
        Charles University in Prague
        (C) 2010
        - updated on 01.06.2014 - Trimming the string before parsing
        - updated on 14.06.2012 - parsing improved. Thanks to Andy!
        - updated on 03.10.2012 - there was a terrible bug in LU, SoLE and Inversion. Thanks to Danilo Neves Cruz for reporting that!
        - updated on 21.01.2014 - multiple changes based on comments -> see git for further info

        This code is distributed under MIT licence.

            Permission is hereby granted, free of charge, to any person
            obtaining a copy of this software and associated documentation
            files (the "Software"), to deal in the Software without
            restriction, including without limitation the rights to use,
            copy, modify, merge, publish, distribute, sublicense, and/or sell
            copies of the Software, and to permit persons to whom the
            Software is furnished to do so, subject to the following
            conditions:

            The above copyright notice and this permission notice shall be
            included in all copies or substantial portions of the Software.

            THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
            EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
            OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
            NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
            HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
            WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
            FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
            OTHER DEALINGS IN THE SOFTWARE.

        See: https://github.com/darkdragon-001/LightweightMatrixCSharp


        NOTE: The Orriginal Code has been modified greatly to add Functionality. I would like to thank the orriginal Author and also the contributors.
    */



    #region Using Statements:



    using System;
    using System.Text;
    using System.Text.RegularExpressions;



    #endregion



    /// <summary>
    /// A Simple Light Weight Matrix: 
    /// </summary>
    public class Matrix
    {


        #region Fields:



        public Matrix L;



        public Matrix U;



        private int[] pi;



        private double detOfP = 1;



        #endregion



        #region Properties:


        /// <summary>
        /// The Number of Columns in the Matrix.
        /// </summary>
        public int Columns { get; set; }


        /// <summary>
        /// Get a Column from the Matrix. Indexes the Matrix Column[i] Where i is the Column Index starting from 0.
        /// </summary>
        public Column Column { get; set; }


        /// <summary>
        /// The Number of Rows in the Matrix.
        /// </summary>
        public int Rows { get; set; }


        /// <summary>
        /// Get a Row from the Matrix. Indexes the Matrix Row[i] Where i is the Row Index starting from 0.
        /// </summary>
        public Row Row { get; set; }


        /// <summary>
        /// The Matrix Base Object.
        /// </summary>
        public double[] mat;



        /// <summary>
        /// Access this matrix as a 2D array
        /// </summary>
        /// <param name="iRow"></param>
        /// <param name="iCol"></param>
        /// <returns></returns>
        public double this[int iRow, int iCol]
        {
            get { return mat[iRow * Columns + iCol]; }
            set { mat[iRow * Columns + iCol] = value; }
        }




        #endregion



        /// <summary>
        /// Matrix Constructor
        /// </summary>
        /// <param name="iRows"></param>
        /// <param name="iCols"></param>
        public Matrix(int rows, int cols)
        {
            Rows = rows;
            Columns = cols;
            mat = new double[Rows * Columns];

            Row = new Row(this);
            Column = new Column(this);
        }



        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public Boolean IsSquare()
        {
            return (Rows == Columns);
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="k"></param>
        /// <returns></returns>
        public Matrix GetCol(int k)
        {
            Matrix m = new Matrix(Rows, 1);
            for (int i = 0; i < Rows; i++) m[i, 0] = this[i, k];
            return m;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="v"></param>
        /// <param name="k"></param>
        public void SetCol(Matrix v, int k)
        {
            for (int i = 0; i < Rows; i++) this[i, k] = v[i, 0];
        }



        /// <summary>
        /// Function for LU decomposition
        /// </summary>
        public void MakeLU()
        {
            if (!IsSquare()) throw new MException("The matrix is not square!");
            L = IdentityMatrix(Rows, Columns);
            U = Duplicate();

            pi = new int[Rows];
            for (int i = 0; i < Rows; i++) pi[i] = i;

            double p = 0;
            double pom2;
            int k0 = 0;
            int pom1 = 0;

            for (int k = 0; k < Columns - 1; k++)
            {
                p = 0;
                for (int i = k; i < Rows; i++)      // find the row with the biggest pivot
                {
                    if (Math.Abs(U[i, k]) > p)
                    {
                        p = Math.Abs(U[i, k]);
                        k0 = i;
                    }
                }
                if (p == 0) // samé nuly ve sloupci
                    throw new MException("The matrix is singular!");

                pom1 = pi[k]; pi[k] = pi[k0]; pi[k0] = pom1;    // switch two rows in permutation matrix

                for (int i = 0; i < k; i++)
                {
                    pom2 = L[k, i]; L[k, i] = L[k0, i]; L[k0, i] = pom2;
                }

                if (k != k0) detOfP *= -1;

                for (int i = 0; i < Columns; i++)                  // Switch rows in U
                {
                    pom2 = U[k, i]; U[k, i] = U[k0, i]; U[k0, i] = pom2;
                }

                for (int i = k + 1; i < Rows; i++)
                {
                    L[i, k] = U[i, k] / U[k, k];
                    for (int j = k; j < Columns; j++)
                        U[i, j] = U[i, j] - L[i, k] * U[k, j];
                }
            }
        }



        /// <summary>
        /// Function solves Ax = v in confirmity with solution vector "v"
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public Matrix SolveWith(Matrix v)
        {
            if (Rows != Columns) throw new MException("The matrix is not square!");
            if (Rows != v.Rows) throw new MException("Wrong number of results in solution vector!");
            if (L == null) MakeLU();

            Matrix b = new Matrix(Rows, 1);
            for (int i = 0; i < Rows; i++) b[i, 0] = v[pi[i], 0];   // switch two items in "v" due to permutation matrix

            Matrix z = SubsForth(L, b);
            Matrix x = SubsBack(U, z);

            return x;
        }



        /// <summary>
        /// Function makes reduced echolon form
        /// TODO check for redundancy with MakeLU() and SolveWith()
        /// </summary>
        public void MakeRref()
        {
            int lead = 0;
            for (int r = 0; r < Rows; r++)
            {
                if (Columns <= lead) break;
                int i = r;
                while (this[i, lead] == 0)
                {
                    i++;
                    if (i == Rows)
                    {
                        i = r;
                        lead++;
                        if (Columns == lead)
                        {
                            lead--;
                            break;
                        }
                    }
                }
                for (int j = 0; j < Columns; j++)
                {
                    double temp = this[r, j];
                    this[r, j] = this[i, j];
                    this[i, j] = temp;
                }
                double div = this[r, lead];
                for (int j = 0; j < Columns; j++) this[r, j] /= div;
                for (int j = 0; j < Rows; j++)
                {
                    if (j != r)
                    {
                        double sub = this[j, lead];
                        for (int k = 0; k < Columns; k++) this[j, k] -= (sub * this[r, k]);
                    }
                }
                lead++;
            }
        }



        /// <summary>
        /// Function returns the inverted matrix
        /// </summary>
        /// <returns></returns>
        public Matrix Invert()
        {
            if (L == null) MakeLU();

            Matrix inv = new Matrix(Rows, Columns);

            for (int i = 0; i < Rows; i++)
            {
                Matrix Ei = Matrix.ZeroMatrix(Rows, 1);
                Ei[i, 0] = 1;
                Matrix col = SolveWith(Ei);
                inv.SetCol(col, i);
            }
            return inv;
        }



        /// <summary>
        /// Function for determinant
        /// </summary>
        /// <returns></returns>
        public double Det()
        {
            if (L == null) MakeLU();
            double det = detOfP;
            for (int i = 0; i < Rows; i++) det *= U[i, i];
            return det;
        }



        /// <summary>
        /// Function returns permutation matrix "P" due to permutation vector "pi"
        /// </summary>
        /// <returns></returns>
        public Matrix GetP()
        {
            if (L == null) MakeLU();

            Matrix matrix = ZeroMatrix(Rows, Columns);
            for (int i = 0; i < Rows; i++) matrix[pi[i], i] = 1;
            return matrix;
        }



        /// <summary>
        /// Function returns the copy of this matrix
        /// </summary>
        /// <returns></returns>
        public Matrix Duplicate()
        {
            Matrix matrix = new Matrix(Rows, Columns);
            for (int i = 0; i < Rows; i++)
                for (int j = 0; j < Columns; j++)
                    matrix[i, j] = this[i, j];
            return matrix;
        }



        /// <summary>
        /// Function solves Ax = b for A as a lower triangular matrix
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Matrix SubsForth(Matrix A, Matrix b)
        {
            if (A.L == null) A.MakeLU();
            int n = A.Rows;
            Matrix x = new Matrix(n, 1);

            for (int i = 0; i < n; i++)
            {
                x[i, 0] = b[i, 0];
                for (int j = 0; j < i; j++) x[i, 0] -= A[i, j] * x[j, 0];
                x[i, 0] = x[i, 0] / A[i, i];
            }
            return x;
        }



        /// <summary>
        /// Function solves Ax = b for A as an upper triangular matrix
        /// </summary>
        /// <param name="A"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static Matrix SubsBack(Matrix A, Matrix b)
        {
            if (A.L == null) A.MakeLU();
            int n = A.Rows;
            Matrix x = new Matrix(n, 1);

            for (int i = n - 1; i > -1; i--)
            {
                x[i, 0] = b[i, 0];
                for (int j = n - 1; j > i; j--) x[i, 0] -= A[i, j] * x[j, 0];
                x[i, 0] = x[i, 0] / A[i, i];
            }
            return x;
        }



        /// <summary>
        /// Function generates the zero matrix
        /// </summary>
        /// <param name="iRows"></param>
        /// <param name="iCols"></param>
        /// <returns></returns>
        public static Matrix ZeroMatrix(int iRows, int iCols)
        {
            Matrix matrix = new Matrix(iRows, iCols);
            for (int i = 0; i < iRows; i++)
                for (int j = 0; j < iCols; j++)
                    matrix[i, j] = 0;
            return matrix;
        }



        /// <summary>
        /// Function generates the identity matrix
        /// </summary>
        /// <param name="iRows"></param>
        /// <param name="iCols"></param>
        /// <returns></returns>
        public static Matrix IdentityMatrix(int iRows, int iCols)
        {
            Matrix matrix = ZeroMatrix(iRows, iCols);
            for (int i = 0; i < Math.Min(iRows, iCols); i++)
                matrix[i, i] = 1;
            return matrix;
        }



        /// <summary>
        /// Function generates the random matrix
        /// </summary>
        /// <param name="iRows"></param>
        /// <param name="iCols"></param>
        /// <param name="dispersion"></param>
        /// <returns></returns>
        public static Matrix RandomMatrix(int iRows, int iCols, int dispersion)
        {
            Random random = new Random();
            Matrix matrix = new Matrix(iRows, iCols);
            for (int i = 0; i < iRows; i++)
                for (int j = 0; j < iCols; j++)
                    matrix[i, j] = random.Next(-dispersion, dispersion);
            return matrix;
        }



        /// <summary>
        /// Function parses the matrix from string
        /// </summary>
        /// <param name="ps"></param>
        /// <returns></returns>
        public static Matrix Parse(string ps)
        {
            string s = NormalizeMatrixString(ps);
            string[] rows = Regex.Split(s, "\r\n");
            string[] nums = rows[0].Split(' ');
            Matrix matrix = new Matrix(rows.Length, nums.Length);
            try
            {
                for (int i = 0; i < rows.Length; i++)
                {
                    nums = rows[i].Split(' ');
                    for (int j = 0; j < nums.Length; j++) matrix[i, j] = double.Parse(nums[j]);
                }
            }
            catch (FormatException) { throw new MException("Wrong input format!"); }
            return matrix;
        }



        /// <summary>
        /// Function returns matrix as a string
        /// </summary>
        /// <returns></returns>
        public override string ToString()
        {
            StringBuilder s = new StringBuilder();
            for (int i = 0; i < Rows; i++)
            {
                for (int j = 0; j < Columns; j++)
                    s.Append(String.Format("{0,5:E2}", this[i, j]) + " ");
                s.AppendLine();
            }
            return s.ToString();
        }



        /// <summary>
        /// Matrix transpose, for any rectangular matrix
        /// </summary>
        /// <param name="m"></param>
        /// <returns></returns>
        public static Matrix Transpose(Matrix m)
        {
            Matrix t = new Matrix(m.Columns, m.Rows);
            for (int i = 0; i < m.Rows; i++)
                for (int j = 0; j < m.Columns; j++)
                    t[j, i] = m[i, j];
            return t;
        }



        /// <summary>
        /// Power matrix to exponent
        /// </summary>
        /// <param name="m"></param>
        /// <param name="pow"></param>
        /// <returns></returns>
        public static Matrix Power(Matrix m, int pow)
        {
            if (pow == 0) return IdentityMatrix(m.Rows, m.Columns);
            if (pow == 1) return m.Duplicate();
            if (pow == -1) return m.Invert();

            Matrix x;
            if (pow < 0) { x = m.Invert(); pow *= -1; }
            else x = m.Duplicate();

            Matrix ret = IdentityMatrix(m.Rows, m.Columns);
            while (pow != 0)
            {
                if ((pow & 1) == 1) ret *= x;
                x *= x;
                pow >>= 1;
            }
            return ret;
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="xa"></param>
        /// <param name="ya"></param>
        /// <param name="B"></param>
        /// <param name="xb"></param>
        /// <param name="yb"></param>
        /// <param name="C"></param>
        /// <param name="size"></param>
        private static void SafeAplusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
        {
            for (int i = 0; i < size; i++)          // rows
                for (int j = 0; j < size; j++)     // cols
                {
                    C[i, j] = 0;
                    if (xa + j < A.Columns && ya + i < A.Rows) C[i, j] += A[ya + i, xa + j];
                    if (xb + j < B.Columns && yb + i < B.Rows) C[i, j] += B[yb + i, xb + j];
                }
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="xa"></param>
        /// <param name="ya"></param>
        /// <param name="B"></param>
        /// <param name="xb"></param>
        /// <param name="yb"></param>
        /// <param name="C"></param>
        /// <param name="size"></param>
        private static void SafeAminusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
        {
            for (int i = 0; i < size; i++)          // rows
                for (int j = 0; j < size; j++)     // cols
                {
                    C[i, j] = 0;
                    if (xa + j < A.Columns && ya + i < A.Rows) C[i, j] += A[ya + i, xa + j];
                    if (xb + j < B.Columns && yb + i < B.Rows) C[i, j] -= B[yb + i, xb + j];
                }
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="xa"></param>
        /// <param name="ya"></param>
        /// <param name="C"></param>
        /// <param name="size"></param>
        private static void SafeACopytoC(Matrix A, int xa, int ya, Matrix C, int size)
        {
            for (int i = 0; i < size; i++)          // rows
                for (int j = 0; j < size; j++)     // cols
                {
                    C[i, j] = 0;
                    if (xa + j < A.Columns && ya + i < A.Rows) C[i, j] += A[ya + i, xa + j];
                }
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="xa"></param>
        /// <param name="ya"></param>
        /// <param name="B"></param>
        /// <param name="xb"></param>
        /// <param name="yb"></param>
        /// <param name="C"></param>
        /// <param name="size"></param>
        private static void AplusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
        {
            for (int i = 0; i < size; i++)          // rows
                for (int j = 0; j < size; j++) C[i, j] = A[ya + i, xa + j] + B[yb + i, xb + j];
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="xa"></param>
        /// <param name="ya"></param>
        /// <param name="B"></param>
        /// <param name="xb"></param>
        /// <param name="yb"></param>
        /// <param name="C"></param>
        /// <param name="size"></param>
        private static void AminusBintoC(Matrix A, int xa, int ya, Matrix B, int xb, int yb, Matrix C, int size)
        {
            for (int i = 0; i < size; i++)          // rows
                for (int j = 0; j < size; j++) C[i, j] = A[ya + i, xa + j] - B[yb + i, xb + j];
        }



        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="xa"></param>
        /// <param name="ya"></param>
        /// <param name="C"></param>
        /// <param name="size"></param>
        private static void ACopytoC(Matrix A, int xa, int ya, Matrix C, int size)
        {
            for (int i = 0; i < size; i++)          // rows
                for (int j = 0; j < size; j++) C[i, j] = A[ya + i, xa + j];
        }



        /// <summary>
        /// Smart matrix multiplication
        /// TODO assume matrix 2^N x 2^N and then directly call StrassenMultiplyRun(A,B,?,1,?)
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        private static Matrix StrassenMultiply(Matrix A, Matrix B)
        {
            if (A.Columns != B.Rows) throw new MException("Wrong dimension of matrix!");

            Matrix R;

            int msize = Math.Max(Math.Max(A.Rows, A.Columns), Math.Max(B.Rows, B.Columns));

            int size = 1; int n = 0;
            while (msize > size) { size *= 2; n++; };
            int h = size / 2;


            Matrix[,] mField = new Matrix[n, 9];

            /*
             *  8x8, 8x8, 8x8, ...
             *  4x4, 4x4, 4x4, ...
             *  2x2, 2x2, 2x2, ...
             *  . . .
             */

            int z;
            for (int i = 0; i < n - 4; i++)          // rows
            {
                z = (int)Math.Pow(2, n - i - 1);
                for (int j = 0; j < 9; j++) mField[i, j] = new Matrix(z, z);
            }

            SafeAplusBintoC(A, 0, 0, A, h, h, mField[0, 0], h);
            SafeAplusBintoC(B, 0, 0, B, h, h, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 1], 1, mField); // (A11 + A22) * (B11 + B22);

            SafeAplusBintoC(A, 0, h, A, h, h, mField[0, 0], h);
            SafeACopytoC(B, 0, 0, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 2], 1, mField); // (A21 + A22) * B11;

            SafeACopytoC(A, 0, 0, mField[0, 0], h);
            SafeAminusBintoC(B, h, 0, B, h, h, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 3], 1, mField); //A11 * (B12 - B22);

            SafeACopytoC(A, h, h, mField[0, 0], h);
            SafeAminusBintoC(B, 0, h, B, 0, 0, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 4], 1, mField); //A22 * (B21 - B11);

            SafeAplusBintoC(A, 0, 0, A, h, 0, mField[0, 0], h);
            SafeACopytoC(B, h, h, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 5], 1, mField); //(A11 + A12) * B22;

            SafeAminusBintoC(A, 0, h, A, 0, 0, mField[0, 0], h);
            SafeAplusBintoC(B, 0, 0, B, h, 0, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 6], 1, mField); //(A21 - A11) * (B11 + B12);

            SafeAminusBintoC(A, h, 0, A, h, h, mField[0, 0], h);
            SafeAplusBintoC(B, 0, h, B, h, h, mField[0, 1], h);
            StrassenMultiplyRun(mField[0, 0], mField[0, 1], mField[0, 1 + 7], 1, mField); // (A12 - A22) * (B21 + B22);

            R = new Matrix(A.Rows, B.Columns);                  // result

            /// C11
            for (int i = 0; i < Math.Min(h, R.Rows); i++)          // rows
                for (int j = 0; j < Math.Min(h, R.Columns); j++)     // cols
                    R[i, j] = mField[0, 1 + 1][i, j] + mField[0, 1 + 4][i, j] - mField[0, 1 + 5][i, j] + mField[0, 1 + 7][i, j];

            /// C12
            for (int i = 0; i < Math.Min(h, R.Rows); i++)          // rows
                for (int j = h; j < Math.Min(2 * h, R.Columns); j++)     // cols
                    R[i, j] = mField[0, 1 + 3][i, j - h] + mField[0, 1 + 5][i, j - h];

            /// C21
            for (int i = h; i < Math.Min(2 * h, R.Rows); i++)          // rows
                for (int j = 0; j < Math.Min(h, R.Columns); j++)     // cols
                    R[i, j] = mField[0, 1 + 2][i - h, j] + mField[0, 1 + 4][i - h, j];

            /// C22
            for (int i = h; i < Math.Min(2 * h, R.Rows); i++)          // rows
                for (int j = h; j < Math.Min(2 * h, R.Columns); j++)     // cols
                    R[i, j] = mField[0, 1 + 1][i - h, j - h] - mField[0, 1 + 2][i - h, j - h] + mField[0, 1 + 3][i - h, j - h] + mField[0, 1 + 6][i - h, j - h];

            return R;
        }



        /// <summary>
        /// A * B into C, level of recursion, matrix field
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="C"></param>
        /// <param name="l"></param>
        /// <param name="f"></param>
        private static void StrassenMultiplyRun(Matrix A, Matrix B, Matrix C, int l, Matrix[,] f)
        {
            int size = A.Rows;
            int h = size / 2;

            AplusBintoC(A, 0, 0, A, h, h, f[l, 0], h);
            AplusBintoC(B, 0, 0, B, h, h, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 1], l + 1, f); // (A11 + A22) * (B11 + B22);

            AplusBintoC(A, 0, h, A, h, h, f[l, 0], h);
            ACopytoC(B, 0, 0, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 2], l + 1, f); // (A21 + A22) * B11;

            ACopytoC(A, 0, 0, f[l, 0], h);
            AminusBintoC(B, h, 0, B, h, h, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 3], l + 1, f); //A11 * (B12 - B22);

            ACopytoC(A, h, h, f[l, 0], h);
            AminusBintoC(B, 0, h, B, 0, 0, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 4], l + 1, f); //A22 * (B21 - B11);

            AplusBintoC(A, 0, 0, A, h, 0, f[l, 0], h);
            ACopytoC(B, h, h, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 5], l + 1, f); //(A11 + A12) * B22;

            AminusBintoC(A, 0, h, A, 0, 0, f[l, 0], h);
            AplusBintoC(B, 0, 0, B, h, 0, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 6], l + 1, f); //(A21 - A11) * (B11 + B12);

            AminusBintoC(A, h, 0, A, h, h, f[l, 0], h);
            AplusBintoC(B, 0, h, B, h, h, f[l, 1], h);
            StrassenMultiplyRun(f[l, 0], f[l, 1], f[l, 1 + 7], l + 1, f); // (A12 - A22) * (B21 + B22);

            /// C11
            for (int i = 0; i < h; i++)          // rows
                for (int j = 0; j < h; j++)     // cols
                    C[i, j] = f[l, 1 + 1][i, j] + f[l, 1 + 4][i, j] - f[l, 1 + 5][i, j] + f[l, 1 + 7][i, j];

            /// C12
            for (int i = 0; i < h; i++)          // rows
                for (int j = h; j < size; j++)     // cols
                    C[i, j] = f[l, 1 + 3][i, j - h] + f[l, 1 + 5][i, j - h];

            /// C21
            for (int i = h; i < size; i++)          // rows
                for (int j = 0; j < h; j++)     // cols
                    C[i, j] = f[l, 1 + 2][i - h, j] + f[l, 1 + 4][i - h, j];

            /// C22
            for (int i = h; i < size; i++)          // rows
                for (int j = h; j < size; j++)     // cols
                    C[i, j] = f[l, 1 + 1][i - h, j - h] - f[l, 1 + 2][i - h, j - h] + f[l, 1 + 3][i - h, j - h] + f[l, 1 + 6][i - h, j - h];
        }



        /// <summary>
        /// Stupid matrix multiplication
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        private static Matrix StupidMultiply(Matrix m1, Matrix m2)
        {
            if (m1.Columns != m2.Rows) throw new MException("Wrong dimensions of matrix!");

            Matrix result = ZeroMatrix(m1.Rows, m2.Columns);
            for (int i = 0; i < result.Rows; i++)
                for (int j = 0; j < result.Columns; j++)
                    for (int k = 0; k < m1.Columns; k++)
                        result[i, j] += m1[i, k] * m2[k, j];
            return result;
        }



        /// <summary>
        /// Matrix multiplication
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        private static Matrix Multiply(Matrix m1, Matrix m2)
        {
            if (m1.Columns != m2.Rows) throw new MException("Wrong dimension of matrix!");
            int msize = Math.Max(Math.Max(m1.Rows, m1.Columns), Math.Max(m2.Rows, m2.Columns));
            // stupid multiplication faster for small matrices
            if (msize < 32)
            {
                return StupidMultiply(m1, m2);
            }
            // stupid multiplication faster for non square matrices
            if (!m1.IsSquare() || !m2.IsSquare())
            {
                return StupidMultiply(m1, m2);
            }
            // Strassen multiplication is faster for large square matrix 2^N x 2^N
            // NOTE because of previous checks msize == m1.cols == m1.rows == m2.cols == m2.cols
            double exponent = Math.Log(msize) / Math.Log(2);
            if (Math.Pow(2, exponent) == msize)
            {
                return StrassenMultiply(m1, m2);
            }
            else
            {
                return StupidMultiply(m1, m2);
            }
        }



        /// <summary>
        /// Multiplication by constant n
        /// </summary>
        /// <param name="n"></param>
        /// <param name="m"></param>
        /// <returns></returns>
        private static Matrix Multiply(double n, Matrix m)
        {
            Matrix r = new Matrix(m.Rows, m.Columns);
            for (int i = 0; i < m.Rows; i++)
                for (int j = 0; j < m.Columns; j++)
                    r[i, j] = m[i, j] * n;
            return r;
        }



        /// <summary>
        /// Sčítání matic
        /// </summary>
        /// <param name="m1"></param>
        /// <param name="m2"></param>
        /// <returns></returns>
        private static Matrix Add(Matrix m1, Matrix m2)
        {
            if (m1.Rows != m2.Rows || m1.Columns != m2.Columns) throw new MException("Matrices must have the same dimensions!");
            Matrix r = new Matrix(m1.Rows, m1.Columns);
            for (int i = 0; i < r.Rows; i++)
                for (int j = 0; j < r.Columns; j++)
                    r[i, j] = m1[i, j] + m2[i, j];
            return r;
        }



        /// <summary>
        /// From Andy - thank you! :)
        /// </summary>
        /// <param name="matStr"></param>
        /// <returns></returns>
        public static string NormalizeMatrixString(string matStr)
        {
            // Remove any multiple spaces
            while (matStr.IndexOf("  ") != -1)
                matStr = matStr.Replace("  ", " ");

            // Remove any spaces before or after newlines
            matStr = matStr.Replace(" \r\n", "\r\n");
            matStr = matStr.Replace("\r\n ", "\r\n");

            // If the data ends in a newline, remove the trailing newline.
            // Make it easier by first replacing \r\n’s with |’s then
            // restore the |’s with \r\n’s
            matStr = matStr.Replace("\r\n", "|");
            while (matStr.LastIndexOf("|") == (matStr.Length - 1))
                matStr = matStr.Substring(0, matStr.Length - 1);

            matStr = matStr.Replace("|", "\r\n");
            return matStr.Trim();
        }



        //   O P E R A T O R S



        public static Matrix operator -(Matrix m)
        { return Matrix.Multiply(-1, m); }



        public static Matrix operator +(Matrix m1, Matrix m2)
        { return Matrix.Add(m1, m2); }



        public static Matrix operator -(Matrix m1, Matrix m2)
        { return Matrix.Add(m1, -m2); }



        public static Matrix operator *(Matrix m1, Matrix m2)
        { return Matrix.Multiply(m1, m2); }



        public static Matrix operator *(double n, Matrix m)
        { return Matrix.Multiply(n, m); }




    }



}
