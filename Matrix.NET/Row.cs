namespace Matrix.NET
{


    #region Using Statements:



    #endregion



    /// <summary>
    /// Indexes the Matrix Row[i] Where i is the Row Index starting from 0.
    /// </summary>
    public class Row
    {



        #region Fields:


        Matrix mat;


        #endregion



        #region Properties:



        /// <summary>
        /// Indexer of the Rows in the Matrix
        /// </summary>
        /// <param name="row"></param>
        /// <returns></returns>
        public double[] this[int row]
        {
            get
            {
                double[] MatrixRow = new double[mat.Columns];

                for (int j = 0; j < mat.Columns; j++)
                    MatrixRow[j] = mat[row, j];

                return MatrixRow;
            }
            set
            {
                for (int j = 0; j < mat.Columns; j++)
                    mat[row, j] = value[j];
            }
        }



        #endregion



        /// <summary>
        /// Row Constructor.
        /// </summary>
        /// <param name="matrix"></param>
        public Row(Matrix matrix)
        {
            mat = matrix;
        }



    }

}
