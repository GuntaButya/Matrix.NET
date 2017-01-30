namespace Matrix.NET
{



    #region Using Statements:



    #endregion



    /// <summary>
    /// Indexes the Matrix Column[i] Where i is the Column Index starting from 0.
    /// </summary>
    public class Column
    {



        #region Fields:


        Matrix mat;


        #endregion



        #region Properties:



        /// <summary>
        /// Indexer of the Columns in the Matrix
        /// </summary>
        /// <param name="row"></param>
        /// <returns></returns>
        public double[] this[int column]
        {
            get
            {
                double[] MatrixColumn = new double[mat.Rows];

                for (int j = 0; j < mat.Rows; j++)
                    MatrixColumn[j] = mat[j, column];

                return MatrixColumn;
            }
            set
            {
                for (int j = 0; j < mat.Rows; j++)
                    mat[j, column] = value[j];
            }
        }



        #endregion



        /// <summary>
        /// Row Constructor.
        /// </summary>
        /// <param name="matrix"></param>
        public Column(Matrix matrix)
        {
            mat = matrix;
        }



    }

}
