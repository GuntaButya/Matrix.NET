# Matrix.NET
A Light Weight Matrix based on the Code of darkdragon-001 
Orriginal code can be viewed here: https://github.com/darkdragon-001/LightweightMatrixCSharp

    // Test on Matrix: 10 Columns in each Row
    Matrix M = new Matrix(3, 10);
    M.Row[0] = new double[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double[] ColZero = M.Column[0];
            
    M.Row[1] = new double[] { 10, 11, 12, 13, 14, 15, 16, 17, 18, 19 };
    M.Row[2] = new double[] { 20, 21, 22, 23, 24, 25, 26, 27, 28, 29 };
            
    double[] RowZero = M.Row[0];
    double[] RowOne = M.Row[1];
    double[] RowTwo = M.Row[2];

    double[] ColNine = M.Column[0];
