//-------------------------------------------------------------------
//
//		Matrix mathematics module
//		(c) maisvendoo, 2014.12 
//
//-------------------------------------------------------------------
module		matrix;

import		std.stdio;

//-------------------------------------------------------------------
//		Matrix (n x m) creation
//-------------------------------------------------------------------
double[][] create_matrix(uint n, uint m)
{
	double[][]	new_matrix = new double[][](n, m); 

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			new_matrix[i][j] = 0;

	return new_matrix;
}

//-------------------------------------------------------------------
//		Matrix multiplication
//-------------------------------------------------------------------
double[][] matrix_mult(double[][] M1, double[][] M2)
{
	ulong n1 = M1.length;
	ulong m1 = M1[0].length;

	ulong n2 = M2.length;
	ulong m2 = M2[0].length;

	if (m1 != n2)
		return null;

	double[][] M3 = new double[][](n1, m2);

	for (ulong j = 0; j < m2; j++)
	{
		for (ulong i = 0; i < n1; i++)
		{
			double c = 0;

			for (ulong k = 0; k < m1; k++)
				c = c + M1[i][k]*M2[k][j];

			M3[i][j] = c;
		}	
	}

	return M3;
}

//-------------------------------------------------------------------
//		Matrix transpose
//-------------------------------------------------------------------
double[][] transpose_matrix(double[][] M)
{
	ulong		n = M.length;
	ulong		m = M[0].length;

	double[][]	Mt = new double[][](m, n);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			Mt[j][i] = M[i][j];

	return Mt;
}

//-------------------------------------------------------------------
//		Matrix addition
//-------------------------------------------------------------------
double[][] sum_matrix(double[][] M1, double[][] M2)
{
	ulong	n1 = M1.length;
	ulong	m1 = M1[0].length;

	ulong	n2 = M2.length;
	ulong	m2 = M2[0].length;

	if ( (n1 != n2) || (m1 != m2) )
		return null;

	double[][] Ms = new double[][](n1, m1);

	for (int i = 0; i < n1; i++)
		for (int j = 0; j < m1; j++)
			Ms[i][j] = M1[i][j] + M2[i][j];

	return Ms;
}

//-------------------------------------------------------------------
//		Number-matrix multiplication
//-------------------------------------------------------------------
double[][] matrix_num_mult(double[][] M, double value)
{
	ulong	n = M.length;
	ulong 	m = M[0].length;

	double[][] ret_M = new double[][](n, m);

	for (int i = 0; i < n; i++)
		for (int j = 0; j < m; j++)
			ret_M[i][j] = value*M[i][j];

	return ret_M;
}

//-------------------------------------------------------------------
//		Print matrix for debug
//-------------------------------------------------------------------
void print_matrix(double[][] M, File file = stdout)
{
	ulong	n = M.length;
	ulong	m = M[0].length;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
			file.writef("%10.3f ", M[i][j]);

		file.writeln();
	}
}
