//-------------------------------------------------------------------
//
//		Matrix mathematics module
//		(c) maisvendoo, 2014/12 
//
//-------------------------------------------------------------------
module		matrix;

import		std.stdio;
import		std.math;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void change(T)(ref T[] array, int p, int q)
{
	T	c;
	
	c = array[p];
	array[p] = array[q];
	array[q] = c;
}

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
			file.writef("%20.5f ", M[i][j]);

		file.writeln();
	}
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
int rows_change(ref double[][] M, int i, int j)
{
	int	n = cast(int) M.length;

	// Indexes are upper than matrix size 
	if ( (i > n-1) || (j > n-1) )
		return -1;

	// Indexes are equate - exit without actions
	if (i == j)
		return 0;

	double[] tmp = M[i];

	M[i] = M[j];
	M[j] = tmp;

	return 0;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
int rows_sum(ref double[][] M, double c_i, double c_j, int i, int j)
{
	int	n = cast(int) M.length;
	int m = cast(int) M[0].length;

	if ( (i > n-1) || (j > n-1) )
		return -1;

	if (i == j)
		return -1;

	for (int k = 0; k < m; k++)
		M[i][k] = c_i*M[i][k] + c_j*M[j][k];

	return 0;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double[][] get_minor_matrix(double[][] M, int p, int q)
{
	int	n = cast(int) M.length;
	int	m = cast(int) M[0].length;

	if ( (p > n-1) || (q > m-1) )
		return null;

	if ( (p <0) || (q < 0) )
		return null;

	double[][] tmp = create_matrix(n-1, m-1);

	for (int i = 0; i < p; i++)
		for (int j = 0; j < q; j++)
			tmp[i][j] = M[i][j];

	for (int i = 0; i < p; i++)
		for (int j = q; j < m-1; j++)
			tmp[i][j] = M[i][j+1];

	for (int i = p; i < n-1; i++)
		for (int j = 0; j < q; j++)
			tmp[i][j] = M[i+1][j];

	for (int i = p; i < n-1; i++)
		for (int j = q; j < m-1; j++)
			tmp[i][j] = M[i+1][j+1];

	return tmp;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double get_absmax_element(double[][] M, ref int p, ref int q)
{
	double	max = 0;
	int		n = cast(int) M.length;
	int		m = cast(int) M[0].length;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (abs(M[i][j]) > max)
			{
				max = M[i][j];
				p = i;
				q = j;
			}
		}
	}

	return M[p][q];
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void row_zero(double[][] M, int p)
{
	int	m = cast(int) M[0].length;

	for (int j = 0; j < m; j++)
		M[p][j] = 0;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void columns_change(double[][] M, int p, int q)
{
	int	n = cast(int) M.length;

	for (int i = 0; i < n; i++)
	{
		double tmp = M[i][p];
		M[i][p] = M[i][q];
		M[i][q] = tmp;
	}
}