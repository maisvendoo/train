//-------------------------------------------------------------------
//
//		Linear equations solving library
//		(c) maisvendoo, 2014/12
//
//-------------------------------------------------------------------
module	LinearEQs;

import	std.math;
import	std.stdio;

import	matrix;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
bool seidel_solver(double[][] A, 
				   double[][] B,
				   ref double[] x, 
				   double eps)
{
	// Check solve error
	bool converge(double[] x1, double[] x)
	{
		ulong	n = x.length;
		double	norm = 0;

		for (int i = 0; i < n; i++)
		{
			norm += (x1[i] - x[i])*(x1[i] - x[i]);
		}

		if (sqrt(norm) >= eps)
			return false;

		return true;
	}	
	
	ulong		n = A.length;		// System dimension
	double[]	p = new double[n];	// Previos iteration solve
	uint		iter_count = 0;		// Iterations count
	uint		max_iters = 1000;	// Max iterations

	// Seider method siquence
	do
	{
		for (int i = 0; i < n; i++)
			p[i] = x[i];
			
		for (int i = 0; i < n; i++)
		{
			double sum = 0;

			for (int j = 0; j < i; j++)
				sum += A[i][j]*x[j];

			for (int j = i+1; j < n; j++)
				sum += A[i][j]*p[j];

			x[i] = (B[i][0] - sum) / A[i][i];
		}

		iter_count++;

	} while ( (!converge(x, p)) && (iter_count <= max_iters) );

	if (iter_count > max_iters)
		return false;
	else
		return true;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void gauss_solver(double[][] A, double[][] b, ref double[] x)
{
	double	c = 0;
	int		n = cast(int) A.length;

	// Forward (A to uptriangle form)
	for (int i = 0; i < n-1; i++)
	{
		for (int j = i; j < n-1; j++)
		{
			c = A[j+1][i] / A[i][i];

			for (int k = 0; k < n; k++)
				A[j+1][k] = A[j+1][k] - c*A[i][k];

			b[j+1][0] = b[j+1][0] - c*b[i][0];
		}
	}

	// Backward (roots found)
	for (int i = n-1; i >=0; i--)
	{
		double sum = b[i][0];

		for (int j = n-1; j >= i+1; j--)
			sum -= A[i][j]*x[j];

		x[i] = sum / A[i][i];
	}
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void gaussLER_solver(double[][] A, double[][] b, ref double[] x)
{
	double	c = 0;
	int		n = cast(int) A.length;
	double	max = 0;
	int 	p = 0;

	for (int i = 0; i < n-1; i++)
	{
		for (int k = i; k < n; k++)
		{
			if (abs(A[k][i]) > max)
			{
				max = A[k][i];
				p = k;
			}
		}

		print_matrix(A);
		writeln();
		
		rows_change(A, i, p);
		rows_change(b, i, p);

		print_matrix(A);
		writeln();

		for (int j = i; j < n-1; j++)
		{
			c = A[j+1][i] / A[i][i];
			
			for (int k = 0; k < n; k++)
				A[j+1][k] = A[j+1][k] - c*A[i][k];
			
			b[j+1][0] = b[j+1][0] - c*b[i][0]; 
		}
	}

	// Backward (roots found)
	for (int i = n-1; i >=0; i--)
	{
		double sum = b[i][0];
		
		for (int j = n-1; j >= i+1; j--)
			sum -= A[i][j]*x[j];
		
		x[i] = sum / A[i][i];
	}
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void gaussLEF_solver(const double[][] A, const double[][] b, ref double[] x)
{
	int			n = cast(int) A.length;
	int			p = 0;
	int 		q = 0;
	double		c = 0;
	double		m = 0;
	int			count = 0;
	int[]		unkw_num = new int[n];
	double		eps_z = 1e-10;

	// Choose leading element
	double get_LE(double[][] M, ref int p, ref int q)
	{
		double max = 0;
		int n = cast(int) M.length;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				if (abs(M[i][j]) > max)
				{
					max = abs(M[i][j]);
					p = i;
					q = j;
				}
			}
		}

		return M[p][q];
	}

	// Tmp matrixes
	double[][] A1 = new double[][](n, n+1); 
	double[][] A2 = new double[][](n, n+1); 
	double[][] A3 = new double[][](n, n+1);

	// Prepare tmp matrixes
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A1[i][j] = A[i][j];
			A2[i][j] = A[i][j];
		}

		A1[i][n] = b[i][0];
		A2[i][n] = b[i][0];
	}

	// Elimination of unknowns with pivoting
	do
	{
		c = get_LE(A1, p, q);

		for (int i = 0; i < n; i++)
		{
			if (i != p)
			{
				m = -(A1[i][q] / A1[p][q]);

				rows_sum(A1, 1, m, i, p);
				rows_sum(A2, 1, m, i, p);
			}		
		}	

		row_zero(A1, p);
		A3[count] = A2[p];

		count++;

	} while (count <= n-1);

	// Transorm matrix to triangle form
	for (int i = 0; i < n; i++)	
		unkw_num[i] = i;

	for (int i = n-1; i >= 0; i--)
	{
		if (abs(A3[i][i]) > eps_z)
			continue;

		for (int k = 0; k < i; k++)
		{
			if (abs(A3[i][k]) > eps_z)
			{
				p = k;
			}
		}

		columns_change(A3, i, p);
		change(unkw_num, i, p);
	}

	// Backward (roots found)
	double[] y = new double[n];

	for (int i = n-1; i >=0; i--)
	{
		double sum = A3[i][n];
		
		for (int j = n-1; j >= i+1; j--)
			sum -= A3[i][j]*y[j];
		
		y[i] = sum / A3[i][i];
		x[unkw_num[i]] = y[i];
	}
}