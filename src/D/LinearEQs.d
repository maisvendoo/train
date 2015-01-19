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