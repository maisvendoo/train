//-------------------------------------------------------------------
//
//		Linear equations solving library
//		(c) maisvendoo, 2014/12
//
//-------------------------------------------------------------------
module	LinearEQs;

import	std.math;

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