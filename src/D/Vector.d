//-----------------------------------------------------------------------
//
//		Work with vectors
//
//-----------------------------------------------------------------------
module	Vector;

import	std.math;

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------
double	Norm(double[] y)
{
	double	sum = 0;
	int		n = cast(int) y.length;

	for (int i = 0; i < n; i++)
	{
		sum += y[i]*y[i];
	}

	return sqrt(sum);
}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------
double[] vAdd(double[] a, double[] b)
{
	if (a.length != b.length)
		return null;

	int n = cast(int) a.length;

	double[] c = new double[n];

	for (int i = 0; i < n; i++)
		c[i] = a[i] + b[i];

	return c;
}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------
double[] vSub(double[] a, double[] b)
{
	if (a.length != b.length)
		return null;
	
	int n = cast(int) a.length;
	
	double[] c = new double[n];
	
	for (int i = 0; i < n; i++)
		c[i] = a[i] - b[i];
	
	return c;
}


//------------------------------------------------------------------------
//
//------------------------------------------------------------------------
double[] vLinearCombine(double[] a, double[] b, double k1, double k2)
{
	if (a.length != b.length)
		return null;

	int n = cast(int) a.length;
	double[] c = new double[n];

	for (int i = 0; i < n; i++)
		c[i] = k1*a[i] + k2*b[i];

	return c;
}

//-----------------------------------------------------------------------
//
//-----------------------------------------------------------------------
double[] vdMult(double[] a, double lambda)
{
	int n = cast(int) a.length;
	
	double[] c = new double[n];
	
	for (int i = 0; i < n; i++)
		c[i] = lambda*a[i];
	
	return c;
}
