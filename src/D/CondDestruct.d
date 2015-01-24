//-------------------------------------------------------------------
//
//		Conditional destruction calculation module
//		(c) maisvendoo, 2015/01/14
//
//-------------------------------------------------------------------
module	CondDestruct;

import	MathFuncs;
import	std.stdio;

private
{
	enum	double		MIN_LONG_FORCE		= 0;
	enum	double		MAX_LONG_FORCE		= 2e6;
	enum	int			INTERVAL_COUNT		= 100;
	enum	double		BASE_FORCE			= MAX_LONG_FORCE;
	enum	int			m					= 6;

	int[INTERVAL_COUNT]	force_count;
	double				dF					= 0;

	int get_interval(double force)
	{
		return cast(int) (abs(force)/dF);
	}

	void cond_destruct_init()
	{
		for (int i = 0; i < INTERVAL_COUNT; i++)
			force_count[i] = 0;
		
		dF = (MAX_LONG_FORCE - MIN_LONG_FORCE) / INTERVAL_COUNT;
	}
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void force_analyze(double force)
{
	int idx = get_interval(force);

	force_count[idx]++;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double get_cond_destruct(double[][] force_array)
{
	double	J = 0;

	cond_destruct_init();

	int		nv = cast(int) force_array.length;
	int		np = cast(int) force_array[0].length;

	for (int i = 0; i < nv; i++)
	{
		for (int j = 0; j < np; j++)
		{
			force_analyze(force_array[i][j]);
		}
	}

	for (int i = 0; i < INTERVAL_COUNT; i++)
	{
		double force = (2*i+1)*dF/2.0;

		J += pow(abs(force)/BASE_FORCE, m)*force_count[i];
	}

	return J;
}