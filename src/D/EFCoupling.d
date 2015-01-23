module	EFCoupling;

import	MathFuncs;

public import	Coupling;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
class	CEFCoupling: CCoupling
{
	private
	{
		double	c1;
		double	c2;
		double	ck;
		double	t0;
	}

	this()
	{
		this.c1 = 2.57e7;
		this.c2 = 2.85e6;
		this.ck = 2.5e8;

		this.t0 = 50.0e3;
		this.T0 = 240e3;
		this.T_prev = this.T0;

		this.c = this.ck;
		this.beta = 100;
	}

	double get_force(double ds, double dv)
	{
		double beta = 0;
		double eps_ds = 1e-3;

		T_cur = abs(T_prev) + c*(abs(ds) - abs(ds_prev)) + beta*dv;

		/*if (ds*ds_prev < 0)
			T_prev = T0;*/		 

		if (ds*dv >= 0)
		{
			c = c1;
			beta = 0;
		}
		else
		{
			if (abs(T_cur) >= F2(abs(ds)))
			{
				c = ck;
				beta = 0;
			}
			else
			{
				c = c2;
				beta = 0;
			}

			if (ds*ds_prev < 0)
				T_prev = T0;
		}

		T_cur = (abs(T_prev) + c*(abs(ds) - abs(ds_prev)))*sign(ds) + beta*dv;

		if (calls_count == 0)
		{
			T_prev = T_cur;
			ds_prev = ds;
		}

		calls_count++;

		return T_cur;
	}

	double F2(double x)
	{
		return t0 + c2*x;
	}
}