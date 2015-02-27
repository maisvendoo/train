//-------------------------------------------------------------------
//
//		Elastic-friction coupling model
//		(c) maisvendoo, 2015/01/19
//
//-------------------------------------------------------------------
module	EFCoupling;

import	MathFuncs;

public 
{
	import	Coupling;
}

//-------------------------------------------------------------------
//		Elastic-friction coupling class
//-------------------------------------------------------------------
class	CEFCoupling: CCoupling
{
	private
	{
		double	c1;
		double	c2;
		double	ck;
		double	t0;

		double	s1;
		double	s2;
		double	s_max;
		double	T_max;
	}

	this()
	{
		this.c1 = 2.57e7;
		this.c2 = 2.85e6;
		this.ck = 2.5e8;

		this.t0 = 50.0e3;
		this.T0 = 240e3;
		this.T_prev = 0;//this.T0;
		this.T_cur = 0;

		this.c = this.ck;
		this.beta = 0;
	}



	//---------------------------------------------------------------
	//		Coupling force calculation
	//---------------------------------------------------------------
	double get_force(double ds, double dv)
	{
		double beta = 0;

		c1 = (T_max - T0)/(s_max - T0/ck);
		c2 = (T0 - t0)*ck/(ck*s_max + T0 - T_max - t0);

		s1 = T0/ck;
		s2 = t0/ck;

		if (ds*dv >= 0)
		{
			if (abs(T_prev) < T0)
			{
				c = ck;
				beta = this.beta;
			}
			else
			{
				c = c1;
				beta = 0;
			}
		}
		else
		{
			if (abs(T_prev) > T0)
			{
				c = ck;
				beta = this.beta;
			}
			else
			{
				if (abs(T_prev) < ck*abs(ds))
				{
					c = c2;
					beta = 0;
				}
				else
				{
					c = ck;
					beta = this.beta;
				}
			}
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




	//---------------------------------------------------------------
	//		Release force calculation
	//---------------------------------------------------------------
	double F2(double x)
	{
		return t0 + c2*x;
	}




	//---------------------------------------------------------------
	//		Set coupling stiffnesses
	//---------------------------------------------------------------
	void set_stiffs(double c1, double c2, double ck)
	{
		this.c1 = c1;
		this.c2 = c2;
		this.ck = ck;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void set_damp_coeff(double beta)
	{
		this.beta = beta;
	}

	void set_s_max(double s_max)
	{
		this.s_max = s_max;
	}

	void set_T_max(double T_max)
	{
		this.T_max = T_max;
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void set_release_init_force(double t0)
	{
		this.t0 = t0;
	}

	void set_prev_force()
	{
		this.T_prev = 0;
	}
}