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

		double	dT;
		double	ds_f;
		double	cm;
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

		this.dT = 0;
		this.ds_f = 0;
		this.cm = 0;
	}



	//---------------------------------------------------------------
	//		Coupling force calculation
	//---------------------------------------------------------------
	double get_force(double ds, double dv)
	{
		double beta = 0;

		//T_cur = (abs(T_prev) + c*(abs(ds) - abs(ds_prev)))*sign(ds) + beta*dv;

		if (ds*dv >= 0)
		{
			if (abs(T_cur) < T0)
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
			if (abs(T_cur) > T0)
			{
				c = ck;
				beta = this.beta;
			}
			else
			{
				c = c2;
				beta = 0;
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