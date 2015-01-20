//-------------------------------------------------------------------
//
//		Base coupling model
//		(c) maisvendoo, 2015/01/18
//
//-------------------------------------------------------------------
module	Coupling;

import	MathFuncs;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
class	CCoupling
{
	protected
	{
		double		c;				// Coupling stiffness
		double		beta;			// Damping coefficient
		double		T0;				// Initial force

		double		ds_prev;		// Previos deformation
		double		T_prev;			// Previois force
		double		T_cur;			// Current force

		int			calls_count;	// Count of get_force method's calls
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	this()
	{
		this.c = 1e7;
		this.beta = 0;
		this.T0 = 240e3;

		this.ds_prev = 0;
		this.T_cur = 0;
		this.T_prev = T0;

		this.calls_count = 0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void set_stiffness(double c)
	{
		this.c = c;
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
	void set_init_force(double T0)
	{
		this.T0 = T0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double get_force(double ds, double dv)
	{
		return (T0 + c*abs(ds))*sign(ds) + beta*dv;
	}

	double get_init_force()
	{
		return T0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void reset()
	{
		calls_count = 0;
	}
}