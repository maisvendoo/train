module	Model;

import	ODEqs;

class CModel
{
	protected	
	{
		// State variables and it's derivatives
		double[]		y;
		double[]		dydt;	

		// Time and integration parameters
		double			t0;
		double			t_end;
		double			t;
		double			dt;
		double			dt_max;
		double			eps;

		// ODE solver method
		ode_solver_t	ode_solver_step;
		// ODE system function 
		ode_system_t	ode_sys;
	}

	this()
	{
		this.ode_solver_step = &rkf5_solver_step;
		this.ode_sys = &this.ode_system;
		this.y = [0];
		this.dydt = [0];

		this.t0 = 0;
		this.t = 0;
		this.t_end = 1.0;
		this.dt = 1e-4;
		this.dt_max = 0.1;
		this.eps = 1e-8;
	}

	//---- Protected methods

	protected void ode_system(double[] Y, ref double[] dYdt, double t)
	{
		dYdt[0] = Y[0];
	}

	//---- Private methods

	protected void step()
	{
		ode_solver_step(y, dydt, t, dt, dt_max, eps, ode_sys);
	}

	//---- Public methods

	// Simulation progress
	void process()
	{
		while (t <= t_end)
		{
			step();
			t += dt;
		}
	}

	// Set state variable
	void set_y(uint idx, double y)
	{
		this.y[idx] = y;
	}

	// Set initial time value
	void set_init_time(double t0)
	{
		this.t0 = t0;
	}

	// Set simulation stop time
	void set_end_time(double t_end)
	{
		this.t_end = t_end;
	}

	// Set time step (initial value for variable step methods)
	void set_time_step(double dt)
	{
		this.dt = dt;
	}

	// Set maximal time step (for variable step methods)
	void set_max_time_step(double dt_max)
	{
		this.dt_max = dt_max;
	}

	// Set local error of method (for variable step methods)
	void set_local_error(double eps)
	{
		this.eps = eps;
	}

	// Set ODE solve method
	void set_integration_method(ode_solver_t ode_solver)
	{
		this.ode_solver_step = ode_solver;
	}

	// Set system dimension
	void set_dimension(ulong dim)
	{
		y.length = dim;
		dydt.length = dim;
	}

	// Get state variable value
	double get_y(uint idx)
	{
		return this.y[idx];
	}

	// Get current time value
	double get_time()
	{
		return this.t;
	}

	// Get current time step
	double get_time_step()
	{
		return this.dt;
	}
}