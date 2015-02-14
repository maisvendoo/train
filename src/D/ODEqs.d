//-------------------------------------------------------------------
//
//		Library for solving of Ordinar Differential Equations (ODE)
//		(c) maisvendoo, 2014/12
//
//-------------------------------------------------------------------
module	ODEqs;

import	std.math;
import	std.stdio;
import	std.c.stdlib;

private
{
	enum			int		MAX_ITER	= 10;

	// For Adams method
	int						rk4_count	= 0;
	enum			int		ORDER		= 4;
	double[][]				f;
	bool					init		= true;
}

//-------------------------------------------------------------------
//		Alias for ODE's system function 
//-------------------------------------------------------------------
alias	ode_system_t = void delegate(double[], ref double[], double);

alias	ode_solver_t = double[] function(ref double[], 
									 ref double[], 
									 double, 
									 ref double,
									 double,
									 double,
									 ode_system_t);


//-------------------------------------------------------------------
//		Euler method (1'st order) step
//-------------------------------------------------------------------
double[] euler_solver_step(ref double[] Y, 
					   ref double[] dYdt, 
					   double t, 
					   ref double dt,
					   double dt_max,
					   double eps,
					   ode_system_t ode_sys)
{
	// Derivatives (right part) calculation 
	ode_sys(Y, dYdt, t);

	// New value of the state vector
	for (int i = 0; i < Y.length; i++)
		Y[i] += dYdt[i]*dt;

	return dYdt;
}

//-------------------------------------------------------------------
//		Runge-Kutta method (4'th order) step
//-------------------------------------------------------------------
double[] rk4_solver_step(ref double[] Y,
					 ref double[] dYdt,
					 double	t,
					 ref double dt,
					 double dt_max,
					 double eps,
					 ode_system_t ode_sys)
{
	ulong n = Y.length;

	double[] k1 = new double[n];
	double[] k2 = new double[n];
	double[] k3 = new double[n];
	double[] k4 = new double[n];
	double[] Y1 = new double[n];

	ode_sys(Y, dYdt, t);
	
	for (int i = 0; i < n; i++)
	{
		k1[i] = dYdt[i];
		Y1[i] = Y[i] + k1[i]*dt/2;
	}

	ode_sys(Y1, dYdt, t + dt/2);

	for (int i = 0; i < n; i++)
	{
		k2[i] = dYdt[i];
		Y1[i] = Y[i] + k2[i]*dt/2;
	}

	ode_sys(Y1, dYdt, t + dt/2);

	for (int i = 0; i < n; i++)
	{
		k3[i] = dYdt[i];
		Y1[i] = Y[i] + k3[i]*dt;
	}

	ode_sys(Y1, dYdt, t + dt);

	for (int i = 0; i < n; i++)
	{
		k4[i] = dYdt[i];
		Y[i] = Y[i] + dt*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;
	}

	return k1;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double[] rkf5_solver_step(ref double[] Y,
					  ref double[] dYdt,
					  double t,
					  ref double dt,
					  double dt_max,
					  double eps,
					  ode_system_t ode_sys)
{
	const	double	c1 = 16.0/135;
	const	double	c3 = 6656.0/12825;
	const	double	c4 = 28561.0/56430;
	const	double	c5 = -9.0/50;
	const	double	c6 = 2.0/55;

	const	double	b21 = 0.25;
	
	const	double	b31 = 3.0/32;
	const	double	b32 = 9.0/32;

	const	double	b41 = 1932.0/2197;
	const	double	b42 = -7200.0/2197;
	const	double	b43 = 7296.0/2197;

	const	double	b51 = 439.0/216;
	const	double	b52 = -8.0;
	const	double	b53 = 3680.0/513;
	const	double	b54 = -845.0/4140;

	const	double	b61 = -8.0/27;
	const	double	b62 = 2.0;
	const	double	b63 = -3544.0/2565;
	const	double	b64 = 1859.0/4104;
	const	double	b65 = -11.0/40;

	const	double	e1 = 1.0/360;
	const	double	e3 = -128.0/4275;
	const	double	e4 = -2197.0/75240;
	const	double	e5 = 1.0/50;
	const	double	e6 = 2.0/55;

	ulong	n = Y.length;

	double[]	k1 = new double[n];
	double[]	k2 = new double[n];
	double[]	k3 = new double[n];
	double[]	k4 = new double[n];
	double[]	k5 = new double[n];
	double[]	k6 = new double[n];

	bool	ready = false;
	int		iter = 0;
	double	delta = 0;

	double[]	Y1 = new double[n];
	double[]	dY1dt = new double[n];
	double[]	eps_y = new double[n];

	do
	{
		delta = 0;
		iter++;

		ode_sys(Y, dYdt, t);

		for (int i = 0; i < n; i++)
		{
			dY1dt = dYdt;
			k1[i] = dt*dYdt[i];
			Y1[i] = Y[i] + b21*k1[i];
		}

		ode_sys(Y1, dYdt, t + dt/4.0);

		for (int i = 0; i < n; i++)
		{
			k2[i] = dt*dYdt[i];
			Y1[i] = Y[i] + b31*k1[i] + b32*k2[i];
		}

		ode_sys(Y1, dYdt, t + 3.0*dt/8.0);

		for (int i = 0; i < n; i++)
		{
			k3[i] = dt*dYdt[i];
			Y1[i] = Y[i] + b41*k1[i] + b42*k2[i] + b43*k3[i];
		}

		ode_sys(Y1, dYdt, t + 12.0*dt/13.0);

		for (int i = 0; i < n; i++)
		{
			k4[i] = dt*dYdt[i];
			Y1[i] = Y[i] + b51*k1[i] + b52*k2[i] + b53*k3[i] + 
					b54*k4[i];
		}

		ode_sys(Y1, dYdt, t + dt);

		for (int i = 0; i < n; i++)
		{
			k5[i] = dt*dYdt[i];
			Y1[i] = Y[i] + b61*k1[i] + b62*k2[i] + b63*k3[i] + 
					b64*k4[i] + b65*k5[i];
		}

		ode_sys(Y1, dYdt, t + dt/2.0);

		for (int i = 0; i < n; i++)
		{
			k6[i] = dt*dYdt[i];
			eps_y[i] = abs(e1*k1[i] + e3*k3[i] + e4*k4[i] +
						   e5*k5[i] + e6*k6[i]);
						   
			if (delta < eps_y[i])
				delta = eps_y[i];						 	
		}

		if (delta >= eps)
		{
			dt = dt/2;
			ready = false;
		}

		if (delta <= eps/32.0)
		{
			dt = 2*dt;

			if (dt > dt_max)
				dt = dt_max;

			ready = true;
		}

		if ( (delta > eps/32.0) && (delta < eps) )
		{
			ready = true;
		}

		
		if (ready)
		{
			for (int i = 0; i < n; i++)
				Y[i] = Y[i] + c1*k1[i] + c3*k3[i] + c4*k4[i] +
							  c5*k5[i] + c6*k6[i];
		}

	} while ( (!ready) && (iter <= MAX_ITER) );

	if (iter > MAX_ITER)
	{
		writeln("FAIL: Iterations limit");
		exit(0);
	}

	return dY1dt;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//		Adams method (4'th order) step
//-------------------------------------------------------------------
double[] adams_solver_step(ref double[] Y,
					   ref double[] dYdt,
					   double	t,
					   ref double dt,
					   double dt_max,
					   double eps,
					   ode_system_t ode_sys)
{
	int n = cast(int) Y.length;

	if (init)
	{
		f = new double[][](ORDER, Y.length);
		init = false;
	}

	if (rk4_count < ORDER)
	{
		f[rk4_count] = rk4_solver_step(Y, dYdt, t, dt, dt_max, eps, ode_sys);

		rk4_count++;
	}
	else
	{
		ode_sys(Y, dYdt, t);
		//f[rk4_count] = dYdt;

		for (int i = 0; i < n; i++)
		{
			Y[i] = Y[i] + (-9*f[rk4_count-3][i] + 55*dYdt[i] + 
				           37*f[rk4_count-2][i] - 59*f[rk4_count-1][i])*dt/24;

			f[ORDER-1] = dYdt;
		}

		for (int i = 0; i < ORDER - 1; i++)
			f[i] = f[i+1];
	}

	return null;
}