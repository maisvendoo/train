//-------------------------------------------------------------------
//
//		Train motion simulation module
//		(c) maisvendoo, 2015/01/13
//
//-------------------------------------------------------------------
module	Train;

import	std.stdio;
import	MathFuncs;
import	LogFile;
import	Model;
import	LuaScript;
import	matrix;
import	EFCoupling;
import	LinearEQs;
import	ConstReacts;

//-------------------------------------------------------------------
//		General train  model class
//-------------------------------------------------------------------
class	CTrainModel: CModel
{
	//---------------------------------------------------------------
	//		Train model parameters
	//---------------------------------------------------------------
	private
	{
		uint		nv;
		uint		nb;
		uint		ode_dim;
		uint		mass_n;

		string		couplig_type;

		double		v0;
		double		railway_coord;

		double[]	m;

		CLogFile	terminal;
		CLogFile	file_log;
		CLuaScript	lua_cfg;

		double[]	F;
		double[]	R1;
		double[]	R2;
		double[]	P1;
		double[]	P2;
		double[]	B;
		double[]	Bmax;

		double[][]	Theta;
		double[][]	Q;

		CEFCoupling[]	fwd_coup;
		CEFCoupling[]	bwd_coup;
	}



	//------------------------------------------------------------------
	//		Constructor
	//------------------------------------------------------------------
	this()
	{
		this.nv = 1;
		this.mass_n = 3;
		this.nb = this.mass_n*this.nv;

		this.couplig_type = "default";

		this.v0 = 0;
		this.railway_coord = 1e6;

		terminal = new CLogFile();
		terminal.init();
		terminal.set_print_func(&this.term_out);

		file_log = new CLogFile();
		file_log.init("log.txt");
		file_log.set_print_func(&this.file_out);

		this.lua_cfg = new CLuaScript();

		set_dimension(2);
	}



	//------------------------------------------------------------------
	//		Motion ODE system
	//------------------------------------------------------------------
	override protected void ode_system(double[] Y,
									   ref double[] dYdt, 
									   double t) 
	{
		for (int i = 0; i < nv; i++)
		{
			int k = mass_n*i; 

			dYdt[k]		= Y[k+nb];
			dYdt[k+1]	= Y[k+1+nb];
			dYdt[k+2]	= Y[k+2+nb];

			double[] a = get_accels(Y, t, i);

			dYdt[k+nb]		= a[0];
			dYdt[k+1+nb]	= a[1];
			dYdt[k+2+nb]	= a[2];
		}
	}



	//------------------------------------------------------------------
	//		Simulation progress
	//------------------------------------------------------------------
	override protected void process() 
	{
		while (t < t_end)
		{
			terminal.print(0.01, dt);
			file_log.print(0.01, dt);
			step();

			for (int i = 0; i < nv; i++)
			{
				fwd_coup[i].reset();
				bwd_coup[i].reset();
			}

			t += dt;
		}
	}



	//------------------------------------------------------------------
	//		Terminal print by simulation process
	//------------------------------------------------------------------
	private void term_out(File term)
	{
		term.writefln("t = %f x0 = %f s10 = %f s20 = %f x1 = %f s11 = %f s21 = %f", t, y[0], y[1], y[2], y[3], y[4], y[5]);
	}

	private void file_out(File file)
	{
		file.writefln("%f %f %f", t, 
			          y[1] - y[4], 
			          P2[0]);
	}



	//------------------------------------------------------------------
	//		Initialization form Lua-script file
	//------------------------------------------------------------------
	int init(string cfg_name)
	{
		int err = 0;

		if (read_lua_config(cfg_name) == -1)
			return -1;

		return err;
	}

	//---------------------------------------------------------------
	//		Read parameters from Lua-config
	//---------------------------------------------------------------
	private int read_lua_config(string cfg_name)
	{
		int	err = LUA_S_OK;

		err = lua_cfg.exec_script(cfg_name);

		if (err == -1)
			return err;

		//--- Set ODE solver parameters
		string method = lua_cfg.get_str_field("solver_params", "method", err);

		// method
		if (err == LUA_S_NOEXIST)
			method = "rkf5";

		if (method == "rkf5")
			set_integration_method(&rkf5_solver_step);

		if (method == "rk4")
			set_integration_method(&rk4_solver_step);

		if (method == "euler")
			set_integration_method(&euler_solver_step);

		// init time
		double init_time = lua_cfg.get_double_field("solver_params", "init_time", err);

		if (err == LUA_S_NOEXIST)
			init_time = 0;

		set_init_time(init_time);

		// stop time
		double stop_time = lua_cfg.get_double_field("solver_params", "stop_time", err);

		if (err == LUA_S_NOEXIST)
			stop_time = 1.0;

		set_stop_time(stop_time);

		// time step
		double time_step = lua_cfg.get_double_field("solver_params", "step", err);

		if (err == LUA_S_NOEXIST)
			time_step = get_time_step();

		set_time_step(time_step);

		// maximal time step
		double max_step = lua_cfg.get_double_field("solver_params", "max_step", err);

		if (err == LUA_S_NOEXIST)
			max_step = 0.1;

		set_max_time_step(max_step);

		// maximal time step
		double local_err = lua_cfg.get_double_field("solver_params", "local_err", err);

		if (err == LUA_S_NOEXIST)
			local_err = 1e-8;

		set_local_error(local_err);

		//----	Set train model parameters

		// number of vehicles
		nv = lua_cfg.get_int_field("train_model", "vehicles_num", err);

		if (err == LUA_S_NOEXIST)
			nv = 1;

		nb = mass_n*nv;

		m = new double[nb];

		ode_dim = 2*nb;

		set_dimension(ode_dim);

		for (int i = 0; i < ode_dim; i++)
		{
			y[i] = dydt[i] = 0;
		}

		R1 = new double[nv];
		R2 = new double[nv];
		B = new double[nv];
		Bmax = new double[nv];
		F = new double[nv];
		P1 = new double[nv];
		P2 = new double[nv];

		// couplig type
		couplig_type = lua_cfg.get_str_field("train_model", "coupling_type", err);

		if (err == LUA_S_NOEXIST)
			couplig_type = "default";

		// railway coordinate
		railway_coord = lua_cfg.get_double_field("train_model", "railway_coord", err);

		if (err == LUA_S_NOEXIST)
			railway_coord = 1e6;

		// initial velocity
		v0 = lua_cfg.get_double_field("train_model", "init_velocity", err);

		if (err == LUA_S_NOEXIST)
			v0 = 0;

		err = mass_init();

		if (err == -1)
			return err;

		err = couplings_init();

		get_accels(y, 0, 0);

		return err;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int mass_init()
	{
		int err = LUA_S_OK;

		double mass_coeff = lua_cfg.get_double("mass_coeff", err);

		if (err == LUA_S_NOEXIST)
			return -1;

		for (int i = 0; i < nv; i++)
		{
			int k = mass_n*i;

			double mass = lua_cfg.get_double_field("vehicle_mass", i, err);

			if (err == LUA_S_NOEXIST)
				return -1;

			m[k] = mass_coeff*mass;
			m[k+2] = mass_coeff*mass;
			m[k+1] = mass - m[k] - m[k+2];
		}

		Theta = create_matrix(mass_n, mass_n);
		Q = create_matrix(mass_n, 1);

		return 0;
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int couplings_init()
	{
		fwd_coup = new CEFCoupling[nv];
		bwd_coup = new CEFCoupling[nv];

		for (int i = 0; i < nv; i++)
		{
			fwd_coup[i] = new CEFCoupling();
			bwd_coup[i] = new CEFCoupling();

			fwd_coup[i].reset();
			bwd_coup[i].reset();
		}

		return 0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double[] get_accels(double[] Y, double t, int idx)
	{
		double	eps_v = 1e-10;
		double	eps_s = 1e-10;

		int k = mass_n*idx;

		double[] x = new double[2*mass_n];
		double[] a = new double[mass_n];

		double[][] A = create_matrix(2*mass_n, 2*mass_n);
		double[][] b = create_matrix(2*mass_n, 1);
	
		A[0][0] = m[k] + m[k+1] + m[k+2];
		A[0][1] = m[k];
		A[0][2] = -m[k+2];
		A[0][3] = 1;

		A[1][0] = m[k];
		A[1][1] = m[k];
		A[1][2] = 0;
		A[1][4] = 1;

		A[2][0] = -m[k+2];
		A[2][1] = 0;
		A[2][2] = m[k+2];
		A[2][5] = 1;

		for (int i = 0; i < mass_n; i++)
			for (int j = 0; j < mass_n; j++)
				Theta[i][j] = A[i][j];

		if (idx == 0)
			R1[idx] = 0;
		else
			R1[idx] = R2[idx-1]; 

		if (idx == nv-1)
			R2[idx] = 0;
		else
			R2[idx] = get_gap_force(y[k] - y[k+4] - (y[k+2] + y[k+3]), 
									y[k+nb] - y[k+4+nb] - (y[k+2+nb] + y[k+3+nb]),
									-0.025,
									0.025);


		F[idx] = 0;
		int err = 0;
		F[0] = lua_cfg.call_func("Traction", [t], err);
		Bmax[idx] = 1000;

		b[0][0] = F[idx] + R1[idx] - R2[idx];
		b[1][0] = R1[idx];
		b[2][0] = R2[idx];

		if ( (abs(y[k+nb]) < eps_v) && (abs(b[0][0]) <= Bmax[idx]) )
		{
			A[3][0] = 1;
			A[3][3] = 0;

			b[3][0] = 0;
		}
		else
		{
			A[3][0] = 0;
			A[3][3] = 1;

			b[3][0] = Bmax[idx]*sign(y[k+nb]);
		}

		//
		if ( (abs(y[k+1]) < eps_s) && (abs(y[k+1+nb]) < eps_v) )
		{
			A[4][1] = 1;
			A[4][4] = 0;

			b[4][0] = 0;
		}
		else
		{
			A[4][1] = 0;
			A[4][4] = 1;

			b[4][0] = fwd_coup[idx].get_force(y[k+1], y[k+1+nb]) + 
				      get_gap_force(y[k+1], y[k+1+nb], -0.09, 0.09);
		}

		//
		if ( (abs(y[k+2]) < eps_s) && (abs(y[k+2+nb]) < eps_v) )
		{
			A[5][2] = 1;
			A[5][5] = 0;
			
			b[5][0] = 0;
		}
		else
		{
			A[5][2] = 0;
			A[5][5] = 1;
			
			b[5][0] = bwd_coup[idx].get_force(y[k+2], y[k+2+nb]) + 
				      get_gap_force(y[k+2], y[k+2+nb], -0.09, 0.09);
		}

		gaussLEF_solver(A, b, x);

		//print_matrix(A);

		//writeln(x);

		if (abs(y[k]) < eps_v)
		{
			if (abs(x[3]) <= Bmax[idx])
			{
				B[idx] = x[3];
			}
			else
			{
				B[idx] = Bmax[idx]*sign(x[3]);
			}
		}
		else
			B[idx] = Bmax[idx]*sign(y[k+nb]);

		if ( (abs(y[k+1]) < eps_s) && (abs(y[k+1+nb]) < eps_v) )
		{
			double T0 = fwd_coup[idx].get_init_force();

			if (abs(x[4]) <= T0)
				P1[idx] = x[4];
			else
				P1[idx] = T0;
		}
		else
			P1[idx] = fwd_coup[idx].get_force(y[k+1], y[k+1+nb]) + 
			          get_gap_force(y[k+1], y[k+1+nb], -0.09, 0.09);

		if ( (abs(y[k+2]) < eps_s) && (abs(y[k+2+nb]) < eps_v) )
		{
			double T0 = bwd_coup[idx].get_init_force();
			
			if (abs(x[5]) <= T0)
				P2[idx] = x[5];
			else
				P2[idx] = T0;
		}
		else
			P2[idx] = bwd_coup[idx].get_force(y[k+2], y[k+2+nb]) + 
					  get_gap_force(y[k+2], y[k+2+nb], -0.09, 0.09);

		Q[0][0] = F[idx] + R1[idx] - R2[idx] - B[idx];
		Q[1][0] = R1[idx] - P1[idx];
		Q[2][0] = R2[idx] - P2[idx];

		gaussLEF_solver(Theta, Q, a);

		//writeln(a);

		return a;
	}
}