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

import	CondDestruct;

import	Brakes;

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
		uint		nv;				// Numbel of vehicles
		uint		nb;				// Number of bodies in train model
		uint		ode_dim;		// ODE system dimencion
		uint		mass_n;			// Number of mass in vahicle model

		string		couplig_type;	// Coupling type

		double		v0;				// Initial velocity
		double		railway_coord;	// Railway coordinate

		double[]	m;				// Vehicles mass array

		CLogFile	terminal;		
		CLogFile	file_log;

		CLuaScript	lua_cfg;

		double[]	F;				// Traction forces
		double[]	R1;				// Forward gap reaction
		double[]	R2;				// Backward gap reaction
		double[]	P1;				// Forward coupling force
		double[]	P2;				// Backward coupling force
		double[]	B;				// Brake force
		double[]	Bmax;			// Maximal brake force

		CEFCoupling[]	fwd_coup;	// Forward couplings
		CEFCoupling[]	bwd_coup;	// Backward couplings

		double			lambda;		// Coupling movement range
		double			delta;		// Coupling gap value

		double[]		dy;			// Initial gap values

		// Calculated data arrays
		double[]		Time;		// Simulation time points
		double[][]		x;			// Vehicle body movement
		double[][]		s1;			// Forward coupling relative movement
		double[][]		s2;			// Backward coupling relative movement
		double[][]		v;			// Velocities
		double[][]		u1;			// Forward coupling relative velocity
		double[][]		u2;			// Backward coupling relative velocity
		double[][]		Trac;		// Traction forces
		double[][]		W;			// Resistens forces
		double[][]		G;			// Gravity forces
		double[][]		FwdCoup;	// Forward coupling forces
		double[][]		BwdCoup;	// Backward coupling forces
		double[][]		FwdGap;		// Forward gap force
		double[][]		BwdGap;		// Backward gap force
		double[][]		Power;		// Power of the all forces


		// Data registator parameters
		bool			first_reg;	// First registration flag
		double			reg_time;	// Registration time count
		double			reg_dtime;	// Registration time intervel

		double			term_dtime;
		double			log_dtime;

		CBrakes			brakes;
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

		this.lambda = 0.09;
		this.delta = 0.01;

		this.first_reg = true;
		this.reg_time = 0;

		this.reg_dtime = 0.01;
		this.term_dtime = 0.01;
		this.log_dtime = 0.01;

		this.brakes = new CBrakes();
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
		double eps_v = 1e-3;

		while ( (t <= t_end) && (abs(y[nb+1]) >= eps_v) )
		{
			save_reg_data(t, reg_dtime, dt);

			// Data output (for debug)
			terminal.print(term_dtime, dt);
			file_log.print(log_dtime, dt);

			// Integration step
			brakes.process(t, dt);
			step();

			// Couplings reset
			for (int i = 0; i < nv; i++)
			{
				fwd_coup[i].reset();
				bwd_coup[i].reset();
			}

			// Simulation time increnent
			t += dt;
		}

		error_estimate();

		double J = get_cond_destruct(BwdCoup)	+ 
			       get_cond_destruct(FwdCoup)	+
				   get_cond_destruct(BwdGap)	+
				   get_cond_destruct(FwdGap);

		stdout.writefln("Conditional destruction: %f", J);
	}



	//------------------------------------------------------------------
	//		Terminal print by simulation process
	//------------------------------------------------------------------
	private void term_out(File term)
	{
		term.writefln("t = %10.4f x0 = %10.4f v0 = %10.4f s10 = %10.4f s20 = %10.4f x1 = %10.4f s11 = %10.4f s21 = %10.4f h = %10.4f ms", t, 
			          y[1],
					  y[nb+1],	
			          y[0], 
			          y[2], 
			          y[4], 
			          y[3], 
			          y[5],
					  dt*1000);
	}

	private void file_out(File file)
	{
		int i = 50;
		int k = i*mass_n;

		file.writefln("%f %f %f %f", t, 
			          y[k+2], 
			          R2[i],
					  y[nb+1]);
	}



	//------------------------------------------------------------------
	//		Initialization form Lua-script file
	//------------------------------------------------------------------
	int init(string cfg_name)
	{
		int err = 0;

		if (read_lua_config(cfg_name) == -1)
			return -1;

		brakes.set_vehicles_num(nv);
		brakes.init();
		brakes.set_valve_pos(SERVICE_BRAKE);

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

		// Set ODE solver parameters
		err = solver_init();

		// Set train model parameters
		err = train_init();

		err = mass_init();

		if (err == -1)
			return err;

		err = couplings_init();

		err = set_initc();

		err = init_reg_data();

		return err;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int solver_init()
	{
		int err = LUA_S_OK;

		// method
		string method = lua_cfg.get_str_field("solver_params", "method", err);

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

		return 0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int train_init()
	{
		int err = 0;

		// number of vehicles
		nv = lua_cfg.get_int_field("train_model", "vehicles_num", err);
		
		if (err == LUA_S_NOEXIST)
			nv = 1;
		
		nb = mass_n*nv;
		
		m = new double[nb];
		
		ode_dim = 2*nb;
		
		set_dimension(ode_dim);
		
		R1 = new double[nv];
		R2 = new double[nv];
		B = new double[nv];
		Bmax = new double[nv];
		F = new double[nv];
		P1 = new double[nv];
		P2 = new double[nv];
		
		for (int i = 0; i < nv; i++)
		{
			R1[i] = R2[i] = 0;
			B[i] = Bmax[i] = 0;
			F[i] = 0;
			P1[i] = P2[i] = 0;
		}
		
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

		return 0;
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

		return 0;
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int couplings_init()
	{
		int err = 0;

		double c_1 = lua_cfg.get_double_field("coupling_params", "c_1", err);

		if (err == LUA_S_NOEXIST)
			c_1 = 2.57e7;

		double c_2 = lua_cfg.get_double_field("coupling_params", "c_2", err);

		if (err == LUA_S_NOEXIST)
			c_2 = 2.85e6;

		double c_k = lua_cfg.get_double_field("coupling_params", "c_k", err);

		if (err == LUA_S_NOEXIST)
			c_k = 2.5e8;

		double beta = lua_cfg.get_double_field("coupling_params", "beta", err);

		if (err == LUA_S_NOEXIST)
			beta = 0;

		double T0 = lua_cfg.get_double_field("coupling_params", "T0", err);

		if (err == LUA_S_NOEXIST)
			T0 = 240e3;

		double t0 = lua_cfg.get_double_field("coupling_params", "t0", err);

		if (err == LUA_S_NOEXIST)
			t0 = 50e3;

		lambda = lua_cfg.get_double_field("coupling_params", "lambda", err);

		if (err == LUA_S_NOEXIST)
			lambda = 0.09;

		delta = lua_cfg.get_double_field("coupling_params", "delta", err);

		if (err == LUA_S_NOEXIST)
			delta = 0.01;

		fwd_coup = new CEFCoupling[nv];
		bwd_coup = new CEFCoupling[nv];

		for (int i = 0; i < nv; i++)
		{
			fwd_coup[i] = new CEFCoupling();
			bwd_coup[i] = new CEFCoupling();

			fwd_coup[i].set_stiffs(c_1, c_2, c_k);
			fwd_coup[i].set_damp_coeff(beta);
			fwd_coup[i].set_init_force(T0);
			fwd_coup[i].set_release_init_force(t0);

			bwd_coup[i].set_stiffs(c_1, c_2, c_k);
			bwd_coup[i].set_damp_coeff(beta);
			bwd_coup[i].set_init_force(T0);
			bwd_coup[i].set_release_init_force(t0);

			fwd_coup[i].reset();
			bwd_coup[i].reset();
		}

		return err;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int set_initc()
	{
		int err = 0;

		dy = new double[nv-1];

		// Set movements
		for (int i = 0; i < ode_dim; i++)
		{
			y[i] = dydt[i] = 0;
		}

		for (int i = 0; i < nv-1; i++)
		{
			dy[i] = lua_cfg.get_double_field("delta_initc", i, err);

			if (err == LUA_S_NOTTABLE)
				dy[i] = 0;
		}

		for (int i = 0; i < nv-1; i++)
		{
			int k = i*mass_n;

			y[k+4] = y[k+1] - dy[i] - (y[k+2] + y[k+3]);
		}

		// Set velocities
		for (int i = 0; i < nv; i++)
		{
			int k = i*mass_n;

			y[k+nb] = 0;
			y[k+1+nb] = v0;
			y[k+2+nb] = 0;
		}

		return err;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double[] get_accels(double[] Y, double t, int idx)
	{
		double	eps_v	= 1e-3;
		double	eps_s	= 1e-4;
		int		err		= 0;

		int k = mass_n*idx;

		double[] a = new double[mass_n];

		// Gap forces calculation
		if (idx == 0)
			R1[idx] = 0;
		else
			R1[idx] = R2[idx-1]; 

		if (idx == nv-1)
			R2[idx] = 0;
		else
			R2[idx] = get_gap_force(y[k+1] - y[k+4] - (y[k+2] + y[k+3]), 
									y[k+1+nb] - y[k+4+nb] - (y[k+2+nb] + y[k+3+nb]),
									-delta/2,
									 delta/2);

		// Other active forces calculation
		F[idx] = 0;

		// Brake forces calculation
		Bmax[idx] = brakes.get_brake_force(idx, y[k+1+nb]);

		// Resistive force precheck
		double Fa = F[idx] + R1[idx] - R2[idx];

		if ( (abs(y[k+1+nb]) < eps_v) && (abs(Fa) <= Bmax[idx]) )
		{
			B[idx] = Fa;
		}
		else
		{
			B[idx] = Bmax[idx]*sign(y[k+1+nb]);
		}

		//  Vehicle body acceleration
		a[1] = (F[idx] + R1[idx] - R2[idx] - B[idx])/(m[k] + m[k+1] + m[k+2]);

		// Forward coupling force
		if (abs(y[k]) < eps_s)
		{
			P1[idx] = R1[idx] - m[k]*a[1];
		}
		else
		{
			P1[idx] = fwd_coup[idx].get_force(y[k], y[k+nb]) + 
				      get_gap_force(y[k], y[k+nb], -lambda, lambda);
		}

		// Backward coupling force
		if (abs(y[k+2]) < eps_s)
		{
			P2[idx] = R2[idx] + m[k+2]*a[1];
		}
		else
		{
			P2[idx] = bwd_coup[idx].get_force(y[k+2], y[k+2+nb]) + 
				      get_gap_force(y[k+2], y[k+2+nb], -lambda, lambda);
		}

		// Brake force calculation
		if (abs(y[k+1]) < eps_v)
		{
			if (abs(B[idx]) > Bmax[idx])
				B[idx] = Bmax[idx]*sign(B[idx]);
		}
		else
			B[idx] = Bmax[idx]*sign(y[k+1+nb]);

		// Coupling forces calculation
		if (abs(y[k]) < eps_s)
		{
			double T0 = fwd_coup[idx].get_init_force();

			if (abs(P1[idx]) > T0)
				P1[idx] = T0;
		}
		else
			P1[idx] = fwd_coup[idx].get_force(y[k], y[k+nb]) + 
			          get_gap_force(y[k], y[k+nb], -lambda, lambda);

		if (abs(y[k+2]) < eps_s)
		{
			double T0 = bwd_coup[idx].get_init_force();
			
			if (abs(P2[idx]) > T0)
				P2[idx] = T0;
		}
		else
			P2[idx] = bwd_coup[idx].get_force(y[k+2], y[k+2+nb]) + 
					  get_gap_force(y[k+2], y[k+2+nb], -lambda, lambda);

		// Final accelleraion calculation
		a[1] = (F[idx] + P1[idx] - P2[idx] - B[idx])/m[k+1];
		a[0] = (R1[idx] - P1[idx] - m[k]*a[1])/m[k];
		a[2] = (R2[idx] - P2[idx] + m[k+2]*a[1])/m[k+2];

		return a;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int init_reg_data()
	{
		x		= new double[][nv];
		s1		= new double[][nv];
		s2		= new double[][nv];
		v		= new double[][nv];
		u1		= new double[][nv];
		u2		= new double[][nv];
		Trac	= new double[][nv];
		W		= new double[][nv];
		G		= new double[][nv];
		FwdCoup	= new double[][nv];
		BwdCoup	= new double[][nv];
		FwdGap	= new double[][nv];
		BwdGap	= new double[][nv];
		Power	= new double[][nv];

		return 0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void save_reg_data(double t, double delta_t, double dt)
	{
		if ( first_reg || (reg_time >= delta_t) )
		{
			first_reg = false;
			reg_time = 0;

			Time ~= t;

			for (int i = 0; i < nv; i++)
			{
				int k = i*mass_n;

				x[i]		~= y[k+1];
				s1[i]		~= y[k];
				s2[i]		~= y[k+2];

				v[i]		~= y[k+1+nb];
				u1[i]		~= y[k+nb];
				u2[i]		~= y[k+2+nb];

				Trac[i]		~= F[i];
				W[i]		~= B[i];
				G[i]		~= 0;
				FwdCoup[i]	~= P1[i];
				BwdCoup[i]	~= P2[i];
				FwdGap[i]	~= R1[i];
				BwdGap[i]	~= R2[i];

				Power[i]	~= F[i]*y[k+1+nb] - B[i]*y[k+1+nb] - P1[i]*y[k+nb] - P2[i]*y[k+2+nb] +
					           R1[i]*y[k+nb] + R2[i]*y[k+2+nb];
			}
		}

		reg_time += dt;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double forces_work()
	{
		double	A = 0;
		int		N = cast(int) Time.length; 

		for (int i = 0; i < nv; i++)
		{
			for (int j = 0; j < N-1; j++)
				A += (Power[i][j] + Power[i][j+1])*(Time[j+1] - Time[j])/2;

		}

		return A;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double kinetic_energy_change()
	{
		double	Ek1 = 0;
		double	Ek2 = 0;
		int		N = cast(int) Time.length;

		for (int i = 0; i < nv; i++)
		{
			int k = i*mass_n;

			Ek1 += (m[k]*pow(v[i][0] + u1[i][0], 2)/2 + m[k+1]*pow(v[i][0], 2)/2 + m[k+2]*pow(v[i][0] - u2[i][0], 2)/2);
			Ek2 += (m[k]*pow(v[i][N-1] + u1[i][N-1], 2)/2 + m[k+1]*pow(v[i][N-1], 2)/2 + m[k+2]*pow(v[i][N-1] - u2[i][N-1], 2)/2);
		}

		return Ek2 - Ek1;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void error_estimate()
	{
		double Af = forces_work();
		double dEk = kinetic_energy_change();
		
		stdout.writefln("Forces work: %.2f J", Af);
		stdout.writefln("Change of kinetic energy: %.2f J", dEk);
		stdout.writefln("Relative integration error: %.2f %c", 
			            abs((Af - dEk)*100/dEk), '%');
	}
}