//-------------------------------------------------------------------
//
//		Train motion simulation module
//		(c) maisvendoo, 2015/01/13
//
//-------------------------------------------------------------------
module	Train;

import	std.stdio;
import	std.getopt;

import	physics;
import	TrainParams;
import	CmdLine;
import	MathFuncs;
import	LogFile;
import	Model;
import	LuaScript;
import	matrix;
import	EFCoupling;
import	LinearEQs;
import	ConstReacts;
import	MainResForces;

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
		double[][]		pc;
		double[]		pM;


		// Data registator parameters
		bool			first_reg;	// First registration flag
		double			reg_time;	// Registration time count
		double			reg_dtime;	// Registration time intervel

		double			term_dtime;
		double			log_dtime;

		CBrakes			brakes;

		double			train_mass;

		TTrainParams	train_params;
		TCmdLine		cmd_line;

		bool			term_out_flag;
		string			log_file_name;

		double			Tmax;
		double			Tmin;
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

		this.lua_cfg = new CLuaScript();

		set_dimension(2);

		this.lambda = 0.09;
		this.delta = 0.01;

		this.first_reg = true;
		this.reg_time = 0;

		this.reg_dtime = 0.01;
		this.term_dtime = 0.001;
		this.log_dtime = 0.01;

		this.brakes = new CBrakes();

		this.train_mass = 0;

		this.term_out_flag = false;
		this.log_file_name = "";

		this.Tmax = 0;
		this.Tmin = 0;
	}



	//------------------------------------------------------------------
	//		Motion ODE system
	//------------------------------------------------------------------
	override protected void ode_system(double[] Y,
									   ref double[] dYdt, 
									   double t) 
	{
		double[] a;

		for (int i = 0; i < nv; i++)
		{
			int k = mass_n*i; 

			dYdt[k]		= Y[k+nb];
			dYdt[k+1]	= Y[k+1+nb];
			dYdt[k+2]	= Y[k+2+nb];

			a = get_accels(Y, t, i);

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
			if (term_out_flag)
				terminal.print(term_dtime, dt);

			/*if (log_file_name != "")
				file_log.print(log_dtime, dt);*/

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

		double error = error_estimate();

		double J = get_cond_destruct(BwdCoup)	+ 
			       get_cond_destruct(FwdCoup)	+
				   get_cond_destruct(BwdGap)	+
				   get_cond_destruct(FwdGap);

		found_max_forces();

		print_log();

		// Program output in stdout
		stdout.writeln(train_mass/TONN);	// 0
		stdout.writeln(J);					// 1
		stdout.writeln(error);				// 2
		stdout.writeln(Tmax/kN);			// 3
		stdout.writeln(Tmin/kN);			// 4
	}



	//------------------------------------------------------------------
	//		Terminal print by simulation process
	//------------------------------------------------------------------
	private void term_out(File term)
	{
		term.writef("Time: %10.4f s | V : %6.2f kmh | pM: %6.3f MPa\r", t, y[1+nb]*kmh, brakes.get_pipe_pressure(0)/1e6);  
	}

	private void file_out(File file)
	{
		file.writef("%f ", t);

		for (int i = 0; i < nv-1; i++)
			file.writef("%f ", R2[i]/1000);

		file.writeln();
	}



	//------------------------------------------------------------------
	//		Initialization form Lua-script file
	//------------------------------------------------------------------
	int init(string[] args)
	{
		int		err = 0;

		// Get config script name from command line args
		with (cmd_line)
		{
			getopt(args, 
			   	   "config", 		&cfg_name,
			   	   "vehicles",		&vehicles_num,
				   "payload_coeff",	&payload_coeff,
				   "railway_coord", &railway_coord,
				   "init_velocity", &init_velocity,
				   "delta",			&delta,
				   "delta_eps",		&delta_eps,
				   "term_out",		&term_out,
				   "log_file",		&log_file
				   );
		}

		if (cmd_line.cfg_name == "")
			cmd_line.cfg_name = "../../cfg/train.lua";

		// Config reading
		if (read_lua_config(cmd_line.cfg_name) == -1)
			return -1;

		cmdline_parse(cmd_line);

		// Set ODE solver parameters
		solver_init();
		train_init();
		mass_init();
		couplings_init();
		set_initc();
		init_reg_data();

		brakes.set_vehicles_num(nv);
		brakes.init();
		brakes.set_valve_pos(SERVICE_BRAKE);

		return err;
	}



	//---------------------------------------------------------------
	//		Set parameters from command line
	//---------------------------------------------------------------
	void cmdline_parse(ref TCmdLine cmd_line)
	{
		with (cmd_line)
		{
			if (vehicles_num > 0)
				train_params.train.vehicles_num = vehicles_num;

			if (payload_coeff >= 0)
			/*{
				train_params.train.payload_coeff = 0.0;
			}
			else*/
			{
				train_params.train.payload_coeff = payload_coeff;
			}

			if (railway_coord > 0)
				train_params.train.railway_coord = railway_coord*km;

			if (init_velocity > 0)
				train_params.train.init_velocity = init_velocity/kmh;

			if (delta > 0)
				train_params.coupling.delta = delta;

			if (abs(delta_eps) < 2)
				train_params.train.delta_eps = delta_eps;

			if (log_file != "")
			{
				log_file_name = log_file;
				file_log.init(log_file_name);
				file_log.set_print_func(&this.file_out);
			}

			if (term_out)
				term_out_flag = true;
		}
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

		// Read solver params
		with (train_params.solver)
		{
			method = lua_cfg.get_str_field("solver_params", "method", err);
			
			if (err == LUA_S_NOEXIST)
				method = "rkf5";

			init_time = lua_cfg.get_double_field("solver_params", "init_time", err);
			
			if (err == LUA_S_NOEXIST)
				init_time = 0;

			stop_time = lua_cfg.get_double_field("solver_params", "stop_time", err);
			
			if (err == LUA_S_NOEXIST)
				stop_time = 1.0;

			step = lua_cfg.get_double_field("solver_params", "step", err);
			
			if (err == LUA_S_NOEXIST)
				step = get_time_step();

			max_step = lua_cfg.get_double_field("solver_params", "max_step", err);
			
			if (err == LUA_S_NOEXIST)
				max_step = 0.1;
									
			local_err = lua_cfg.get_double_field("solver_params", "local_err", err);
			
			if (err == LUA_S_NOEXIST)
				local_err = 1e-8;
		}

		with (train_params.train)
		{
			vehicles_num = lua_cfg.get_int_field("train_model", "vehicles_num", err);

			if (err == LUA_S_NOEXIST)
				vehicles_num = 1;

			// couplig type
			couplig_type = lua_cfg.get_str_field("train_model", "coupling_type", err);
			
			if (err == LUA_S_NOEXIST)
				couplig_type = "default";
			
			// railway coordinate
			railway_coord = lua_cfg.get_double_field("train_model", "railway_coord", err);
			
			if (err == LUA_S_NOEXIST)
				railway_coord = 1e6;
			
			// initial velocity
			init_velocity = lua_cfg.get_double_field("train_model", "init_velocity", err);

			if (err == LUA_S_NOEXIST)
				init_velocity = 0;

			mass_coeff = lua_cfg.get_double_field("train_model", "mass_coeff", err);
			mass_coup = lua_cfg.get_double_field("train_model", "mass_coup", err); 
			payload_mass = lua_cfg.get_double_field("train_model", "payload_mass", err);
			payload_coeff = lua_cfg.get_double_field("train_model", "payload_coeff", err);
			empty_mass = lua_cfg.get_double_field("train_model", "empty_mass", err);
			loco_section_mass = lua_cfg.get_double_field("train_model", "loco_section_mass", err);
			loco_sections_num = lua_cfg.get_int_field("train_model", "loco_sections_num", err);

			delta_eps = lua_cfg.get_double_field("train_model", "delta_eps", err); 
		}

		with (train_params.coupling)
		{
			c_1 = lua_cfg.get_double_field("coupling_params", "c_1", err);
			
			if (err == LUA_S_NOEXIST)
				c_1 = 2.57e7;
			
			c_2 = lua_cfg.get_double_field("coupling_params", "c_2", err);
			
			if (err == LUA_S_NOEXIST)
				c_2 = 2.85e6;
			
			c_k = lua_cfg.get_double_field("coupling_params", "c_k", err);
			
			if (err == LUA_S_NOEXIST)
				c_k = 2.5e8;
			
			beta = lua_cfg.get_double_field("coupling_params", "beta", err);
			
			if (err == LUA_S_NOEXIST)
				beta = 0;
			
			T0 = lua_cfg.get_double_field("coupling_params", "T0", err);
			
			if (err == LUA_S_NOEXIST)
				T0 = 240e3;
			
			t0 = lua_cfg.get_double_field("coupling_params", "t0", err);
			
			if (err == LUA_S_NOEXIST)
				t0 = 50e3;
			
			lambda = lua_cfg.get_double_field("coupling_params", "lambda", err);
			
			if (err == LUA_S_NOEXIST)
				lambda = 0.09;
			
			delta = lua_cfg.get_double_field("coupling_params", "delta", err);
			
			if (err == LUA_S_NOEXIST)
				delta = 0.01;

			T_max = lua_cfg.get_double_field("coupling_params", "T_max", err);

			if (err == LUA_S_NOEXIST)
				T_max = 2.45e6;
		}

		return err;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int solver_init()
	{
		int err = LUA_S_OK;

		// method
		with (train_params.solver)
		{
			if (method == "rkf5")
				set_integration_method(&rkf5_solver_step);
		
			if (method == "rk4")
				set_integration_method(&rk4_solver_step);
		
			if (method == "euler")
				set_integration_method(&euler_solver_step);

			if (method == "adams")
				set_integration_method(&adams_solver_step);

			if (method == "impeuler")
				set_integration_method(&impeuler_solver_step);

			if (method == "adams-multhon")
				set_integration_method(&AM_solver_step);

			if (method == "adams-multhon5")
				set_integration_method(&am5_solver_step);
		
			set_init_time(init_time);
			set_stop_time(stop_time);
			set_time_step(step);
			set_max_time_step(max_step);
			set_local_error(local_err);
		}

		return 0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int train_init()
	{
		int err = 0;

		// number of vehicles
		nv = train_params.train.vehicles_num; 
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
		couplig_type = train_params.train.coupling_type;
			
		// railway coordinate
		railway_coord = train_params.train.railway_coord;
			
		// initial velocity
		v0 = train_params.train.init_velocity; 

		return 0;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int mass_init()
	{
		int err = LUA_S_OK;

		double mass_coeff = train_params.train.mass_coeff;
		double mass_coup = train_params.train.mass_coup;

		for (int i = 0; i < nv; i++)
		{
			int k = mass_n*i;

			double mass = 0;

			if ( (i >= 0) && (i < train_params.train.loco_sections_num) )
			{
				mass = train_params.train.loco_section_mass;
			}
			else
			{
				with (train_params.train)
				{
					mass = empty_mass + payload_coeff*payload_mass;
				}
			}

			train_mass += mass;

			/*if (err == LUA_S_NOEXIST)
				return -1;*/

			//m[k] = mass_coeff*mass;
			m[k] = mass_coup;
			//m[k+2] = mass_coeff*mass;
			m[k+2] = mass_coup;
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

		lambda = train_params.coupling.lambda;
		delta = train_params.coupling.delta;

		fwd_coup = new CEFCoupling[nv];
		bwd_coup = new CEFCoupling[nv];

		with (train_params.coupling)
		{
			for (int i = 0; i < nv; i++)
			{
				fwd_coup[i] = new CEFCoupling();
				bwd_coup[i] = new CEFCoupling();

				fwd_coup[i].set_stiffs(c_1, c_2, c_k);
				fwd_coup[i].set_damp_coeff(beta);
				fwd_coup[i].set_init_force(T0);
				fwd_coup[i].set_release_init_force(t0);
				fwd_coup[i].set_s_max(lambda);
				fwd_coup[i].set_T_max(T_max);

				bwd_coup[i].set_stiffs(c_1, c_2, c_k);
				bwd_coup[i].set_damp_coeff(beta);
				bwd_coup[i].set_init_force(T0);
				bwd_coup[i].set_release_init_force(t0);
				bwd_coup[i].set_s_max(lambda);
				bwd_coup[i].set_T_max(T_max);

				fwd_coup[i].reset();
				bwd_coup[i].reset();
			}
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
			dy[i] = train_params.train.delta_eps*delta/2; 
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

		double Ft = lua_cfg.call_func("traction", [t, y[1+nb]*kmh], err);

		for (int i = 0; i < train_params.train.loco_sections_num; i++)
		{
			F[i] = Ft;
		}

		double dpM = brakes.get_dpM(0);

		int valve_pos = cast(int) lua_cfg.call_func("valve_pos", [t, y[1+nb]*kmh, dpM], err);

		brakes.set_valve_pos(valve_pos);

		// Brake forces calculation
		double mv = m[k] + m[k+1] + m[k+2];

		Bmax[idx] = brakes.get_brake_force(idx, y[k+1+nb]) + get_main_res(y[k+1+nb], mv);

		P1[idx] = fwd_coup[idx].get_force(y[k], y[k+nb]) + get_gap_force(y[k], y[k+nb], -lambda, lambda);
		P2[idx] = bwd_coup[idx].get_force(y[k+2], y[k+2+nb]) + get_gap_force(y[k+2], y[k+2+nb], -lambda, lambda);

		// Resistive force precheck
		double Fa = F[idx] + P1[idx] - P2[idx];

		// Brake force calculation
		if (abs(y[k+1]) < ZERO)
		{
			if (abs(Fa) >= Bmax[idx])
				B[idx] = Bmax[idx]*sign(Fa);
			else
				B[idx] = Fa;
		}
		else
			B[idx] = Bmax[idx]*sign(y[k+1+nb]);

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
		pc 		= new double[][nv];

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
			pM ~= brakes.get_pipe_pressure(0);

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

				pc[i]		~= brakes.get_cylinder_pressure(i);
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
	double error_estimate()
	{
		double error = 0;

		double Af = forces_work();
		double dEk = kinetic_energy_change();
		
		if (abs(dEk) < ZERO)
		{
			error = 0;
		}
		else
		{
			error = (1 - abs(Af)/abs(dEk))*100;
		}

		return error;
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void found_max_forces()
	{
		int		N = cast(int) Time.length;

		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < nv; j++)
			{
				if (Tmax < BwdGap[j][i])
					Tmax = BwdGap[j][i];

				if (Tmin > BwdGap[j][i])
					Tmin = BwdGap[j][i];
			}
		}
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void print_log()
	{
		if (log_file_name == "")
			return;

		int N = cast(int) Time.length;

		File log = File(log_file_name, "wt");

		for (int i = 0; i < N; i++)
		{
			log.writef("%f ", Time[i]);

			for (int j = 0; j < nv-1; j++)
			{
				log.writef("%f ", BwdGap[j][i]/kN/g);
			}

			for (int j = 0; j < nv-1; j++)
			{
				//log.writef("%f ", (x[j][i] - x[j+1][i] + s1[j+1][i] + s2[j][i])*1000.0);
				log.writef("%f ", s2[j][i]*1000.0);
			}

			log.writeln();
		}

		log.writeln();
		log.writeln();

		for (int i = 0; i < N; i++)
		{
			log.writef("%f ", Time[i]);

			for (int j = 0; j < nv; j++)
			{
				log.writef("%f ", W[j][i]/kN);
			}

			log.writeln();
		}

		log.writeln();
		log.writeln();

		for (int i = 0; i < N; i++)
		{
			log.writef("%f %f ", Time[i], pM[i]/1e5);

			for (int j = 0; j < nv; j++)
			{
				log.writef("%f ", pc[j][i]/1e5);
			}

			log.writeln();
		}

		log.close();
	}
}