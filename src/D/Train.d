//-------------------------------------------------------------------
//
//		Train motion simulation module
//		(c) maisvendoo, 2015.01.13
//
//-------------------------------------------------------------------
module	Train;

import	std.stdio;
import	LogFile;
import	Model;
import	LuaScript;

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
		uint		ode_dim;
		string		couplig_type;

		double		v0;
		double		railway_coord;

		double[]	m;

		CLogFile	terminal;
		CLuaScript	lua_cfg;
	}



	//------------------------------------------------------------------
	//		Constructor
	//------------------------------------------------------------------
	this()
	{
		this.nv = 1;
		this.couplig_type = "default";

		this.v0 = 0;
		this.railway_coord = 1e6;

		terminal = new CLogFile();
		terminal.init();
		terminal.set_print_func(&this.term_out);

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
		dYdt[0] = Y[1];
		dYdt[1] = 2;
	}



	//------------------------------------------------------------------
	//		Simulation progress
	//------------------------------------------------------------------
	override void process() 
	{
		while (t < t_end)
		{
			terminal.print(0.01, dt);
			step();
			t += dt;
		}
	}



	//------------------------------------------------------------------
	//		Terminal print by simulation process
	//------------------------------------------------------------------
	private void term_out(File term)
	{
		term.writefln("t = %f y = %f", t, y[0]);
	}



	//------------------------------------------------------------------
	//		Initialization form Lua-script file
	//------------------------------------------------------------------
	int init(string cfg_name)
	{
		int err = LUA_S_OK;

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

		return err;
	}
}