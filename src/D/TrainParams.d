//-------------------------------------------------------------------
//
//		Structures for model initialization
//		(c) maisvendoo, 2015/01/27
//
//-------------------------------------------------------------------
module TrainParams;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
struct	TTrainParams
{
	TSolverParams	solver;
	TTrainModel		train;
	TCoupling		coupling;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
struct TSolverParams
{
	string		method;
	double		init_time;
	double		stop_time;
	double		step;
	double		max_step;
	double		local_err;

	this(this)
	{
		this.method		= "rkf5";
		this.init_time	= 0;
		this.stop_time	= 1.0;
		this.step		= 1e-4;
		this.max_step	= 0.1;
		this.local_err	= 1e-8;
	}
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
struct	TTrainModel
{
	int			vehicles_num;
	string		coupling_type;
	double		railway_coord;
	double		init_velocity;
	double		mass_coeff;
	double		mass_coup;
	double		payload_mass;
	double		empty_mass;
	double		payload_coeff;
	double		loco_section_mass;
	int			loco_sections_num;
	double		delta_eps;

	this(this)
	{
		this.vehicles_num		= 1;
		this.coupling_type		= "default";
		this.railway_coord		= 1000;
		this.init_velocity		= 0;
		this.mass_coeff			= 0.02;
		this.payload_mass		= 60e3;
		this.empty_mass			= 25e3;
		this.payload_coeff		= 1;
		this.loco_section_mass	= 96e3;
		this.loco_sections_num	= 1;
		this.delta_eps			= -1;
	}
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
struct	TCoupling
{
	double		c_1;
	double		c_2;
	double		c_k;
	double		beta;
	double		T0;
	double		t0;
	double		lambda;
	double		delta;
	double		T_max;

	this(this)
	{
		this.c_1	= 2.57e7;
		this.c_2	= 2.85e6;
		this.c_k	= 2.5e8;
		this.beta	= 100;
		this.T0		= 240e3;
		this.t0		= 50e3;
		this.lambda	= 0.09;
		this.delta	= 0.05;
		this.T_max	= 2.45e6;
	}
}