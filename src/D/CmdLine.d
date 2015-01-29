//-------------------------------------------------------------------
//
//		Structure for reading command line parameters
//		(c) maisvendoo, 2015/01/27
//
//-------------------------------------------------------------------
module CmdLine;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
struct TCmdLine
{
	string	cfg_name;
	int		vehicles_num;
	double	payload_coeff;

	double	railway_coord;
	double	init_velocity;

	double	delta;
	double	delta_eps;

	bool	term_out;
	string	log_file;

	this(this)
	{
		this.cfg_name		= "";
		this.vehicles_num	= -1;
		this.payload_coeff	= -1;
		this.railway_coord	= -1;
		this.init_velocity	= -1;
		this.delta			= -1;
		this.delta_eps		= 2;
		this.term_out		= false;
		this.log_file		= "";
	}
}