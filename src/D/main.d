//-------------------------------------------------------------------
//
//		Main application module
//		(c) maisvendoo, 2015/01/13
//
//-------------------------------------------------------------------
module	main;

import	std.stdio;
import	std.getopt;

import	Train;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
struct TTrainParams
{
	int		nv;
	double	railway_coord;
	double	v0;
	double	delta_eps;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void main(string[] args)
{
	// Get config Lua script name from command line
	string	cfg_script;

	// Read options from command line
	getopt(args, "cfg", &cfg_script);

	// Model creation
	CTrainModel model = new CTrainModel();

	// Model initialization
	int err = model.init(cfg_script);

	if (err == -1)
	{
		writeln("FAIL: train model is't inizialized");
		return;
	}

	model.start();
}