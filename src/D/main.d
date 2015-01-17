module	main;

import	Train;
import	std.stdio;
import	core.thread;

CTrainModel	model;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void main()
{
	// Model creation
	model = new CTrainModel();

	// Model initialization
	int err = model.init("../../cfg/train.lua");

	if (err == -1)
		return;

	// Starting simulation thread
	Thread sim_tread = new Thread(&model.process);

	sim_tread.start();
}