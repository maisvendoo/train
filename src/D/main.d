//-------------------------------------------------------------------
//
//		Main application module
//		(c) maisvendoo, 2015/01/13
//
//-------------------------------------------------------------------
module	main;

import	std.stdio;

import	Train;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void main(string[] args)
{
	// Model creation
	CTrainModel model = new CTrainModel();

	// Model initialization
	int err = model.init(args);

	if (err == -1)
	{
		writeln("FAIL: train model is't inizialized");
		return;
	}

	model.start();
}