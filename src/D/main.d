//-------------------------------------------------------------------
//
//		Main application module
//		(c) maisvendoo, 2015/01/13
//
//-------------------------------------------------------------------
module	main;

import	Train;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void main()
{
	// Model creation
	CTrainModel model = new CTrainModel();

	// Model initialization
	int err = model.init("../../cfg/train.lua");

	if (err == -1)
		return;

	model.start();
}