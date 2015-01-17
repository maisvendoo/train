module	main;

import	Train;
import	std.stdio;

CTrainModel	model = null;

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
}
