module	main;

import	Train;
import	std.stdio;
import	std.conv;
import	LuaScript;

CTrainModel	model = null;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
void main()
{
	model = new CTrainModel();

	int err = model.init("../../cfg/train.lua");

	if (err == -1)
		return;
}
