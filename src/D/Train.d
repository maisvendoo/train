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

//-------------------------------------------------------------------
//		General train  model class
//-------------------------------------------------------------------
class	CTrainModel: CModel
{
	//---- Train model parameters
	private
	{
		uint		nv;
		uint		ode_dim;

		double[]	m;

		CLogFile	terminal;
	}

	//---- Constructor
	this()
	{
		this.nv = 1;

		terminal = new CLogFile();
		terminal.init();
		terminal.set_print_func(&this.term_out);

		set_dimension(2);
	}

	// Motion ODE system
	override protected void ode_system(double[] Y,
									   ref double[] dYdt, 
									   double t) 
	{
		dYdt[0] = Y[1];
		dYdt[1] = 2;
	}

	// Simulation progress
	override void process() 
	{
		while (t < t_end)
		{
			terminal.print(0.01, dt);
			step();
			t += dt;
		}
	}

	// Terminal print by simulation process
	private void term_out(File term)
	{
		term.writefln("t = %f y = %f", t, y[0]);
	}
}