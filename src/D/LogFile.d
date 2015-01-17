//-------------------------------------------------------------------
//
//		Module for data output
//		(c) maisvendoo, 2015.01.13
//
//-------------------------------------------------------------------
module	LogFile;

import	std.stdio;

alias	print_func_t = void delegate(File);

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
class	CLogFile
{
	this()
	{
		this.first_print = true;
		this.print_time = 0;
		this.log_file = stdout;
		this.log_file_name = "";
		this.print_func = null;
	}

	bool			first_print;	// First printing flag
	double			print_time;		// Printing time counter

	File			log_file;		// Log file object

	string			log_file_name;	// Log file name

	print_func_t	print_func;		// Printing callback function

	// Log initialization
	void init(string file_name = "")
	{
		if (file_name == "")
		{
			log_file = stdout;
		}
		else
		{
			log_file_name = file_name;
			log_file = File(log_file_name, "wt");
			log_file.close();
		}
	}

	// Print data to log
	void print(double print_dt, double dt)
	{
		if ( first_print || (print_time >= print_dt) )
		{
			if (print_func != null)
			{
				first_print = false;
				print_time = 0;

				if (log_file == stdout)
				{
					print_func(log_file);
				}
				else
				{
					log_file = File(log_file_name, "a+");
					print_func(log_file);
					log_file.close();
				}
			}
			else
			{
				stdout.writefln("FAIL: Print function is't defined");
			}
		}

		print_time += dt;
	}

	// Set print function
	void set_print_func(print_func_t print_func)
	{
		this.print_func = print_func;
	}
}