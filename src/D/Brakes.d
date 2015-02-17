//-------------------------------------------------------------------
//
//		Simplex brakes model
//		(c) maisvendoo, 2015/01/24
//
//-------------------------------------------------------------------
module	Brakes;

import	physics;
import	MathFuncs;
import	std.stdio;

//-------------------------------------------------------------------
//		Braking modes
//-------------------------------------------------------------------
enum	int		FAST_RELEASE		= 0;
enum	int		RELEASE				= 1;
enum	int		SLOW_BRAKE			= 2;
enum	int		SERVICE_BRAKE		= 3;
enum	int		HOLD				= 4;
enum	int		EMERGENCY_BRAKE		= 5;

enum	double	VEHICLE_LENGTH		= 20.0;

//-------------------------------------------------------------------
//		Braking wave speed
//-------------------------------------------------------------------
enum	double	SERVICE_WAVE		= 280.0;
enum	double	EMERGENCY_WAVE		= 300.0;
enum	double	RELEASE_WAVE		= 280.0;

enum	double	MAX_CYL_PRESS		= 4.0e5;
enum	double	MIN_CYL_PRESS		= 0;
enum	double	LADEN_PRESS_STEP	= 1.2e5;
enum	double	MIDDLE_PRESS_STEP	= 0.8e5;
enum	double	EMPTY_PRESS_STEP	= 0.8e5;

enum	int		LADEN_MODE			= 0;
enum	int		MIDDLE_MODE			= 1;
enum	int		EMPTY_MODE			= 2;

enum	double	dpdt_BRAKE			= 0.3e5;
enum	double	dpdt_RELEASE		= 0.05e5;

enum	double	SHOE_FORCE			= 3.80;
enum	int		AXIS				= 4;
enum	int		SHOES				= 2;

enum	double	TONN_2_FORCE		= 9810.0;
//enum	double	kmh					= 3.6;

enum	double	pM_max				= 5.5e5;
enum	double	pM_min				= 0.0;

enum	double	dpMdt_BRAKE			= 0.05e5;
enum	double	dpMdt_RELEASE		= 0.05e5;
enum	double	dpMdt_HOLD			= 0.0;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
class	CBrakes
{
	protected
	{
		int			train_valve_pos;
		int			brake_mode;

		int			nv;

		double[]	pc;
		int[]		vehicle_brake_mode;
		double[]	brake_force;
		double[]	vehicle_coord;

		double		pM;
		double		dpMdt;

		double		mode_time;
		double		wave_speed;

		double[]	dpdt;

		double		max_cyl_press;
		double		min_cyl_press;

		double		press_step;
		double		laden_press_step;
		double		middle_press_step;
		double		empty_press_step;

		int			valve_mode;
	}

	this()
	{
		this.nv = 0;
		this.mode_time = 0;
		this.wave_speed = SERVICE_WAVE;

		this.train_valve_pos = brake_mode = RELEASE;

		this.max_cyl_press = MAX_CYL_PRESS;
		this.min_cyl_press = MIN_CYL_PRESS;
		this.press_step = LADEN_PRESS_STEP;

		this.pM = pM_max;
		this.dpMdt = dpMdt_HOLD;
	}

	protected
	{
		void loco_eq_process(double t, double dt)
		{
			if (brake_mode != train_valve_pos)
			{
				brake_mode = train_valve_pos;
				mode_time = 0;
			}

			switch (brake_mode)
			{
				case RELEASE:
				{
					wave_speed = RELEASE_WAVE;
					break;
				}

				case SERVICE_BRAKE:
				{
					wave_speed = SERVICE_WAVE;
					break;
				}

				case HOLD:
				{
					wave_speed = SERVICE_WAVE;
					break;
				}
			}
		}


		void brake_pipe_process(double t, double dt)
		{
			for (int i = 0; i < nv; i++)
			{
				if (mode_time >= vehicle_coord[i]/wave_speed)
					vehicle_brake_mode[i] = brake_mode;

				switch (vehicle_brake_mode[i])
				{
					case RELEASE:
					{
						dpdt[i] = -dpdt_RELEASE;
						break;
					}

					case SERVICE_BRAKE:
					{
						dpdt[i] = dpdt_BRAKE;

						if (pc[i] < 0.95*press_step)
							pc[i] = press_step;

						break;
					}

					case HOLD:
					{
						dpdt[i] = 0;
						break;
					}
				}
			}
		}

		void brake_cylinder_process(double t, double dt)
		{
			for (int i = 0; i < nv; i++)
			{
				pc[i] += dpdt[i]*dt;

				if (pc[i] > max_cyl_press)
					pc[i] = max_cyl_press;

				if (pc[i] < min_cyl_press)
					pc[i] = min_cyl_press;
			}
		}
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	int init()
	{
		int	err = 0;

		if (nv > 0)
		{
			pc 					= new double[nv];
			brake_force			= new double[nv];
			vehicle_coord		= new double[nv];
			dpdt				= new double[nv];
			vehicle_brake_mode	= new int[nv];

			err = 0;

			for (int i = 0; i < nv; i++)
			{
				pc[i] = 0;
				brake_force[i] = 0;
				dpdt[i] = 0;

				vehicle_coord[i] = VEHICLE_LENGTH*(cast(double) i + 0.5);
				vehicle_brake_mode[i] = RELEASE;
			}
		}
		else
		{
			err = -1;
		}

		return err;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void process(double t, double dt)
	{
		loco_eq_process(t, dt);
		brake_pipe_process(t, dt);
		brake_cylinder_process(t, dt);

		mode_time += dt;
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double get_brake_force(int idx, double v)
	{
		double	k1 = SHOE_FORCE / max_cyl_press;
		double	K = k1*pc[idx];

		double	v_kmh = abs(v)*kmh;

		double phi =  0.6*(16*K + 100)*(v_kmh + 100)/(80*K + 100)/(5*v_kmh + 100);

		double shoe_brake = K*phi*TONN_2_FORCE;

		brake_force[idx] = shoe_brake*AXIS*SHOES;

		return brake_force[idx];
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void set_vehicles_num(int nv)
	{
		this.nv = nv;
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	void set_valve_pos(int train_valve_pos)
	{
		this.train_valve_pos = train_valve_pos;
	}
}