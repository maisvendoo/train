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
import	ODEqs;

//-------------------------------------------------------------------
//		Braking modes
//-------------------------------------------------------------------
enum	int		FAST_RELEASE		= 0;
enum	int		RELEASE				= 1;
enum	int		SLOW_BRAKE			= 2;
enum	int		SERVICE_BRAKE		= 3;
enum	int		HOLD				= 4;
enum	int		EMERGENCY_BRAKE		= 5;

enum	double	VEHICLE_LENGTH		= 17.5;

//-------------------------------------------------------------------
//		Braking wave speed
//-------------------------------------------------------------------
enum	double	SERVICE_WAVE		= 280.0;
enum	double	EMERGENCY_WAVE		= 300.0;
enum	double	RELEASE_WAVE		= 70.0;

enum	double	MAX_CYL_PRESS		= 4.5e5;
enum	double	MIN_CYL_PRESS		= 0;
enum	double	LADEN_PRESS_STEP	= 1.2e5;
enum	double	MIDDLE_PRESS_STEP	= 0.8e5;
enum	double	EMPTY_PRESS_STEP	= 0.8e5;

enum	int		LADEN_MODE			= 0;
enum	int		MIDDLE_MODE			= 1;
enum	int		EMPTY_MODE			= 2;

enum	double	dpdt_BRAKE			= 0.13e5;
enum	double	dpdt_RELEASE		= 0.05e5;

enum	double	SHOE_FORCE			= 3.8;
enum	int		AXIS				= 4;
enum	int		SHOES				= 2;

enum	double	TONN_2_FORCE		= 9810.0;

enum	double	pM_max				= 5.5e5;
enum	double	pM_min				= 0.0;

enum	double	dpMdt_BRAKE			= 0.3e5;
enum	double	dpMdt_RELEASE		= 0.2e5;
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
		double[]	pc_ref;
		int[]		vehicle_brake_mode;
		double[]	brake_force;
		double[]	vehicle_coord;

		double[]	pM;
		double[]	pM_ref;
		double[]	dpMdt;
		double[]	dpM;

		double		mode_time;
		double		wave_speed;

		double[]	dpdt;
		double[]	gamma;

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
		this.press_step = EMPTY_PRESS_STEP;
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
					gamma[0] = 0.04;
					pM_ref[0] = pM_max;
					break;
				}

				case SERVICE_BRAKE:
				{
					wave_speed = SERVICE_WAVE;
					gamma[0] = 0.02;
					pM_ref[0] = 0;
					break;
				}

				case HOLD:
				{
					wave_speed = SERVICE_WAVE;
					gamma[0] = 0;
					pM_ref[0] = 0;
					break;
				}
			}
		}


		void brake_pipe_process(double t, double dt)
		{
			for (int i = 1; i < nv; i++)
			{
				double tau = vehicle_coord[i]/wave_speed;

				if (mode_time >= tau)
				{
					vehicle_brake_mode[i] = brake_mode;
					pM[i] = pM[0];
				}
			}

			rk4_solver_step(pM, dpMdt, t, dt, 0.1, 1e-8, &this.pipe_eqs);
		}

		void brake_cylinder_process(double t, double dt)
		{
			rk4_solver_step(pc, dpdt, t, dt, 0.1, 1e-8, &this.brake_cyl_eqs);
		}

		void brake_cyl_eqs(double[] pc, ref double[] dpcdt, double t)
		{
			for (int i = 0; i < nv; i++)
			{
				double	ref_p = get_ref_pressure(get_dpM(i));
				double	k = 10e-6*(ref_p - pc[i]);
				double	gamma_b = 0.08;
				double	gamma_r = 0.15;

				if (k >=0 )
				{
					pc_ref[i] = max_cyl_press;

					if (k > gamma_b)
						k = gamma_b;
				}
				else
				{
					pc_ref[i] = 0;

					k = abs(k);

					if (k > gamma_r)
						k = gamma_r;
				}		

				dpcdt[i] = k*(pc_ref[i] - pc[i]);
			}
		}

		void pipe_eqs(double[] pM, ref double[] dpMdt, double t)
		{
			//for (int i = 0; i < nv; i++)
			//{
				dpMdt[0] = gamma[0]*(pM_ref[0] - pM[0]);
			//}
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
			pc_ref				= new double[nv];
			brake_force			= new double[nv];
			vehicle_coord		= new double[nv];
			dpdt				= new double[nv];
			vehicle_brake_mode	= new int[nv];
			gamma				= new double[nv];
			dpM					= new double[nv];

			pM					= new double[nv];
			pM_ref				= new double[nv];
			dpMdt				= new double[nv];

			err = 0;

			for (int i = 0; i < nv; i++)
			{
				pc[i] = 0;
				pc_ref[i] = 0;
				brake_force[i] = 0;
				dpdt[i] = 0;
				dpM[i] = 0;

				pM[i] = pM_max;
				dpMdt[i] = 0;
				gamma[i] = 0;

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
		double	p0 = 4e4;
		double	k1 = SHOE_FORCE / (max_cyl_press - p0);
		double	K = 0;

		if (pc[idx] > p0)
			K = k1*(pc[idx] - p0);
		else
			K = 0;

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

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double get_pipe_pressure(int idx)
	{
		return pM[idx];
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double get_dpM(int idx)
	{
		return pM_max - pM[idx];
	}

	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double get_cylinder_pressure(int idx)
	{
		return pc[idx];
	}



	//---------------------------------------------------------------
	//
	//---------------------------------------------------------------
	double get_ref_pressure(double dp)
	{
		double ref_p = 0;

		if (dp < 0.4e5)
			ref_p = 0;

		if ( (dp >= 0.4e5) && (dp < 0.5e5) )
			ref_p = 1.25e5;

		if ( (dp >= 0.5e5) && (dp < 1.2e5) )
			ref_p = 1.25e5 + 4.21*(dp - 0.5e5);

		if (dp > 1.2e5)
			ref_p = 4.2e5;

		return ref_p;
	}
}