//-------------------------------------------------------------------
//
//		Main resistence forces calculation
//		(c) maisvendoo, 2015/01/25
//
//-------------------------------------------------------------------
module	MainResForces;

import	MathFuncs;
import	physics;

enum	int		AXIS		= 4;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double get_main_res(double v, double m)
{
	double res_force = 0;
	double v_kmh = abs(v)*kmh;

	double q0 = m/AXIS/TONN;

	double w = 0;

	if (q0 <= 6)
	{
		w = 1.0 + 0.044*v_kmh + 0.00024*v_kmh*v_kmh;
	}
	else
	{
		w = 0.7 + (3 + 0.1*v_kmh + 0.0025*v_kmh*v_kmh)/q0;
	}

	res_force = (m/TONN)*w*g;

	return res_force;
}