module	ConstReacts;

import	std.math;

enum	double	STIFF		= 1e9;
enum	double	DAMP_COEFF	= 1e7;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double get_kf_force(double ds, double dv)
{
	double	c = STIFF;
	double	b = DAMP_COEFF;

	return c*ds + b*dv;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double get_mg_force(double ds, double dv)
{
	double	c = STIFF;
	double	b = DAMP_COEFF;
	double	force = 0;

	if (ds <= 0)
		force = 0;
	else
		force = c*pow(ds, 1.5) + b*pow(ds, 0.25)*dv;

	return force;
}

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
double get_gap_force(double ds, double dv, double ds_min, double ds_max)
{
	double force = 0;

	if ( (ds >= ds_min) && (ds <= ds_max) )
		force = 0;

	if (ds < ds_min)
	{
		force = -get_mg_force(ds_min - ds, dv);
	}

	if (ds > ds_max)
	{
		force = get_mg_force(ds - ds_max, dv);
	}

	return force;
}