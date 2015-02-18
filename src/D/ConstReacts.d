//-------------------------------------------------------------------
//
//		Contact and strike reaction calculation module
//		(c) maisvendoo, 2015/01/19
//
//-------------------------------------------------------------------
module	ConstReacts;

import	std.math;

enum	double	STIFF		= 1e8;
enum	double	DAMP_COEFF	= 1e5;

//-------------------------------------------------------------------
//		Kelvin-Focht force
//-------------------------------------------------------------------
double get_kf_force(double ds, double dv)
{
	double	c = STIFF;
	double	b = DAMP_COEFF;

	return c*ds + b*dv;
}

//-------------------------------------------------------------------
//		Modifed Herz contact-strike force
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

	if (force < 0)
		force = 0;

	return force;
}

//-------------------------------------------------------------------
//		Reaction in mechanism's gaps
//-------------------------------------------------------------------
double get_gap_force(double ds, double dv, double ds_min, double ds_max)
{
	double force = 0;
	double c = STIFF;
	double b = DAMP_COEFF;

	if ( (ds >= ds_min) && (ds <= ds_max) )
		force = 0;

	if (ds < ds_min)
	{
		//force = -c*pow(ds_min - ds, 1.5) + b*pow(ds_min - ds, 0.25)*dv;
		force = -get_mg_force(ds_min - ds, -dv);
	}

	if (ds > ds_max)
	{
		//force = c*pow(ds - ds_max, 1.5) + b*pow(ds - ds_max, 0.25)*dv;
		force = get_mg_force(ds - ds_max, dv);
	}

	return force;
}