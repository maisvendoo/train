---------------------------------------------------------------------
--
--	Train simulation config 
--
---------------------------------------------------------------------

local	kmh		= 3.6
local	km		= 1000.0

---------------------------------------------------------------------
--	ODE solver tunning
---------------------------------------------------------------------
solver_params = 
{
	method		= "rkf5",	-- integration method
	init_time	= 0,		-- initial time
	stop_time	= 180.0,		-- stop simulation time
	step		= 1e-4,		-- time step
	max_step	= 0.1,		-- maximal time step
	local_err	= 1e-9		-- loacal solver error
}

---------------------------------------------------------------------
--	Train model parameters
---------------------------------------------------------------------
train_model = 
{
	vehicles_num 	= 70,
	coupling_type	= "default",
	railway_coord	= 1500,			-- initial railway coordinate (km)
	init_velocity	= 0			-- initial velocity (km/h)
}

train_model.railway_coord = train_model.railway_coord*km
train_model.init_velocity = train_model.init_velocity / kmh

vehicle_mass = {}

mass_coeff = 0.02

for i = 1, train_model.vehicles_num do

	vehicle_mass[i] = 85e3

end

Traction = function(t)

	Fmax = 600e3
	dFdt = 20e3
	tmax = Fmax /dFdt
	force = 0

	if (t <= tmax) then
		force = dFdt*t
	else
		force = Fmax
	end

	return force

end