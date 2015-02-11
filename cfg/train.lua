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
	method		= "adams",	-- integration method
	init_time	= 0,		-- initial time
	stop_time	= 1000.0,	-- stop simulation time
	step		= 1e-5,		-- time step
	max_step	= 1e-4,		-- maximal time step
	local_err	= 1e-10		-- local solver error
}

---------------------------------------------------------------------
--	Train model parameters
---------------------------------------------------------------------
train_model = 
{
	vehicles_num	= 70,
	coupling_type	= "default",
	railway_coord	= 1500,			-- initial railway coordinate (km)
	init_velocity	= 10,			-- initial velocity (km/h)
	mass_coeff	= 0.02,
	payload_mass	= 60e3,
	empty_mass	= 25e3,
	payload_coeff	= 0.0,
	loco_section_mass	= 96e3,
	loco_sections_num	= 2,
	delta_eps	= -1.0
}

train_model.railway_coord = train_model.railway_coord*km
train_model.init_velocity = train_model.init_velocity / kmh

--[[vehicle_mass = {}

mass_coeff = 0.02

payload_coeff = 0.0
local payload_mass = 60e3
local empty_mass = 25e3

for i = 1, train_model.vehicles_num do

	if (i == 1) or (i == 2) then
		vehicle_mass[i] = 96e3
	else
		vehicle_mass[i] = empty_mass + payload_coeff*payload_mass
	end

end]]--

---------------------------------------------------------------------
--		Coupling parameters
---------------------------------------------------------------------
coupling_params = 
{
	c_1 = 2.57e7,
	c_2 = 2.85e6,
	c_k = 2.5e8,
	beta = 100,
	T0 = 240e3,
	t0 = 50e3,
	lambda = 0.09,
	delta = 0.05
}

---------------------------------------------------------------------
--	Initial conditions (coupling's gaps)
---------------------------------------------------------------------

--[[local delta = coupling_params.delta
local delta_eps = -1.0

delta_initc = {}

for i = 1, train_model.vehicles_num - 1 do

	delta_initc[i] = delta_eps*delta/2

end]]--

---------------------------------------------------------------------
--		Brakes program
---------------------------------------------------------------------
valve_pos = function(t, v)

	v_pos = 3

	return v_pos

end

---------------------------------------------------------------------
--    Traction program
---------------------------------------------------------------------
traction = function(t, v)
  
  force = 0
  dFdt = 1e4
  Fmax = 300e3
  
  --[[if (math.abs(v) < 78) then
    force = dFdt*t
    
    if (force > Fmax) then
      force = Fmax
    end
    
  else
    
    force = 0
    
  end ]]--
  
  
  return force
  
end

res_file = "v_90_kmh.txt"