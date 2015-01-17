---------------------------------------------------------------------
--
--	Train simulation config 
--
---------------------------------------------------------------------

solver_params = 
{
	method = "rkf5",	-- integration method
	init_time = 0,		-- initial time
	stop_time = 2.0,	-- stop simulation time
	step = 1e-4,		-- time step
	max_step = 0.1,		-- maximal time step
	local_err = 1e-8	-- loacal solver error
}