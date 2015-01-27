---------------------------------------------------------------------
--
--    Test Lua simulation scenario
--    (c) maisvendoo, 2015/01/27
--
---------------------------------------------------------------------
local lfs = require "lfs"

---------------------------------------------------------------------
--    Path variables
---------------------------------------------------------------------
local train_path = "/home/maisvendoo/work/mono-d/train/bin/Release/"
local train_exec = "train"
local train = train_path .. train_exec

local cfg_path = "/home/maisvendoo/work/mono-d/train/cfg/"
local cfg_name = "train.lua"

local result_path = "/home/maisvendoo/work/conddest/"

train = train .. " --config=" .. cfg_path .. cfg_name

---------------------------------------------------------------------
--    Variables
---------------------------------------------------------------------
local output
local result = 0

local vehicles_num = 10
local v0 = 0
local init_velocity = v0

local vehicles_step = 5
local velocity_step = 1.0

--os.execute("mkdir " .. result_path .. results_dir)
local results_dir = "half_vehicles_num_v0/"

lfs.mkdir(result_path .. results_dir);

while (vehicles_num <= 100) do

  local results_file = "vehicles_num-" 
  
  results_file = results_file .. vehicles_num .. ".txt"
    
  init_velocity = v0  
    
  while (init_velocity <= 120) do
    
    local cmd_line = " --vehicles=" .. vehicles_num .. " --init_velocity=" .. init_velocity .. " --payload_coeff=0.5"
    
    print("Start with nv = " .. vehicles_num .. " v0 = " .. init_velocity .. " km/h\n")
    
    output = io.popen(train .. cmd_line)
    
    result = output:read("*a")
    
    output:close()
    
    local out_file = io.open(result_path .. results_dir .. results_file, "a")
    
    io.output(out_file)
    io.write(init_velocity .. " " .. result .. "\n")
    
    io.close(out_file)
    
    init_velocity = init_velocity + velocity_step
    
  end
  
  vehicles_num = vehicles_num + vehicles_step

end

