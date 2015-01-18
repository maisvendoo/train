//-------------------------------------------------------------------
//
//		Module for work with Lua configuration scripts
//		(c) maisvendoo, 2015/01/16
//
//-------------------------------------------------------------------
module	LuaScript;

import	std.stdio;

import	lua;
import	lualib;
import	lauxlib;

enum	int		LUA_S_OK			= 0;
enum	int		LUA_S_NOTTABLE		= 1;
enum	int		LUA_S_NONUMBER		= 2;
enum	int		LUA_S_NOEXIST		= 3;

//-------------------------------------------------------------------
//
//-------------------------------------------------------------------
class	CLuaScript
{
	// Private parameters
	private
	{
		lua_State		*L;		// Lua VM state structure pointer
		int				top;	// Lua stack top
	}

	// Constructor
	this()
	{
		this.L = luaL_newstate();
		luaL_openlibs(this.L);

		this.top = lua_gettop(this.L);
	}

	// Destructor
	~this()
	{
		lua_close(L);
	}

	private
	{
		void restore_stack()
		{
			while (top - lua_gettop(L))
				lua_pop(L, 1);
		}
	}

	public
	{
		// Execute Lua-script
		int exec_script(string filename)
		{
			if (!L)
			{
				writeln("FAIL: Lua VM is't initialized");
				return -1;
			}

			if (luaL_dofile(L, cast(char*) filename))
			{
				string err = lua_tostring(L, -1);
				writeln(err);

				return -1;
			}

			return 0;
		}

		// Get double value by name
		double get_double(string param_name)
		{
			top = lua_gettop(L);

			lua_getglobal(L, cast(char*) param_name);

			double ret = lua_tonumber(L, -1);

			restore_stack();

			return ret;
		}

		// Get integer value by name
		int get_int(string param_name)
		{
			top = lua_gettop(L);

			lua_getglobal(L, cast(char*) param_name);

			int ret = cast(int) lua_tointeger(L, -1);

			restore_stack();

			return ret;
		}

		// Get string value by name
		string get_string(string param_name)
		{
			top = lua_gettop(L);

			lua_getglobal(L, cast(char*) param_name);

			string ret = lua_tostring(L, -1);

			restore_stack();

			return ret;
		}

		// Check is variable exist
		bool is_exist(string var)
		{
			bool ret = false;

			lua_getglobal(L, cast(char*) var);

			if (lua_isnil(L, -1))
				ret = false;
			else
				ret = true;

			restore_stack();

			return true;
		}

		// Check table field is exist
		bool is_field_exist(string table, string field)
		{
			bool ret = false;

			top = lua_gettop(L);

			if (lua_istable(L, -1))
			{
				lua_getfield(L, -1, cast(char*) field);

				if (lua_isnil(L, -1))
					ret = false;
				else
					ret = true;
			}
			else
			{
				ret = false;
			}

			restore_stack();

			return ret;
		}

		// Get table double field
		double get_double_field(string table, string field, ref int err)
		{
			double ret = 0;

			top = lua_gettop(L);

			lua_getglobal(L, cast(char*) table);

			if (lua_istable(L, -1))
			{
				if (is_field_exist(table, field))
					lua_getfield(L, -1, cast(char*) field);
				else
				{
					err = LUA_S_NOEXIST;
					restore_stack();
					return ret;
				}

				if (lua_isnumber(L, -1))
				{
					ret = cast(double) lua_tonumber(L, -1);
					err = LUA_S_OK;
				}
				else
				{
					err = LUA_S_NONUMBER;
					ret = 0;
				}
			}
			else
			{
				err = LUA_S_NOTTABLE;
				ret = 0;
			}

			restore_stack();

			return ret;
		}

		// Get double field by index
		double get_double_field(string table, int idx, int err)
		{
			double ret = 0;

			top = lua_gettop(L);

			lua_getglobal(L, cast(char*) table);

			if (lua_istable(L, -1))
			{
				int index = lua_gettop(L);
				int count = 0;
				lua_pushnil(L);

				while ( (lua_next(L, index) != 0) && (count <= idx) )
				{
					ret = lua_tonumber(L, -1);
					lua_pop(L, 1);
					count++;
				}

				err = LUA_S_OK;
			}
			else
			{
				err = LUA_S_NOTTABLE;
				ret = 0;
			}

			restore_stack();

			return ret;
		}

		// Get table integer field
		int get_int_field(string table, string field, ref int err)
		{
			int ret = 0;
			
			top = lua_gettop(L);
			
			lua_getglobal(L, cast(char*) table);
			
			if (lua_istable(L, -1))
			{
				if (is_field_exist(table, field))
					lua_getfield(L, -1, cast(char*) field);
				else
				{
					err = LUA_S_NOEXIST;
					restore_stack();
					return ret;
				}
				
				if (lua_isnumber(L, -1))
				{
					ret = cast(int) lua_tonumber(L, -1);
					err = LUA_S_OK;
				}
				else
				{
					err = LUA_S_NONUMBER;
					ret = 0;
				}
			}
			else
			{
				err = LUA_S_NOTTABLE;
				ret = 0;
			}
			
			restore_stack();
			
			return ret;
		}

		// Get table string field
		string get_str_field(string table, string field, ref int err)
		{
			string ret = "";
			
			top = lua_gettop(L);
			
			lua_getglobal(L, cast(char*) table);
			
			if (lua_istable(L, -1))
			{
				if (is_field_exist(table, field))
					lua_getfield(L, -1, cast(char*) field);
				else
				{
					err = LUA_S_NOEXIST;
					restore_stack();
					return ret;
				}
				
				if (lua_isstring(L, -1))
				{
					ret = lua_tostring(L, -1);
					err = LUA_S_OK;
				}
				else
				{
					err = LUA_S_NONUMBER;
					ret = "";
				}
			}
			else
			{
				err = LUA_S_NOTTABLE;
				ret = "";
			}
			
			restore_stack();
			
			return ret;
		}

		// Call Lua function
		double call_func(string func, double[] args, ref int err)
		{
			int nargs = cast(int) args.length;

			top = lua_gettop(L);

			lua_getglobal(L, cast(char*) func);

			for (int i = 0; i < nargs; i++)
			{
				lua_pushnumber(L, args[i]);
			}

			err = lua_pcall(L, nargs, 1, 0);

			if (err)
			{
				writeln("FAIL: Lua run-time error");
				return 0;
			}

			double ret = lua_tonumber(L, -1);

			restore_stack();

			return ret;
		}

	
	}
}