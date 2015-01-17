#--------------------------------------------------------------------
#		Project globals
#--------------------------------------------------------------------
source_dir = 'src/'
d_include_path = 'src/D/:src/D/lua_d_api/'

release_target_dir = 'bin/release/'
debug_target_dir = 'bin/debug/'

target_name = 'train'

#--------------------------------------------------------------------
#	Release build configuration
#--------------------------------------------------------------------
release_env = Environment(

	CC='gcc',
	CXX='g++',
	DMD='dmd',
	DPATH=d_include_path,
	LINK='gcc',
	
	CPPFLAGS='-O3',
	DFLAGS='-O'
	)

release_env.VariantDir(release_target_dir,
					   source_dir,
					   duplicate=0)

d_sources = Glob(release_target_dir + 'D/*.d')
lua_d_sources = Glob(release_target_dir + 'D/lua_d_api/*.d')

c_sources = Glob(release_target_dir + 'C/*.c')

#c_obj = release_env.Object(c_sources)
d_obj = release_env.Object(d_sources)
lua_d_obj = release_env.Object(lua_d_sources)



release_env.Program(release_target_dir + target_name, d_obj + lua_d_obj, LIBS=['phobos2', 'lua'])


#--------------------------------------------------------------------
#	Debug build configuration
#--------------------------------------------------------------------
debug_env = Environment(

	CC='gcc',
	CXX='g++',
	DMD='dmd',
	DPATH=d_include_path,
	LINK='gcc',
	
	CPPFLAGS='-g3',
	DFLAGS='-g'
	)

debug_env.VariantDir(debug_target_dir,
					 source_dir,
					 duplicate=0)

d_sources = Glob(debug_target_dir + 'D/*.d')
lua_d_sources = Glob(debug_target_dir + 'D/lua_d_api/*.d')

c_sources = Glob(debug_target_dir + 'C/*.c')

#c_obj = debug_env.Object(c_sources)
d_obj = debug_env.Object(d_sources)
lua_d_obj = debug_env.Object(lua_d_sources)

debug_env.Program(debug_target_dir + target_name, d_obj + lua_d_obj, LIBS=['phobos2', 'lua'])
