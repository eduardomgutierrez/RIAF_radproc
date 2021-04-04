#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

headers = [ 'Particle.h' ]
sources = [ 'Particle.cpp' ]

includes = include_directories( ['./', '../'] )

fparticle = static_library( 'inout',
                        sources : sources + headers,
                        #dependencies : _dep,
                        include_directories : includes )

fparticle_dep = declare_dependency( include_directories : './',
                                link_with : fparticle )
