#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

headers = [ 'ioutil.h' ]
sources = [ 'ioutil.cpp' ]

includes = include_directories( ['./', '../'] )

inout = static_library( 'inout',
                        sources : sources + headers,
                        #dependencies : _dep,
                        include_directories : includes )

inout_dep = declare_dependency( include_directories : './',
                                link_with : inout )
