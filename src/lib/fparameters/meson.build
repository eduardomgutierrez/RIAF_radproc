#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

cc = meson.get_compiler('cpp')
#fmath_dep = cc.find_library( 'fmath',
#                             dirs : '/home/emgutierrez/codes/new_radproc/src/lib/fmath/build/' )

headers = [ 'Dimension.h',
            'DimensionCoord.h',
            'DimensionIterator.h',
            'parameters.h',
            'ParamSpace.h',
            'ParamSpaceValues.h',
            'SpaceCoord.h',
            'SpaceIterator.h' ]

sources = [ 'Dimension.cpp',       
            'DimensionIterator.cpp',
            'parameters.cpp',
            'ParamSpace.cpp',
            'ParamSpaceValues.cpp',
            'SpaceCoord.cpp',
            'SpaceIterator.cpp' ]


fparameters = static_library( 'fparameters',
                              sources : sources + headers,
                              dependencies : fmath_dep,
                              include_directories : include_directories('../') )

fparameters_dep = declare_dependency( include_directories : include_directories('.'),
                                      link_with : fparameters )
