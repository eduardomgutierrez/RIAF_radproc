#project( 'fmath',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

headers = [ 'configure.h',
            'fbisection.h',
            'mathematics.h',
            'physics.h',
            'constants.h',
            'matrixInit.h',
            'tridiagonal.h',
            'interpolation.h',
            'mathFunctions.h',
            'RungeKutta.h' ]

sources = [ 'fbisection.cpp',       
            'matrixInit.cpp',
            'tridiagonal.cpp',
            'interpolation.cpp',
            'mathFunctions.cpp',
            'RungeKutta.cpp' ]


fmath = static_library( 'fmath',
                        sources : sources + headers )

fmath_dep = declare_dependency( include_directories : include_directories('.'),
                                link_with : fmath )
