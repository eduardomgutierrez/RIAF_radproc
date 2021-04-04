#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

headers = [ 'bessel.h',
            'bisection.h',
            'brent.h',
            'integrators.h',
            'laguerre.h',
            'mathematics.h',
            'nr.h',
            'nrutil.h',
            'random.h',
            'ran1.h' ]

sources = [ 'bessel.cpp',
            'bisection.cpp',
            'brent.cpp',
            'interpol.cpp',
            'laguerre.cpp',
            'broydn.cpp',
            'erff.cpp',
            'integrators.cpp',
            'lnsrch.cpp',
            'lubksb.cpp',
            'ludcmp.cpp',
            'newt.cpp',
            'nrutil.cpp',
            'ran1.cpp',
            'zbrac.cpp',
            'zbrent.cpp' ]


includes = include_directories( ['./', '../'] )

nrMath = static_library( 'nrMath',
                         sources : sources + headers,
                         dependencies : gsl_dep,
                         include_directories : includes )

nrMath_dep = declare_dependency( include_directories : include_directories('.'),
                                 link_with : nrMath )
