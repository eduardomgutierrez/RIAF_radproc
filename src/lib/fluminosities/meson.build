#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

headers = [ 'blackBody.h',
            'compton.h',
            'luminosityHadronic.h',
            'luminosityIC.h',
            'luminosityNTHadronic.h',
            'luminosityPairAnnihilation.h',
            'luminosityPhotoHadronic.h',
            'luminositySynchrotron.h',
            'opticalDepthSSA.h',
            'probaprox.h',
            'probexact.h',
            'reflection.h',
            'thermalBremss.h',
            'thermalSync.h' ]

sources = [ 'blackBody.cpp',
            'compton.cpp',
            'luminosityHadronic.cpp',
            'luminosityIC.cpp',
            'luminosityNTHadronic.cpp',
            'luminosityPairAnnihilation.cpp',
            'luminosityPhotoHadronic.cpp',
            'luminositySynchrotron.cpp',
            'opticalDepthSSA.cpp', 
            'probaprox.cpp', 
            'probexact.cpp', 
            'reflection.cpp', 
            'thermalBremss.cpp', 
            'thermalSync.cpp' ]

includes = include_directories( ['./', '../'] )

fluminosities = static_library( 'fluminosities',
                                sources : sources + headers,
                                dependencies : fparticle_dep,
                                include_directories : includes )

fluminosities_dep = declare_dependency( include_directories : include_directories('.'),
                                        link_with : fluminosities )
