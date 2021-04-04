#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

headers = [ 'crossSectionInel.h',
            'lossesBrem.h',
            'lossesHadronics.h',
            'lossesIC.h',
            'lossesPhotoHadronic.h',
            'lossesSyn.h',
            'nonThermalLosses.h' ]

sources = [ 'crossSectionInel.cpp',       
            'lossesBrem.cpp',
            'lossesHadronics.cpp',
            'lossesIC.cpp',
            'lossesPhotoHadronic.cpp',
            'lossesSyn.cpp',
            'nonThermalLosses.cpp' ]

includes = include_directories( ['./', '../'] )


flosses = static_library( 'flosses',
                          sources : sources + headers,
                          include_directories : includes )

flosses_dep = declare_dependency( include_directories : include_directories('.'),
                                      link_with : flosses )
