#project( 'fparameters',
#         'cpp',
#         default_options : ['cpp_std=c++17'],
#         version : '0.1' )

boost_dep = dependency('boost')

headers = [ 'muonInj.h',
            'neutrinoInj.h',
            'neutronDecay.h',
            'neutronPp.h',
            'pairBH.h',
            'pairMuonDecay.h',
            'ppPionInj.h',
            'pairgg.h',
            'pgammaPionInj.h',
            'neutronPgamma.h',
            'pairInjectionExact.h' ]

sources = [ 'muonInj.cpp',
            'neutrinoInj.cpp',
            'neutronDecay.cpp',
            'neutronPp.cpp',
            'pairBH.cpp',
            'pairMuonDecay.cpp',
            'ppPionInj.cpp', 
            'pairgg.cpp', 
            'pgammaPionInj.cpp', 
            'neutronPgamma.cpp', 
            'pairInjectionExact.cpp' ]

includes = include_directories( ['./', '../'] )

finjection = static_library( 'finjection',
                             sources : sources + headers,
                             include_directories : includes,
                             dependencies : boost_dep )

finjection_dep = declare_dependency( include_directories : include_directories('.'),
                                     link_with : finjection )
