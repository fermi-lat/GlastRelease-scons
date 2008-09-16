# $Header$
def generate(env, **kw):
     if not kw.get('depsOnly', 0):
          env.Tool('addLibrary', library = ['AnalysisNtuple'])
     env.Tool('TkrUtilLib')
     env.Tool('ntupleWriterSvcLib')
     env.Tool('OnboardFilterLib')
     env.Tool('GlastSvcLib')
     env.Tool('FluxSvcLib')
     env.Tool('G4PropagatorLib')
def exists(env):
    return 1;
