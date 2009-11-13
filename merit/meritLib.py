# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['merit'])

    env.Tool('EventLib')
    env.Tool('LdfEventLib')
    env.Tool('facilitiesLib')
    env.Tool('AnalysisNtupleLib')
    env.Tool('GuiSvcLib')
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['cppunitLibs'])
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('addLibrary', library = env['cfitsioLibs'])
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootLibs'])

def exists(env):
    return 1;
