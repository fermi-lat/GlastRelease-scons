# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['calibGenTRG'])
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('commonRootDataLib')
    env.Tool('digiRootDataLib')
    env.Tool('reconRootDataLib')
    env.Tool('mcRootDataLib')
    env.Tool('CalUtilLib')

def exists(env):
    return 1;
