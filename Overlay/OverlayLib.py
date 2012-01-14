#$Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library=['Overlay'])

    env.Tool('RootConvertLib')
    env.Tool('rootUtilLib')
    env.Tool('OverlayEventLib')
    env.Tool('overlayRootDataLib')
    env.Tool('astroLib')
    ##env.Tool('configDataLib')   # doesn't help
    env.Tool('xmlUtilLib')
    env.Tool('AcdUtilLib')
    env.Tool('CalUtilLib')
    ###env.Tool('addLibrary', library = env['geant4Libs']) # neither does this

    env.Tool('addLibrary', library=env['gaudiLibs'])
    #env.Tool('addLibrary', library=env['obfLibs'])

def exists(env):
    return 1;
