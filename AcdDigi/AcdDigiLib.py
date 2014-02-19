# $Header
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AcdDigi'])
    env.Tool('AcdUtilLib')
    env.Tool('EventLib')
    env.Tool('OverlayEventLib')
    env.Tool('GlastSvcLib')
    env.Tool('CalibDataLib')
    env.Tool('xmlBaseLib')
    env.Tool('addLibrary', library=env['gaudiLibs'])
    env.Tool('addLibrary', library=env['clhepLibs'])
def exists(env):
    return 1;
