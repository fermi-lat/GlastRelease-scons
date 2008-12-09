# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrRecon'])
    env.Tool('GuiSvcLib')
    env.Tool('CalUtilLib')
    env.Tool('CalibDataLib')
    env.Tool('GlastSvcLib')
    env.Tool('guiLib')
    env.Tool('EventLib')
    env.Tool('TkrUtilLib')
    env.Tool('LdfEventLib')
    env.Tool('geometryLib')
    env.Tool('RootIoLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
def exists(env):
    return 1;
