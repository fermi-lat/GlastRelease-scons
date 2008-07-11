# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['GuiSvc'])
    env.Tool('guiLib')

def exists(env):
    return 1;
