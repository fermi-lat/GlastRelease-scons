# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['TkrDigi'])
    env.Tool('facilitiesLib')
    env.Tool('EventLib')
    env.Tool('GuiSvcLib')
    env.Tool('TkrUtilLib')
def exists(env):
    return 1;
