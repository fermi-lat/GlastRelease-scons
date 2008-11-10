# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['HepRepXml'])
    env.Tool('addLibrary', library = env['zlibLibs'])
    env.Tool('HepRepSvcLib')	
def exists(env):
    return 1;
