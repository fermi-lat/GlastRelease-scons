# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['HepRepCorba'])
    env.Tool('HepRepSvcLib')	
    if env['PLATFORM'] == 'win32':
	env.Tool('addLibrary', library = ['oldnames'])
    
def exists(env):
    return 1;
