# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['LdfConverter'])
    env.Tool('ldfReaderLib')	
    env.Tool('GlastSvcLib')
    env.Tool('EventLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])

    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'AdfEvent') 
        env.Tool('findPkgPath', package = 'LdfEvent') 
        env.Tool('findPkgPath', package = 'EbfWriter') 
def exists(env):
    return 1;
