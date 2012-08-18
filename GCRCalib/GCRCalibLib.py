# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['GCRCalib'])
    env.Tool('addLibrary', library = env['rootLibs'])
    env.Tool('addLibrary', library = env['rootGuiLibs'])
    env.Tool('EventLib')
    env.Tool('TkrUtilLib')
    env.Tool('OnboardFilterTdsLib')
    if env['PLATFORM']=='win32' and env.get('CONTAINERNAME','')=='GlastRelease':
        env.Tool('findPkgPath', package = 'GlastSvc') 
        env.Tool('findPkgPath', package = 'lsfData') 
def exists(env):
    return 1;
