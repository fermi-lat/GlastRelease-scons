# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['detModel'])
        if env['PLATFORM'] == "win32" and env.get('CONTAINERNAME','') == 'GlastRelease':
            env.Tool('findPkgPath', package = 'detModel') 
            env.Tool('findPkgPath', package = 'xmlBase') 
            env.Tool('findPkgPath', package = 'facilities') 

    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('xmlUtilLib')
    env.Tool('xmlBaseLib')
    env.Tool('identsLib')
    env.Tool('addLibrary', library = env['xercesLibs'])

    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'facilities')        
        env.Tool('findPkgPath', package = 'xmlBase')        
        env.Tool('findPkgPath', package = 'xmlUtil')        
def exists(env):
    return 1;
