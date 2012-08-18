# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AncillaryDataUtil'])
        if env['PLATFORM'] == "win32" and env.get('CONTAINERNAME','') == 'GlastRelease':
            env.Tool('findPkgPath', package = 'AncillaryDataUtil') 
            env.Tool('findPkgPath', package = 'AncillaryDataEvent') 
            env.Tool('findPkgPath', package = 'AdfEvent') 
            env.Tool('findPkgPath', package = 'facilities') 
    env.Tool('facilitiesLib')
    if kw.get('incsOnly', 0) == 1: 
        env.Tool('findPkgPath', package = 'AncillaryDataEvent') 
        env.Tool('findPkgPath', package = 'AdfEvent') 
        env.Tool('findPkgPath', package = 'facilities') 
def exists(env):
    return 1;
