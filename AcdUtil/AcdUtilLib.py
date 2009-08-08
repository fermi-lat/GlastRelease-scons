# $Header$
def generate(env, **kw):
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['AcdUtil'])
    env.Tool('addLibrary', library = ['AcdUtilCommon'])
    env.Tool('CalibDataLib')
    env.Tool('geometryLib')
    env.Tool('xmlBaseLib')
    env.Tool('rdbModelLib')
    env.Tool('calibUtilLib')
    env.Tool('EventLib')
    env.Tool('GlastSvcLib')
    env.Tool('mootCoreLib')
    env.Tool('facilitiesLib')
    env.Tool('addLibrary', library = env['gaudiLibs'])
    env.Tool('addLibrary', library = env['mysqlLibs'])
    env.Tool('addLibrary', library = env['clhepLibs'])
    env.Tool('addLibrary', library = env['cppunitLibs'])
    env.Tool('addLibrary', library = env['xercesLibs'])

def exists(env):
    return 1;
