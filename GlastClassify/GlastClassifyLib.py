# $Header$
def generate(env, **kw):
#    if not kw.get('noGaudi', 0):
#        env.Tool('addLibrary', library = env['gaudiLibs'])
    if not kw.get('depsOnly', 0):
        env.Tool('addLibrary', library = ['GlastClassify'])
    env.Tool('classifierLib')
    env.Tool('ntupleWriterSvcLib')
    env.Tool('xmlBaseLib')
    env.Tool('facilitiesLib')
    env.Tool('addLibrary', library = env['xercesLibs'])
    env.Tool('addLibrary', library = env['TMineLibs'])
def exists(env):
    return 1;
