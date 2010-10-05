import os
import SCons.Action
import SCons.Builder
from   fermidebug import fdebug

### template for .sh script
shellScript = '''#!/bin/sh
# Autogenerated by SCons; do not edit!

${REPLACE-WRAPPER-SCRIPT}

export BASE_DIR=${REPLACE-BASEDIR}

export LD_LIBRARY_PATH=${REPLACE-LIBDIRS}:$LD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=${REPLACE-LIBDIRS}:$DYLD_LIBRARY_PATH
export PATH=${REPLACE-PATHS}:$PATH

${REPLACE-PFILES_SETUP}

export PYTHONPATH=${REPLACE-PYTHONPATHS}:$PYTHONPATH
export ROOTSYS=${REPLACE-ROOTSYS}

${REPLACE-WRAPPER-EXECUTE}
'''
# end template for .sh script

### template for .csh script
cshellScript = '''#!/bin/csh
# Autogenerated by SCons; do not edit!

${REPLACE-WRAPPER-SCRIPT}

setenv BASE_DIR ${REPLACE-BASEDIR}

if ($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH ${REPLACE-LIBDIRS}:$LD_LIBRARY_PATH
else
  setenv LD_LIBRARY_PATH ${REPLACE-LIBDIRS}
endif

if ($?DYLD_LIBRARY_PATH) then
  setenv DYLD_LIBRARY_PATH ${REPLACE-LIBDIRS}:$DYLD_LIBRARY_PATH
else
  setenv DYLD_LIBRARY_PATH ${REPLACE-LIBDIRS}
endif

if ($?PATH) then
  setenv PATH ${REPLACE-PATHS}:$PATH
else
  setenv PATH ${REPLACE-PATHS}
endif

${REPLACE-PFILES_SETUP}

if ($?PYTHONPATH) then
  setenv PYTHONPATH ${REPLACE-PYTHONPATHS}:${PYTHONPATH}
else
  setenv PYTHONPATH ${REPLACE-PYTHONPATHS}
endif

setenv ROOTSYS ${REPLACE-ROOTSYS}

${REPLACE-WRAPPER-EXECUTE}
'''
# end template .csh script

### template for bat (Windows)
batScript = '''//Autogenerated by SCons; do not edit!
${REPLACE-WRAPPER-SCRIPT}

set BASE_DIR=${REPLACE-BASEDIR}

set PATH=${REPLACE-PATHS};%PATH%
set PATH=%PATH%;${REPLACE-LIBDIRS}

set ROOTSYS=${REPLACE-ROOTSYS}
set PYTHONPATH=${REPLACE-PYTHONPATHS}:%PYTHONPATH%

${REPLACE-PFILES_SETUP}

${REPLACE-WRAPPER-EXECUTE}
'''

# Helper functions to determine relative paths
# Obtained from http://code.activestate.com/recipes/208993/
def commonpath(l1, l2, common=[]):
    if len(l1) < 1: return (common, l1, l2)
    if len(l2) < 1: return (common, l1, l2)
    if l1[0] != l2[0]: return (common, l1, l2)
    return commonpath(l1[1:], l2[1:], common+[l1[0]])

def relpath(p1, p2):


    (common,l1,l2) = commonpath(p1.split(os.path.sep), p2.split(os.path.sep))
    p = []
    if len(l1) > 0:
        p = [ '../' * len(l1) ]
    p = p + l2
    return os.path.join( *p )

## Another helper to restore /nfs/slac... type paths
def resolve_nfs_path(path):
    tokens = path.split(":")
    for i in range(len(tokens)):
        if tokens[i].find('g.glast.'):
            tokens[i] = os.path.join('/nfs/farm/g/glast',
                                     tokens[i].split('g.glast.')[-1])
    return ":".join(tokens)

## Fill contents of wrapper scripts and setup script for an SCons installation
def fillScript(scriptFile, env, wrapper, script, executable):
    finalScript = script.get_contents()

    inst = ''
    if env['PLATFORM'] == 'win32':
        separator= ';'
        bs = '%BASE_DIR%'
        inst = '%INST_DIR%'
    else:
        separator = ':'
        bs='$BASE_DIR'
        inst='$INST_DIR'

    scriptdir = os.path.join(inst, 'bin', env['VARIANT'])

    if env.GetOption('supersede') != '.':
	basedirAbs = env.Dir('.').abspath
	if env['PLATFORM'] == "posix":  # might be nfs path
            basedirAbs = resolve_nfs_path(basedirAbs)
	finalScript = finalScript.replace('${REPLACE-BASEDIR}', '"' + basedirAbs+ '"')
    else:
        print "inst is ", inst
        finalScript = finalScript.replace('${REPLACE-BASEDIR}', inst)
        
    # Handle pfiles setup
    if 'usePfiles' in env:
        if env['PLATFORM'] == 'win32':
            finalScript = finalScript.replace('${REPLACE-PFILES_SETUP}', 
                                              'call ${REPLACE-SCRIPTDIR}\\pfiles_setup.bat')
        else:
            filename = scriptFile.name
            if wrapper > 0:
                ext = 'sh'
            else:
                ext = filename[filename.rfind('.') + 1:]
            finalScript = finalScript.replace('${REPLACE-PFILES_SETUP}', 
                                              'source ${REPLACE-SCRIPTDIR}/pfiles_setup.' + ext)
    else:
        ##print "usePfiles not in env"
        finalScript = finalScript.replace('${REPLACE-PFILES_SETUP}', '')

    # Handle SCRIPTDIR references
    finalScript = finalScript.replace('${REPLACE-SCRIPTDIR}', scriptdir)

    #Set up LD_LIBRARY_PATH and DYLD_LIBRARY_PATH
    ldLibraryPath = [ os.path.join(inst, relpath(env.Dir(env.GetOption('supersede')).abspath, env['LIBDIR'].abspath))]
    ldLibraryPath.append(os.path.join(bs, 'lib', env['VARIANT']))

    ldLibraryPath.extend(env['WRAPPERLIBS'])

    finalScript = finalScript.replace('${REPLACE-LIBDIRS}', separator.join(ldLibraryPath))
    #Setup PATH
    path = [os.path.join(inst, relpath(env.Dir(env.GetOption('supersede')).abspath, env['BINDIR'].abspath))]
    path.append(os.path.join(bs, 'exe', env['VARIANT']))
                
    path.extend(env['WRAPPERBINS'])

    finalScript = finalScript.replace('${REPLACE-PATHS}', separator.join(path))

    #Setup ROOTSYS
    rootSys = env['ROOTSYS']
    finalScript = finalScript.replace('${REPLACE-ROOTSYS}', rootSys)


    #Setup PYTHONPATH
    pythonPath = [os.path.join(inst,'python')]
    pythonPath.append(os.path.join(inst, relpath(env.Dir(env.GetOption('supersede')).abspath, env['LIBDIR'].abspath)))
    pythonPath.append(os.path.join(bs, 'python'))
    pythonPath.append(os.path.join(bs, 'lib', env['VARIANT']))
    pythonPath.append(os.path.join(env['ROOTSYS'], 'lib'))
    finalScript = finalScript.replace('${REPLACE-PYTHONPATHS}', separator.join(pythonPath))
    
    if env['PLATFORM'] == 'win32':
        finalScript = finalScript.replace('$GLAST_EXT', '%GLAST_EXT%')

    if wrapper > 0:
        if env['PLATFORM'] != 'win32':
            finalScript = finalScript.replace('${REPLACE-WRAPPER-SCRIPT}', 'export INST_DIR=`dirname $0`\nexport INST_DIR=`cd $INST_DIR/../../; pwd`\n')
            finalScript = finalScript.replace('${REPLACE-WRAPPER-EXECUTE}', os.path.join('$INST_DIR', relpath(env.Dir(env.GetOption('supersede')).abspath, executable.abspath)+' "$@"\n'))
        else:
            finalScript = finalScript.replace('${REPLACE-WRAPPER-SCRIPT}', 
                                              'set INST_DIR=%~p0\\..\\..')
            finalScript = finalScript.replace('${REPLACE-WRAPPER-EXECUTE}', os.path.join('%INST_DIR%', relpath(env.Dir(env.GetOption('supersede')).abspath, executable.abspath)+' "%*" \n'))
    else:
        if env['PLATFORM'] == 'win32':
            finalScript = finalScript.replace('${REPLACE-WRAPPER-EXECUTE}', '')
            finalScript = finalScript.replace('${REPLACE-WRAPPER-SCRIPT}', 'if not defined INST_DIR (\necho First set environment variable INST_DIR\n goto EOF:\n)\n')
        else:
            finalScript = finalScript.replace('${REPLACE-WRAPPER-EXECUTE}', '')
            finalScript = finalScript.replace('${REPLACE-WRAPPER-SCRIPT}', '')
    scriptFile.write(finalScript)

def generate(env):
    def createWrapper(target = None, source = None, env = None):
        scriptContents = source.pop(0)
        for trgt in target:
            scriptFile = open(str(trgt), 'w')
            fillScript(scriptFile, env, 1, scriptContents, source.pop(0))
            scriptFile.close()
        return 0

    def createWrapperEmitter(target, source, env):
        target = []
        wrapperSuffix = ''
        if env['PLATFORM'] == 'win32':
            wrapperSuffix = '.bat'
        for src in source:
            target.append(env['SCRIPTDIR'].File(os.path.splitext(os.path.basename(src.abspath))[0] + wrapperSuffix))
        if env["PLATFORM"] != "win32":
            source = [ env.Value(shellScript) ] + source
        else:
            source = [ env.Value(batScript) ] + source
        return (target, source)

    def createWrapperGenerator(source, target, env, for_signature):
        actions = [env.Action(createWrapper, "Creating wrapper for '$TARGETS'")]
        for trgt in target:
            if 'setupTarget' in env:
                env.Depends(trgt, env['setupTarget'])
            actions.append(SCons.Defaults.Chmod(trgt, 0755))
        return actions



    def createSetup(target = None, source = None, env = None):
        for trgt in target:
            scriptFile = open(str(trgt), 'w')
            fillScript(scriptFile, env, 0, source.pop(0), '')
            scriptFile.close()
        return 0

    def createSetupEmitter(target, source, env):
        if env['PLATFORM'] == 'win32':
            fdebug("Generating windows setup scripts")
            source = [env.Value(batScript)]
        else:
            target.append(str(target[0]).replace('.sh', '') + '.csh')
            source = [env.Value(shellScript), env.Value(cshellScript)]
        return (target, source)

    def createSetupGenerator(source, target, env, for_signature):
        actions = [env.Action(createSetup, "Creating setup scripts")]
        for trgt in target:
            actions.append(SCons.Defaults.Chmod(trgt, 0755))
        return actions

    if env['PLATFORM'] == 'win32': suf = '.bat'
    else: suf = '.sh'
    env.Append(BUILDERS = {'GenerateWrapperScript' : env.Builder(generator = createWrapperGenerator,
                                                                 emitter = createWrapperEmitter, suffix=suf)})

    env.Append(BUILDERS = {'GenerateSetupScript' : env.Builder(generator = createSetupGenerator,
                                                               emitter = createSetupEmitter,
                                                               suffix=suf)})

def exists(env):
    return 1;

