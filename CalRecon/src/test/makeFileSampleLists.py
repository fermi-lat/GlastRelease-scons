import ROOT
import os

MC_FILE_PATH = "/data/users/pesce/v18r8p5/work/output"

DATA_FILE_PATH = "/data47b/glast/flight/CalibSets"

TRAIN_FILE_PATH_DICT = {}
CLASS_FILE_PATH_DICT = {}

MC_REF_DICT = {'gam':'all_gamma_180GeV',
	       'had':'CrProtonMix',
	       'mip':'high_e_surface_muons'}
#	       'ghost':'PTSkim'}



DATA_REF_DICT = {'gam': 'EarthLimb_LEO200_Diffuse_Skim',
		 'had': 'HighLatBkg_Skim',
		 'mip': 'Mips_Skim_Dgn'}
#		 'ghost': 'PTSkim'}


def setupTrainDict():
     #   TrainDict = {}
        for (topology,fileBaseName) in REF_DICT.items():
            #FILE_LIST = []
           # fileName = os.path.basename(filePath)
           # fileName = fileName.split('-')[0]
            print topology,fileBaseName
            (CLASS_FILE_LIST,TRAIN_FILE_LIST) = getTrainFilePath(fileBaseName)
            print topology,len(CLASS_FILE_LIST)
            TRAIN_FILE_PATH_DICT[topology] = TRAIN_FILE_LIST
            CLASS_FILE_PATH_DICT[topology] = CLASS_FILE_LIST
        return TRAIN_FILE_PATH_DICT,CLASS_FILE_PATH_DICT


def getTrainFilePath(fileName):
    #fileName = os.path.basename(filePath)
    #fileName = fileName.split('-')[0]
    #print fileName
    list = []
    TrainFileList = []
    ClassFileList = []  
    for file in os.listdir(MAIN_FILE_PATH):
        if fileName in file:
            if 'root' in file:
                list.append(file)
    numFiles = len(list)
    
    for i,trainfile in enumerate(list):
        if i<=numFiles/2.:
            fullPath = os.path.join(MAIN_FILE_PATH,trainfile)
            TrainFileList.append(fullPath)
        else:
            fullPath = os.path.join(MAIN_FILE_PATH,trainfile)
            ClassFileList.append(fullPath)
        i+1
        
   # print "TrainFileList is",trainFileList
#    print len(trainFileList)
    #raw_input()
    return TrainFileList,ClassFileList


def getGRversion():
        GRversion = MAIN_FILE_PATH.split('/')[4]
        GRversion = GRversion.split('/')[0]
        return GRversion


def getDataSampleList():

	grVersion = 'v18r8p5'
	for (topology,filepath) in DATA_REF_DICT.items():
		ClassFileList = []
		fullfilePath = "%s/%s/%s/merit/"%\
			       (DATA_FILE_PATH,filepath,grVersion)
		for filename in os.listdir(fullfilePath):
			if 'merit' in filename:
				fullpath = os.path.join(fullfilePath,filename)
				ClassFileList.append(fullpath)

		CLASS_FILE_PATH_DICT[topology] = ClassFileList
		print "Data files:\n %s (tot %s)"%(topology,len(ClassFileList))
	return CLASS_FILE_PATH_DICT


def getMCSampleList():
	for (topology,basename) in MC_REF_DICT.items():
		TrainFileList = []

		for filename in os.listdir(MC_FILE_PATH):
			if basename in filename:
				fullpath = os.path.join(MC_FILE_PATH,filename)
				TrainFileList.append(fullpath)

		TRAIN_FILE_PATH_DICT[topology] = TrainFileList
		print "Mc files:\n %s (tot %s)"%(topology,len(TrainFileList))
	return TRAIN_FILE_PATH_DICT



#filePath = '/data47b/glast/flight/CalibSets/HighLatBkg_Skim/v18r8p5/merit'
label = 'ClassData'
#topology = 'had'
#GR_VERSION = getGRversion()

CLASS_FILE_PATH_DICT = getDataSampleList()
TRAIN_FILE_PATH_DICT = getMCSampleList()
outputFile = file("FileSampleDicts_%s.py"%label,"w")
#outputFile.writelines("GR_VERSION = '%s'\n"%GR_VERSION)
outputFile.writelines("TRAIN_FILE_PATH_DICT = %s\n"%TRAIN_FILE_PATH_DICT)
outputFile.writelines("CLASS_FILE_PATH_DICT = %s\n"%CLASS_FILE_PATH_DICT)
outputFile.close()



"""
(TRAIN_FILE_PATH_DICT,CLASS_FILE_PATH_DICT) = setupTrainDict()
GR_VERSION = getGRversion()

outputFile = file("FileSampleDicts_NoGhosts.py","w")
outputFile.writelines("GR_VERSION = %s\n"%GR_VERSION)
outputFile.writelines("TRAIN_FILE_PATH_DICT = %s\n"%TRAIN_FILE_PATH_DICT)
outputFile.writelines("CLASS_FILE_PATH_DICT = %s\n"%CLASS_FILE_PATH_DICT)
outputFile.close()
"""
#print TRAIN_FILE_PATH_DICT
#print
#print 
#print CLASS_FILE_PATH_DICT
