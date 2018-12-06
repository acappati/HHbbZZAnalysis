
#DATA_TAG = "ReReco"   # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2018    # current default = 2018 
ELECORRTYPE = "None"   # "None" to switch off
ELEREGRESSION = "None" # "None" to switch off
APPLYMUCORR = False    # Switch off muon scale corrections
APPLYJEC = False       #
APPLYJER = False       #
RECORRECTMET = False   #
KINREFIT = False       # control KinZFitter (very slow)
PROCESS_CR = False     # Uncomment to run CR paths and trees
ADDLOOSEELE = False    # Run paths for loose electrons
APPLYTRIG = False      # Skip events failing required triggers. They are stored with sel<0 if set to False 
#KEEPLOOSECOMB = True  # Do not skip loose lepton ZZ combinations (for debugging)

PD = ""
MCFILTER = ""

#For DATA: 
#IsMC = False
#PD = ""

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/HHbbZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

process.source.fileNames = cms.untracked.vstring(

### HH SM MC
     '/eos/cms/store/user/covarell/HH/SM/testMINIAOD_HHSM_VV2l_0.root'
    )



process.maxEvents.input = -1
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:8670")

# Debug
# process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
#      dumpTrigger = cms.untracked.bool(True),
#      muonSrcs = cms.PSet(
# #       slimmedMuons = cms.InputTag("slimmedMuons"),
#         muons = cms.InputTag("appendPhotons:muons"),
#      ),
#      electronSrcs = cms.PSet(
# #       slimmedElectron = cms.InputTag("slimmedElectrons"),
#         electrons = cms.InputTag("appendPhotons:electrons"),
# #        RSE = cms.InputTag("appendPhotons:looseElectrons"),
# #        TLE = cms.InputTag("appendPhotons:electronstle"), #These are actually photons, should add a photonSrcs section for them.
#      ),
#      candidateSrcs = cms.PSet(
#         Z     = cms.InputTag("ZCand"),
# #        ZRSE     = cms.InputTag("ZCandlooseEle"),
# #        ZTLE     = cms.InputTag("ZCandtle"),
#         ZZ  = cms.InputTag("ZZCand"),
# #        ZZRSE     = cms.InputTag("ZZCandlooseEle"),
# #        ZZTLE     = cms.InputTag("ZZCandtle"),
# #        ZLL  = cms.InputTag("ZLLCand"),
# #        ZL  = cms.InputTag("ZlCand"),
#      ),
#      jetSrc = cms.InputTag("cleanJets"),
# )

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
