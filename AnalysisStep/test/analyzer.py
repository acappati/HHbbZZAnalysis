from HHbbZZAnalysis.AnalysisStep.defaults import *
from HHbbZZAnalysis.AnalysisStep.miscenums import *

### ----------------------------------------------------------------------
###
### Example analyzer
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # for MC,  for data
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())
declareDefault("GENXSEC", 1, globals())
declareDefault("GENBR", 1, globals())
declareDefault("PROCESS_CR", False, globals())


# K factors
# declareDefault("APPLY_K_NNLOQCD_ZZGG", 0, globals()) # 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
# declareDefault("APPLY_K_NNLOQCD_ZZQQB", False, globals())
# declareDefault("APPLY_K_NLOEW_ZZQQB", False, globals())

#failed events
declareDefault("SKIP_EMPTY_EVENTS", True, globals())
declareDefault("FAILED_TREE_LEVEL", 0, globals())

#ggF uncertainties for HTXS
# declareDefault("APPLY_QCD_GGF_UNCERT", False, globals() )

if FAILED_TREE_LEVEL and not SKIP_EMPTY_EVENTS:
    raise ValueError(
                     "Inconsistent options: FAILED_TREE_LEVEL={}, SKIP_EMPTY_EVENTS={}\n"
                     "If you want to write a failed tree, set SKIP_EMPTY_EVENTS=True"
                     .format(FAILED_TREE_LEVEL, SKIP_EMPTY_EVENTS)
                    )

# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/HHbbZZAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "MasterPy/HHbbZZAnalysis.py")


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------
process.maxEvents.input = -1
#process.options.wantSummary = False


### ----------------------------------------------------------------------
### Output root file (monitoring histograms)
### ----------------------------------------------------------------------
process.TFileService=cms.Service('TFileService',
                                fileName=cms.string('HHbbZZAnalysis.root')
                                )

### ----------------------------------------------------------------------
### Analyzers for Plots
### ----------------------------------------------------------------------

# All events together
process.PlotsHH    = cms.EDAnalyzer("HHbbZZAnalyzer",
                                    # channel = cms.untracked.string('ZZ'),
                                    # candCollection = cms.untracked.string('ZZCand'),
                                    # isMC = cms.untracked.bool(IsMC),
                                    # sampleType = cms.int32(SAMPLE_TYPE),
                                    # setup = cms.int32(LEPTON_SETUP),
                                    # skimPaths = cms.vstring(SkimPaths),
                                    # PD = cms.string(PD),
                                    # MCFilterPath = cms.string(MCFILTER),
                                    # sampleName = cms.string(SAMPLENAME),
                                    # dumpForSync = cms.untracked.bool(False),
                                    )




### ----------------------------------------------------------------------
### Analyzer for Trees
### ----------------------------------------------------------------------

TreeSetup = cms.EDAnalyzer("HHbbZZNtupleMaker",
                           # channel = cms.untracked.string('aChannel'),
                           # CandCollection = cms.untracked.string('ZZCand'),
                           # fileName = cms.untracked.string('candTree'),
                           # isMC = cms.untracked.bool(IsMC),
                           # sampleType = cms.int32(SAMPLE_TYPE),
                           # setup = cms.int32(LEPTON_SETUP),
                           # skimPaths = cms.vstring(SkimPaths),
                           # PD = cms.string(PD),
                           # MCFilterPath = cms.string(MCFILTER),
                           # metSrc = metTag,
                           # applyTrigger = cms.bool(APPLYTRIG), #Skip events failing required triggers. They are stored with sel<0 if set to false
                           # applyTrigEff = cms.bool(False), #Add trigger efficiency as a weight, for samples where the trigger cannot be applied (obsoltete)
                           # skipEmptyEvents = cms.bool(SKIP_EMPTY_EVENTS),
                           # failedTreeLevel = cms.int32(FAILED_TREE_LEVEL),
                           # sampleName = cms.string(SAMPLENAME),
			   #      					GenXSEC = cms.double(GENXSEC),
			   #      					GenBR = cms.double(GENBR),
									
                           # # MELA parameters
                           # superMelaMass = cms.double(SUPERMELA_MASS),

                           # # Reco MEs to pick from the candidate
                           # recoProbabilities = cms.vstring(),

                           # # LHE info. parameters
                           # lheProbabilities = cms.vstring(),
                           # xsec = cms.double(XSEC),
                           # VVMode = cms.int32(VVMODE),
                           # VVDecayMode = cms.int32(VVDECAYMODE),
                           # AddLHEKinematics = cms.bool(ADDLHEKINEMATICS),
                           # Apply_K_NNLOQCD_ZZGG = cms.int32(APPLY_K_NNLOQCD_ZZGG),
                           # Apply_K_NNLOQCD_ZZQQB = cms.bool(APPLY_K_NNLOQCD_ZZQQB),
                           # Apply_K_NLOEW_ZZQQB = cms.bool(APPLY_K_NLOEW_ZZQQB),
			   #      					Apply_QCD_GGF_UNCERT = cms.bool(APPLY_QCD_GGF_UNCERT),
                           )

### Signal region
process.HHTree = TreeSetup.clone()
process.HHTree.channel = 'HH'

### Trees for control regions




# Debug
#Define candidates to be dumped
