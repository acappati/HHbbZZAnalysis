from HHbbZZAnalysis.AnalysisStep.defaults import *
from HHbbZZAnalysis.AnalysisStep.miscenums import *

### ----------------------------------------------------------------------
###
### Example analyzer
###
###----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)
declareDefault("PD", "", globals()) # "" for MC, "DoubleEle", "DoubleMu", or "MuEG" for data
declareDefault("MCFILTER", "", globals())
declareDefault("XSEC", 1, globals())
declareDefault("GENXSEC", 1, globals())
declareDefault("GENBR", 1, globals())
declareDefault("PROCESS_CR", False, globals())
declareDefault("ADDZTREE", False, globals())

# LHE info                 -> FIXME: remove this comment
#  VVDECAYMODE\VVMODE  / ZZ==1 / WW==0  / Yukawa==2 / Zgam=3 / gamgam=4 / Z+nj=5
#                     0: 4l    / lnulnu / 2l        / 2l     / gam      / 2l
#                     1: 4q    / 4q     / 2q        / 2q     / -        / 2q
#                     2: 2l2q  / lnu2q  / -         / -      / -        / -
#                     3: 2l2nu / -      / -         / -      / -        / -
#                     4: 2q2nu / -      / -         / -      / -        / -
#                     5: 4nu   / -      / -         / 2nu    / -        / 2nu
#                    -1: [ Any                                                 ]
#                    -2: [ 2l2X         ]
#                    -3: [ 2nu2X        ]
#                    -4: [ 2q2X         ]
declareDefault("VVMODE", 1, globals())                   # FIXME: remove?
declareDefault("VVDECAYMODE", 0, globals())              # FIXME: remove?
declareDefault("ADDLHEKINEMATICS", False, globals())

# K factors
declareDefault("APPLY_K_NNLOQCD_ZZGG", 0, globals()) # 0: Do not; 1: NNLO/LO; 2: NNLO/NLO; 3: NLO/LO
declareDefault("APPLY_K_NNLOQCD_ZZQQB", False, globals())
declareDefault("APPLY_K_NLOEW_ZZQQB", False, globals())

#failed events
declareDefault("SKIP_EMPTY_EVENTS", True, globals())
declareDefault("FAILED_TREE_LEVEL", 0, globals())

#ggF uncertainties for HTXS
declareDefault("APPLY_QCD_GGF_UNCERT", False, globals() )

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
# process.PlotsbbZZ  = cms.EDAnalyzer("HHbbZZAnalyzer",                                #FIXME
#                                     channel = cms.untracked.string('bbZZ'),
#                                     candCollection = cms.untracked.string('ZZCand'),
#                                     isMC = cms.untracked.bool(IsMC),
#                                     sampleType = cms.int32(SAMPLE_TYPE),
#                                     setup = cms.int32(LEPTON_SETUP),
#                                     skimPaths = cms.vstring(SkimPaths),
#                                     PD = cms.string(PD),
#                                     MCFilterPath = cms.string(MCFILTER),
#                                     sampleName = cms.string(SAMPLENAME),
#                                     dumpForSync = cms.untracked.bool(False),
#                                     )

### Debug 


#Count events with at least 1 Z
# process.ZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
#     src = cms.InputTag("ZCand"),
#     cut = cms.string("userFloat('GoodLeptons')")
# )
# process.sStep4 = cms.EDFilter("CandViewCountFilter",
#                               src = cms.InputTag("ZFiltered"),
#                               minNumber = cms.uint32(1)
#                               )
#process.step4 = cms.Path(process.SkimSequence + process.ZFiltered + process.sStep4 )



### ----------------------------------------------------------------------
### Analyzer for Trees
### ----------------------------------------------------------------------

# TreeSetup = cms.EDAnalyzer("HHbbZZNtupleMaker",                                       #FIXME
#                            channel = cms.untracked.string('aChannel'),
#                            CandCollection = cms.untracked.string('ZZCand'),
#                            fileName = cms.untracked.string('candTree'),
#                            isMC = cms.untracked.bool(IsMC),
#                            sampleType = cms.int32(SAMPLE_TYPE),
#                            setup = cms.int32(LEPTON_SETUP),
#                            skimPaths = cms.vstring(SkimPaths),
#                            PD = cms.string(PD),
#                            MCFilterPath = cms.string(MCFILTER),
#                            metSrc = metTag,
#                            applyTrigger = cms.bool(APPLYTRIG), #Skip events failing required triggers. They are stored with sel<0 if set to false
#                            applyTrigEff = cms.bool(False), #Add trigger efficiency as a weight, for samples where the trigger cannot be applied (obsoltete)
#                            skipEmptyEvents = cms.bool(SKIP_EMPTY_EVENTS),
#                            failedTreeLevel = cms.int32(FAILED_TREE_LEVEL),
#                            sampleName = cms.string(SAMPLENAME),
# 			   GenXSEC = cms.double(GENXSEC),
# 			   GenBR = cms.double(GENBR),
									
#                            # MELA parameters
#                            superMelaMass = cms.double(SUPERMELA_MASS),

#                            # Reco MEs to pick from the candidate
#                            recoProbabilities = cms.vstring(),

#                            # LHE info. parameters
#                            lheProbabilities = cms.vstring(),
#                            xsec = cms.double(XSEC),
#                            VVMode = cms.int32(VVMODE),
#                            VVDecayMode = cms.int32(VVDECAYMODE),
#                            AddLHEKinematics = cms.bool(ADDLHEKINEMATICS),
#                            Apply_K_NNLOQCD_ZZGG = cms.int32(APPLY_K_NNLOQCD_ZZGG),
#                            Apply_K_NNLOQCD_ZZQQB = cms.bool(APPLY_K_NNLOQCD_ZZQQB),
#                            Apply_K_NLOEW_ZZQQB = cms.bool(APPLY_K_NLOEW_ZZQQB),
# 			   Apply_QCD_GGF_UNCERT = cms.bool(APPLY_QCD_GGF_UNCERT),
#                            )

### Signal region
process.bbZZTree = TreeSetup.clone()
process.bbZZTree.channel = 'bbZZ'




# Debug
#Define candidates to be dumped
process.bbZZFiltered = cms.EDFilter("PATCompositeCandidateRefSelector",
                                    src = cms.InputTag("bbZZCand"),
                                    cut = cms.string("userFloat('isBestCand')")
                                   )
### Select only events with one such candidate
process.bbZZSelection= cms.EDFilter("CandViewCountFilter",
                                    src = cms.InputTag("bbZZFiltered"),
                                    minNumber = cms.uint32(1)
                                   )



process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
                                        dumpTrigger = cms.untracked.bool(True),
                                        muonSrcs =  cms.PSet(
                                            muons = cms.InputTag("appendPhotons:muons"),
                                        ),
                                        electronSrcs = cms.PSet(
                                            electrons = cms.InputTag("appendPhotons:electrons"),
                                        ),
                                        candidateSrcs = cms.PSet(
                                            ZZ = cms.InputTag("bbZZCand"),
                                        ),
                                        jetSrc = cms.InputTag("cleanJets"),
)


    
process.trees = cms.EndPath(process.bbZZTree)
#process.plots = cms.EndPath(process.PlotsbbZZ)



