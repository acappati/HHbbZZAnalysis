import FWCore.ParameterSet.Config as cms
from HHbbZZAnalysis.AnalysisStep.defaults import *
import os, sys

process = cms.Process("bbZZ")

### ----------------------------------------------------------------------
### Flags that need to be set
### ----------------------------------------------------------------------

# Set defaults for variables used in this file (in case they are not defined by a caller script)

declareDefault("IsMC", True, globals())

# Set of effective areas, rho corrections, etc. (can be 2011, 2012, 2015 or 2016)
declareDefault("LEPTON_SETUP", 2016, globals())

# Flag that reflects the actual sqrts of the sample (can be 2011, 2012, 2015 or 2016)
# Can differ from SAMPLE_TYPE for samples that are rescaled to a different sqrts.
declareDefault("SAMPLE_TYPE", LEPTON_SETUP, globals())

# Control global tag to be used (for data or rereco)
declareDefault("DATA_TAG", "ReReco", globals())

#Optional name of the sample/dataset being analyzed
declareDefault("SAMPLENAME", "", globals())

#Type of electron scale correction/smearing: "None", "RunII"
declareDefault("ELECORRTYPE", "RunII", globals())

#Apply electron escale regression: "None", "Moriond17v1"
declareDefault("ELEREGRESSION", "Moriond17v1", globals())

#Apply muon scale correction
declareDefault("APPLYMUCORR", True, globals())

#Reapply JEC
declareDefault("APPLYJEC", True, globals())

#Apply JER
declareDefault("APPLYJER", True, globals())

#Recorrect MET
declareDefault("RECORRECTMET", True, globals())

#FSR mode
declareDefault("FSRMODE", "RunII", globals())

#Bunch spacing (can be 25 or 50)
declareDefault("BUNCH_SPACING", 25, globals())

#Mass used for SuperMELA
declareDefault("SUPERMELA_MASS", 125, globals())

#Selection flow strategy
declareDefault("SELSETUP", "allCutsAtOncePlusSmart", globals())

#Best candidate comparator (see interface/Comparators.h)
declareDefault("BESTCANDCOMPARATOR", "byBestKD", globals())

# Set to True to make candidates with the full combinatorial of loose leptons (for debug; much slower)
declareDefault("KEEPLOOSECOMB", False, globals())

# Activate the Z kinematic refit (very slow)
declareDefault("KINREFIT", False, globals())

# Activate paths for loose electron categories
declareDefault("ADDLOOSEELE", False, globals())

# Activate trigger paths in MC; note that for 2016, only reHLT samples have the correct triggers!!!
declareDefault("APPLYTRIG", True, globals())

# CMSSW version
CMSSW_VERSION = os.environ['CMSSW_VERSION']
CMSSWVERSION = int(CMSSW_VERSION.split("_")[1])




# The isolation cuts for electrons and muons. FIXME: there is an hardcoded instance of these values in src/LeptonIsoHelper.cc !!
ELEISOCUT = 0.35
MUISOCUT = 0.35


### ----------------------------------------------------------------------
### Set the GT
### ----------------------------------------------------------------------
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag


elif (SAMPLE_TYPE == 2017): 
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, '94X_mc2017_realistic_v14', '') #For JEC
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_v6', '')

elif (SAMPLE_TYPE == 2018):
    if IsMC:
        process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '') 
    else:
        process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9', '')


print '\t',process.GlobalTag.globaltag

### ----------------------------------------------------------------------
### Standard stuff
### ----------------------------------------------------------------------
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.GeometryDB_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


### ----------------------------------------------------------------------
### Source
### ----------------------------------------------------------------------
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
      '/eos/cms/store/user/covarell/HH/SM/testMINIAOD_HHSM_VV2l_0.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


### ----------------------------------------------------------------------
### Trigger bit Requests
### ----------------------------------------------------------------------
import HLTrigger.HLTfilters.hltHighLevel_cfi

process.hltFilterDiMu  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterMuEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterSingleEle = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterSingleMu = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.hltFilterDiMu.TriggerResultsTag  = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterMuEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterSingleEle.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterSingleMu.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.hltFilterDiMu.throw  = cms.bool(False) #FIXME: beware of this!
process.hltFilterDiEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterMuEle.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterSingleMu.throw = cms.bool(False) #FIXME: beware of this!
process.hltFilterSingleEle.throw = cms.bool(False) #FIXME: beware of this!


### triggers - FIXME
if (LEPTON_SETUP == 2018): 
    process.hltFilterDiEle.HLTPaths = ["HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                                       "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"]
    process.hltFilterDiMu.HLTPaths = ["HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
                                      "HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*"]
    process.hltFilterMuEle.HLTPaths = ["HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_v*",
                                       "HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ_v*",
                                       "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_v*",
                                       "HLT_Mu23_TrkIsoVVL_Ele8_CaloIdL_TrackIdL_IsoVL_DZ_v*"]
    process.hltFilterSingleEle.HLTPaths = ["HLT_Ele25_eta2p1_WPTight_Gsf_v*",
                                           "HLT_Ele27_WPTight_Gsf_v*"]
    process.hltFilterSingleMu.HLTPaths = ["HLT_IsoMu24_v*","HLT_IsoTkMu24_v*"]

process.triggerSingleEle = cms.Path(process.hltFilterSingleEle)
process.triggerSingleMu = cms.Path(process.hltFilterSingleMu)
process.triggerMuEle = cms.Path(process.hltFilterMuEle)
process.triggerDiEle = cms.Path(process.hltFilterDiEle)
process.triggerDiMu = cms.Path(process.hltFilterDiMu)
   

### ----------------------------------------------------------------------
### MET FILTERS
### ----------------------------------------------------------------------
process.METFilters  = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()
process.METFilters.TriggerResultsTag  = cms.InputTag("TriggerResults","","RECO")
if (IsMC):
	process.METFilters.TriggerResultsTag  = cms.InputTag("TriggerResults","","PAT")

# if (LEPTON_SETUP == 2017):#MET Filters available in miniAOD as described here https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
# 	if (IsMC):
# 		process.METFilters.HLTPaths = ["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadChargedCandidateFilter"]
# 	else:
# 		process.METFilters.HLTPaths = ["Flag_goodVertices","Flag_globalSuperTightHalo2016Filter","Flag_HBHENoiseFilter","Flag_HBHENoiseIsoFilter","Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_BadPFMuonFilter","Flag_BadChargedCandidateFilter","Flag_eeBadScFilter"]

process.triggerMETFilters = cms.Path(process.METFilters)

### ----------------------------------------------------------------------
### MC Filters and tools
### ----------------------------------------------------------------------

process.heavyflavorfilter = cms.EDFilter('HeavyFlavorFilter2',
#                                 src= cms.InputTag("genParticles"), # genParticles available only in PAT
                                 src= cms.InputTag("prunedGenParticles"),
                                 status2 = cms.bool(True),
                                 status3 = cms.bool(False),
                                 hDaughterVeto = cms.bool(False),
                                 zDaughterVeto = cms.bool(True),
                                 ptcut=cms.double(0)
                                 )


process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.drawTree = cms.EDAnalyzer("ParticleTreeDrawer",
                                   src = cms.InputTag("prunedGenParticles"),
                                   printP4 = cms.untracked.bool(False),
                                   printPtEtaPhi = cms.untracked.bool(False),
                                   printVertex = cms.untracked.bool(False),
                                   printStatus = cms.untracked.bool(True),
                                   printIndex = cms.untracked.bool(False) )


process.printTree = cms.EDAnalyzer("ParticleListDrawer",
                                   maxEventsToPrint = cms.untracked.int32(-1),
                                   printVertex = cms.untracked.bool(False),
                                   src = cms.InputTag("prunedGenParticles")
                                   )


process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                                                   calibratedPatElectrons = cms.PSet(
                                                   initialSeed = cms.untracked.uint32(1),
                                                   engineName = cms.untracked.string('TRandom3')
                                                   ),
                                                   )

# FIXME Add total kinematics filter for MC

process.goodPrimaryVertices = cms.EDFilter("VertexSelector",
  src = cms.InputTag("offlineSlimmedPrimaryVertices"),
  cut = cms.string("!isFake && ndof > 4 && abs(z) <= 24 && position.Rho <= 2"),
  filter = cms.bool(True),
)

# Filter for Zbb enrichment

process.zJetsFilter = cms.EDFilter("ZJetsFilterMerger",
  theLHESrc = cms.InputTag('externalLHEProducer'),
  theGenSrc = cms.InputTag('prunedGenParticles'),
  option = cms.untracked.int32(0),
)


### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### Loose lepton selection + cleaning + embeddding of user data
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------

SIP =  "userFloat('SIP')<4"
GOODLEPTON = "userFloat('ID') && " + SIP  # Lepton passing tight ID + SIP [ISO is asked AFTER FSR!!!]

if (LEPTON_SETUP >= 2016) : # Run II
    TIGHTMUON = "userFloat('isPFMuon') || (userFloat('isTrackerHighPtMuon') && pt>200)"
else:  
    TIGHTMUON = "userFloat('isPFMuon')" # Run I



#------- MUONS -------

#--- Mu e-scale corrections (Rochester)
# process.calibratedMuons =  cms.EDProducer("RochesterPATMuonCorrector",
#                                                src = cms.InputTag("patMuonsWithTrigger")
#                                               )

#--- Mu e-scale corrections (MuScleFit)
#process.calibratedMuons = cms.EDProducer("MuScleFitPATMuonCorrector",
#                         src = cms.InputTag("slimmedMuons"),
#                         debug = cms.bool(False),
#                         identifier = cms.string("Summer12_DR53X_smearReReco"),
#                         applySmearing = cms.bool(IsMC),
#                         fakeSmearing = cms.bool(False)
#                         )

#--- Mu e-scale corrections (KalmanMuonCalibrator, 2015)
process.calibratedMuons = cms.EDProducer("KalmanPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string(""),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )

#--- Set correct identifier for muon corrections

# elif LEPTON_SETUP == 2016: # (KalmanMuonCalibrator, ICHEP 2016) FIXME: still using the version from ICHEP16
#      if IsMC:
#          process.calibratedMuons.identifier = cms.string("MC_80X_13TeV")
#      else:
#          process.calibratedMuons.identifier = cms.string("DATA_80X_13TeV")
elif LEPTON_SETUP == 2017:# (Rochester corrections, Moriond 2018)
     process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string("RoccoR2017v0"),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )
elif LEPTON_SETUP == 2018:# FIXME: still 2017 version
     process.calibratedMuons = cms.EDProducer("RochesterPATMuonCorrector",
                                         src = cms.InputTag("slimmedMuons"),
                                         identifier = cms.string("RoccoR2017v0"),
                                         isMC = cms.bool(IsMC),
                                         isSynchronization = cms.bool(False),
                                         )

else:
    if APPLYMUCORR:
        print "APPLYMUCORR not configured for LEPTON_SETUP =", LEPTON_SETUP
        sys.exit()

#--- Mu Ghost cleaning
process.cleanedMu = cms.EDProducer("PATMuonCleanerBySegments",
                                   src = cms.InputTag("calibratedMuons"),
                                   preselection = cms.string("track.isNonnull"),
                                   passthrough = cms.string("isGlobalMuon && numberOfMatches >= 2"),
                                   fractionOfSharedSegments = cms.double(0.499))


process.bareSoftMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("cleanedMu"),
    cut = cms.string("pt>5 && abs(eta)<2.4 && (isGlobalMuon || (isTrackerMuon && numberOfMatches>0)) && muonBestTrackType!=2")
)


# MC matching. As the genParticles are no more available in cmg, we re-match with genParticlesPruned.
process.muonMatch = cms.EDProducer("MCMatcher", # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                   src     = cms.InputTag("softMuons"), # RECO objects to match
                                   matched = cms.InputTag("prunedGenParticles"),   # mc-truth particle collection
                                   mcPdgId     = cms.vint32(13), # one or more PDG ID (13 = muon); absolute values (see below)
                                   checkCharge = cms.bool(True), # True = require RECO and MC objects to have the same charge
                                   mcStatus = cms.vint32(1),     # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                   maxDeltaR = cms.double(0.5),  # Minimum deltaR for the match
                                   maxDPtRel = cms.double(0.5),  # Minimum deltaPt/Pt for the match
                                   resolveAmbiguities = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                   resolveByMatchQuality = cms.bool(False), # False = just match input in order; True = pick lowest deltaR pair first
                                   )

process.softMuons = cms.EDProducer("MuFiller",
    src = cms.InputTag("bareSoftMuons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1."),
    flags = cms.PSet(
        ID = cms.string(TIGHTMUON), # tight muon ID
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<" + str(MUISOCUT)),
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
    )
)


if APPLYMUCORR :
    process.muons =  cms.Sequence(process.calibratedMuons + process.cleanedMu + process.bareSoftMuons + process.softMuons)
else:
    process.cleanedMu.src = cms.InputTag("slimmedMuons")
    process.muons =  cms.Sequence(process.cleanedMu + process.bareSoftMuons + process.softMuons)




#------- ELECTRONS -------

#--- Run1 types of corrections below here, just as a reference
#
##--- Electron regression+calibrarion must be applied after BDT is recomputed
## NOTE patElectronsWithRegression->eleRegressionEnergy;  calibratedElectrons-> calibratedPatElectrons
## Default: NEW ECAL regression + NEW calibration + NEW combination
#process.load('EgammaAnalysis.ElectronTools.electronRegressionEnergyProducer_cfi')
#process.eleRegressionEnergy.inputElectronsTag = cms.InputTag('slimmedElectrons')
#process.eleRegressionEnergy.energyRegressionType = 2 ## 1: ECAL regression w/o subclusters 2 (default): ECAL regression w/ subclusters)
##process.eleRegressionEnergy.vertexCollection = cms.InputTag('goodPrimaryVertices')

#process.load("EgammaAnalysis.ElectronTools.calibratedPatElectrons_cfi")
#process.calibratedPatElectrons.correctionsType = 2 # 1 = old regression, 2 = new regression, 3 = no regression, 0 = nothing
#process.calibratedPatElectrons.combinationType = 3
#process.calibratedPatElectrons.lumiRatio = cms.double(1.0)
#process.calibratedPatElectrons.isMC    = IsMC
#process.calibratedPatElectrons.synchronization = cms.bool(False)



#--- Run2 electron momentum scale and resolution corrections

process.selectedSlimmedElectrons = cms.EDFilter("PATElectronSelector",
    ## this protects against a crash in electron calibration
    ## due to electrons with eta > 2.5
    src = cms.InputTag("slimmedElectrons"),
    cut = cms.string("pt>5 && abs(eta)<2.5 && abs(-log(tan(superClusterPosition.theta/2.)))<2.5")
)


## Preliminary Moriond 18 corrections
if (LEPTON_SETUP == 2017):
    process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducerRun2",
        electrons = cms.InputTag('selectedSlimmedElectrons'),
        gbrForestName = cms.vstring('electron_eb_ECALTRK_lowpt', 'electron_eb_ECALTRK',
                                    'electron_ee_ECALTRK_lowpt', 'electron_ee_ECALTRK',
                                    'electron_eb_ECALTRK_lowpt_var', 'electron_eb_ECALTRK_var',
                                    'electron_ee_ECALTRK_lowpt_var', 'electron_ee_ECALTRK_var'),
        isMC = cms.bool(IsMC),
        isSynchronization = cms.bool(False),
        correctionFile = cms.string("EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc"),
        recHitCollectionEB = cms.InputTag('reducedEgamma:reducedEBRecHits'),
        recHitCollectionEE = cms.InputTag('reducedEgamma:reducedEERecHits')
   )


if (BUNCH_SPACING == 50):
    process.calibratedPatElectrons.grbForestName = cms.string("gedelectron_p4combination_50ns")


if (LEPTON_SETUP <= 2017) :
    #--- Set up electron ID (VID framework)
    from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
    # turn on VID producer, indicate data format to be DataFormat.MiniAOD, as appropriate
    dataFormat = DataFormat.MiniAOD
    switchOnVIDElectronIdProducer(process, dataFormat)
    # define which IDs we want to produce
    my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_HZZ_V1_cff']
    # add them to the VID producer
    for idmod in my_id_modules:
        setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)
    # and don't forget to add the producer 'process.egmGsfElectronIDSequence' to the path, i.e. process.electrons


process.bareSoftElectrons = cms.EDFilter("PATElectronRefSelector",
   src = cms.InputTag("calibratedPatElectrons"),
   cut = cms.string("pt>7 && abs(eta)<2.5")
   )

process.softElectrons = cms.EDProducer("EleFiller",
   src    = cms.InputTag("bareSoftElectrons"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("userFloat('dxy')<0.5 && userFloat('dz')<1"),
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+str(ELEISOCUT))
#       Note: passCombRelIsoPFFSRCorr is currently set in LeptonPhotonMatcher for new FSR strategy; in ZZCandidateFiller for the old one
        ),
        mvaValuesMap = cms.InputTag(""),
   	correctionFile = cms.string(""),
   )

if (LEPTON_SETUP < 2017):
   process.softElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Spring16HZZV1Values")
#94X BDT with ID and Isolation
if (LEPTON_SETUP == 2017):
#	process.softElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV1Values")
   process.softElectrons.mvaValuesMap = cms.InputTag("electronMVAValueMapProducer:ElectronMVAEstimatorRun2Fall17IsoV2Values")

if (LEPTON_SETUP < 2018):
   process.softElectrons.correctionFile = process.calibratedPatElectrons.correctionFile


#process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons) # (use this version when running VID)
#process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons) # (use this version without VID)

# Handle special cases
if ELEREGRESSION == "None" and (ELECORRTYPE == "None" or BUNCH_SPACING == 50) :   # No correction at all. Skip correction modules.
#	 process.bareSoftElectrons.src = cms.InputTag('slimmedElectrons')
#	 process.electrons = cms.Sequence(process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons)
    process.bareSoftElectrons.src = cms.InputTag('selectedSlimmedElectrons')
    process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.bareSoftElectrons + process.softElectrons)
elif ELECORRTYPE == "RunII" :
	 process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag("calibratedPatElectrons")
	 process.electronMVAValueMapProducer.srcMiniAOD=cms.InputTag("calibratedPatElectrons")
if (ELEREGRESSION == "Moriond17v1" and LEPTON_SETUP == 2016):
	from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
	process = regressionWeights(process)
	process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

	process.selectedSlimmedElectrons.src = cms.InputTag("slimmedElectrons")
	process.electrons = cms.Sequence(process.regressionApplication + process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons)

if (LEPTON_SETUP == 2017): #For the moment regresion is applied on RECO level so no additional procedure is needed https://twiki.cern.ch/twiki/bin/view/CMS/Egamma2017DataRecommendations#Overview_of_E_gamma_Energy_Corre
	process.selectedSlimmedElectrons.src = cms.InputTag("slimmedElectrons")
	process.electrons = cms.Sequence(process.selectedSlimmedElectrons + process.calibratedPatElectrons + process.egmGsfElectronIDSequence + process.bareSoftElectrons + process.softElectrons)


#elif ELEREGRESSION == "None" and ELECORRTYPE == "RunII" :
#    process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('calibratedPatElectrons') # (when running VID)
#    process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('calibratedPatElectrons') # (when running VID)

#elif ELEREGRESSION == "Moriond" and ELECORRTYPE == "Moriond" : # Moriond corrections: OLD ECAL regression + OLD calibration + OLD combination
#    process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)
#    if (LEPTON_SETUP == 2011):
#        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2011Weights_V1.root")
#    else :
#        process.eleRegressionEnergy.regressionInputFile = cms.string("EgammaAnalysis/ElectronTools/data/eleEnergyReg2012Weights_V1.root")
#    process.eleRegressionEnergy.energyRegressionType = 1
#    process.calibratedPatElectrons.correctionsType   = 1
#    process.calibratedPatElectrons.combinationType   = 1
#
#elif ELEREGRESSION == "Paper" and ELECORRTYPE == "None" : # NEW ECAL regression + NO calibration + NO combination
#    process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)
#    process.eleRegressionEnergy.energyRegressionType = 2
#    process.calibratedPatElectrons.correctionsType   = 0
#    process.calibratedPatElectrons.combinationType   = 0
#
#elif ELEREGRESSION == "Paper" and ELECORRTYPE == "PaperNoComb" : # NEW ECAL regression + NEW calibration + NO combination
#    process.electrons = cms.Sequence(process.eleRegressionEnergy + process.calibratedPatElectrons + process.bareSoftElectrons + process.softElectrons)
#    process.eleRegressionEnergy.energyRegressionType = 2
#    process.calibratedPatElectrons.correctionsType   = 1
#    process.calibratedPatElectrons.combinationType   = 0


#--- TrackLess Electrons
process.bareSoftPhotons = cms.EDFilter("PATPhotonRefSelector",
   src = cms.InputTag("slimmedPhotons"),
   cut = cms.string("pt>7 && abs(eta)<2.5")
   )

process.softPhotons = cms.EDProducer("Philler",
   src    = cms.InputTag("bareSoftPhotons"),
   srcElectron = cms.InputTag("softElectrons"),
   mvaValuesMap = cms.InputTag("photonMVAValueMapProducer:TLEMVAEstimatorRun2Fall15V1Values"),
   mvaValuesMap2 = cms.InputTag("photonMVAValueMapProducer:PhotonMVAEstimatorRun2Spring15NonTrig25nsV2Values"),
   sampleType = cms.int32(SAMPLE_TYPE),
   setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
   cut = cms.string("1"), # dxy, dz not applied
   flags = cms.PSet(
        ID = cms.string("userFloat('isBDT')"),
        isSIP = cms.string(SIP),
        isGood = cms.string(GOODLEPTON),
        pass_lepton_ID = cms.string("userFloat('isBDT')"),
        pass_lepton_SIP = cms.string(SIP),
#        isIsoFSRUncorr  = cms.string("userFloat('combRelIsoPF')<"+ELEISOCUT), #TLE isolation is not corrected for FSR gammas.
        ),
   )


process.electronMatch = cms.EDProducer("MCMatcher",       # cut on deltaR, deltaPt/Pt; pick best by deltaR
                                       src         = cms.InputTag("bareSoftElectrons"), # RECO objects to match
                                       matched     = cms.InputTag("prunedGenParticles"), # mc-truth particle collection
                                       mcPdgId     = cms.vint32(11),               # one or more PDG ID (11 = electron); absolute values (see below)
                                       checkCharge = cms.bool(True),               # True = require RECO and MC objects to have the same charge
                                       mcStatus    = cms.vint32(1),                # PYTHIA status code (1 = stable, 2 = shower, 3 = hard scattering)
                                       maxDeltaR   = cms.double(0.5),              # Minimum deltaR for the match
                                       maxDPtRel   = cms.double(0.5),              # Minimum deltaPt/Pt for the match
                                       resolveAmbiguities    = cms.bool(True),     # Forbid two RECO objects to match to the same GEN object
                                       resolveByMatchQuality = cms.bool(False),    # False = just match input in order; True = pick lowest deltaR pair first
                                       )


### ----------------------------------------------------------------------
### Lepton Cleaning (clean electrons collection from muons)
### ----------------------------------------------------------------------

process.cleanSoftElectrons = cms.EDProducer("PATElectronCleaner",
    # pat electron input source
    src = cms.InputTag("softElectrons"),
    # preselection (any string-based cut for pat::Electron)
    preselection = cms.string(''),
    # overlap checking configurables
    checkOverlaps = cms.PSet(
        muons = cms.PSet(
           src       = cms.InputTag("softMuons"), # Start from loose lepton def
           algorithm = cms.string("byDeltaR"),
           preselection        = cms.string("userFloat('isGood')"),
           deltaR              = cms.double(0.05),
           checkRecoComponents = cms.bool(False), # don't check if they share some AOD object ref
           pairCut             = cms.string(""),
           requireNoOverlaps   = cms.bool(True), # overlaps don't cause the electron to be discared
        )
    ),
    # finalCut (any string-based cut for pat::Electron)
    finalCut = cms.string(''),
)



### ----------------------------------------------------------------------
### Search for FSR candidates
### ----------------------------------------------------------------------

# Create a photon collection; cfg extracted from "UFHZZAnalysisRun2.FSRPhotons.fsrPhotons_cff"
process.fsrPhotons = cms.EDProducer("PhotonFiller",
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE)  # "skip", "passThrough", "Legacy", "RunII"
)

import PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi
process.boostedFsrPhotons = PhysicsTools.PatAlgos.producersLayer1.pfParticleProducer_cfi.patPFParticles.clone(
    pfCandidateSource = 'fsrPhotons'
)

process.appendPhotons = cms.EDProducer("LeptonPhotonMatcher",
    muonSrc = cms.InputTag("softMuons"),
    electronSrc = cms.InputTag("cleanSoftElectrons"),
    photonSrc = cms.InputTag("boostedFsrPhotons"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    photonSel = cms.string(FSRMODE),  # "skip", "passThrough", "Legacy", "RunII"
    muon_iso_cut = cms.double(MUISOCUT),
    electron_iso_cut = cms.double(ELEISOCUT),
    debug = cms.untracked.bool(False),
    )

if ADDLOOSEELE:
    process.appendPhotons.looseElectronSrc = cms.InputTag("cleanSoftLooseElectrons")
#    process.appendPhotons.tleSrc = cms.InputTag("softPhotons")
#    process.appendPhotons.TLEMinPt = cms.double(25.)

# All leptons, any F/C.
# CAVEAT: merging creates copies of the objects, so that CandViewShallowCloneCombiner is not able to find
# overlaps between merged collections and the original ones.
process.softLeptons = cms.EDProducer("CandViewMerger",
#    src = cms.VInputTag(cms.InputTag("softMuons"), cms.InputTag("cleanSoftElectrons"))
    src = cms.VInputTag(cms.InputTag("appendPhotons:muons"), cms.InputTag("appendPhotons:electrons"))
)




### ----------------------------------------------------------------------
### ----------------------------------------------------------------------
### BUILD CANDIDATES
### ----------------------------------------------------------------------
### ----------------------------------------------------------------------



### ----------------------------------------------------------------------
### Dileptons: combine/merge leptons into intermediate (bare) collections;
###            Embed additional user variables into final collections
### ----------------------------------------------------------------------
TWOGOODLEPTONS = ("userFloat('d0.isGood') && userFloat('d1.isGood')") # Z made of 2 good leptons (ISO not yet applied)

### NOTE: Isolation cut has been moved to ZZ candidates as we now correct for FSR of all four photons.
### Because if this, isBestZ flags are no longer correct; BESTZ_AMONG is set to "" for safety
# ZISO           = ("( (abs(daughter(0).pdgId)==11 && userFloat('d0.combRelIsoPFFSRCorr')<0.5) || (abs(daughter(0).pdgId)==13 && userFloat('d0.combRelIsoPFFSRCorr')<0.4) ) && ( (abs(daughter(1).pdgId)==11 && userFloat('d1.combRelIsoPFFSRCorr')<0.5) || (abs(daughter(1).pdgId)==13 && userFloat('d1.combRelIsoPFFSRCorr')<0.4) )") #ISO after FSR
# ZLEPTONSEL     = TWOGOODLEPTONS + "&&" + ZISO
# BESTZ_AMONG = ( ZLEPTONSEL ) # "Best Z" chosen among those with 2 leptons with ID, SIP, ISO
# BESTZ_AMONG = ("")

ZLEPTONSEL     = TWOGOODLEPTONS # Note: this is without ISO

Z1PRESEL    = (ZLEPTONSEL + " && mass > 40") # Note: this is without ISO

BESTZ_AMONG = ( Z1PRESEL + "&& userFloat('d0.passCombRelIsoPFFSRCorr') && userFloat('d1.passCombRelIsoPFFSRCorr')" )

TWOGOODISOLEPTONS = ( TWOGOODLEPTONS + "&& userFloat('d0.passCombRelIsoPFFSRCorr') && userFloat('d1.passCombRelIsoPFFSRCorr')" )

# Cut to filter out unneeded ll combinations as upstream as possible
if KEEPLOOSECOMB:
    KEEPLOOSECOMB_CUT = 'mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())' # Propagate also combinations of loose leptons (for debugging); just require same-flavour
else:
    if FSRMODE == "RunII" : # Just keep combinations of tight leptons (passing ID, SIP and ISO)
        KEEPLOOSECOMB_CUT = "mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood') && daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr') &&  daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')"
    else :
        print "KEEPLOOSECOMB == False && FSRMODE =! RunII", FSRMODE, "is no longer supported"
        sys.exit()

### ----------------------------------------------------------------------
### Dileptons (Z->ee, Z->mm)
### ----------------------------------------------------------------------

# l+l- (SFOS, both e and mu)
process.bareZCand = cms.EDProducer("PATCandViewShallowCloneCombiner",
    decay = cms.string('softLeptons@+ softLeptons@-'),
    cut = cms.string(KEEPLOOSECOMB_CUT), # see below
    checkCharge = cms.bool(True)
)


if KEEPLOOSECOMB:
    process.bareZCand.cut = cms.string('mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId())') # Propagate also combinations of loose leptons (for debugging)
else:
    if FSRMODE == "RunII" : # Just keep combinations of tight leptons (passing ID, SIP and ISO)
        process.bareZCand.cut = cms.string("mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood') && daughter(0).masterClone.userFloat('passCombRelIsoPFFSRCorr') &&  daughter(1).masterClone.userFloat('passCombRelIsoPFFSRCorr')")
    else : # Just keep combinations of tight leptons (passing ID and SIP; iso cannot be required at this point, with the legacy FSR logic)
        process.bareZCand.cut = cms.string("mass > 0 && abs(daughter(0).pdgId())==abs(daughter(1).pdgId()) && daughter(0).masterClone.userFloat('isGood') && daughter(1).masterClone.userFloat('isGood')")


process.ZCand = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZCand"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ_AMONG),
    FSRMode = cms.string(FSRMODE), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string(ZLEPTONSEL),
        GoodIsoLeptons = cms.string(TWOGOODISOLEPTONS),
        Z1Presel = cms.string(Z1PRESEL),
    )
)




### ----------------------------------------------------------------------
### Jets
### ----------------------------------------------------------------------

from RecoJets.JetProducers.PileupJetIDParams_cfi import full_80x_chs
from RecoJets.JetProducers.PileupJetIDCutParams_cfi import full_80x_chs_wp

process.load("CondCore.CondDB.CondDB_cfi")



process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone(
    jets=cms.InputTag("slimmedJets"),
    inputIsCorrected=True,
    applyJec=True,
    vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
    )

### q/g likelihood
qgDatabaseVersion = 'cmssw8020_v2'
process.QGPoolDBESSource = cms.ESSource("PoolDBESSource",
      DBParameters = cms.PSet(messageLevel = cms.untracked.int32(1)),
      timetype = cms.string('runnumber'),
      toGet = cms.VPSet(
        cms.PSet(
            record = cms.string('QGLikelihoodRcd'),
            tag    = cms.string('QGLikelihoodObject_'+qgDatabaseVersion+'_AK4PFchs'),
            label  = cms.untracked.string('QGL_AK4PFchs')
        ),
      ),
      connect = cms.string('sqlite_fip:HHbbZZAnalysis/AnalysisStep/data/QGTagging/QGL_'+qgDatabaseVersion+'.db')
)
process.es_prefer_qg = cms.ESPrefer('PoolDBESSource','QGPoolDBESSource')
process.load('RecoJets.JetProducers.QGTagger_cfi')
process.QGTagger.srcJets = cms.InputTag( 'slimmedJets' )
process.QGTagger.jetsLabel = cms.string('QGL_AK4PFchs')

process.dressedJets = cms.EDProducer("JetFiller",
    src = cms.InputTag("slimmedJets"),
    sampleType = cms.int32(SAMPLE_TYPE),
    setup = cms.int32(LEPTON_SETUP),
    cut = cms.string("pt>20 && abs(eta)<4.7 && userFloat('looseJetID') && userFloat('PUjetID')"),
    isMC = cms.bool(IsMC),
    bTaggerName = cms.string("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
    bTaggerThreshold = cms.double(0.8484), #CSVv2M, from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    jecType = cms.string("AK4PFchs"),
    applyJER = cms.bool(APPLYJER),
    jerType = cms.string("AK4PFchs"),
    bTagSFFile = cms.string("HHbbZZAnalysis/AnalysisStep/data/BTagging/CSVv2Moriond17_2017_1_26_BtoH.csv"), #Preliminary Moriond17 SF, from https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
    bTagMCEffFile = cms.string("HHbbZZAnalysis/AnalysisStep/data/BTagging/bTagEfficiencies_80X_ICHEP.root"),
    flags = cms.PSet()
    )
if (LEPTON_SETUP == 2017):
    process.dressedJets.bTaggerName = cms.string("pfDeepCSVJetTags:probb") #Moving to Moriond18 new recommended DeepCSV btagger
    process.dressedJets.bTaggerThreshold = cms.double(0.4941) #https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
    process.dressedJets.bTagSFFile =  cms.string("HHbbZZAnalysis/AnalysisStep/data/BTagging/DeepCSV_94XSF_V1_B_F.csv")
    process.dressedJets.bTagMCEffFile  = cms.string("HHbbZZAnalysis/AnalysisStep/data/BTagging/bTagEfficiencies_94X_Moriond18_v1.root") #FIXME!!! Update to Moriond18 file


### Load JEC
if (APPLYJEC and SAMPLE_TYPE == 2017):
    if IsMC:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Fall17_17Nov2017_V8_MC_AK4PFchs'), #for 94X/MET fix 
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
              connect = cms.string('sqlite_fip:HHbbZZAnalysis/AnalysisStep/data/JEC/Fall17_17Nov2017_V8_MC.db'), #for 94X/MET fix
            )
    else:
        process.jec = cms.ESSource("PoolDBESSource",
            DBParameters = cms.PSet(
                messageLevel = cms.untracked.int32(1)
                ),
            timetype = cms.string('runnumber'),
            toGet = cms.VPSet(
                cms.PSet(
                    record = cms.string('JetCorrectionsRecord'),
                    tag    = cms.string('JetCorrectorParametersCollection_Fall17_17Nov2017BCDEF_V6_DATA_AK4PFchs'), #for 94X/Moriond18
                    label  = cms.untracked.string('AK4PFchs')
                    ),
                ),
            connect = cms.string('sqlite_fip:HHbbZZAnalysis/AnalysisStep/data/JEC/Fall17_17Nov2017BCDEF_V6_DATA.db'),
            )

    ## add an es_prefer statement to resolve a possible conflict from simultaneous connection to a global tag
    process.es_prefer_jec = cms.ESPrefer('PoolDBESSource','jec')

    ### reapply JEC
    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJetCorrFactors
    process.patJetCorrFactorsReapplyJEC = updatedPatJetCorrFactors.clone(
        src = cms.InputTag("slimmedJets"),
        levels = ['L1FastJet','L2Relative','L3Absolute'],
        payload = 'AK4PFchs' )
    if not IsMC:
        process.patJetCorrFactorsReapplyJEC.levels = ['L1FastJet','L2Relative','L3Absolute','L2L3Residual']

    from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cff import updatedPatJets
    process.patJetsReapplyJEC = updatedPatJets.clone(
        jetSource = cms.InputTag("slimmedJets"),
        jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJEC"))
        )

    ### add pileup id and discriminant to patJetsReapplyJEC
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']

    ### Replace inputs in QGTagger and dressedJets
    process.QGTagger.srcJets = cms.InputTag( 'patJetsReapplyJEC')
    process.dressedJets.src = cms.InputTag('patJetsReapplyJEC')



### Clean jets wrt. good (preFSR-)isolated leptons
process.cleanJets = cms.EDProducer("JetsWithLeptonsRemover",
                                   Jets      = cms.InputTag("dressedJets"),
                                   Muons     = cms.InputTag("appendPhotons:muons"),
                                   Electrons = cms.InputTag("appendPhotons:electrons"),
                                   Diboson   = cms.InputTag(""),
                                   JetPreselection      = cms.string(""),
                                   MuonPreselection     = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   ElectronPreselection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   DiBosonPreselection  = cms.string(""),
                                   MatchingType = cms.string("byDeltaR"),
                                   cleanFSRFromLeptons = cms.bool(True),
                                   DebugPlots = cms.untracked.bool(False),
                                   DebugPrintOuts = cms.untracked.bool(False)
                                   )

if FSRMODE=="Legacy" :
    process.cleanJets.MuonPreselection     = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJets.ElectronPreselection = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJets.cleanFSRFromLeptons = False


### ----------------------------------------------------------------------
### Fat jets
### ----------------------------------------------------------------------

### reapply JEC
from PhysicsTools.SelectorUtils.pfJetIDSelector_cfi import pfJetIDSelector

process.patJetCorrFactorsReapplyJECFat = updatedPatJetCorrFactors.clone(
                                    src     = cms.InputTag("slimmedJetsAK8"),
                                    levels  = ['L1FastJet','L2Relative','L3Absolute'],
                                    payload = 'AK8PFchs')

process.patJetsReapplyJECFat = updatedPatJets.clone(
                                    jetSource = cms.InputTag("slimmedJetsAK8"),
                                    jetCorrFactorsSource = cms.VInputTag(cms.InputTag("patJetCorrFactorsReapplyJECFat") ))

### FIXME: should we be dressing the fat jets as well?
process.goodJetsFat = cms.EDFilter("PFJetIDSelectionFunctorFilter",
                         filterParams = pfJetIDSelector.clone(),
                         src = cms.InputTag("patJetsReapplyJECFat"),
                         filter = cms.bool(True) )

# Clean jets wrt. good (preFSR-)isolated leptons
process.cleanJetsFat = cms.EDProducer("JetsWithLeptonsRemover",
                                   Jets      = cms.InputTag("goodJetsFat"),
                                   Muons     = cms.InputTag("appendPhotons:muons"),
                                   Electrons = cms.InputTag("appendPhotons:electrons"),
                                   Diboson   = cms.InputTag(""),
                                   JetPreselection      = cms.string(""),
                                   MuonPreselection     = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   ElectronPreselection = cms.string("userFloat('isGood') && userFloat('passCombRelIsoPFFSRCorr')"),
                                   DiBosonPreselection  = cms.string(""),
                                   MatchingType = cms.string("byDeltaR"),
                                   DeltaRCut = cms.untracked.double(0.8),  
                                   cleanFSRFromLeptons = cms.bool(True),
                                   DebugPlots = cms.untracked.bool(False)
                                   )

process.corrJetsProducer = cms.EDProducer ( "CorrJetsProducer",
                                    jets    = cms.InputTag( "cleanJetsFat" ),
                                    vertex  = cms.InputTag( "goodPrimaryVertices" ), 
                                    rho     = cms.InputTag( "fixedGridRhoFastjetAll"   ),
                                    payload = cms.string  ( "AK8PFchs" ),
                                    isData  = cms.bool    (  False ))

if FSRMODE=="Legacy" :
    process.cleanJetsFat.MuonPreselection     = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJetsFat.ElectronPreselection = "userFloat('isGood') && userFloat('isIsoFSRUncorr')"
    process.cleanJetsFat.cleanFSRFromLeptons = False


### ----------------------------------------------------------------------
### 2l2j cand
### ----------------------------------------------------------------------

process.bareZjjCand = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string('cleanJets cleanJets'),
    cut = cms.string('pt>100 && mass>40 && mass <180'), # protect against ghosts
    checkCharge = cms.bool(False)
)
process.ZjjCand = cms.EDProducer("ZCandidateFiller",
    src = cms.InputTag("bareZjjCand"),
    sampleType = cms.int32(SAMPLE_TYPE),  
    embedDaughterFloats = cms.untracked.bool(True),                   
    setup = cms.int32(LEPTON_SETUP), # define the set of effective areas, rho corrections, etc.
    bestZAmong = cms.string(BESTZ2_AMONG),
    FSRMode = cms.string("skip"), # "skip", "Legacy", "RunII"
    flags = cms.PSet(
        GoodLeptons = cms.string('0'),
        Z1Presel = cms.string(Z2PRESEL),
    )
)




### ----------------------------------------------------------------------
### Missing ET
### ----------------------------------------------------------------------

metTag = cms.InputTag("slimmedMETsMuEGClean")


### Recorrect MET, cf. https://indico.cern.ch/event/759372/contributions/3149378/attachments/1721436/2779341/metreport.pdf slide 10
###                and https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETUncertaintyPrescription#Instructions_for_9_4_X_X_9_for_2
if (RECORRECTMET and SAMPLE_TYPE == 2017):

    from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
   		
    runMetCorAndUncFromMiniAOD(process,
                               isData=(not IsMC),
                               fixEE2017 = True,
                               fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
                               postfix = "ModifiedMET"
                               )
    metTag = cms.InputTag("slimmedMETsModifiedMET","","bbZZ")
	 
    ### somehow MET recorrection gets this lost again...
    process.patJetsReapplyJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']
    process.patJetsReapplyJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId']


### ----------------------------------------------------------------------
### Paths
### ----------------------------------------------------------------------

process.preSkimCounter = cms.EDProducer("EventCountProducer")
process.PVfilter =  cms.Path(process.preSkimCounter+process.goodPrimaryVertices)

if APPLYJEC:
    process.Jets = cms.Path(process.pileupJetIdUpdated + process.patJetCorrFactorsReapplyJEC + process.patJetsReapplyJEC + process.QGTagger + process.dressedJets )
else:
    process.Jets = cms.Path( process.QGTagger + process.dressedJets )

if (RECORRECTMET and SAMPLE_TYPE == 2017):
    if IsMC:
        process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)
    else:
        process.MET = cms.Path(process.fullPatMetSequenceModifiedMET)

### ----------------------------------------------------------------------
### Filters
### ----------------------------------------------------------------------
### Create filter for events with one candidate in the SR
process.ZZCandSR = cms.EDFilter("PATCompositeCandidateRefSelector",
    src = cms.InputTag("bbZZCand"),
    cut = cms.string(SR) # That is, all candidates with mZ1>70 ... #FIXME: to be defined
 )

process.bbZZCandFilter = cms.EDFilter("CandViewCountFilter",
                                src = cms.InputTag("bbZZCandSR"),
                                minNumber = cms.uint32(1)
                            )

# Prepare lepton collections
process.Candidates = cms.Path(
       process.muons             +
       process.electrons         + process.cleanSoftElectrons +
       process.fsrPhotons        + process.boostedFsrPhotons +
       process.appendPhotons     +
       process.softLeptons       +
       process.cleanJets         +
       process.cleanJetsFat      + process.corrJetsProducer + 
# Build 4-lepton candidates
       process.bareZCand         + process.ZCand     +
       process.bareZjjCand       + process.ZjjCand   +  
       process.bareZZCand        + process.ZZCand    +
       process.bareZZCandFat     + process.ZZCandFat
    )



### Skim, triggers and MC filters (Only store filter result, no filter is applied)

# 2012 skim, Reimplementation by Giovanni
#process.load("HHbbZZAnalysis.AnalysisStep.Skim2012_cfg")
#process.SkimSequence = cms.Sequence(process.HZZSkim2012)
#process.Skim = cms.Path(process.SkimSequence)

#SkimPaths = cms.vstring('Skim')
SkimPaths = cms.vstring('PVfilter') #Do not apply skim, just require a good PV

# process.HF = cms.Path(process.heavyflavorfilter)

# FIXME total kin filter?


if (ADDLOOSEELE) :
    import os
    execfile(os.environ['CMSSW_BASE'] + "/src/HHbbZZAnalysis/AnalysisStep/test/MasterPy/LooseEle.py")
#    execfile(os.environ['CMSSW_BASE'] + "/src/HHbbZZAnalysis/AnalysisStep/test/MasterPy/TracklessEle.py")


process.GlobalTag.toGet = cms.VPSet( (cms.PSet(
        connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
        label = cms.untracked.string('electron_eb_ECALonly'),
        record = cms.string('GBRDWrapperRcd'),
        tag = cms.string('GEDelectron_EBCorrection_80X_EGM_v4')
    ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_eb_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('electron_ee_ECALTRK_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDelectron_track_lowpt_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBCorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_eb_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EBUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EECorrection_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_EEUncertainty_80X_EGM_v4')
        ), 
        cms.PSet(
            connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
            label = cms.untracked.string('photon_ee_ECALonly_lowpt_var'),
            record = cms.string('GBRDWrapperRcd'),
            tag = cms.string('GEDphoton_lowpt_EEUncertainty_80X_EGM_v4')
        ),
    )
)
