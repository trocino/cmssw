# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --no_exec --conditions auto:upgradePLS3 -n 10 --eventcontent FEVTDEBUGHLT,DQM -s RAW2DIGI,L1Reco,RECO,VALIDATION,DQM --datatier GEN-SIM-RECO,DQM --customise SLHCUpgradeSimulations/Configuration/combinedCustoms.cust_2023Muon --geometry Extended2023Muon,Extended2023MuonReco --magField 38T_PostLS1 --filein file:step2.root --fileout file:step3.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('ME0REFIT')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2023MuonReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.Validation_cff')
#process.load('DQMOffline.Configuration.DQMOfflineMC_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    #input = cms.untracked.int32(3398)
)

# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/tmp/dnash/step3_TryWithoutMET.root'),
#     fileNames = cms.untracked.vstring(#'file:/tmp/dnash/step3.root',
#                                       #'file:/tmp/dnash/step3_Post130.root',
#                                       #'file:/tmp/dnash/step3_Post262.root',
#                                       #'file:/tmp/dnash/step3_Post600.root'
#                                       )
                            #skipEvents = cms.untracked.uint32(3399)
                            skipEvents = cms.untracked.uint32(4109)
                            )

process.options = cms.untracked.PSet(

)

# # Production Info
# process.configurationMetadata = cms.untracked.PSet(
#     version = cms.untracked.string('$Revision: 1.20 $'),
#     annotation = cms.untracked.string('step3 nevts:10'),
#     name = cms.untracked.string('Applications')
# )

# Output definition

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring('keep *_*_*_*'),
    #outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('file:testME0refit_Post262.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RECO')
    )
)
#process.FEVTDEBUGHLToutput.outputCommands.append('keep *_*_*_ME0REFIT')

# process.DQMoutput = cms.OutputModule("PoolOutputModule",
#     splitLevel = cms.untracked.int32(0),
#     outputCommands = process.DQMEventContent.outputCommands,
#     fileName = cms.untracked.string('file:step3_inDQM.root'),
#     dataset = cms.untracked.PSet(
#         filterName = cms.untracked.string(''),
#         dataTier = cms.untracked.string('DQM')
#     )
# )

# Additional output definition

# Other statements
#process.mix.playback = True
#process.mix.digitizers = cms.PSet()
#for a in process.aliases: delattr(process, a)
#process.RandomNumberGeneratorService.restoreStateLabel=cms.untracked.string("randomEngineStateProducer")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.load("TrackingTools.TrackRefitter.ME0TracksToTrajectories_cff")
process.me0Muons = cms.EDProducer("ME0TracksToTrajectories",
                                  Tracks = cms.InputTag("me0SegmentMatcher"),
                                  Type = cms.string("ME0Muons"),
                                  TrackTransformer = cms.PSet(
    Fitter = cms.string('KFFitterForRefitInsideOut'),
    TrackerRecHitBuilder = cms.string('WithTrackAngle'),
    Smoother = cms.string('KFSmootherForRefitInsideOut'),
    MuonRecHitBuilder = cms.string('MuonRecHitBuilder'),
    RefitDirection = cms.string('insideOut'),
    RefitRPCHits = cms.bool(True),
    DoPredictionsOnly = cms.bool(False), 
    Propagator = cms.string('SmartPropagatorAnyRK'),
    MuonPropagator = cms.string('SteppingHelixPropagatorAny')
    )
                                  )

process.testMe0muons = cms.Path(process.me0Muons)

# Path and EndPath definitions
#process.raw2digi_step = cms.Path(process.RawToDigi)
#process.L1Reco_step = cms.Path(process.L1Reco)
#process.reconstruction_step = cms.Path(process.reconstruction)
#process.prevalidation_step = cms.Path(process.prevalidation)
#process.dqmoffline_step = cms.Path(process.DQMOffline)
#process.validation_step = cms.EndPath(process.validation)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.prevalidation_step,process.validation_step,process.dqmoffline_step,process.endjob_step,process.FEVTDEBUGHLToutput_step,process.DQMoutput_step)
process.schedule = cms.Schedule(process.testMe0muons, process.endjob_step, process.FEVTDEBUGHLToutput_step)

# customisation of the process.

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023Muon 

#call to customisation function cust_2023Muon imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
process = cust_2023Muon(process)

# Automatic addition of the customisation function from SimGeneral.MixingModule.fullMixCustomize_cff
#from SimGeneral.MixingModule.fullMixCustomize_cff import setCrossingFrameOn 

#call to customisation function setCrossingFrameOn imported from SimGeneral.MixingModule.fullMixCustomize_cff
#process = setCrossingFrameOn(process)

# End of customisation functions
