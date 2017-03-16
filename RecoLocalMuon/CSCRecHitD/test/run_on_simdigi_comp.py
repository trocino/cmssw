## Dump  10  events in CSC rechit builder - Tim Cox - 07.11.2012
## This version runs in 6_0_1_PostLS1 on a simulated data DIGI sample.

import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.StandardSequences.EndOfProcess_cff")

process.load("RecoLocalMuon.CSCRecHitD.cscRecHitFromComparator_cfi")


# --- MATCH GT TO RELEASE AND DATA SAMPLE

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')


# --- NUMBER OF EVENTS

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#process.options   = cms.untracked.PSet( SkipEvent = cms.untracked.vstring("ProductNotFound") )
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.source    = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
         "/store/relval/CMSSW_8_0_16/RelValZpMM_13/GEN-SIM-DIGI-RAW-HLTDEBUG/80X_mcRun2_asymptotic_v16_gs7120p2-v1/10000/12D7D3ED-6557-E611-9149-0CC47A4D767C.root"
    )
)

# ME1/1A is unganged Post-LS1

process.CSCGeometryESModule.useGangedStripsInME1a = False
##process.CSCGeometryESModule.debugV = True
##process.idealForDigiCSCGeometry.useGangedStripsInME1a = False

# Turn off some flags for CSCRecHitD that are turned ON in default config

process.csc2DRecHitsFromComparator.readBadChannels = cms.bool(False)
##process.csc2DRecHitsFromComparator.CSCUseTimingCorrections = cms.bool(False)
process.csc2DRecHitsFromComparator.CSCUseGasGainCorrection = cms.bool(False)

## Switch input for CSCRecHitD to simulated digis (reco with Comparator digis) 
process.csc2DRecHitsFromComparator.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
process.csc2DRecHitsFromComparator.compDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCComparatorDigi")

## Switch input for CSCRecHitD to simulated digis (standard reco) 
process.csc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
process.csc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")


## TO ACTIVATE LogTrace IN CSCRecHitD, NEED TO COMPILE IT WITH 
##  scram b -j8 USER_CXXFLAGS="-DEDM_ML_DEBUG" 
## (If the release is already compiled, clean it up with 
##  scram b clean; scram b vclean 
## first) 
## LogTrace output goes to cout; all other output to "junk.log"

process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load("FWCore.MessageLogger.MessageLogger_cfi")
##process.MessageLogger.categories.append("CSCGeometry")
process.MessageLogger.categories.append("CSCRecHit")  ## COMMENT OUT THIS LINE to turn off the verbose printout 
##process.MessageLogger.categories.append("CSCRecHitDBuilder")
##process.MessageLogger.categories.append("CSCMake2DRecHit")
## process.MessageLogger.categories.append("CSCHitFromStripOnly")
## process.MessageLogger.categories.append("CSCRecoConditions")
# module label is something like "muonCSCDigis"...
process.MessageLogger.debugModules = cms.untracked.vstring("*")
process.MessageLogger.destinations = cms.untracked.vstring("cout","junk")
process.MessageLogger.cout = cms.untracked.PSet(
    threshold = cms.untracked.string("DEBUG"),
    default   = cms.untracked.PSet( limit = cms.untracked.int32(0)  ),
    FwkReport = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#   , CSCGeometry = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
    , CSCRecHit = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#   , CSCRecHitDBuilder = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#   , CSCMake2DRecHit = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#   , CSCHitFromStripOnly = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
#   , CSCRecoConditions = cms.untracked.PSet( limit = cms.untracked.int32(-1) )
)


process.FEVTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('file:reco_csc_from_comp.root'),
    outputCommands = process.FEVTDEBUGEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)
## Watch it: if you run the full CMS reco, the following line will save LOTS of collections! 
## (You may want to restrict the collections you keep) 
process.FEVTDEBUGoutput.outputCommands.append('keep *_*_*_TEST')

## MODIFIER: replace standard CSC hit reco with new CSC Comparator reco 
## Comment out the next four lines to use the standard CSC hit reco 
modifier = cms.Modifier()
modifier._setChosen()
if hasattr(process, "csc2DRecHits"):
    modifier.toReplaceWith(process.csc2DRecHits, process.csc2DRecHitsFromComparator)

## Path and EndPath def
#process.reco = cms.Path(process.csc2DRecHits)                ## CSC hit reco only (standard) 
#process.reco = cms.Path(process.csc2DRecHitsFromComparator)  ## CSC hit reco only (using Comparator digis) 
process.reco = cms.Path(process.csclocalreco)                ## CSC hit + segment reco (standard or with Comparator, depending on the "modifier" above) 
#process.reco = cms.Path(process.reconstruction)              ## Full CMS reconstruction (standard or with Comparator, depending on the "modifier" above) 

#process.endjob = cms.EndPath(process.endOfProcess)
process.FEVTDEBUGoutput_step = cms.EndPath(process.FEVTDEBUGoutput)

# Schedule definition
#process.schedule = cms.Schedule(process.reco, process.endjob)
process.schedule = cms.Schedule(process.reco, process.FEVTDEBUGoutput_step)
