import FWCore.ParameterSet.Config as cms
#process = cms.Process("ME0SegmentMatchingLocalTest")

from Configuration.StandardSequences.Eras import eras
process = cms.Process("ME0SegmentMatchingLocalTest", eras.Phase2C1)


## Standard sequence
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

# process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDevReco_cff')
# process.load('Configuration.Geometry.GeometryExtended2015MuonGEMDev_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4_cff')     # 1 eta partition 
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff') # 1 eta partition 

#process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
#process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## TrackingComponentsRecord required for matchers
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi')
process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi')

## global tag for 2019 upgrade studies
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgrade2019', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '') 



# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.combinedCustoms
#from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023SHCal 
#from SLHCUpgradeSimulations.Configuration.combinedCustoms import cust_2023Muon 

#call to customisation function cust_2023SHCal imported from SLHCUpgradeSimulations.Configuration.combinedCustoms
#process = cust_2023SHCal(process)
#process = cust_2023Muon(process)



# # Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.me0Customs
# from SLHCUpgradeSimulations.Configuration.me0Customs import customise 
# process = customise(process)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:///somewhere/simevent.root') ##/somewhere/simevent.root" 
                            )
process.PoolSource.fileNames = [    
    'file:/afs/cern.ch/user/p/piet/public/ForDaniele/TenMuGun_820p1_D4_100Evts/step3.root'
    ]

# #process.load('RecoLocalMuon.GEMRecHit.me0RecHits_cfi')
# #process.load('RecoLocalMuon.GEMSegments.me0Segments_cfi')
# process.load('RecoMuon.MuonIdentification.me0MuonReco_cff')

# #process.p = cms.Path(process.me0RecHits*process.me0Segments*process.me0MuonReco)
# process.p = cms.Path(process.me0MuonReco)
# #process.p = cms.Path(process.me0RecHits*process.me0Segments)

# setattr(process, "newMe0SegmentMatching", process.me0SegmentMatching.clone())
# process.pathProva = cms.Path( getattr(process, "newMe0SegmentMatching") ) 

process.newMe0SegmentMatching = process.me0SegmentMatching.clone() 
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useGEM)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useME0)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useMuon)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useCalo)
process.newMe0SegmentMatching.TrackAssociatorParameters.useGEM = cms.bool(False)
process.newMe0SegmentMatching.TrackAssociatorParameters.useME0 = cms.bool(True)
process.newMe0SegmentMatching.TrackAssociatorParameters.useEcal = cms.bool(False)
process.newMe0SegmentMatching.TrackAssociatorParameters.useHO = cms.bool(False)
process.newMe0SegmentMatching.TrackAssociatorParameters.useMuon = cms.bool(True)
process.newMe0SegmentMatching.TrackAssociatorParameters.useCalo = cms.bool(False)
process.newMe0SegmentMatching.TrackAssociatorParameters.usePreshower = cms.bool(False)
process.newMe0SegmentMatching.TrackAssociatorParameters.useHcal = cms.bool(False)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useGEM)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useME0)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useMuon)
print str(process.newMe0SegmentMatching.TrackAssociatorParameters.useCalo)

process.pathProva = cms.Path(process.newMe0SegmentMatching) 

process.o1 = cms.OutputModule("PoolOutputModule",
                              outputCommands = cms.untracked.vstring(
        'keep *',
        #'drop *_newMe0SegmentMatching_*_*'
        #'drop *_me0SegmentMatcher_*_*'
        #'drop *',
        ##'keep *_me0SegmentMatching_*_*',
        #'keep *_me0MuonConverting_*_*',
        ),
                              fileName = cms.untracked.string('out_me0Reco.root')
                              )
process.outpath = cms.EndPath(process.o1)
