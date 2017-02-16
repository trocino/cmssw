import FWCore.ParameterSet.Config as cms

from TrackingTools.TrackAssociator.default_cfi import *

me0SegmentMatching = cms.EDProducer("ME0SegmentMatcher",
                                    TrackAssociatorParameterBlock,
                                    maxPullX = cms.double (3.0),
                                    maxDiffX = cms.double (4.0),
                                    maxPullY = cms.double (-1.0),
                                    maxDiffY = cms.double (-1.0),
                                    maxDiffPhiDirection = cms.double (-1.0),
                                    me0SegmentTag = cms.InputTag("me0Segments"),
                                    tracksTag = cms.InputTag("generalTracks")
                                    )
me0SegmentMatching.TrackAssociatorParameters.useGEM = cms.bool(True)
me0SegmentMatching.TrackAssociatorParameters.useME0 = cms.bool(True)
me0SegmentMatching.TrackAssociatorParameters.useEcal = cms.bool(False)
me0SegmentMatching.TrackAssociatorParameters.useHO = cms.bool(False)
me0SegmentMatching.TrackAssociatorParameters.useMuon = cms.bool(True)
me0SegmentMatching.TrackAssociatorParameters.useCalo = cms.bool(False)
me0SegmentMatching.TrackAssociatorParameters.usePreshower = cms.bool(False)
me0SegmentMatching.TrackAssociatorParameters.useHcal = cms.bool(False)
 
