#ifndef TrackingTools_TrackRefitter_ME0TrackTransformer_H
#define TrackingTools_TrackRefitter_ME0TrackTransformer_H

/** \class ME0TrackTransformer
 *  This class takes a reco::Track and refits the rechits inside it.
 *  The final result is a Trajectory refitted and smoothed.
 *  To make the refitting (and the smoothing) the usual KF tools are used.
 *
 *  CAVEAT: till now (it will be changed in the near future) the class stores the
 *  pointers to the services, therefore EACH event the setServices(const edm::EventSetup&)
 *  method MUST be called in the code in which the ME0TrackTransformer is used.
 *
 *  $Date: 2013/05/28 17:59:53 $
 *  $Revision: 1.16 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */

#include "TrackingTools/TrackRefitter/interface/TrackTransformerBase.h"

#include "TrackingTools/TrackRefitter/interface/RefitDirection.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"

//#include "ME0Reconstruction/ME0Segment/interface/ME0MuonCollection.h"
#include "DataFormats/MuonReco/interface/ME0MuonCollection.h"

// Only for this one time!!!!
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

namespace edm {class ParameterSet; class EventSetup;}
namespace reco {class TransientTrack; class ME0Muon;}

class TrajectoryFitter;
class TrajectorySmoother;
class Propagator;
class TransientTrackingRecHitBuilder;
class Trajectory;
//class ME0SegmentRef; 

class ME0TrackTransformer: public TrackTransformerBase{

public:

  /// Constructor
  ME0TrackTransformer(const edm::ParameterSet&);

  /// Destructor
  virtual ~ME0TrackTransformer();
  
  // Operations

  /// Convert a reco::Track into Trajectory
  virtual std::vector<Trajectory> transform(const reco::Track&) const; // this is not used, but must be re-implemented (pure virtual in TrackTransformerBase) 
  //virtual std::vector<Trajectory> transform(const reco::ME0Muon&) const;
  virtual std::vector<Trajectory> transform(const reco::ME0Muon&, const edm::Handle<reco::GenParticleCollection>&) const;

  /// Convert a reco::TrackRef into Trajectory
  //std::vector<Trajectory> transform(const reco::TrackRef&) const;

  /// Convert a reco::TrackRef into Trajectory, refit with a new set of hits
//   std::vector<Trajectory> transform(const reco::TransientTrack&,
//                                     const TransientTrackingRecHit::ConstRecHitContainer&, 
// 				    const EmulatedME0SegmentRef) const;
  std::vector<Trajectory> transform(const reco::TransientTrack&,
                                    const TransientTrackingRecHit::ConstRecHitContainer&, 
				    const EmulatedME0SegmentRef,
				    const edm::Handle<reco::GenParticleCollection>&) const;

  /// the magnetic field
  const MagneticField* magneticField() const {return &*theMGField;}
  
  /// the tracking geometry
  edm::ESHandle<GlobalTrackingGeometry> trackingGeometry() const {return theTrackingGeometry;}

  /// set the services needed by the ME0TrackTransformer
  virtual void setServices(const edm::EventSetup&);

  /// the refitter used to refit the reco::Track
  edm::ESHandle<TrajectoryFitter> refitter() const {return theFitter;}
  
  /// the smoother used to smooth the trajectory which came from the refitting step
  edm::ESHandle<TrajectorySmoother> smoother() const {return theSmoother;}

  TransientTrackingRecHit::ConstRecHitContainer
    getTransientRecHits(const reco::TransientTrack& track) const;
  
 protected:
  
 private:

  std::string thePropagatorName;
  edm::ESHandle<Propagator> propagator() const {return thePropagator;}
  edm::ESHandle<Propagator> thePropagator;
  
  std::string theMuonPropagatorName;
  edm::ESHandle<Propagator> muonPropagator() const {return theMuonPropagator;}
  edm::ESHandle<Propagator> theMuonPropagator;
  
  unsigned long long theCacheId_TC;
  unsigned long long theCacheId_GTG;
  unsigned long long theCacheId_MG;
  unsigned long long theCacheId_TRH;
  
  bool theRPCInTheFit;

  bool theDoPredictionsOnly;
  RefitDirection theRefitDirection;

  edm::ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
  edm::ESHandle<MagneticField> theMGField;
  
  std::string theFitterName;
  edm::ESHandle<TrajectoryFitter> theFitter;
  
  std::string theSmootherName;
  edm::ESHandle<TrajectorySmoother> theSmoother;
 
  RefitDirection::GeometricalDirection
    checkRecHitsOrdering(TransientTrackingRecHit::ConstRecHitContainer&) const;

  //  void reorder(TransientTrackingRecHit::ConstRecHitContainer& recHits) const;

  std::string theTrackerRecHitBuilderName;
  edm::ESHandle<TransientTrackingRecHitBuilder> theTrackerRecHitBuilder;
  
  std::string theMuonRecHitBuilderName;
  edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;
};
#endif

