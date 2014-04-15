#include "TrackingTools/TrackRefitter/interface/ME0TrackTransformer.h"

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/DetId/interface/DetId.h"

//#include "ME0Reconstruction/ME0Segment/interface/ME0MuonCollection.h"
#include "DataFormats/Math/interface/invertPosDefMatrix.h"

// Only for this one time!!!!
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;
using namespace edm;

/// Constructor
ME0TrackTransformer::ME0TrackTransformer(const ParameterSet& parameterSet){
  
  // Refit direction
  string refitDirectionName = parameterSet.getParameter<string>("RefitDirection");
  theRefitDirection = RefitDirection(refitDirectionName);
  
  theFitterName = parameterSet.getParameter<string>("Fitter");  
  theSmootherName = parameterSet.getParameter<string>("Smoother");  
  thePropagatorName = parameterSet.getParameter<string>("Propagator");
  theMuonPropagatorName = parameterSet.getParameter<string>("MuonPropagator");

  theTrackerRecHitBuilderName = parameterSet.getParameter<string>("TrackerRecHitBuilder");
  theMuonRecHitBuilderName = parameterSet.getParameter<string>("MuonRecHitBuilder");

  theRPCInTheFit = parameterSet.getParameter<bool>("RefitRPCHits");
  theDoPredictionsOnly = parameterSet.getParameter<bool>("DoPredictionsOnly");

  theCacheId_TC = theCacheId_GTG = theCacheId_MG = theCacheId_TRH = 0;
}

/// Destructor
ME0TrackTransformer::~ME0TrackTransformer(){}


void ME0TrackTransformer::setServices(const EventSetup& setup){
  
  const std::string metname = "Reco|TrackingTools|ME0TrackTransformer";

  setup.get<TrajectoryFitter::Record>().get(theFitterName,theFitter);
  setup.get<TrajectoryFitter::Record>().get(theSmootherName,theSmoother);

  
  unsigned long long newCacheId_TC = setup.get<TrackingComponentsRecord>().cacheIdentifier();

  if ( newCacheId_TC != theCacheId_TC ){
    LogTrace(metname) << "Tracking Component changed!";
    theCacheId_TC = newCacheId_TC;
    setup.get<TrackingComponentsRecord>().get(thePropagatorName,thePropagator);
    setup.get<TrackingComponentsRecord>().get(theMuonPropagatorName,theMuonPropagator); 
  }

  // Global Tracking Geometry
  unsigned long long newCacheId_GTG = setup.get<GlobalTrackingGeometryRecord>().cacheIdentifier();
  if ( newCacheId_GTG != theCacheId_GTG ) {
    LogTrace(metname) << "GlobalTrackingGeometry changed!";
    theCacheId_GTG = newCacheId_GTG;
    setup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry); 
  }
  
  // Magfield Field
  unsigned long long newCacheId_MG = setup.get<IdealMagneticFieldRecord>().cacheIdentifier();
  if ( newCacheId_MG != theCacheId_MG ) {
    LogTrace(metname) << "Magnetic Field changed!";
    theCacheId_MG = newCacheId_MG;
    setup.get<IdealMagneticFieldRecord>().get(theMGField);
  }
  
  // Transient Rechit Builders
  unsigned long long newCacheId_TRH = setup.get<TransientRecHitRecord>().cacheIdentifier();
  if ( newCacheId_TRH != theCacheId_TRH ) {
    theCacheId_TRH = newCacheId_TRH;
    LogTrace(metname) << "TransientRecHitRecord changed!";
    setup.get<TransientRecHitRecord>().get(theTrackerRecHitBuilderName,theTrackerRecHitBuilder);
    setup.get<TransientRecHitRecord>().get(theMuonRecHitBuilderName,theMuonRecHitBuilder);
  }
}


// vector<Trajectory> ME0TrackTransformer::transform(const reco::TrackRef& track) const {
//   return transform(*track);
// }


TransientTrackingRecHit::ConstRecHitContainer
ME0TrackTransformer::getTransientRecHits(const reco::TransientTrack& track) const {

  TransientTrackingRecHit::ConstRecHitContainer result;
  
  for (trackingRecHit_iterator hit = track.recHitsBegin(); hit != track.recHitsEnd(); ++hit) {
    if((*hit)->isValid()) {
      if ( (*hit)->geographicalId().det() == DetId::Tracker ) {
	result.push_back(theTrackerRecHitBuilder->build(&**hit));
      } else if ( (*hit)->geographicalId().det() == DetId::Muon ){
	if( (*hit)->geographicalId().subdetId() == 3 && !theRPCInTheFit){
	  LogTrace("Reco|TrackingTools|ME0TrackTransformer") << "RPC Rec Hit discarged"; 
	  continue;
	}
	result.push_back(theMuonRecHitBuilder->build(&**hit));
      }
    }
  }
  
  return result;
}

// FIXME: check this method!
RefitDirection::GeometricalDirection
ME0TrackTransformer::checkRecHitsOrdering(TransientTrackingRecHit::ConstRecHitContainer& recHits) const {
  
  if (!recHits.empty()){
    GlobalPoint first = trackingGeometry()->idToDet(recHits.front()->geographicalId())->position();
    GlobalPoint last = trackingGeometry()->idToDet(recHits.back()->geographicalId())->position();
    
    double rFirst = first.mag();
    double rLast  = last.mag();
    if(rFirst < rLast) return RefitDirection::insideOut;
    else if(rFirst > rLast) return RefitDirection::outsideIn;
    else{
      LogDebug("Reco|TrackingTools|ME0TrackTransformer") << "Impossible to determine the rechits order" <<endl;
      return RefitDirection::undetermined;
    }
  }
  else{
    LogDebug("Reco|TrackingTools|ME0TrackTransformer") << "Impossible to determine the rechits order" <<endl;
    return RefitDirection::undetermined;
  }
}

// void reorder(TransientTrackingRecHit::ConstRecHitContainer& recHits, RefitDirection::GeometricalDirection recHitsOrder) const{

//   if(theRefitDirection.geometricalDirection() != recHitsOrder) reverse(recHits.begin(),recHits.end());

//   if(theRefitDirection.geometricalDirection() == RefitDirection::insideOut &&recHitsOrder){}
//   else if(theRefitDirection.geometricalDirection() == RefitDirection::outsideIn){} 
//   else LogWarning("Reco|TrackingTools|ME0TrackTransformer") << "Impossible to determine the rechits order" <<endl;
// }


/// Convert Tracks into Trajectories
vector<Trajectory> ME0TrackTransformer::transform(const reco::Track& newTrack) const {
  LogWarning("Reco|TrackingTools|ME0TrackTransformer") << "Cannot use function "; 
  LogWarning("Reco|TrackingTools|ME0TrackTransformer") << "  vector<Trajectory> ME0TrackTransformer::transform(const reco::Track& newTrack) const "; 
  LogWarning("Reco|TrackingTools|ME0TrackTransformer") << "Re-implemented just because it's pure virtual in TrackTransformerBase. Skipping track. "; 
  return vector<Trajectory>(); 
}

/// Convert MEMuons into Trajectories
//vector<Trajectory> ME0TrackTransformer::transform(const reco::ME0Muon& oldMe0Muon) const {
vector<Trajectory> ME0TrackTransformer::transform(const reco::ME0Muon& oldMe0Muon, const Handle<reco::GenParticleCollection>& genmuoncoll) const {

  const reco::Track& oldTrack = *(oldMe0Muon.innerTrack()); 
  const std::string metname = "Reco|TrackingTools|ME0TrackTransformer";

  reco::TransientTrack track(oldTrack,magneticField(),trackingGeometry());   
  const EmulatedME0SegmentRef me0Segment = oldMe0Muon.me0segment(); 

  // Build the transient Rechits
  TransientTrackingRecHit::ConstRecHitContainer recHitsForReFit = getTransientRecHits(track);

  //return transform(track, recHitsForReFit, me0Segment);
  return transform(track, recHitsForReFit, me0Segment, genmuoncoll); 
}


/// Convert Tracks into Trajectories with a given set of hits
// vector<Trajectory> ME0TrackTransformer::transform(const reco::TransientTrack& track,
// 						  const TransientTrackingRecHit::ConstRecHitContainer& _recHitsForReFit,
// 						  const EmulatedME0SegmentRef me0Segm) const {
vector<Trajectory> ME0TrackTransformer::transform(const reco::TransientTrack& track,
						  const TransientTrackingRecHit::ConstRecHitContainer& _recHitsForReFit,
						  const EmulatedME0SegmentRef me0Segm,
						  const Handle<reco::GenParticleCollection>& genmuons) const {
  
  ////////////////////////////////
  // Just for this one time!!!! //
  ////////////////////////////////

//   Handle<reco::GenParticleCollection> genmuons;
//   event.getByLabel("genParticles", genmuons);
  
  unsigned int gensize = genmuons->size();

  unsigned int genidx = 999; 
  float minDeltaR = 999.;
  for(unsigned int i=0; i<gensize; ++i) {
    const reco::GenParticle& gp = (*genmuons)[i];
    if( gp.status()==1 && abs(gp.pdgId())==13 ) {
      float tmpDeltaR = deltaR(gp, track.track()); 
      if(tmpDeltaR<0.3 && tmpDeltaR<minDeltaR) { 
	minDeltaR = tmpDeltaR; 
	genidx = i; 
      }
//       std::cout << " ###### GenMuon " << i << ": " 
// 		<< gp.eta() << " / " 
// 		<< gp.phi() << " / " 
// 		<< gp.charge() << " / " 
// 		<< gp.pt() << std::endl; 
    }
  }

  TransientTrackingRecHit::ConstRecHitContainer recHitsForReFit =  _recHitsForReFit;
  const std::string metname = "Reco|TrackingTools|ME0TrackTransformer";

  if(recHitsForReFit.size() < 2) return vector<Trajectory>();
 
  // 8 cases are foreseen: 
  // [RH = rec hit order, P = momentum dir, FD = fit direction. IO/OI = inside-out/outside-in, AM/OM = along momentum/opposite to momentum]
  // (1) RH IO | P IO | FD AM  ---> Start from IN
  // (2) RH IO | P IO | FD OM  ---> Reverse RH and start from OUT
  // (3) RH IO | P OI | FD AM  ---> Reverse RH and start from IN
  // (4) RH IO | P OI | FD OM  ---> Start from OUT
  // (5) RH OI | P IO | FD AM  ---> Reverse RH and start from IN
  // (6) RH OI | P IO | FD OM  ---> Start from OUT
  // (7) RH OI | P OI | FD AM  ---> Start from IN
  // (8) RH OI | P OI | FD OM  ---> Reverse RH and start from OUT
  //
  // *** Rules: ***
  // -A- If RH-FD agree (IO-AM,OI-OM) do not reverse the RH
  // -B- If FD along momentum start from innermost state, otherwise use outermost
  
  // Other special cases can be handled:
  // (1 bis) RH IO | P IO | GFD IO => FD AM  ---> Start from IN
  // (2 bis) RH IO | P IO | GFD OI => FD OM  ---> Reverse RH and start from OUT
  // (3 bis) RH IO | P OI | GFD OI => FD AM  ---> Reverse RH and start from OUT
  // (4 bis) RH IO | P OI | GFD IO => FD OM  ---> Start from IN
  // (5 bis) RH OI | P IO | GFD IO => FD AM  ---> Reverse RH and start from IN
  // (6 bis) RH OI | P IO | GFD OI => FD OM  ---> Start from OUT
  // (7 bis) RH OI | P OI | GFD OI => FD AM  ---> Start from OUT
  // (8 bis) RH OI | P OI | GFD IO => FD OM  ---> Reverse RH and start from IN
  // 
  // *** Additional rule: ***
  // -A0- If P and GFD agree, then FD is AM otherwise is OM
  // -A00- rechit must be ordered as GFD in order to handle the case of cosmics
  // -B0- The starting state is decided by GFD

  // Determine the RH order
  RefitDirection::GeometricalDirection recHitsOrder = checkRecHitsOrdering(recHitsForReFit); // FIXME change nome of the *type*  --> RecHit order!
  LogTrace(metname) << "RH order (0-insideOut, 1-outsideIn): " << recHitsOrder;

  PropagationDirection propagationDirection = theRefitDirection.propagationDirection();

  // Apply rule -A0-
  if(propagationDirection == anyDirection){
    GlobalVector momentum = track.innermostMeasurementState().globalMomentum();
    GlobalVector position = track.innermostMeasurementState().globalPosition() - GlobalPoint(0,0,0);
    RefitDirection::GeometricalDirection p = (momentum.x()*position.x() > 0 || momentum.y()*position.y() > 0) ? RefitDirection::insideOut : RefitDirection::outsideIn;

    propagationDirection = p == theRefitDirection.geometricalDirection() ? alongMomentum : oppositeToMomentum;
    LogTrace(metname) << "P  (0-insideOut, 1-outsideIn): " << p;
    LogTrace(metname) << "FD (0-OM, 1-AM, 2-ANY): " << propagationDirection;
  }
  // -A0-

  // Apply rule -A-
  if(theRefitDirection.propagationDirection() != anyDirection){
    if((recHitsOrder == RefitDirection::insideOut && propagationDirection == oppositeToMomentum) ||
       (recHitsOrder == RefitDirection::outsideIn && propagationDirection == alongMomentum) ) 
      reverse(recHitsForReFit.begin(),recHitsForReFit.end());}
  // -A-
  // Apply rule -A00-
  else{
    // reorder the rechit as defined in theRefitDirection.geometricalDirection(); 
    if(theRefitDirection.geometricalDirection() != recHitsOrder) reverse(recHitsForReFit.begin(),recHitsForReFit.end()); 
  }
  // -A00-
    
  // Apply rule -B-
  TrajectoryStateOnSurface firstTSOS = track.innermostMeasurementState();
  unsigned int innerId = track.track().innerDetId();
  if(theRefitDirection.propagationDirection() != anyDirection){
    if(propagationDirection == oppositeToMomentum){
      innerId   = track.track().outerDetId();
      firstTSOS = track.outermostMeasurementState();
    }
  }
  else { // if(theRefitDirection.propagationDirection() == anyDirection)
    // Apply rule -B0-
    if(theRefitDirection.geometricalDirection() == RefitDirection::outsideIn){
      innerId   = track.track().outerDetId();
      firstTSOS = track.outermostMeasurementState();
    }
    // -B0-
  }
  // -B-

  if(!firstTSOS.isValid()){
    LogTrace(metname)<<"Error wrong initial state!"<<endl;
    return vector<Trajectory>();
  }

  TrajectorySeed seed(PTrajectoryStateOnDet(),TrajectorySeed::recHitContainer(),propagationDirection);

  if(recHitsForReFit.front()->geographicalId() != DetId(innerId)){
    LogTrace(metname)<<"Propagation occured"<<endl;
    firstTSOS = propagator()->propagate(firstTSOS, recHitsForReFit.front()->det()->surface());
    if(!firstTSOS.isValid()){
      LogTrace(metname)<<"Propagation error!"<<endl;
      return vector<Trajectory>();
    }
  }

  if(theDoPredictionsOnly){
    Trajectory aTraj(seed,propagationDirection);
    TrajectoryStateOnSurface predTSOS = firstTSOS;
    for(TransientTrackingRecHit::ConstRecHitContainer::const_iterator ihit = recHitsForReFit.begin(); 
	ihit != recHitsForReFit.end(); ++ihit ) {
      predTSOS = propagator()->propagate(predTSOS, (*ihit)->det()->surface());
      if (predTSOS.isValid()) aTraj.push(TrajectoryMeasurement(predTSOS, *ihit));
    }
    return vector<Trajectory>(1, aTraj);
  }


  vector<Trajectory> trajectories = theFitter->fit(seed,recHitsForReFit,firstTSOS);
  
  if(trajectories.empty()){
    LogTrace(metname)<<"No Track refitted!"<<endl;
    return vector<Trajectory>();
  }

  // -----------------------------------------------------------------------------------------
  //  Shot at refitting with ME0 segment (dissecting KFUpdator, adding a touch of dark magic) 
  // -----------------------------------------------------------------------------------------

  // First, propagate track to ME0 (taken from SegmentMatcher producer) 
  Trajectory trajectoryTrkOnly = trajectories.front(); 
  FreeTrajectoryState *ftsTrkOnly = trajectoryTrkOnly.lastMeasurement().updatedState().freeState(); 
  GlobalVector momTrkOnly = ftsTrkOnly->momentum(); 
  GlobalVector posTrkOnly( ftsTrkOnly->parameters().position().x(), ftsTrkOnly->parameters().position().y(), ftsTrkOnly->parameters().position().z() ); 
  float zSign  = momTrkOnly.z()/fabs(momTrkOnly.z());
  float zValue = 560. * zSign;
  Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());
  // Getting the initial variables for propagation
  int chargeReco = ftsTrkOnly->charge(); 
  //AlgebraicSymMatrix66 covReco = ftsTrkOnly->hasError() ? ftsTrkOnly->cartesianError().matrix() : AlgebraicSymMatrix66(); 

  // Now we propagate and get the propagated variables from the propagated state
  SteppingHelixStateInfo startrecostate(*ftsTrkOnly); 
  const SteppingHelixPropagator* theShPropagator = dynamic_cast<const SteppingHelixPropagator*>(&*theMuonPropagator);
  if(theShPropagator==0) { 
    LogWarning("Reco|TrackingTools|ME0TrackTransformer") << "No SteppingHelixPropagator! Skipping track. "; 
    return vector<Trajectory>(); 
  } 
  SteppingHelixStateInfo lastrecostate = theShPropagator->propagate(startrecostate, *plane); 
  FreeTrajectoryState ftsTrkOnlyAtMu; 
  lastrecostate.getFreeState(ftsTrkOnlyAtMu); 
  AlgebraicSymMatrix66 covFinalReco = ftsTrkOnlyAtMu.hasError() ? ftsTrkOnlyAtMu.cartesianError().matrix() : AlgebraicSymMatrix66(); 
  double x_glb = ftsTrkOnlyAtMu.position().x(); 
  double y_glb = ftsTrkOnlyAtMu.position().y(); 
  double z_glb = ftsTrkOnlyAtMu.position().z(); 
  //double r_glb = sqrt(x_glb*x_glb + y_glb*y_glb + z_glb*z_glb); 
  GlobalVector posTrkOnlyAtMu(x_glb, y_glb, z_glb); 
  double px_glb = ftsTrkOnlyAtMu.momentum().x(); 
  double py_glb = ftsTrkOnlyAtMu.momentum().y(); 
  double pz_glb = ftsTrkOnlyAtMu.momentum().z(); 
  double p_glb  = sqrt(px_glb*px_glb + py_glb*py_glb + pz_glb*pz_glb); 
  GlobalVector momTrkOnlyAtMu(px_glb, py_glb, pz_glb); 
  chargeReco = ftsTrkOnlyAtMu.charge(); 

  ////////////////////////////////
  // Just for this one time!!!! //
  ////////////////////////////////

  GlobalPoint genPoint((*genmuons)[genidx].vx(), (*genmuons)[genidx].vy(), (*genmuons)[genidx].vz()); 
  GlobalVector genVector((*genmuons)[genidx].px(), (*genmuons)[genidx].py(), (*genmuons)[genidx].pz()); 
  FreeTrajectoryState genfts(genPoint, genVector, (*genmuons)[genidx].charge(), magneticField()); 
  SteppingHelixStateInfo startgenstate(genfts); 
  SteppingHelixStateInfo lastgenstate = theShPropagator->propagate(startgenstate, *plane); 
  FreeTrajectoryState genftsAtMu; 
  lastgenstate.getFreeState(genftsAtMu); 
  GlobalVector genMomAtMu = genftsAtMu.momentum(); 

  ////////////////////////////////
  //            END             //
  ////////////////////////////////

  // Here comes the fun
  // Global  coordinates (6 coord.):  (x, y, x, px, py, pz)
  // "Local" coordinates (5 coord.):  (q/p, dx/dz, dy/dz, x, y)
  // Jacobian must transform covariance from global to "local":  6x6 -> 5x5
  //  so Jacobian must be 5x6:  C(loc) =  J  * C(glb) * J^T
  //                              5x5  = 5x6 *  6x6   * 6x5
  // 
  // 0: q/p,  1: dx/dz,  2: dy/dz,  3: x,    4: y
  // 0: x,    1: y,      2: z,      3: px,   4: py,   5: pz
  //AlgebraicMatrix jacobGlbToLoc(5, 6, 0); 
  AlgebraicMatrix56 jacobGlbToLoc; 
  jacobGlbToLoc[0][0] = 0.0;                                               // d(q/p) / d(x)
  jacobGlbToLoc[0][1] = 0.0;                                               // d(q/p) / d(y)
  jacobGlbToLoc[0][2] = 0.0;                                               // d(q/p) / d(z)
  jacobGlbToLoc[0][3] = (-1.) * chargeReco * px_glb / (p_glb*p_glb*p_glb); // d(q/p) / d(px)
  jacobGlbToLoc[0][4] = (-1.) * chargeReco * py_glb / (p_glb*p_glb*p_glb); // d(q/p) / d(py)
  jacobGlbToLoc[0][5] = (-1.) * chargeReco * pz_glb / (p_glb*p_glb*p_glb); // d(q/p) / d(pz)

  jacobGlbToLoc[1][0] = 0.0;                                               // d(px/pz) / d(x)
  jacobGlbToLoc[1][1] = 0.0;                                               // d(px/pz) / d(y)
  jacobGlbToLoc[1][2] = 0.0;                                               // d(px/pz) / d(z)
  jacobGlbToLoc[1][3] = 1. / pz_glb;                                       // d(px/pz) / d(px)
  jacobGlbToLoc[1][4] = 0.0;                                               // d(px/pz) / d(py)
  jacobGlbToLoc[1][5] = (-1.) * px_glb / (pz_glb*pz_glb);                  // d(px/pz) / d(pz)

  jacobGlbToLoc[2][0] = 0.0;                                               // d(py/pz) / d(x)
  jacobGlbToLoc[2][1] = 0.0;                                               // d(py/pz) / d(y)
  jacobGlbToLoc[2][2] = 0.0;                                               // d(py/pz) / d(z)
  jacobGlbToLoc[2][3] = 0.0;                                               // d(py/pz) / d(px)
  jacobGlbToLoc[2][4] = 1. / pz_glb;                                       // d(py/pz) / d(py)
  jacobGlbToLoc[2][5] = (-1.) * px_glb / (pz_glb*pz_glb);                  // d(py/pz) / d(pz)

  jacobGlbToLoc[3][0] = 1.0;                                               // d(x) / d(x)
  jacobGlbToLoc[3][1] = 0.0;                                               // d(x) / d(y)
  jacobGlbToLoc[3][2] = 0.0;                                               // d(x) / d(z)
  jacobGlbToLoc[3][3] = 0.0;                                               // d(x) / d(px)
  jacobGlbToLoc[3][4] = 0.0;                                               // d(x) / d(py)
  jacobGlbToLoc[3][5] = 0.0;                                               // d(x) / d(pz)

  jacobGlbToLoc[4][0] = 0.0;                                               // d(y) / d(x)
  jacobGlbToLoc[4][1] = 1.0;                                               // d(y) / d(y)
  jacobGlbToLoc[4][2] = 0.0;                                               // d(y) / d(z)
  jacobGlbToLoc[4][3] = 0.0;                                               // d(y) / d(px)
  jacobGlbToLoc[4][4] = 0.0;                                               // d(y) / d(py)
  jacobGlbToLoc[4][5] = 0.0;                                               // d(y) / d(pz)

  // Predicted state in "local" coordinates 
  AlgebraicVector5 x( chargeReco/p_glb,  px_glb/pz_glb,  py_glb/pz_glb, x_glb, y_glb ); 
  // Predicted state covariance matrix in "local" coordinates 
  //AlgebraicMatrix Ctmp = (jacobGlbToLoc * covFinalReco); 
  //AlgebraicSymMatrix55 C =  Ctmp * jacobGlbToLocTransp; 
  //AlgebraicSymMatrix55 C =  (jacobGlbToLoc * covFinalReco) * ROOT::Math::Transpose(jacobGlbToLoc); 
  //AlgebraicMatrix56 Ctmp = jacobGlbToLoc * covFinalReco; 
  //AlgebraicSymMatrix55 C = Ctmp * ROOT::Math::Transpose(jacobGlbToLoc); 

  AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * covFinalReco) * ROOT::Math::Transpose(jacobGlbToLoc); 
  AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
  for(int i=0; i<5; ++i) {
    for(int j=0; j<5; ++j) {
      C[i][j] = Ctmp[i][j]; 
    }
  }  

  // Projection matrix: 4x5
  AlgebraicMatrix Htmp = me0Segm->projectionMatrix(); 
  AlgebraicMatrix45 H; // I couldn't find any other way, so I resort to the brute force
  for(int i=0; i<Htmp.num_row(); ++i) {
    for(int j=0; j<Htmp.num_col(); ++j) {
      H[i][j] = Htmp[i][j]; 
    }
  }

  // Predicted state in "local" coordinates,  dim = 4 
  AlgebraicVector4 r = H * x; 
  // Predicted state covariance matrix in "local" coordinates,  dim = 4x4 
  //AlgebraicSymMatrix44 V = (H * C) * ROOT::Math::Transpose(H); 
  AlgebraicMatrix44 Vtmp = (H * C) * ROOT::Math::Transpose(H); 
  AlgebraicSymMatrix44 V;  // I couldn't find any other way, so I resort to the brute force
  for(int i=0; i<4; ++i) {
    for(int j=0; j<4; ++j) {
      V[i][j] = Vtmp[i][j]; 
    }
  }  

  // Measured state in "local" coordinates,  dim = 4 
  AlgebraicVector r_meas_tmp = me0Segm->parameters(); 
  AlgebraicVector4 r_meas; 
  // Measured state covariance matrix in "local" coordinates,  dim = 4x4 
  AlgebraicSymMatrix V_meas_tmp = me0Segm->parametersError(); 
  //AlgebraicSymMatrix44 V_meas; 
  AlgebraicSymMatrix44 V_meas; 

  for(int i=0; i<V_meas_tmp.num_row(); ++i) {
    r_meas[i] = r_meas_tmp[i]; 
    for(int j=0; j<V_meas_tmp.num_col(); ++j) {
      V_meas[i][j] = V_meas_tmp[i][j]; 
    }
  }

  // Math taken from https://github.com/cms-sw/cmssw/blob/CMSSW_6_2_X_SLHC/TrackingTools/KalmanUpdators/src/KFUpdator.cc
  r -= r_meas;
  AlgebraicSymMatrix44 R = V + V_meas;  // sum of covariance matrices, still symmetric and invertible 
  //AlgebraicMatrix44 R = V + V_meas;  // sum of covariance matrices, still symmetric and invertible 
  bool ok = invertPosDefMatrix(R); 
  if(!ok) {
    LogWarning("Reco|TrackingTools|ME0TrackTransformer") << "Inversion of total covariance matrix (predicted state + ME0 segment)! Skipping track. "; 
    return vector<Trajectory>(); 
  }

  // Kalman gain matrix
  AlgebraicMatrix54 K;
  K = (C * ROOT::Math::Transpose(H)) * R; 

  AlgebraicSymMatrix55 M = AlgebraicMatrixID(); 
  //M -= K * H; 

  // Filtered state in "local" coordinates 
  AlgebraicVector5 fsv = x + K * r; 
  // Filtered state covariance matrix in "local" coordinates 
  //AlgebraicSymMatrix55 fse = ROOT::Math::Similarity(M, C) + ROOT::Math::Similarity(K, V); 
  AlgebraicMatrix55 fse = ROOT::Math::Similarity(M, C) + ROOT::Math::Similarity(K, V); 

  //const double pi = 3.14159265358979312; 
  double qbp = fsv[0]; 
  double dxdz = fsv[1]; 
  double dydz = fsv[2]; 
  double p_upd = fabs(1./qbp); 
  double charge_upd = qbp/fabs(qbp); 
  double pz_upd = zSign * p_upd / sqrt(dxdz*dxdz + dydz*dydz + 1); 
  double px_upd = pz_upd * dxdz; 
  double py_upd = pz_upd * dydz; 

  GlobalPoint posUpdatedAtMu(fsv[3], fsv[4], z_glb); 
  GlobalVector momUpdatedAtMu(px_upd, py_upd, pz_upd); 

  double q = charge_upd*zSign;
  double sqr = sqrt(dxdz*dxdz + dydz*dydz + 1);
  double den = -q/(sqr*sqr*sqr*qbp);

  // 0: x,    1: y,      2: z,      3: px,   4: py,   5: pz
  // 0: q/p,  1: dx/dz,  2: dy/dz,  3: x,    4: y
  AlgebraicMatrix65 jacobLocToGlb; 
  jacobLocToGlb[0][0] = 0.0;                            // d(x) / d(q/p)
  jacobLocToGlb[0][1] = 0.0;                            // d(x) / d(dx/dz)
  jacobLocToGlb[0][2] = 0.0;                            // d(x) / d(dy/dz)
  jacobLocToGlb[0][3] = 1.0;                            // d(x) / d(x)
  jacobLocToGlb[0][4] = 0.0;                            // d(x) / d(y)    

  jacobLocToGlb[1][0] = 0.0;                            // d(y) / d(q/p)	
  jacobLocToGlb[1][1] = 0.0;                            // d(y) / d(dx/dz)
  jacobLocToGlb[1][2] = 0.0;                            // d(y) / d(dy/dz)
  jacobLocToGlb[1][3] = 0.0;                            // d(y) / d(x)	
  jacobLocToGlb[1][4] = 1.0;                            // d(y) / d(y)    

  jacobLocToGlb[2][0] = 0.0;                            // d(z) / d(q/p)	
  jacobLocToGlb[2][1] = 0.0;                            // d(z) / d(dx/dz)
  jacobLocToGlb[2][2] = 0.0;                            // d(z) / d(dy/dz)
  jacobLocToGlb[2][3] = 0.0;                            // d(z) / d(x)	
  jacobLocToGlb[2][4] = 0.0;                            // d(z) / d(y)    

  jacobLocToGlb[3][0] = dxdz*(-q/(sqr*qbp*qbp));        // d(px) / d(q/p)	
  jacobLocToGlb[3][1] = q/(sqr*qbp) + (den*dxdz*dxdz);  // d(px) / d(dx/dz)
  jacobLocToGlb[3][2] = den*dxdz*dydz;                  // d(px) / d(dy/dz)
  jacobLocToGlb[3][3] = 0.0;                            // d(px) / d(x)	
  jacobLocToGlb[3][4] = 0.0;                            // d(px) / d(y)    

  jacobLocToGlb[4][0] = dydz*(-q/(sqr*qbp*qbp));        // d(py) / d(q/p)	
  jacobLocToGlb[4][1] = den*dxdz*dydz;                  // d(py) / d(dx/dz)
  jacobLocToGlb[4][2] = q/(sqr*qbp) + (den*dydz*dydz);  // d(py) / d(dy/dz)
  jacobLocToGlb[4][3] = 0.0;                            // d(py) / d(x)	
  jacobLocToGlb[4][4] = 0.0;                            // d(py) / d(y)    

  jacobLocToGlb[5][0] = -q/(sqr*qbp*qbp);               // d(pz) / d(q/p)	
  jacobLocToGlb[5][1] = den*dxdz;                       // d(pz) / d(dx/dz)
  jacobLocToGlb[5][2] = den*dydz;                       // d(pz) / d(dy/dz)
  jacobLocToGlb[5][3] = 0.0;                            // d(pz) / d(x)	
  jacobLocToGlb[5][4] = 0.0;                            // d(pz) / d(y)    

  //AlgebraicSymMatrix66 covRecoUpdated = (jacobLocToGlb * fse) * ROOT::Math::Transpose(jacobLocToGlb); 
  //AlgebraicMatrix65 covRecoUpdatedTmp = (jacobLocToGlb * fse); 
  //AlgebraicSymMatrix66 covRecoUpdated = covRecoUpdatedTmp * ROOT::Math::Transpose(jacobLocToGlb); 

  AlgebraicMatrix66 covRecoUpdatedTmp = (jacobLocToGlb * fse) * ROOT::Math::Transpose(jacobLocToGlb); 
  AlgebraicSymMatrix66 covRecoUpdated; 
  for(int i=0; i<6; ++i) {
    r_meas[i] = r_meas_tmp[i]; 
    for(int j=0; j<6; ++j) {
      covRecoUpdated[i][j] = covRecoUpdatedTmp[i][j]; 
    }
  }  

  GlobalTrajectoryParameters updPars(posUpdatedAtMu, momUpdatedAtMu, charge_upd, magneticField());
  CartesianTrajectoryError updCov(covRecoUpdated);
  FreeTrajectoryState ftsUpdatedAtMu(updPars, updCov); 
  SteppingHelixStateInfo murecorecostateUpdated(ftsUpdatedAtMu); 
  //SteppingHelixStateInfo lastrecostateuUpdated = theMuonPropagator->propagate(murecorecostateUpdated, trajectoryTrkOnly.lastMeasurement().recHit()->surface()); 
  SteppingHelixStateInfo lastrecostateuUpdated = theShPropagator->propagate(murecorecostateUpdated, *(trajectoryTrkOnly.lastMeasurement().recHit()->surface())); 
  TrajectoryStateOnSurface tsosTrkPlusMuUpdated = lastrecostateuUpdated.getStateOnSurface(*(trajectoryTrkOnly.lastMeasurement().recHit()->surface())); 
  TrajectoryMeasurement newTM(tsosTrkPlusMuUpdated, tsosTrkPlusMuUpdated, trajectoryTrkOnly.lastMeasurement().recHit(), trajectoryTrkOnly.lastMeasurement().estimate()); 

  Trajectory trajectoryBW = trajectoryTrkOnly; 
  trajectoryBW.pop(); 
  trajectoryBW.push(newTM); 

  // ------------------------------------------------------------------------------
  //  End of embarassing "refit" with ME0 segment -- you can now relax and breathe
  // ------------------------------------------------------------------------------
  
  ////////////////////////////////
  // Just for this one time!!!! //
  ////////////////////////////////

  // Comparing genMomAtMu, momTrkOnlyAtMu, momUpdatedAtMu 
  std::cout << " ###### Generated Muon " << genidx << "/" << gensize 
	    << ": GEN/" 
	    << genftsAtMu.charge()      << "/" << genMomAtMu.eta()     << "/" << genMomAtMu.phi()     << "/" << genMomAtMu.mag() 
	    << " PRD/" 
	    << ftsTrkOnlyAtMu.charge()  << "/" << momTrkOnlyAtMu.eta() << "/" << momTrkOnlyAtMu.phi() << "/" << momTrkOnlyAtMu.mag() 
	    << " UPD/" 
	    << charge_upd               << "/" << momUpdatedAtMu.eta() << "/" << momUpdatedAtMu.phi() << "/" << momUpdatedAtMu.mag() 
	    << std::endl; 

  ////////////////////////////////
  //            END             //
  ////////////////////////////////


  //Trajectory trajectoryBW = trajectories.front();
    
  vector<Trajectory> trajectoriesSM = theSmoother->trajectories(trajectoryBW);

  if(trajectoriesSM.empty()){
    LogTrace(metname)<<"No Track smoothed!"<<endl;
    return vector<Trajectory>();
  }
  
  return trajectoriesSM;

}


