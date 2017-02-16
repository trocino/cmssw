/** \file ME0SegmentMatcher.cc
 *
 * \author David Nash
 */


#include <FWCore/PluginManager/interface/ModuleDef.h>
#include <FWCore/Framework/interface/MakerMacros.h>

#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h> 

#include <DataFormats/MuonReco/interface/ME0Muon.h>

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "DataFormats/GeometrySurface/interface/Plane.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixStateInfo.h"
#include "DataFormats/Math/interface/deltaR.h"


#include "DataFormats/GeometrySurface/interface/LocalError.h"


#include "TLorentzVector.h"

#include "DataFormats/TrajectoryState/interface/LocalTrajectoryParameters.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianCartesianToLocal.h"
#include "TrackingTools/AnalyticalJacobians/interface/JacobianLocalToCartesian.h"

#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include <DataFormats/GeometrySurface/interface/SimpleDiskBounds.h>


/** \class ME0SegmentMatcher 
 *
 * \author David Nash
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/GeometryVector/interface/GlobalVector.h"

#include "DataFormats/Math/interface/AlgebraicROOTObjects.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "Geometry/GEMGeometry/interface/ME0Geometry.h"
#include <Geometry/GEMGeometry/interface/ME0EtaPartition.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <DataFormats/MuonDetId/interface/ME0DetId.h>

#include "FWCore/ServiceRegistry/interface/Service.h"

#include <DataFormats/GEMRecHit/interface/ME0SegmentCollection.h>
#include <DataFormats/CSCRecHit/interface/CSCSegmentCollection.h>

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"




class FreeTrajectoryState;
class MagneticField;
class ME0SegmentMatcher : public edm::stream::EDProducer<> {
public:
  /// Constructor
  explicit ME0SegmentMatcher(const edm::ParameterSet&);
  /// Destructor
  ~ME0SegmentMatcher();
  /// Produce the ME0Segment collection
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

    
  virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;



  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
			     int , const AlgebraicSymMatrix66& ,
			     const MagneticField* );

  FreeTrajectoryState getFTS(const GlobalVector& , const GlobalVector& , 
			     int , const AlgebraicSymMatrix55& ,
			     const MagneticField* );

  void getFromFTS(const FreeTrajectoryState& ,
		  GlobalVector& , GlobalVector& , 
		  int& , AlgebraicSymMatrix66& );

private:



  double theX_RESIDUAL_CUT, theX_PULL_CUT, theY_RESIDUAL_CUT, theY_PULL_CUT, thePHIDIR_RESIDUAL_CUT;
  edm::InputTag ourSegmentsTag, generalTracksTag;
  edm::EDGetTokenT<ME0SegmentCollection> ourSegmentsToken_;
  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;

  TrackDetectorAssociator trackAssociator_; 
  TrackAssociatorParameters parameters_; 
};


ME0SegmentMatcher::ME0SegmentMatcher(const edm::ParameterSet& pas) {
  produces<std::vector<reco::ME0Muon> >();  
  theX_PULL_CUT   = pas.getParameter<double>("maxPullX");
  theX_RESIDUAL_CUT   = pas.getParameter<double>("maxDiffX");
  theY_PULL_CUT   = pas.getParameter<double>("maxPullY");
  theY_RESIDUAL_CUT   = pas.getParameter<double>("maxDiffY");
  thePHIDIR_RESIDUAL_CUT   = pas.getParameter<double>("maxDiffPhiDirection");
  //Might need to replace "ourSegments" with an edm::InputTag of "ourSegments"
  ourSegmentsTag = pas.getParameter<edm::InputTag>("me0SegmentTag");
  generalTracksTag = pas.getParameter<edm::InputTag>("tracksTag");
  ourSegmentsToken_ = consumes<ME0SegmentCollection>(ourSegmentsTag);
  generalTracksToken_ = consumes<reco::TrackCollection>(generalTracksTag);

   // Load TrackDetectorAssociator parameters
   edm::ParameterSet parameters = pas.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
   edm::ConsumesCollector iC = consumesCollector();
   parameters_.loadParameters( parameters, iC );

}

ME0SegmentMatcher::~ME0SegmentMatcher() {}

void ME0SegmentMatcher::produce(edm::Event& ev, const edm::EventSetup& setup) {

    //Getting the objects we'll need
    using namespace edm;

    ESHandle<ME0Geometry> me0Geom;
    setup.get<MuonGeometryRecord>().get(me0Geom);
    ESHandle<MagneticField> bField;
    setup.get<IdealMagneticFieldRecord>().get(bField);
    ESHandle<Propagator> thisShProp;
    setup.get<TrackingComponentsRecord>().get("SteppingHelixPropagatorAny", thisShProp);
    trackAssociator_.setPropagator(thisShProp.product());

    using namespace reco;

    Handle<ME0SegmentCollection> ourSegments;
    //ev.getByLabel("me0Segments","",ourSegments);
    ev.getByToken(ourSegmentsToken_,ourSegments);

    auto oc = std::make_unique<std::vector<ME0Muon>>(); 
    std::vector<ME0Muon> TempStore; 

    Handle <TrackCollection > generalTracks;
    //ev.getByLabel <TrackCollection> ("generalTracks", generalTracks);
    ev.getByToken(generalTracksToken_,generalTracks);


    int TrackNumber = 0;
    
    for (std::vector<Track>::const_iterator thisTrack = generalTracks->begin();
	 thisTrack != generalTracks->end(); ++thisTrack,++TrackNumber) {
      //Initializing our plane

      //Remove later
      if (std::abs(thisTrack->eta()) < 1.8) continue;

      float zSign = thisTrack->pz() > 0 ? 1.0f : -1.0f;

      // const float zValue = 526.75 * zSign;
      // Plane *plane = new Plane(Surface::PositionType(0,0,zValue),Surface::RotationType());

      // Track-detector association as in RecoMuon/MuonIdentification/plugins/MuonIdProducer 
      TrackDetMatchInfo info = trackAssociator_.associate(ev, setup, *thisTrack, parameters_, 
							  //TrackDetectorAssociator::InsideOut); 
							  TrackDetectorAssociator::Any); 

      //Getting the initial variables for propagation

      int chargeReco = thisTrack->charge(); 
      // GlobalVector p3reco, r3reco;

      // p3reco = GlobalVector(thisTrack->outerPx(), thisTrack->outerPy(), thisTrack->outerPz());
      // r3reco = GlobalVector(thisTrack->outerX(), thisTrack->outerY(), thisTrack->outerZ());

      // AlgebraicSymMatrix66 covReco;
      // //This is to fill the cov matrix correctly
      // AlgebraicSymMatrix55 covReco_curv;
      // covReco_curv = thisTrack->outerStateCovariance();
      // FreeTrajectoryState initrecostate = getFTS(p3reco, r3reco, chargeReco, covReco_curv, &*bField);
      // getFromFTS(initrecostate, p3reco, r3reco, chargeReco, covReco);

      // std::cout << "[TMP:DT] Track " << TrackNumber << " (eta, phi, pt) = (" 
      // 		<< thisTrack->eta() << ", " << thisTrack->phi() << ", " 
      // 		<< thisTrack->pt() << ")" << std::endl; 

      for (const auto& achamber : info.chambers) { 
	//std::cout << "[TMP:DT]   New Chamber: " << achamber.id.subdetId() << std::endl; 
	//if(achamber.id.subdetId()!=5) continue; 
	if(achamber.id.subdetId()!=MuonSubdetId::ME0) continue; 
	//const GeomDet* geomDet = muonDetIdAssociator_->getGeomDet((*matchedChamber).id);
	
	//std::cout << "[TMP:DT]     New Chamber: ME0" << std::endl; 

	//Now we propagate and get the propagated variables from the propagated state
	//SteppingHelixStateInfo startrecostate(initrecostate);
	//SteppingHelixStateInfo lastrecostate;
	//TrajectoryStateOnSurface lastrecostate;
	TrajectoryStateOnSurface lastrecostate = achamber.tState;

	//const SteppingHelixPropagator* thisShProp = 
	//dynamic_cast<const SteppingHelixPropagator*>(&*shProp);

	//lastrecostate = thisShProp->propagate(startrecostate, *plane);
	//lastrecostate = thisShProp->propagateWithPath(startrecostate, *plane);
	//thisShProp->propagate(startrecostate, *plane,lastrecostate);
	// lastrecostate = thisShProp->propagate(initrecostate,*plane);
	if (!lastrecostate.isValid()) continue;
	//std::cout << "[TMP:DT]     TSOS valid" << std::endl; 

	FreeTrajectoryState finalrecostate(*lastrecostate.freeTrajectoryState());
	//lastrecostate.getFreeState(finalrecostate);
	//finalrecostate = lastrecostate.freeTrajectoryState();

	AlgebraicSymMatrix66 covFinalReco;
	GlobalVector p3FinalReco_glob, r3FinalReco_globv;
	getFromFTS(finalrecostate, p3FinalReco_glob, r3FinalReco_globv, chargeReco, covFinalReco);


	//To transform the global propagated track to local coordinates
	//int SegmentNumber = 0;

	reco::ME0Muon MuonCandidate;
	double ClosestDelR2 = 999.;

	//for (auto thisSegment = ourSegments->begin(); thisSegment != ourSegments->end(); 
	//     ++thisSegment,++SegmentNumber) { 
	for(const auto& asegment : achamber.segments) { 
	  //ME0SegmentRef thisSegment = asegment.me0SegmentRef; 
	  //std::cout << "[TMP:DT]     New Segment" << std::endl; 
	  const ME0Segment *thisSegment = asegment.me0SegmentRef.get(); 
	  //ME0DetId id(achamber.id);
	  auto chamber = me0Geom->chamber(ME0DetId(achamber.id)); 
	  if(zSign*chamber->toGlobal(thisSegment->localPosition()).z()<0) continue; 
	  //std::cout << "[TMP:DT]     Correct Sign/Direction" << std::endl; 

	  GlobalPoint r3FinalReco_glob(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z()); 
	  LocalPoint r3FinalReco = chamber->toLocal(r3FinalReco_glob);
	  LocalVector p3FinalReco=chamber->toLocal(p3FinalReco_glob);

	  LocalPoint thisPosition(thisSegment->localPosition());
	  LocalVector thisDirection(thisSegment->localDirection().x(),thisSegment->localDirection().y(),thisSegment->localDirection().z());  //FIXME

	  //The same goes for the error
	  AlgebraicMatrix thisCov(4,4,0);   
	  for (int i = 1; i <=4; i++){
	    for (int j = 1; j <=4; j++){
	      thisCov(i,j) = thisSegment->parametersError()(i,j);
	    }
	  }

	  /////////////////////////////////////////////////////////////////////////////////////////


	  LocalTrajectoryParameters ltp(r3FinalReco,p3FinalReco,chargeReco);
	  JacobianCartesianToLocal jctl(chamber->surface(),ltp);
	  AlgebraicMatrix56 jacobGlbToLoc = jctl.jacobian(); 

	  AlgebraicMatrix55 Ctmp =  (jacobGlbToLoc * covFinalReco) * ROOT::Math::Transpose(jacobGlbToLoc); 
	  AlgebraicSymMatrix55 C;  // I couldn't find any other way, so I resort to the brute force
	  for(int i=0; i<5; ++i) {
	    for(int j=0; j<5; ++j) {
	      C[i][j] = Ctmp[i][j]; 
	    }
	  } 

	  Double_t sigmax = sqrt(C[3][3]+thisSegment->localPositionError().xx()); 
	  Double_t sigmay = sqrt(C[4][4]+thisSegment->localPositionError().yy()); 
	  bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false; 
	
	  // std::cout << "[TMP:DT]     X: " << std::abs(thisPosition.x()-r3FinalReco.x()) 
	  // 	    << " >< " << theX_PULL_CUT << " * " << sigmax << " = " << theX_PULL_CUT*sigmax 
	  // 	    << " || " << theX_RESIDUAL_CUT << std::endl; 
	  // std::cout << "[TMP:DT]     Y: " << std::abs(thisPosition.x()-r3FinalReco.y()) 
	  // 	    << " >< " << theY_PULL_CUT << " * " << sigmay << " = " << theY_PULL_CUT*sigmay 
	  // 	    << " || " << theY_RESIDUAL_CUT << std::endl; 
	  // std::cout << "[TMP:DT]     Phi: " << reco::deltaPhi(p3FinalReco_glob.barePhi(), chamber->toGlobal(thisSegment->localDirection()).barePhi()) 
	  // 	    << " || " << thePHIDIR_RESIDUAL_CUT << std::endl; 
	  if( (std::abs(thisPosition.x()-r3FinalReco.x()) < theX_PULL_CUT*sigmax) || 
	      (std::abs(thisPosition.x()-r3FinalReco.x()) < theX_RESIDUAL_CUT   )    ) X_MatchFound = true; 

	  if(theY_PULL_CUT<0. && theY_RESIDUAL_CUT<0.) Y_MatchFound = true;
	  else if( (std::abs(thisPosition.y()-r3FinalReco.y()) < theY_PULL_CUT*sigmay) || 
		   (std::abs(thisPosition.y()-r3FinalReco.y()) < theY_RESIDUAL_CUT   )    ) Y_MatchFound = true; 

	  if(thePHIDIR_RESIDUAL_CUT<0.) Dir_MatchFound = true;
	  else if( std::abs(reco::deltaPhi(p3FinalReco_glob.barePhi(), chamber->toGlobal(thisSegment->localDirection()).barePhi()))<thePHIDIR_RESIDUAL_CUT ) 
	    Dir_MatchFound = true; 

	  //Check for a Match, and if there is a match, check the delR from the segment, keeping only the closest in MuonCandidate
	  if (X_MatchFound && Y_MatchFound && Dir_MatchFound) { 
	    TrackRef thisTrackRef(generalTracks,TrackNumber); 
	    GlobalPoint SegPos(chamber->toGlobal(thisSegment->localPosition())); 
	    GlobalPoint TkPos(r3FinalReco_globv.x(),r3FinalReco_globv.y(),r3FinalReco_globv.z()); 

	    double thisDelR2 = reco::deltaR2(SegPos,TkPos);
	    if(thisDelR2<ClosestDelR2){
	      //std::cout << "[TMP:DT]     Good segment match: DeltaR = " << thisDelR2 << std::endl; 
	      ClosestDelR2 = thisDelR2;
	      //unsigned int segmentNumber = std::find(ourSegments->begin(), ourSegments->end(), *(thisSegment.get())) - ourSegments->begin(); 
	      unsigned int segmentNumber = 0; 
	      for(auto tmpsgt=ourSegments->begin(); tmpsgt!=ourSegments->end(); ++tmpsgt, ++segmentNumber) { 
		if( (&(*tmpsgt)) == thisSegment ) break; 
	      } 
	      //std::cout << "[TMP:DT]       SegmentNumber >< " << ourSegments->size() << std::endl; 
	      if(segmentNumber>=ourSegments->size()) 
		throw cms::Exception("ME0SegmentCollection/find") 
		  << "Associated segment not found in ME0SegmentCollection.\n"; 
	      MuonCandidate = reco::ME0Muon(thisTrackRef, (*thisSegment), segmentNumber, chargeReco); 
	      MuonCandidate.setGlobalTrackPosAtSurface(r3FinalReco_glob);
	      MuonCandidate.setGlobalTrackMomAtSurface(p3FinalReco_glob);
	      MuonCandidate.setLocalTrackPosAtSurface(r3FinalReco);
	      MuonCandidate.setLocalTrackMomAtSurface(p3FinalReco);
	      MuonCandidate.setGlobalTrackCov(covFinalReco);
	      MuonCandidate.setLocalTrackCov(C);
	    }
	  }
	} // End loop for(const auto& thisSegment : chamber.segments) 

	// As long as the delR of the MuonCandidate is sensible, store the track-segment pair
	if (ClosestDelR2 < 500.) { 
	  oc->push_back(MuonCandidate);
	} 
      } // End for (const auto& chamber : info.chambers)
    } 
    
    // Put collection in event 
    ev.put(std::move(oc));
} 

FreeTrajectoryState
ME0SegmentMatcher::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix55& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CurvilinearTrajectoryError tCov(cov);
  
  return cov.kRows == 5 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

FreeTrajectoryState
ME0SegmentMatcher::getFTS(const GlobalVector& p3, const GlobalVector& r3, 
			   int charge, const AlgebraicSymMatrix66& cov,
			   const MagneticField* field){

  GlobalVector p3GV(p3.x(), p3.y(), p3.z());
  GlobalPoint r3GP(r3.x(), r3.y(), r3.z());
  GlobalTrajectoryParameters tPars(r3GP, p3GV, charge, field);

  CartesianTrajectoryError tCov(cov);
  
  return cov.kRows == 6 ? FreeTrajectoryState(tPars, tCov) : FreeTrajectoryState(tPars) ;
}

void ME0SegmentMatcher::getFromFTS(const FreeTrajectoryState& fts,
				    GlobalVector& p3, GlobalVector& r3, 
				    int& charge, AlgebraicSymMatrix66& cov){
  GlobalVector p3GV = fts.momentum();
  GlobalPoint r3GP = fts.position();

  GlobalVector p3T(p3GV.x(), p3GV.y(), p3GV.z());
  GlobalVector r3T(r3GP.x(), r3GP.y(), r3GP.z());
  p3 = p3T;
  r3 = r3T;  
  // p3.set(p3GV.x(), p3GV.y(), p3GV.z());
  // r3.set(r3GP.x(), r3GP.y(), r3GP.z());
  
  charge = fts.charge();
  cov = fts.hasError() ? fts.cartesianError().matrix() : AlgebraicSymMatrix66();

}


void ME0SegmentMatcher::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{

}


 DEFINE_FWK_MODULE(ME0SegmentMatcher);
