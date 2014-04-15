#include "TrackingTools/TrackRefitter/plugins/ME0TracksToTrajectories.h"

#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
//#include "TrackingTools/TrackRefitter/interface/TrackTransformerForGlobalCosmicMuons.h"
//#include "TrackingTools/TrackRefitter/interface/TrackTransformerForCosmicMuons.h"
#include "TrackingTools/TrackRefitter/interface/ME0TrackTransformer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"

// Only for this one time!!!!
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace std;
using namespace edm;

/// Constructor
ME0TracksToTrajectories::ME0TracksToTrajectories(const ParameterSet& parameterSet):theTrackTransformer(0),
									     theNTracks(0),theNFailures(0){

  theME0TracksLabel = parameterSet.getParameter<InputTag>("Tracks");

  ParameterSet trackTransformerParam = parameterSet.getParameter<ParameterSet>("TrackTransformer");

  string type = parameterSet.getParameter<string>("Type");

  if(type == "ME0Muons") theTrackTransformer = new ME0TrackTransformer(trackTransformerParam);
  //else if(type == "GlobalCosmicMuonsForAlignment") theTrackTransformer = new TrackTransformerForGlobalCosmicMuons(trackTransformerParam);
  //else if(type == "CosmicMuonsForAlignment") theTrackTransformer = new TrackTransformerForCosmicMuons(trackTransformerParam);
  else{
    throw cms::Exception("ME0TracksToTrajectories") 
      <<"The selected algorithm does not exist"
      << "\n"
      << "The only possible choice is ME0Muons"
      << "\n"
      << "Flexible, uh?!";
  }

  produces<vector<Trajectory> >("Refitted");
  //produces<TrajTrackAssociationCollection>("Refitted");
}
 

/// Destructor
ME0TracksToTrajectories::~ME0TracksToTrajectories(){
  if(theTrackTransformer) delete theTrackTransformer;
}

void ME0TracksToTrajectories::endJob(){
  const string metname = "Reco|TrackingTools|ME0TracksToTrajectories";
  
  if(theNFailures!=0)
    LogWarning(metname) << "During the refit there were " 
			<< theNFailures << " out of " << theNTracks << " tracks, i.e. failure rate is: " << double(theNFailures)/theNTracks;
  else{
    LogTrace(metname) << "Refit of the tracks done without any failure";
  }
}


/// Convert Tracks into Trajectories
void ME0TracksToTrajectories::produce(Event& event, const EventSetup& setup){

  const string metname = "Reco|TrackingTools|ME0TracksToTrajectories";

  theTrackTransformer->setServices(setup);
  
  // Collection of Trajectory
  auto_ptr<vector<Trajectory> > trajectoryCollection(new vector<Trajectory>);
  
  // Get the reference
  RefProd<vector<Trajectory> > trajectoryCollectionRefProd 
    = event.getRefBeforePut<vector<Trajectory> >("Refitted");
  
  // Association map between Trajectory and Track
  //auto_ptr<TrajTrackAssociationCollection> trajTrackMap(new TrajTrackAssociationCollection);
 
 
  ////////////////////////////////
  // Just for this one time!!!! //
  ////////////////////////////////

  Handle<reco::GenParticleCollection> genmuons;
  event.getByLabel("genParticles", genmuons);
  
//   unsigned int gensize = genmuons->size();

//   for(unsigned int i=0; i<gensize; ++i) {
//     const reco::GenParticle& gp = (*genmuons)[i];
//     if( gp.status()==1 && abs(gp.pdgId())==13 ) {
//       std::cout << " ###### GenMuon " << i << ": " 
// 		<< gp.eta() << " / " 
// 		<< gp.phi() << " / " 
// 		<< gp.charge() << " / " 
// 		<< gp.pt() << std::endl; 
//     }
//   }

  ////////////////////////////////
  //            END             //
  ////////////////////////////////


 // Get the RecTrack collection from the event
  // Handle<reco::TrackCollection> tracks;
  // event.getByLabel(theTracksLabel,tracks);
  Handle<ME0MuonCollection> muons;
  event.getByLabel(theME0TracksLabel,muons);
  
  Ref<vector<Trajectory> >::key_type trajectoryIndex = 0;
  reco::TrackRef::key_type trackIndex = 0;

  // Loop over the Rec tracks
  for (ME0MuonCollection::const_iterator newMuon = muons->begin(); 
       newMuon != muons->end(); ++newMuon) {
    
    ++theNTracks;

//     std::cout << " ###### ME0Muon " << theNTracks << ": " 
// 	      << newMuon->eta() << " / " 
// 	      << newMuon->phi() << " / " 
// 	      << newMuon->charge() << " / " 
// 	      << newMuon->pt() << std::endl; 

    //vector<Trajectory> trajectoriesSM = theTrackTransformer->transform(*newMuon);
    vector<Trajectory> trajectoriesSM = theTrackTransformer->transform(*newMuon, genmuons);
    
    if(!trajectoriesSM.empty()){
      // Load the trajectory in the Trajectory Container
      trajectoryCollection->push_back(trajectoriesSM.front());
      //std::cout << " ME0 trajectory size: " << trajectoriesSM.size() << std::endl; 

      // Make the association between the Trajectory and the original Track
      //trajTrackMap->insert(Ref<vector<Trajectory> >(trajectoryCollectionRefProd,trajectoryIndex++),
      //		   reco::TrackRef(tracks,trackIndex++));

//       std::cout << " ###### RefMuon " << theNTracks << ": " 
// 		<< trajectoriesSM.front().firstMeasurement().updatedState().globalMomentum().eta() << " / " 
// 		<< trajectoriesSM.front().firstMeasurement().updatedState().globalMomentum().phi() << " / " 
// 		<< trajectoriesSM.front().firstMeasurement().updatedState().charge() << " / " 
// 		<< trajectoriesSM.front().firstMeasurement().updatedState().globalMomentum().perp() << std::endl; 
    }
    else{
      LogTrace(metname) << "Error in the Track refitting. This should not happen";
      ++theNFailures;
      //std::cout << " Error in the Track refitting. This should not happen. N. failures: " << theNFailures << std::endl; 
    }
  }
  LogTrace(metname)<<"Load the Trajectory Collection";
  event.put(trajectoryCollection,"Refitted");
  //event.put(trajTrackMap,"Refitted");
}
