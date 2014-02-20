#ifndef TrackingTools_TrackRefitter_ME0TracksToTrajectories_H
#define TrackingTools_TrackRefitter_ME0TracksToTrajectories_H

/** \class ME0TracksToTrajectories
 *  This class, which is a EDProducer, takes a reco::TrackCollection from the Event and refits the rechits 
 *  strored in the reco::Tracks. The final result is a std::vector of Trajectories (objs of the type "Trajectory"), 
 *  which is loaded into the Event in a transient way
 *
 *  $Date: 2010/02/11 00:15:17 $
 *  $Revision: 1.4 $
 *  \author R. Bellan - INFN Torino <riccardo.bellan@cern.ch>
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Utilities/interface/InputTag.h"

namespace edm {class ParameterSet; class Event; class EventSetup;}
class TrackTransformerBase;
class ME0TrackTransformer;

class ME0TracksToTrajectories: public edm::EDProducer{
public:

  /// Constructor
  ME0TracksToTrajectories(const edm::ParameterSet&);

  /// Destructor
  virtual ~ME0TracksToTrajectories();
  
  // Operations
  void endJob();

  /// Convert a reco::TrackCollection into std::vector<Trajectory>
  virtual void produce(edm::Event&, const edm::EventSetup&);
  
 protected:
  
 private:
  
  edm::InputTag theME0TracksLabel;
  ME0TrackTransformer *theTrackTransformer;

  int theNTracks;
  int theNFailures;
};
#endif

