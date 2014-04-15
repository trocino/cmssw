#include "FWCore/PluginManager/interface/ModuleDef.h"

#include "FWCore/Framework/interface/MakerMacros.h"

#include "TrackingTools/TrackRefitter/plugins/TracksToTrajectories.h"
#include "TrackingTools/TrackRefitter/plugins/ME0TracksToTrajectories.h"


DEFINE_FWK_MODULE(TracksToTrajectories);
DEFINE_FWK_MODULE(ME0TracksToTrajectories);
