#include <RecoLocalMuon/CSCRecHitD/src/CSCRecHitFromComparatorProducer.h>
#include <RecoLocalMuon/CSCRecHitD/src/CSCRecHitFromComparatorBuilder.h>
//#include <RecoLocalMuon/CSCRecHitD/src/CSCRecHitDBuilder.h>
#include <RecoLocalMuon/CSCRecHitD/src/CSCRecoConditions.h>

#include <FWCore/Framework/interface/Frameworkfwd.h>
#include <FWCore/Framework/interface/EDProducer.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/MakerMacros.h>
#include <DataFormats/Common/interface/Handle.h>
#include <FWCore/Framework/interface/ESHandle.h>
#include <FWCore/ParameterSet/interface/ParameterSet.h>
#include <FWCore/Utilities/interface/Exception.h>
#include <FWCore/MessageLogger/interface/MessageLogger.h>

#include <Geometry/Records/interface/MuonGeometryRecord.h>

#include <DataFormats/CSCRecHit/interface/CSCRecHit2DCollection.h>

CSCRecHitFromComparatorProducer::CSCRecHitFromComparatorProducer( const edm::ParameterSet& ps ) : 
  iRun( 0 ),   
  useCalib( ps.getParameter<bool>("CSCUseCalibrations") ),
  useStaticPedestals( ps.getParameter<bool>("CSCUseStaticPedestals") ),
  useTimingCorrections(ps.getParameter<bool>("CSCUseTimingCorrections") ),
  useGasGainCorrections(ps.getParameter<bool>("CSCUseGasGainCorrections") )

{
  c_token = consumes<CSCComparatorDigiCollection>( ps.getParameter<edm::InputTag>("compDigiTag") );
  //s_token = consumes<CSCStripDigiCollection>( ps.getParameter<edm::InputTag>("stripDigiTag") );
  w_token = consumes<CSCWireDigiCollection>( ps.getParameter<edm::InputTag>("wireDigiTag") );

  recHitBuilder_     = new CSCRecHitFromComparatorBuilder( ps ); // pass on the parameter sets
  //recHitBuilder_     = new CSCRecHitDBuilder( ps ); // pass on the parameter sets
  recoConditions_    = new CSCRecoConditions( ps ); // access to conditions data

  recHitBuilder_->setConditions( recoConditions_ ); // pass down to who needs access

  // register what this produces
  produces<CSCRecHit2DCollection>();

}

CSCRecHitFromComparatorProducer::~CSCRecHitFromComparatorProducer()
{
  delete recHitBuilder_;
  delete recoConditions_;
}


void  CSCRecHitFromComparatorProducer::produce( edm::Event& ev, const edm::EventSetup& setup )
{
  // Dumps the message TWICE if both categories are set!
  //  LogTrace("CSCRecHitFromComparatorProducer|CSCRecHit")<< "[CSCRecHitFromComparatorProducer] starting event " << ev.id().event() << " of run " << ev.id().run();
  LogTrace("CSCRecHit")<< "[CSCRecHitFromComparatorProducer] starting event " << ev.id().event() << " of run " << ev.id().run();
  // find the geometry for this event & cache it in the builder
  edm::ESHandle<CSCGeometry> h;
  setup.get<MuonGeometryRecord>().get( h );
  const CSCGeometry* pgeom = &*h;
  recHitBuilder_->setGeometry( pgeom );

  // access conditions data for this event 
  if ( useCalib || useStaticPedestals || useTimingCorrections || useGasGainCorrections) {  
    recoConditions_->initializeEvent( setup ); 
  }
	
  // Get the collections of strip & wire digis from event
  edm::Handle<CSCComparatorDigiCollection> compDigis;
  //edm::Handle<CSCStripDigiCollection> stripDigis;
  edm::Handle<CSCWireDigiCollection> wireDigis;

  ev.getByToken( c_token, compDigis);
  //ev.getByToken( s_token, stripDigis);
  ev.getByToken( w_token, wireDigis);

  // Create empty collection of rechits  
  auto oc = std::make_unique<CSCRecHit2DCollection>();

  // Fill the CSCRecHit2DCollection
  recHitBuilder_->build( compDigis.product(), wireDigis.product(), *oc);
  //recHitBuilder_->build( stripDigis.product(), wireDigis.product(), *oc);

  // Put collection in event
  LogTrace("CSCRecHit")<< "[CSCRecHitFromComparatorProducer] putting collection of " << oc->size() << " rechits into event.";
  ev.put(std::move(oc));

}

//define this as a plug-in
DEFINE_FWK_MODULE(CSCRecHitFromComparatorProducer);

