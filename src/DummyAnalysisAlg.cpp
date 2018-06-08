// $Id: DummyAnalysisAlg.cpp,v 1.1 2009-10-02 07:24:53 cattanem Exp $
// Include files 

// from Gaudi
#include "GaudiKernel/AlgFactory.h" 

// local
#include "DummyAnalysisAlg.h"

//-----------------------------------------------------------------------------
// Implementation file for class : DummyAnalysisAlg
//
// 2009-10-02 : Marco Cattaneo
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_ALGORITHM_FACTORY( DummyAnalysisAlg )

#include <vector>
#include "LoKi/Trees.h"
#include "LoKi/PrintMCDecay.h"
#include "Event/ODIN.h"
#include "Event/CaloHypo.h"
#include "CaloInterfaces/IGammaPi0SeparationTool.h"


using namespace std;
using namespace LHCb;
using namespace LoKi::Types;

typedef LHCb::CaloHypo::Container    Hypos;
typedef LHCb::CaloHypo               Hypo;


// internal stuff
namespace {

  unsigned long eventCounter = 0;

  const double zEcal = 12693.;
  const double dEcal = 121.2;

  CaloCellID direction2Ecal (double dxdz, double dydz) {
    double xCell = dxdz * zEcal / dEcal;
    double yCell = dydz * zEcal / dEcal;
    if (fabs(xCell) <= 8 and fabs(yCell) <= 6) { // # inner region 16x12 modules
      double xCellLocal = xCell*3;
      double yCellLocal = yCell*3;
      if (fabs(xCellLocal) <= 8 and fabs(yCellLocal) <= 6) { // # beam hole 16x12 inner cells
	return CaloCellID (0, 0, 0, 0); // # missed Ecal
      }
      return CaloCellID (2, 2, int(floor (yCellLocal))+32, int(floor (xCellLocal))+32);
    }
    else if (fabs(xCell) <= 16 and fabs(yCell) <= 10) { // # middle region 32x20 modules
      double xCellLocal = xCell*2;
      double yCellLocal = yCell*2;
      return CaloCellID (2, 1, int(floor (yCellLocal))+32, int(floor (xCellLocal))+32);
    }
    else if (fabs(xCell) <= 32 and fabs(yCell) <= 26) { // # outer region 64x52 modules
      return CaloCellID (2, 0, int(floor (yCell))+32, int(floor (xCell))+32);
    }
    return CaloCellID (0, 0, 0, 0); // # missed ECAL
  }

  const Hypo* matchHypoCluster (const Hypos* hypos, const MCParticle* mcp) {
    double dxdzMC = mcp->momentum().Px()/mcp->momentum().Pz();
    double dydzMC = mcp->momentum().Py()/mcp->momentum().Pz();
    const Hypo* bestHypo = 0;
    double bestDr = 1e6;
    double xseedBest = 1e6;
    double yseedBest = 1e6;
    for (auto ihypo = hypos->begin(); ihypo != hypos->end(); ++ihypo) {
      Hypo* hypo = *ihypo;
      const CaloCellID& seed = (*(hypo->clusters().begin()))->seed();
      double xseed = (float(seed.col())-31.5) * dEcal;
      double yseed = (float(seed.row())-31.5) * dEcal;
      if (seed.area() == 1) {
	xseed /= 2.;
	yseed /= 2.;
      }
      else if (seed.area() == 2) {
	xseed /= 3.;
	yseed /= 3.;
      }
      double dxdzSeed = xseed / zEcal;
      double dydzSeed = yseed / zEcal;

      double dr2 = (dxdzMC - dxdzSeed) * (dxdzMC - dxdzSeed) + (dydzMC - dydzSeed) * (dydzMC - dydzSeed);
      if (dr2 < bestDr) {
	bestDr = dr2;
	bestHypo = hypo;
	xseedBest = xseed;
	yseedBest = yseed;
      }
    }
    return bestHypo;
  }
  
}



//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
DummyAnalysisAlg::DummyAnalysisAlg( const std::string& name, ISvcLocator* pSvcLocator)
  :
  DaVinciTupleAlgorithm ( name , pSvcLocator ),
  m_finder_BKPiGamma (Decays::Trees::Invalid_<const LHCb::MCParticle*> ()),
  m_finder_BKPiPi0 (Decays::Trees::Invalid_<const LHCb::MCParticle*> ()),
  m_GammaPi0 (0),
  m_GammaPi0_XGB (0)
{}


//=============================================================================
// Destructor
//=============================================================================
DummyAnalysisAlg::~DummyAnalysisAlg() {} 

//=============================================================================
// Initialization
//=============================================================================
StatusCode DummyAnalysisAlg::initialize() 
{
  StatusCode sc = DaVinciTupleAlgorithm::initialize(); 
  if ( sc.isFailure() ) return sc;
  
  SmartIF<Decays::IMCDecay> decay = tool<Decays::IMCDecay> ( "LoKi::MCDecay" , this ) ;
  if ( !decay ) { return Error ( "Unable to locate MC-decay finder" ) ; }
  // try to parse/decode the decay descriptor
  // B->KPiGamma
  string bkpigamma = "[(B0 => K+ pi- ^gamma)]CC";
  Decays::IMCDecay::Tree tree = decay->tree ( bkpigamma ) ;
  if ( !tree  ) { return Error ( "Unable to build valid decay tree: '" 
				 + bkpigamma + "'" ) ; }
  cout << "DummyAnalysisAlg::initialize-> the decay tree is: " << tree << endl ;
  // initialize the decay finder 
  m_finder_BKPiGamma = Decays::IMCDecay::Finder ( tree ) ;
  if ( !m_finder_BKPiGamma  ) { return Error ( "Unable to build valid decay finder: '" 
					       + bkpigamma + "'" ) ; }
  // B->KPiPi0
  string bkpipi0 = "[(B0 => K+ pi- ^pi0) || (B0 => K+ K- ^pi0)]CC";
  //  string bkpipi0 = "[(B0 ==> K+ pi- (^pi0 -> gamma gamma))]CC";
  //string bkpipi0 = "[(B0 ==> K+ pi- ^gamma)]CC";
  tree = decay->tree ( bkpipi0 ) ;
  if ( !tree  ) { return Error ( "Unable to build valid decay tree: '" 
				  + bkpipi0 + "'" ) ; }
  cout << "DummyAnalysisAlg::initialize-> the decay tree is: " << tree << endl ;
   // initialize the decay finder 
   m_finder_BKPiPi0 = Decays::IMCDecay::Finder ( tree ) ;
   if ( !m_finder_BKPiPi0  ) { return Error ( "Unable to build valid decay finder: '" 
   					       + bkpipi0 + "'" ) ; }

  m_GammaPi0  = tool<IGammaPi0SeparationTool>("GammaPi0SeparationTool" , "GammaPi0SeparationTool", this);
  m_GammaPi0_XGB  = tool<IGammaPi0SeparationTool>("GammaPi0XGBoostTool" , "GammaPi0XGBoostTool", this);
  
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Initialize" << endmsg;
  
  return StatusCode::SUCCESS;
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode DummyAnalysisAlg::execute()
{
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;
  
  // code goes here
  eventCounter = eventCounter + 1;
  //print "====> processing event ",  self.totalEvents
  //        event = {}
  //  cout << "DummyAnalysisAlg::execute() ->  " << endl;
  // Load the ODIN
  const LHCb::ODIN* odin = getIfExists<ODIN>(evtSvc(),LHCb::ODINLocation::Default);
  if ( !odin ) { odin = getIfExists<ODIN>(evtSvc(),LHCb::ODINLocation::Default,false); }
  if ( !odin )
    {
      // should always be available ...
      return Error( "Cannot load the ODIN data object", StatusCode::SUCCESS );
    }
  auto eventID= odin->eventNumber();
  auto runID= odin->runNumber();
  unsigned long longEventId (runID);
  longEventId = longEventId << 32;
  longEventId  |= (eventID & 0xffffffff);
  if (eventID > 0xffffffff || runID > 0xffffffff) {
    cout << "eventid overflow:" << eventID<< ':'<< runID<<endl;
  }
  Tuple t_EvId = nTuple("ev_id");
  t_EvId->column ("eventID", longEventId); 
  t_EvId->write();
  cout<< "now processing run "<<runID<<", event: "<< eventID<<endl;
  
  // get all MC-particles:
  const LHCb::MCParticle::Container* mcparticles = 
    get<LHCb::MCParticle::Container>( LHCb::MCParticleLocation::Default ) ;

  // for (auto imc = mcparticles->begin(); imc !=  mcparticles->end(); ++imc) {
  //   if ((*imc)->particleID().abspid() == 511) {
  //     cout << "B-mepson decay chain:"<<endl; 
  //     LoKi::PrintMC::printDecay ( *imc , cout ) ;
  //     cout << **imc << endl;
  //     const MCVertex* vx = *((*imc)->endVertices().begin());
  //     cout << *vx << endl;
  //     for (auto ida = vx->products().begin(); ida != vx->products().end(); ++ida) {
  // 	cout << **ida << endl;
  //     }
  //   }
  // }

  
  // BKPIGAMMA
  LHCb::MCParticle::ConstVector  mcBGamma ;
  m_finder_BKPiGamma.findDecay
    ( mcparticles -> begin() ,   // begin of input 
      mcparticles -> end  () ,   // end of the input 
      mcBGamma                 ) ; // the output container 
  
  //cout << " bkpigamma found #" << mcBGamma.size() << " decays" ;
  
  const LHCb::MCParticle* bGamma = 0;
  double xGamma = 0;
  double yGamma = 0;
  double eGamma = 0;
  CaloCellID gammaCell (0, 0, 0, 0);
  Tuple t_mcBGamma = nTuple("mcBGamma");
  if (!mcBGamma.empty()) {
    std::cout<<"common bgamma"<<std::endl;
    const LHCb::MCParticle* gamma = *(mcBGamma.begin());
    //    cout << "gamma->mother() " << *(gamma->mother()) << endl;
    plot1D (log10(gamma->momentum().P()), "Energy_Gamma_B", "; log10(p_{#gamma})", 0, 10, 200);
    CaloCellID gCell = direction2Ecal (gamma->momentum().Px() / gamma->momentum().Pz(), gamma->momentum().Py() / gamma->momentum().Pz());
    //    cout << "found photon cell-> " << gCell << endl;
    if (gCell.calo() > 0) {
      plot1D (log10(gamma->momentum().P()), "Energy_Gamma_B_acceptance", "; log10(p_{#gamma})", 0, 10, 200);
      auto ivx = gamma->endVertices().begin();
      if (ivx != gamma->endVertices().end()) {
	const LHCb::MCVertex* vx = *ivx;
	plot1D (vx->position().Z(), "Gamma_B_endpoint", "; photon decay point Z", 0, 15000, 3000);
	if (vx->position().Z() > 11750) { // hit ECAL directly
	  bGamma = gamma;
	  t_mcBGamma->column ("eventID", longEventId);
	  t_mcBGamma->column ("px", gamma->momentum().Px()); 
	  t_mcBGamma->column ("py", gamma->momentum().Py()); 
	  t_mcBGamma->column ("pz", gamma->momentum().Pz()); 
	  t_mcBGamma->write();
	  xGamma = gamma->momentum().Px() / gamma->momentum().Pz() * zEcal;
	  yGamma = gamma->momentum().Py() / gamma->momentum().Pz() * zEcal;
	  eGamma = gamma->momentum().P();
	  gammaCell = gCell;
	  //	  std::cout << "written photon: " << *gamma << endl;
	}
      }
    }
  }

// BKPIPI0
  LHCb::MCParticle::ConstVector  mcBPi0 ;
  LHCb::MCParticle::ConstVector  mcPi0 ;
  m_finder_BKPiPi0.findDecay
    ( mcparticles -> begin() ,   // begin of input 
      mcparticles -> end  () ,   // end of the input 
      mcBPi0                 ) ; // the output container 
  
  //cout << " bkpipi0 found #" << mcBPi0.size() << " decays" << endl;

  const LHCb::MCParticle* bPi0 = 0;
  double xPi0 = 0;
  double yPi0 = 0;
  double ePi0 = 0;
  CaloCellID pi0Cell0 (0, 0, 0, 0);
  CaloCellID pi0Cell1 (0, 0, 0, 0);
  CaloCellID pi0Cell2 (0, 0, 0, 0);


  Tuple t_mcPi0 = nTuple("mcPi0");
  if (!mcPi0.empty()) {
    std::cout<<"common pi0"<<std::endl;
    std::cout<<"event id pi0"<<longEventId<<std::endl;
  }

  Tuple t_mcBPi0 = nTuple("mcBPi0");
  if (!mcBPi0.empty()) {
    bool writing = false;
    std::cout<<"common bpi0"<<std::endl;
    std::cout<<"event id common bpi0"<<longEventId<<std::endl;
    if (longEventId == 13153543503003945){
        writing = true;
        std::cout<<"writing... "<<std::endl;
    }
    const LHCb::MCParticle* pi0 = *(mcBPi0.begin());
    //    cout << "pi0->mother() " << *(pi0->mother()) << endl;
    plot1D (log10(pi0->momentum().P()), "Energy_Pi0_B", "; log10(p_{#pi^{0}})", 0, 10, 200);
    CaloCellID gCell0 = direction2Ecal (pi0->momentum().Px() / pi0->momentum().Pz(), pi0->momentum().Py() / pi0->momentum().Pz());
    //    cout << "found pi0 cell-> " << gCell0 << endl;
    if (gCell0.calo() > 0) {
      plot1D (log10(pi0->momentum().P()), "Energy_Pi0_B_acceptance", "; log10(p_{#pi^{0}})", 0, 10, 200);
    }
    const MCVertex* vx = *(pi0->endVertices().begin());
    const MCParticle* gamma1 = *(vx->products().begin());
    const MCParticle* gamma2 = *(++(vx->products().begin()));
    pi0Cell1 = direction2Ecal (gamma1->momentum().Px() / gamma1->momentum().Pz(), gamma1->momentum().Py() / gamma1->momentum().Pz());
    pi0Cell2 = direction2Ecal (gamma2->momentum().Px() / gamma2->momentum().Pz(), gamma2->momentum().Py() / gamma2->momentum().Pz());
    if (writing) {
      std::cout<<"point 0"<<std::endl;
    }
    if (pi0Cell1.calo() > 0 && pi0Cell2.calo() > 0) {
      if (writing) {
        std::cout<<"point 1"<<std::endl;
      }
      if (pi0Cell1.area() == pi0Cell2.area()) {
      	if (abs (int(pi0Cell1.col()) - int(pi0Cell2.col())) < 3 && abs (int(pi0Cell1.row()) - int(pi0Cell2.row())) < 3) {
          //std::cout<<"bpi0 after ifs"<<std::endl;
      	  plot1D (log10(pi0->momentum().P()), "Energy_Pi0_B_acceptance", "; log10(p_{#pi^{0}})", 0, 10, 200);
      	  auto ivx1 = gamma1->endVertices().begin();
      	  auto ivx2 = gamma2->endVertices().begin();
      	  if (ivx1 != gamma1->endVertices().end() && ivx2 != gamma2->endVertices().end()) {
      	    double z1 = (*ivx1)->position().Z();
      	    double z2 = (*ivx2)->position().Z();
      	    plot1D (z1, "Gamma1_B_endpoint", "; photon1 decay point Z", 0, 15000, 3000);
      	    plot1D (z2, "Gamma2_B_endpoint", "; photon2 decay point Z", 0, 15000, 3000);
      	    if ((*ivx1)->position().Z() > 11750 && (*ivx2)->position().Z() > 11750) { // hit ECAL directly
              if (writing) {
                std::cout<<"point 5"<<std::endl;
              }
              std::cout<<"write bpi0"<<std::endl;
              std::cout<<"event id bpi0"<<longEventId<<std::endl;
      	      bPi0 = pi0;
      	      t_mcBPi0->column ("eventID", longEventId);
      	      t_mcBPi0->column ("px", pi0->momentum().Px()); 
      	      t_mcBPi0->column ("py", pi0->momentum().Py()); 
      	      t_mcBPi0->column ("pz", pi0->momentum().Pz()); 
      	      t_mcBPi0->write();
      	      xPi0 = pi0->momentum().Px() / pi0->momentum().Pz() * zEcal;
      	      yPi0 = pi0->momentum().Py() / pi0->momentum().Pz() * zEcal;
      	      ePi0 = pi0->momentum().P();
      	      pi0Cell0 = gCell0;
      	      //	      std::cout << "written photon: " << *pi0 << endl;
      	    }
      	  }
      	}
      }
    }
  }

  if ((gammaCell.calo() > 0) || (pi0Cell0.calo() > 0)) {
    // hypo
    if (exist<Hypos>( "Rec/Calo/Photons" )){
      Hypos* photons = get<Hypos>( "Rec/Calo/Photons" );
      //      cout << "Rec/Calo/Photons size = " << photons->size() << endl;
      
      if (bGamma) {
	const Hypo* hGamma = matchHypoCluster (photons, bGamma);
	if (hGamma) { 
	  double isPhoton = m_GammaPi0->isPhoton (hGamma);
	  double isPhotonBDT = m_GammaPi0_XGB->isPhoton (hGamma);
	  Tuple t_Gamma = nTuple("Gamma");
	  t_Gamma->column ("eventID", longEventId); 
	  t_Gamma->column ("px", bGamma->momentum().Px()); 
	  t_Gamma->column ("py", bGamma->momentum().Py()); 
	  t_Gamma->column ("pz", bGamma->momentum().Pz());
	  auto seed = (*(hGamma->clusters().begin()))->seed();
	  t_Gamma->column ("area", seed.area());
	  t_Gamma->column ("row", seed.row());
	  t_Gamma->column ("col", seed.col());
	  t_Gamma->column ("isPhoton", isPhoton);
	  t_Gamma->column ("isPhotonBDT", isPhotonBDT);
	  t_Gamma->write();
	}
      }
      
      if (bPi0) {
	const Hypo* hPi0 = matchHypoCluster (photons, bPi0);
	if (hPi0) { 
	  double isPhoton = m_GammaPi0->isPhoton (hPi0);
	  double isPhotonBDT = m_GammaPi0_XGB->isPhoton (hPi0);
	  Tuple t_Pi0 = nTuple("Pi0");
	  t_Pi0->column ("eventID", longEventId); 
	  t_Pi0->column ("px", bPi0->momentum().Px()); 
	  t_Pi0->column ("py", bPi0->momentum().Py()); 
	  t_Pi0->column ("pz", bPi0->momentum().Pz());
	  auto seed = (*(hPi0->clusters().begin()))->seed();
	  t_Pi0->column ("area", seed.area());
	  t_Pi0->column ("row", seed.row());
	  t_Pi0->column ("col", seed.col());
	  t_Pi0->column ("isPhoton", isPhoton);
	  t_Pi0->column ("isPhotonBDT", isPhotonBDT);
	  t_Pi0->write();
	}
      }
    }
  }
  setFilterPassed(true);  // Mandatory. Set to true if event is accepted. 
  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode DummyAnalysisAlg::finalize()
{
  if ( msgLevel(MSG::DEBUG) ) debug() << "==> Finalize" << endmsg;

  return DaVinciTupleAlgorithm::finalize();
}
//=============================================================================
