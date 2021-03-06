// Include files

#include <cmath>

#include "CaloUtils/CaloMomentum.h"
#include "CaloUtils/CaloAlgUtils.h"
// from Event
#include "Event/RecHeader.h"
#include "Event/CaloHypo.h"
#include "Event/CaloCluster.h"

// local
#include "GammaPi0SeparationTool.h"

using namespace LHCb;
using namespace Gaudi::Units;
//-----------------------------------------------------------------------------
// Implementation file for class : GammaPi0SeparationTool
//
// 2010-03-24 : Miriam Calvo Gomez
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_COMPONENT( GammaPi0SeparationTool )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
GammaPi0SeparationTool::GammaPi0SeparationTool( const std::string& type,
                                                const std::string& name,
                                                const IInterface* parent )
: base_class ( type, name , parent )
{
  declareInterface<IGammaPi0SeparationTool>(this);
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode GammaPi0SeparationTool::initialize() {

  StatusCode sc = base_class::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Initialize" << endmsg;

  /// Retrieve geometry of detector
  m_ecal = getDetIfExists<DeCalorimeter>( DeCalorimeterLocation::Ecal );

  // TMVA discriminant
  static const std::vector<std::string> inputVars = {
           "fr2", "fr2r4", "abs(asym)",
           "kappa", "Eseed/Ecl", "(E2+Eseed)/Ecl",
           "eSumPS>0?eMaxPS/eSumPS:0",
           "eSumPS>0?e2ndPS/eSumPS:0",
           "r2PS", "abs(asymPS)",
           "multiPS", "multiPS15",
           "multiPS30", "multiPS45" } ;

  m_reader0 = std::make_unique<ReadMLPOuter>(inputVars);
  m_reader1 = std::make_unique<ReadMLPMiddle>(inputVars);
  m_reader2 = std::make_unique<ReadMLPInner>(inputVars);

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode GammaPi0SeparationTool::finalize() {

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Finalize" << endmsg;

  m_reader0.reset();
  m_reader1.reset();
  m_reader2.reset();

  return base_class::finalize(); // must be executed last
}

//=============================================================================
// Main execution
//=============================================================================


double GammaPi0SeparationTool::isPhoton(const LHCb::CaloHypo* hypo){

  // clear all data
  m_data.clear();
  m_prsdata.clear();

  if ( !m_ecal ) return m_def;
  if ( LHCb::CaloMomentum(hypo).pt() < m_minPt) return m_def;

  double fr2 = 0;
  double fasym = 0;
  double fkappa = 0;
  double fr2r4 = 0;
  double Eseed = 0;
  double E2 = 0;
  double Ecl = 0;
  int area =0;

  double r2PS = 0.;
  double r2r4PS = 0;
  double asymPS = 0;
  double kappaPS = 0;
  double ePrs = 0;
  double eMaxPS = 0;
  double e2ndPS = 0;
  double ecornerPS = 0;
  double eSumPS = 0.;
  int multiPS = 0;
  int multiPS15 = 0;
  int multiPS30 = 0;
  int multiPS45 = 0;

  // evaluate the NN inputs
  bool ecalV = ClusterVariables(hypo, fr2, fasym, fkappa, fr2r4, Ecl, Eseed, E2, area);
  bool prsV  = PrsVariables(hypo, r2PS, asymPS, kappaPS, r2r4PS, eSumPS, ePrs, eMaxPS, e2ndPS, ecornerPS,
                            multiPS, multiPS15, multiPS30, multiPS45);
  if( !ecalV || !prsV )return m_def;
  // return NN output
  return photonDiscriminant(area, fr2, fr2r4, fasym, fkappa, Eseed, E2,
                            r2PS, asymPS, eMaxPS, e2ndPS, multiPS, multiPS15, multiPS30, multiPS45);
}


bool GammaPi0SeparationTool::ClusterVariables(const LHCb::CaloHypo* hypo,
                                              double& fr2, double& fasym, double& fkappa, double& fr2r4, double& etot,
                                              double& Eseed, double& E2, int& area) {
  m_data.clear();

  if ( !hypo ) return false;
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );   // OD 2014/05 - change to Split Or Main  cluster
  if ( !cluster ) return false;

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) )debug()<<"Inside ClusterVariables ------"<<endmsg;

  const auto cxx = cluster->position().spread()(0,0);
  const auto cyy = cluster->position().spread()(1,1);
  const auto cxy = cluster->position().spread()(0,1);
  fr2 = cxx + cyy ;
  fasym = ( cxx > 0 && cyy > 0 ? cxy / std::sqrt( cxx*cyy ) : 0 );
  const auto arg = ( fr2 > 0. ? 1. - 4.* ( cxx*cyy - cxy*cxy ) / ( fr2 * fr2 ) : 0. );
  fkappa = std::sqrt( arg > 0. ? arg : 0. );

  const auto xmean = cluster->position().x();
  const auto ymean = cluster->position().y();

  // OD : WARNING cluster->e() is cluster-shape dependent (3x3 / 2x2 ...) RE-EVALUATE E3x3 instead for back. compatibility
  // etot = cluster->e(); //same as position.e
  etot = 0.;

  const Gaudi::XYZPoint position ( xmean, ymean, cluster->position().z() );

  double r4 = 0.;
  area = -1;
  int ncells = 0;
  double secondE = 0.;

  const auto & entries = cluster->entries() ;
  for ( auto entry = entries.begin() ; entries.end() != entry ; ++entry ){
    const LHCb::CaloDigit* digit = entry->digit()  ;
    if ( !digit ) { continue ; }
    const auto fraction = entry->fraction();
    const auto energy   = digit->e() * fraction ;

    if( abs( (int)digit->cellID().col() - (int)cluster->seed().col() ) <= 1 &&
        abs( (int)digit->cellID().row() - (int)cluster->seed().row() ) <= 1 &&
        digit->cellID().area() == cluster->seed().area() ) { etot += energy; }

    if ( entry->status() & LHCb::CaloDigitStatus::SeedCell ){
      area = digit->cellID().area();
      Eseed = energy;
    }else{
      if(energy>secondE) {
        secondE = energy;
      }
    }

    const auto& pos =  m_ecal->cellCenter( digit->cellID() );
    const auto x =  pos.x() ;
    const auto y =  pos.y() ;

    if ( energy <= 0 ) { continue ; }
    const double weight = energy > 0.0 ? energy : 0.0 ;

    auto rr = std::pow(x-xmean,2) + std::pow(y-ymean,2);
    if ( entries.size() <= 1 || rr < 1.e-10 ) { rr=0; } // to avoid huge unphysical value due to machine precision
    r4 += weight * rr*rr;

    ncells++;
  }//loop cluster cells

  if( etot > 0. ){
    r4 /= etot;
    fr2r4 = (r4 !=0) ? (r4 - fr2*fr2)/r4 : 0.;
    E2 = (secondE+Eseed)/etot;
    Eseed = Eseed/etot;
  }else{
    // should never happen
    r4 = 0;
    fr2r4 = 0.;
    E2 = 0.;
    Eseed = 0.;
  }

  m_data["Ecl"]   = etot;
  m_data["Fr2"]   = fr2;
  m_data["Fr2r4"] = fr2r4;
  m_data["Asym"] = fasym;
  m_data["Kappa"] = fkappa;
  m_data["Eseed"] = Eseed;
  m_data["E2"]    = E2;
  return true;
}

bool GammaPi0SeparationTool::PrsVariables(const LHCb::CaloHypo* hypo,
                                          double& r2PS, double& asymPS, double& kappaPS, double& r2r4PS,
                                          double& eSumPS, double& ePrs, double& eMaxPS, double& e2ndPS, double& ecornerPS,
                                          int& multiPS, int& multiPS15, int& multiPS30, int& multiPS45) {

  // clear prs data
  m_prsdata.clear();
  e2ndPS = 0.;
  ePrs = 0.;
  eSumPS=0.;
  eMaxPS=0.;
  multiPS=0;
  multiPS15=0;
  multiPS30=0;
  multiPS45=0;
  asymPS=0.;
  r2PS=0.;
  kappaPS=0.;
  r2r4PS=0.;
  ecornerPS=0.;

  if( !hypo)return false;
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );
  if( !cluster)return false;


  LHCb::CaloCellID cellPrs = cluster->seed();  // use cluster seed projection (mostly the same as line projection + faster access)
  cellPrs.setCalo(1);

  if( LHCb::CaloCellID() == cellPrs  )return true; // no valid Prs info

  // compute energy sum / max / ...
  const LHCb::CaloHypo::Digits& digits = hypo->digits();
  double xPs = 0.;
  double yPs = 0.;
  for( auto d = digits.begin() ; digits.end() != d ; ++d ){
    const LHCb::CaloDigit* digit = *d;
    if     ( !digit ) continue           ;
    if ((int) digit->cellID().calo() != 1 )continue;  // select Prs digits
    int dcol = int(digit->cellID().col()) - int(cellPrs.col());
    int drow = int(digit->cellID().row()) - int(cellPrs.row());
    if( abs(dcol) > 1 || abs(drow)> 1 )continue; // keep only neighbors

    double e = digit->e();

    // energy
    eSumPS  += e ;
    if( cellPrs == digit->cellID() )ePrs = e ;

    // barycenter
    xPs += double(dcol) * digit->e();
    yPs += double(drow) * digit->e();

    // multiplicities
    if( e > eMaxPS ) eMaxPS = e;
    if( e > 0.0    ) multiPS++;
    if( e > 15.0   ) multiPS15++;
    if( e > 30.0   ) multiPS30++;
    if( e > 45.0   ) multiPS45++;
  }


  if(  eSumPS <= 0. )return true;
  xPs = xPs/eSumPS;
  yPs = yPs/eSumPS;

  // spread and related shape variables
  double cxxPS = 0.;
  double cyyPS = 0.;
  double cxyPS = 0.;
  double r4PS    = 0.;
  double c1 = 0.0;
  double c2 = 0.0;
  double c3 = 0.0;
  double c4 = 0.0;
  for (const auto& digit : digits) {
    if ( !digit ) continue           ;
    if ((int) digit->cellID().calo() != 1 )continue;  // select Prs digits
    int dcol = ( int(digit->cellID().col()) - int(cellPrs.col()) );
    int drow = ( int(digit->cellID().row()) - int(cellPrs.row()) );
    if( abs(dcol) > 1 || abs(drow)> 1 )continue; // keep only neighbors

    // 2Nd max
    double e = digit->e();
    if( e > e2ndPS && e < eMaxPS )e2ndPS = e;

    // spread
    double dxPS = (xPs - double(dcol));
    double dyPS = (yPs - double(drow));
    cxxPS += dxPS * dxPS * e;
    cyyPS += dyPS * dyPS * e;
    cxyPS += dxPS * dyPS * e;

    // shape variables
    r4PS += e * ( dxPS * dxPS + dyPS * dyPS) * ( dxPS * dxPS + dyPS * dyPS);

    // corner energy
    if(dcol==-1 && drow== 1 ){ c1 += e;}
    if(dcol==-1 && drow== 0 ){ c1 += e; c3 += e;}
    if(dcol==-1 && drow==-1 ){ c3 += e;}
    if(dcol== 0 && drow== 1 ){ c1 += e; c2 += e;}
    if(dcol== 0 && drow==-1 ){ c3 += e; c4 += e;}
    if(dcol== 1 && drow== 1 ){ c2 += e;}
    if(dcol== 1 && drow== 0 ){ c2 += e; c4 += e;}
    if(dcol== 1 && drow==-1 ){ c4 += e;}
  }
  cxxPS = cxxPS/eSumPS;
  cxyPS = cxyPS/eSumPS;
  cyyPS = cyyPS/eSumPS;
  r2PS    = cxxPS + cyyPS;
  asymPS  = ( cxxPS > 0.0 && cyyPS > 0.0) ? cxyPS / std::sqrt( cxxPS * cyyPS) : 0.;
  const auto arg = ( r2PS > 0.0 ? 1.0 - 4.0 *  (cxxPS * cyyPS - cxyPS * cxyPS) / ( r2PS * r2PS ) : 0. );
  kappaPS = std::sqrt( arg > 0. ? arg : 0. );
  r4PS = r4PS/eSumPS;
  if( r4PS!=0.0) r2r4PS = ( r4PS - r2PS * r2PS) /r4PS;

  ecornerPS = std::max( { c1, c2, c3, c4 } );
  eMaxPS = eMaxPS/eSumPS;
  e2ndPS = e2ndPS/eSumPS;


  // === store the data
  m_prsdata["PrsFr2"]     = r2PS;
  m_prsdata["PrsAsym"]    = asymPS;
  m_prsdata["PrsM"]       = multiPS;
  m_prsdata["PrsM15"]     = multiPS15;
  m_prsdata["PrsM30"]     = multiPS30;
  m_prsdata["PrsM45"]     = multiPS45;
  m_prsdata["PrsEmax"]    = eMaxPS;
  m_prsdata["PrsE2"]      = e2ndPS;
  return true;
}


double GammaPi0SeparationTool::photonDiscriminant(int area,
                                                  double r2, double r2r4, double asym,
                                                  double kappa, double Eseed, double E2,
                                                  double r2PS, double asymPS, double eMaxPS, double e2ndPS,
                                                  int multiPS, int multiPS15, int multiPS30, int multiPS45)
{
        std::vector<double> input = {
                            r2,
                            r2r4,
                            std::abs(asym),
                            kappa,
                            Eseed, //already divided by Ecl
                            E2,    //means (e2+eseed)/ecl
                            eMaxPS,//divided by Esum
                            e2ndPS,//divided by Esum
                            r2PS,
                            std::abs(asymPS),
                            double(multiPS),
                            double(multiPS15),
                            double(multiPS30),
                            double(multiPS45) };

        double value = ( area == 0 ? m_reader0->GetMvaValue(input)
                       : area == 1 ? m_reader1->GetMvaValue(input)
                       : area == 2 ? m_reader2->GetMvaValue(input)
                       : -1e10 );
        //info() << " INPUT TO GAMMA/PI0 : NN[" << input << "]= " << value << " (area="<<area<<")"<< endmsg;
        return value;
}
