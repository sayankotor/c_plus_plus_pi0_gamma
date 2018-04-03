// Include files

#include "CaloUtils/CaloMomentum.h"
#include "CaloUtils/CaloAlgUtils.h"
// from Event
#include "Event/RecHeader.h"
#include "Event/CaloHypo.h"
#include "Event/CaloCluster.h"
#include <iostream>
#include <iomanip>
// local
#include "GammaPi0XGBoostTool.h"

using namespace LHCb;
using namespace Gaudi::Units;
//-----------------------------------------------------------------------------
// Implementation file for class : GammaPi0XGBoostTool
//
// 2018-03-24 : author @sayankotor
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( GammaPi0XGBoostTool )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
GammaPi0XGBoostTool::GammaPi0XGBoostTool( const std::string& type,
                                                const std::string& name,
                                                const IInterface* parent )
: base_class ( type, name , parent )
{
  declareInterface<IGammaPi0SeparationTool>(this);
  declareProperty("MinPt", m_minPt );
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode GammaPi0XGBoostTool::initialize() {

  std::cout<<"Initialize()"<<std::endl;

  StatusCode sc = base_class::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Initialize" << endmsg;

  /// Retrieve geometry of detector
  m_ecal = getDetIfExists<DeCalorimeter>( DeCalorimeterLocation::Ecal );
  

  // TMVA discriminant

  m_xgb0 = std::make_unique<XGBClassifierPhPi0>("/afs/cern.ch/user/v/vchekali/public/DaVinci2/DaVinciDev/models/classes/0_simpl.model");
  m_xgb1 = std::make_unique<XGBClassifierPhPi0>("/afs/cern.ch/user/v/vchekali/public/DaVinci2/DaVinciDev/models/classes/1_simpl.model");
  m_xgb2 = std::make_unique<XGBClassifierPhPi0>("/afs/cern.ch/user/v/vchekali/public/DaVinci2/DaVinciDev/models/classes/2_simpl.model");

  return StatusCode::SUCCESS;
}

//=============================================================================
//  Finalize
//=============================================================================
StatusCode GammaPi0XGBoostTool::finalize() {

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Finalize" << endmsg;

  m_xgb0.reset();
  m_xgb1.reset();
  m_xgb2.reset();


  return base_class::finalize(); // must be executed last
}

//=============================================================================
// Main execution
//=============================================================================


double GammaPi0XGBoostTool::isPhoton(const LHCb::CaloHypo* hypo){
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

  std::vector<double> rawEnergyVector;

  // evaluate the NN inputs
  bool ecalV = ClusterVariables(hypo, fr2, fasym, fkappa, fr2r4, Ecl, Eseed, E2, area);
  bool rawEnergy = GetRawEnergy(hypo, rawEnergyVector);

  if(!rawEnergy) return m_def;
  // return NN output
  std::cout<<"IS photon XGBoost new"<<std::endl;
  if (rawEnergyVector.size() != 25)
  {
    std::cout<<"else: raw energy size: "<<rawEnergyVector.size()<<std::endl;
  }
  double prediction_xgb = XGBDiscriminant(area, rawEnergyVector);
  std::cout<<"prediction_xgb: "<<prediction_xgb<<std::endl;
  //return photonDiscriminant(area, fr2, fr2r4, fasym, fkappa, Eseed, E2,
                            //r2PS, asymPS, eMaxPS, e2ndPS, multiPS, multiPS15, multiPS30, multiPS45);
  return prediction_xgb;
}

double GammaPi0XGBoostTool::XGBDiscriminant(int area, std::vector<double>& row_energies)
{

        double value = ( area == 0 ? m_xgb0->getClassifierValues(row_energies)
                       : area == 1 ? m_xgb1->getClassifierValues(row_energies)
                       : area == 2 ? m_xgb2->getClassifierValues(row_energies)
                       : -1e10 );
        //info() << " INPUT TO GAMMA/PI0 : NN[" << input << "]= " << value << " (area="<<area<<")"<< endmsg;
        return value;
}

bool GammaPi0XGBoostTool::GetRawEnergy(const LHCb::CaloHypo* hypo, std::vector<double>& rowEnergy){
  if( NULL == hypo)return false;
  LHCb::CaloDigits * digits_full = getIfExists<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Ecal);
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );   // OD 2014/05 - change to Split Or Main  cluster
  if( NULL == cluster)return false;

  LHCb::CaloCellID centerID = cluster->seed();
  
  std::cout<<"count neibours: "<<std::endl;
  CaloNeighbors n_vector = m_ecal->zsupNeighborCells(centerID);
  LHCb::CaloCellID::Set n_set = LHCb::CaloCellID::Set(n_vector.begin(), n_vector.end());
  for ( CaloNeighbors::const_iterator neighbor =  n_vector.begin(); n_vector.end() != neighbor ; ++neighbor ){
      CaloNeighbors local_vector = m_ecal->zsupNeighborCells(*neighbor);
      LHCb::CaloCellID::Set new_set = LHCb::CaloCellID::Set(local_vector.begin(), local_vector.end());
      n_set.insert(new_set.begin(), new_set.end());
  }

  std::cout<<n_set.size()<<std::endl;
  std::cout<<std::endl;
  if (n_set.size() < 25){
    std::cout<<"less than 25"<<std::endl;
    return false;
  }
  std::vector<LHCb::CaloDigit> digit_v_sort;
  for( const auto& one_digit : *digits_full ){
        if( abs( (int)one_digit->cellID().col() - (int)centerID.col() ) <= 2 &&
        abs( (int)one_digit->cellID().row() - (int)centerID.row() ) <= 2 &&
        one_digit->cellID().area() == cluster->seed().area() ){
            digit_v_sort.push_back(*one_digit);
            //rowEnergy.push_back(one_digit->e());
            //std::cout<<"col, row, area, en: "<<one_digit->cellID().col()<<" "<<one_digit->cellID().row()<<" "<<one_digit->cellID().area();
            //std::cout<<std::endl;
        }
  }
  std::sort(digit_v_sort.begin(), digit_v_sort.end(), comparer);

  for( const auto& one_digit : digit_v_sort ){
    rowEnergy.push_back(one_digit.e());
  }
  return true;
}

bool GammaPi0XGBoostTool::ClusterVariables(const LHCb::CaloHypo* hypo,
                                              double& fr2, double& fasym, double& fkappa, double& fr2r4, double& etot,
                                              double& Eseed, double& E2, int& area) {
  m_data.clear();

  if( NULL == hypo)return false;
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );   // OD 2014/05 - change to Split Or Main  cluster
  if( NULL == cluster)return false;

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) )debug()<<"Inside ClusterVariables ------"<<endmsg;

  double cxx = cluster->position().spread()(0,0);
  double cyy = cluster->position().spread()(1,1);
  double cxy = cluster->position().spread()(0,1);
  fr2 = cxx + cyy ;
  fasym = (cxx > 0 && cyy > 0) ? cxy /std::sqrt( cxx*cyy ) : 0;
  fkappa = (fr2 > 0. ) ? std::sqrt( 1. - 4.* ( cxx*cyy - cxy*cxy ) / fr2 / fr2  ) : 0;

  double xmean = cluster->position().x();
  double ymean = cluster->position().y();

  // OD : WARNING cluster->e() is cluster-shape dependent (3x3 / 2x2 ...) RE-EVALUATE E3x3 instead for back. compatibility
  // etot = cluster->e(); //same as position.e
  etot = 0.;

  const Gaudi::XYZPoint position (xmean, ymean,cluster->position().z());

  double r4 = 0.;
  area = -1;
  int ncells = 0;
  double secondE = 0.;

  const LHCb::CaloCluster::Entries& entries = cluster->entries() ;
  for ( auto entry = entries.begin() ; entries.end() != entry ; ++entry ){
    const LHCb::CaloDigit* digit = entry->digit()  ;
    if ( 0 == digit ) { continue ; }
    const double fraction = entry->fraction();
    const double energy   = digit->e() * fraction ;

    if( abs( (int)digit->cellID().col() - (int)cluster->seed().col() ) <= 1 &&
        abs( (int)digit->cellID().row() - (int)cluster->seed().row() ) <= 1 &&
        digit->cellID().area() == cluster->seed().area() )etot += energy;


    const Gaudi::XYZPoint& pos =  m_ecal->cellCenter( digit->cellID() );
    const double x =  pos.x() ;
    const double y =  pos.y() ;

    if ( entry->status() & LHCb::CaloDigitStatus::SeedCell ){
      area = digit->cellID().area();
      Eseed = energy;
    }else{
      if(energy>secondE) {
        secondE = energy;
      }
    }

    if ( energy <= 0 ) { continue ; }
    const double weight = energy > 0.0 ? energy : 0.0 ;

    double rr = (x-xmean)*(x-xmean) + (y-ymean)*(y-ymean);
    if( entries.size() <= 1 || rr < 1.e-10 )rr=0; // to avoid huge unphysical value due to machine precision
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