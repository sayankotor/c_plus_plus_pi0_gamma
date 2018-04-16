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

  if (!m_cgeom.init()){
    std::cout<<"init do not succes"<<std::endl;
    return 0;
  }
  

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

  std::vector<double> rawEnergyVector(25, 0.0);

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

//bool recreate(std::vector<std::vector<<double>>& initial_vector, ){

//}

double GammaPi0XGBoostTool::XGBDiscriminant(int area, std::vector<double>& row_energies)
{

        double value = ( area == 0 ? m_xgb0->getClassifierValues(row_energies)
                       : area == 1 ? m_xgb1->getClassifierValues(row_energies)
                       : area == 2 ? m_xgb2->getClassifierValues(row_energies)
                       : -1e10 );
        //info() << " INPUT TO GAMMA/PI0 : NN[" << input << "]= " << value << " (area="<<area<<")"<< endmsg;
        return value;
}

std::vector<double> GammaPi0XGBoostTool::RestoreOne(int home_area, int n_area, std::vector<LHCb::CaloDigit>& additional_elems){
  std::vector<std::vector<double>> home_energy = std::vector<std::vector<double>>(5*6, std::vector<double>(6, 0.0));
  std::vector<double> answ (5, 0.0);
  if (n_area == 0){
    for (int i = 0; i< additional_elems.size(); i++)
      for (int j = 0; j < 6; j++)
        for (int k = 0; k < 6; k++){
         home_energy[i*6 + k][j] = additional_elems[i].e()/36.0;
        }

  }

  if (n_area == 1){
    for (int i = 0; i< additional_elems.size(); i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++){
         home_energy[i*3 + k][j] = additional_elems[i].e()/9.0;
        }

  }

  if (n_area == 2){
    for (int i = 0; i< additional_elems.size(); i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++){
         home_energy[i*2 + k][j] = additional_elems[i].e()/2.0;
        }

  }
  
  if (home_area == 0){
    for (int i = 0; i< 5; i++)
      for (int j = 0; j < 6; j++)
        for (int k = 0; k < 6; k++){
         answ [i] += home_energy[i*6 + k][j];
        }
  }

  if (home_area == 1){
    for (int i = 0; i< 5; i++)
      for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++){
         answ [i] += home_energy[i*3 + k][j];
        }
  }

  if (home_area == 2){
    for (int i = 0; i< 5; i++)
      for (int j = 0; j < 2; j++)
        for (int k = 0; k < 2; k++){
         answ [i] += home_energy[i*2 + k][j];
        }
  }

  return answ;

}

std::vector <std::vector<std::pair<int,int> > > GammaPi0XGBoostTool::CheckVector(std::vector<std::vector<double> >& vector_cells){
  std::vector <std::vector<std::pair<int,int> > > answ;
  answ.reserve(25);
  for (int i = 0; i < 5; i++){
    if (vector_cells[i][0] == 0.0 && vector_cells[i][1] == 0.0 && vector_cells[i][2] == 0.0 && vector_cells[i][3] == 0.0 && vector_cells[i][4] == 0.0){
      std::vector<std::pair<int,int> > answ_local = {std::make_pair(i, 0), std::make_pair(i, 1), std::make_pair(i, 2), std::make_pair(i, 3), std::make_pair(i, 4)};
      //std::cout<<std::endl<<"check vector"<<std::endl;
      //std::cout<<"i "<<i<<std::endl;
      //std::cout<<"vector_cells[i][0] == vector_cells[i][1]"<<std::endl;
      answ.emplace_back(answ_local);
    } 
    if (vector_cells[0][i] == 0.0 && vector_cells[1][i] == 0.0 && vector_cells[2][i] == 0.0 && vector_cells[3][i] == 0.0 && vector_cells[4][i] == 0.0){
      std::vector<std::pair<int,int> > answ_local = {std::make_pair(0, i), std::make_pair(1, i), std::make_pair(2, i), std::make_pair(3, i), std::make_pair(4, i)};
      //std::cout<<std::endl<<"check vector"<<std::endl;
      //std::cout<<"i "<<i<<std::endl;
      //std::cout<<"vector_cells[0][i] == vector_cells[1][i]"<<std::endl;
      answ.emplace_back(answ_local);
    } 
  }
  return answ;
}

bool GammaPi0XGBoostTool::GetRawEnergy(const LHCb::CaloHypo* hypo, std::vector<double>& rowEnergy){
  if( NULL == hypo)return false;
  LHCb::CaloDigits * digits_full = getIfExists<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Ecal);
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );   // OD 2014/05 - change to Split Or Main  cluster
  if( NULL == cluster)return false;

  LHCb::CaloCellID centerID = cluster->seed();
  std::cout<<std::endl;
  std::cout<<"calo area row col "<<centerID.calo()<<" "<<centerID.area()<<" "<<centerID.row()<<" "<<centerID.col()<<std::endl;
  std::cout<<"get x "<<m_cgeom.get_x(centerID.area(), centerID.row())<<std::endl;
  std::cout<<"get y "<<m_cgeom.get_y(centerID.area(), centerID.col())<<std::endl;
  
  
  CaloNeighbors n_vector = m_ecal->zsupNeighborCells(centerID);
  LHCb::CaloCellID::Set n_set = LHCb::CaloCellID::Set(n_vector.begin(), n_vector.end());
  for ( CaloNeighbors::const_iterator neighbor =  n_vector.begin(); n_vector.end() != neighbor ; ++neighbor ){
      CaloNeighbors local_vector = m_ecal->zsupNeighborCells(*neighbor);
      LHCb::CaloCellID::Set new_set = LHCb::CaloCellID::Set(local_vector.begin(), local_vector.end());
      n_set.insert(new_set.begin(), new_set.end());
  }

  CaloNeighbors n_vector1 = m_ecal->neighborCells(centerID);
  LHCb::CaloCellID::Set n_set1 = LHCb::CaloCellID::Set(n_vector1.begin(), n_vector1.end());
  for ( CaloNeighbors::const_iterator neighbor =  n_vector1.begin(); n_vector1.end() != neighbor ; ++neighbor ){
      CaloNeighbors local_vector1 = m_ecal->neighborCells(*neighbor);
      LHCb::CaloCellID::Set new_set1 = LHCb::CaloCellID::Set(local_vector1.begin(), local_vector1.end());
      n_set1.insert(new_set1.begin(), new_set1.end());
  }

  CaloNeighbors additional_neib;
  additional_neib.reserve(n_set1.size() - n_set.size() + 1);
  for (const auto & elem : n_set1){
    if (n_set.find(elem) == n_set.end()){
      additional_neib.push_back(elem);
    }
  }

  std::cout<<n_set.size()<<std::endl;
  std::cout<<std::endl;

  std::vector<std::vector<double>> vector_cells (5, std::vector<double>(5, 0.0));
  
  std::vector<int> col_numbers = {(int)centerID.col() - 2, (int)centerID.col() - 1, (int)centerID.col(), (int)centerID.col() + 1, (int)centerID.col() + 2};
  std::vector<int> row_numbers = {(int)centerID.row() - 2, (int)centerID.row() - 1, (int)centerID.row(), (int)centerID.row() + 1, (int)centerID.row() + 2};

  for (auto& col_number: col_numbers){
    for (auto& row_number: row_numbers){
      const auto id_ = LHCb::CaloCellID(centerID.calo(), centerID.area(), row_number, col_number);
      auto * test = digits_full->object(id_);
      if (test) {
         vector_cells[col_number - (int)centerID.col() + 2][row_number - (int)centerID.row() + 2] = test->e();
      } else {
        if (n_set.size() < 25){
          std::cout<<std::endl<<"missed in bound"<<std::endl;
        } else {
          std::cout<<std::endl<<"missed simple"<<std::endl;
        }
        continue;
      }
    }
  }

  std::cout<<"vector_cells"<<std::endl;

    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
          std::cout<<i<<" "<<j<<" "<<vector_cells[i][j];
          std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }

  if (n_set.size() < 20){
    return false;
    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
          std::cout<<i<<" "<<j<<" "<<vector_cells[i][j];
          std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
  

    //makestd::vector<double> restored = restore_one(centerID.area(), additional_neib.begin()->area(), additional_neib);
    std::vector<std::vector<std::pair<int,int> > > missed = CheckVector(vector_cells);
    std::cout<<"missed size"<<missed.size()<<std::endl;
    std::cout<<"missed < 20 0"<<std::endl;
    for (auto & elem : missed[0]){
      std::cout<<elem.first<<" "<<elem.second<<std::endl;
    }

    std::cout<<"missed < 20 1"<<std::endl;
    for (auto & elem : missed[1]){
      std::cout<<elem.first<<" "<<elem.second<<std::endl;
    }
    
    
    std::cout<<std::endl;
  }


  if (n_set.size() < 25){
    std::cout<<n_set.size()<<std::endl;
    std::cout<<"full set "<<std::endl;
    for (auto & elem : n_set1){
      std::cout<<"r c a "<<elem.row()<<" "<<elem.col()<<" "<<elem.area()<<std::endl;
    }
    std::cout<<std::endl;

    for (auto & elem : n_set){
      std::cout<<"r c a "<<elem.row()<<" "<<elem.col()<<" "<<elem.area()<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"additional set "<<std::endl;
    for (auto & elem : additional_neib){
      std::cout<<"r c a "<<elem.row()<<" "<<elem.col()<<" "<<elem.area()<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"vector_cells1"<<std::endl;

    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
          std::cout<<i<<" "<<j<<" "<<vector_cells[i][j];
          std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
  
    std::vector<LHCb::CaloDigit> neighbor_for_restore;
    neighbor_for_restore.reserve(additional_neib.size()+1);
    auto prev_test = LHCb::CaloDigit(LHCb::CaloCellID(2,0,0,0), 0.0);
    for (const auto& id_: additional_neib){
        auto test = digits_full->object(id_);
        if (test){
          prev_test = *test;
          neighbor_for_restore.emplace_back(*test);
        }
        else {
          std::cout<<std::endl<<"missed in bound"<<std::endl;
          neighbor_for_restore.emplace_back(prev_test);
        }
    }

    std::cout<<"neighbor_for_restore"<<std::endl;
    for (auto & elem : neighbor_for_restore){
      std::cout<<elem<<std::endl;
    }
    
    std::cout<<std::endl;

    std::vector<double> restored = RestoreOne(centerID.area(), additional_neib.begin()->area(), neighbor_for_restore);
    std::vector<std::pair<int,int> > missed = CheckVector(vector_cells)[0];
    std::cout<<"missed size "<<missed.size()<<std::endl;
    std::cout<<"missed"<<std::endl;
    for (auto & elem : missed){
      std::cout<<elem.first<<" "<<elem.second<<std::endl;
    }
    std::cout<<std::endl;

    std::cout<<"restored"<<std::endl;
    for (auto & elem : restored){
      std::cout<<elem<<std::endl;
    }
    
    std::cout<<std::endl;

    std::cout<<"vector_cells3"<<std::endl;

    for (int i = 0; i < 5; i++){
        for (int j = 0; j < 5; j++){
          std::cout<<i<<" "<<j<<" "<<vector_cells[i][j];
          std::cout<<std::endl;
        }
        std::cout<<std::endl;
    }
  
    return true;
  } 

  
  for (int i = 0; i < 5; i++){
    for (int j = 0; j < 5; j++){
      rowEnergy[i*5 + j] = vector_cells[i][j];
    }
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