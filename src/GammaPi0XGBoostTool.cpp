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
  m_xgb0_b = std::make_unique<XGBClassifierPhPi0>("/afs/cern.ch/user/v/vchekali/public/DaVinci2/DaVinciDev/models/classes/0_bound.model");
  m_xgb1_b = std::make_unique<XGBClassifierPhPi0>("/afs/cern.ch/user/v/vchekali/public/DaVinci2/DaVinciDev/models/classes/1_bound.model");
  m_xgb2_b = std::make_unique<XGBClassifierPhPi0>("/afs/cern.ch/user/v/vchekali/public/DaVinci2/DaVinciDev/models/classes/2_bound.model");

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
  bool isBorder = false;
  bool rawEnergy = GetRawEnergy(hypo, isBorder, rawEnergyVector);

  if(!rawEnergy) return m_def;
  // return NN output
  std::cout<<"IS photon XGBoost new"<<std::endl;
  if (rawEnergyVector.size() != 25)
  {
    std::cout<<"else: raw energy size: "<<rawEnergyVector.size()<<std::endl;
  }
  double prediction_xgb = XGBDiscriminant(area, isBorder, rawEnergyVector);
  std::cout<<"prediction_xgb: "<<prediction_xgb<<std::endl;
  //return photonDiscriminant(area, fr2, fr2r4, fasym, fkappa, Eseed, E2,
                            //r2PS, asymPS, eMaxPS, e2ndPS, multiPS, multiPS15, multiPS30, multiPS45);
  return prediction_xgb;
}

//bool recreate(std::vector<std::vector<<double>>& initial_vector, ){

//}

double GammaPi0XGBoostTool::XGBDiscriminant(int area, bool isBorder, std::vector<double>& row_energies)
{

        double value = ( (area == 0 && isBorder == false) ? m_xgb0->getClassifierValues(row_energies)
                       : (area == 1 && isBorder == false) ? m_xgb1->getClassifierValues(row_energies)
                       : (area == 2 && isBorder == false) ? m_xgb2->getClassifierValues(row_energies)
                       : (area == 0 && isBorder == true) ? m_xgb0_b->getClassifierValues(row_energies)
                       : (area == 1 && isBorder == true) ? m_xgb1_b->getClassifierValues(row_energies)
                       : (area == 2 && isBorder == true) ? m_xgb2_b->getClassifierValues(row_energies)
                       : -1e10 );
        //info() << " INPUT TO GAMMA/PI0 : NN[" << input << "]= " << value << " (area="<<area<<")"<< endmsg;
        return value;
}

bool GammaPi0XGBoostTool::ClusterVariables(const LHCb::CaloHypo* hypo,
                                              double& fr2, double& fasym, double& fkappa, double& fr2r4, double& etot,
double& Eseed, double& E2, int& area) {
  return true;
}

bool GammaPi0XGBoostTool::GetRawEnergy(const LHCb::CaloHypo* hypo, bool isBorder, std::vector<double>& rowEnergy){
  if( NULL == hypo)return false;
  LHCb::CaloDigits * digits_full = getIfExists<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Ecal);
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );   // OD 2014/05 - change to Split Or Main  cluster
  if( NULL == cluster)return false;

  LHCb::CaloCellID centerID = cluster->seed();
  
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
  std::vector<std::vector<double>> vector_cells1 (5, std::vector<double>(5, 0.0));
  vector_cells1 = GetCluster(centerID, digits_full);
  
  std::vector<int> col_numbers = {(int)centerID.col() - 2, (int)centerID.col() - 1, (int)centerID.col(), (int)centerID.col() + 1, (int)centerID.col() + 2};
  std::vector<int> row_numbers = {(int)centerID.row() - 2, (int)centerID.row() - 1, (int)centerID.row(), (int)centerID.row() + 1, (int)centerID.row() + 2};


  if (n_set.size() < 25){
    return false;
  }

  else {
        for (auto& col_number: col_numbers){
            for (auto& row_number: row_numbers){
                const auto id_ = LHCb::CaloCellID(centerID.calo(), centerID.area(), row_number, col_number);
                auto * test = digits_full->object(id_);
                if (test) {
                    vector_cells[col_number - (int)centerID.col() + 2][row_number - (int)centerID.row() + 2] = test->e();
                } else {
                    continue;
                }
            }
          }
        }

  

  for (int i = 0; i < 5; i++){
      for (int j = 0; j < 5; j++){
        rowEnergy[i*5 + j] = vector_cells[i][j];
        std::cout<<i<<" "<<j<<" vector_cells "<<vector_cells[i][j]<<std::endl;;
        std::cout<<i<<" "<<j<<" vector_cells1 "<<vector_cells1[i][j]<<std::endl;
        if (abs(vector_cells[i][j] - vector_cells1[i][j]) > 0.001){
            std::cout<<"NOT EQUAL"<<std::endl;
     //std::cout<<"calo area row col "<<centerID.calo()<<" "<<centerID.area()<<" "<<centerID.row()<<" "<<centerID.col()<<std::endl;
     //std::cout<<"get x "<<m_cgeom.get_r(centerID.area(), centerID.row())<<std::endl;
     //std::cout<<"get y "<<m_cgeom.get_r(centerID.area(), centerID.col())<<std::endl;

            //std::cout<<i<<" "<<j<<" vector_cells "<<vector_cells[i][j]<<std::endl;;
            //std::cout<<i<<" "<<j<<" vector_cells1 "<<vector_cells1[i][j]<<std::endl;      
        }
        else {
            //std::cout<<"EQUAL"<<std::endl;
        }
      rowEnergy[i*5 + j] = vector_cells[i][j];
    }
  }
   
  //double a = vector_cells[200][200];
  //std::cout<<a;
  return true;
}

std::vector<std::vector<double>> GammaPi0XGBoostTool::GetCluster(LHCb::CaloCellID centerID, LHCb::CaloDigits * digits_full){
      int start_x = m_cgeom.get_r(centerID.area(), centerID.col());
      int start_y = m_cgeom.get_r(centerID.area(), centerID.row());
      int expected_area = centerID.area();
      int shift = m_cgeom.cell_size[expected_area];
      std::vector<std::vector<double>> vector_cells (5, std::vector<double>(5, 0.0));
      for (int i = -2; i < 3; i++){
        for (int j = -2; j < 3; j++){
          int local_x = start_x + i*shift;
          int local_y = start_y + j*shift;
          if (local_x < 0 || local_y < 0 || local_x > 383 || local_y > 383){
            continue;
          }
          for (int k = local_x; k < local_x + shift; k++){
            for (int l = local_y; l < local_y + shift; l++){
              int local_area = m_cgeom.c_geometry[k][l].first;
              int local_col = m_cgeom.get_R(centerID.area(), k);
            
              int local_row = m_cgeom.get_R(centerID.area(), l);
              const auto id_ = LHCb::CaloCellID(centerID.calo(), local_area, local_row, local_col);
              //std::cout<<"id_ calo area row col 2 "<<id_.calo()<<" "<<id_.area()<<" "<<local_row<<" "<<local_col<<std::endl;
              auto * test = digits_full->object(id_);
              if (test){
                vector_cells[i+2][j+2] +=  double(test->e())/m_cgeom.c_geometry[k][l].second;
              }
            }
          }
        }
      }
      return vector_cells;
  }
