// Include files

#include "CaloUtils/CaloMomentum.h"
#include "CaloUtils/CaloAlgUtils.h"
// from Event
#include "Event/RecHeader.h"
#include "Event/CaloHypo.h"
#include "Event/CaloCluster.h"
#include <unordered_map>
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>

// local
#include "GammaPi0XGBoostTool.h"
#include <fstream>

using namespace LHCb;
using namespace Gaudi::Units;
//-----------------------------------------------------------------------------
// Implementation file for class : GammaPi0XGBoostTool
//
// 2018-03-24 : author @sayankotor
//-----------------------------------------------------------------------------

//Declare nessesary help functionality
namespace {

const int CaloNCol[4] = {64, 32, 16, 16};
const int CaloNRow[4] = {52, 20, 12, 12};
const unsigned Granularity[3] = {1, 2, 3};

size_t getClusterType (const CaloCellID& id) {
  const unsigned ClusterSize = 5;
  int type (id.area());
  int xOffsetOut = std::min (int(id.col() - (32 - CaloNCol[type]*Granularity[type]/2)), // left edge
                 int(31 + CaloNCol[type]*Granularity[type]/2 - id.col())); // right edge
  int yOffsetOut = std::min (int(id.row() - (32 - CaloNRow[type]*Granularity[type]/2)),
                 int(31 + CaloNRow[type]*Granularity[type]/2 - id.row()));
  int innerWidth = CaloNCol[type+1] * (type != 2 ? Granularity[type] : 1); // process inner hole specially
  int innerHeight = CaloNRow[type+1] * (type != 2 ? Granularity[type] : 1); // process inner hole specially

  int xOffsetIn = std::min (int(id.col() - (31 - innerWidth/2)),
                int(32 + innerWidth/2 - id.col()));
  int yOffsetIn = std::min (int(id.row() - (31 - innerHeight/2)),
                int(32 + innerHeight/2 - id.row()));
  const int margin = (ClusterSize-1)/2;
  bool outerBorder = xOffsetOut < margin || yOffsetOut < margin;
  bool innerBorder = xOffsetIn > -margin && yOffsetIn > -margin;
  if (innerBorder) return type+3;
  else if (outerBorder) return type+6;
  return type;
}



struct functor_cell { 
  bool operator()(LHCb::CaloDigit& cell_a, LHCb::CaloDigit& cell_b) {
    if ((int)cell_a.cellID().area() == (int)cell_b.cellID().area()) {
      if ((int)cell_a.cellID().col() == (int)cell_b.cellID().col()) {
        return (int)cell_a.cellID().row() < (int)cell_b.cellID().row();
      }
      return ((int)cell_a.cellID().col() < (int)cell_b.cellID().col());
    }
    return ((int)cell_a.cellID().area() < (int)cell_b.cellID().area());
  }
};

struct cellID_hash {
  std::size_t operator () (const LHCb::CaloCellID &cell_id) const {
      auto h1 = std::hash<int>{}(cell_id.calo());
      auto h2 = std::hash<int>{}(cell_id.area());
      auto h3 = std::hash<int>{}(cell_id.row());
      auto h4 = std::hash<int>{}(cell_id.col());

      // Mainly for demonstration purposes, i.e. works but is overly simple
      // In the real world, use sth. like boost.hash_combine
      return h1 ^ h2 ^ h3 ^ h4;  
  }
};

struct calorimeter_geometry {
  int row_size = 384;
  int col_size = 384;
  std::vector<int> cell_size = {6, 3, 2};
  std::map<int,int> cell_area = {
                { 6, 0 },
                { 3, 1 },
                { 2, 2 } };

  std::vector<std::pair<double,double>> help = std::vector<std::pair<double,double>>(col_size, std::make_pair(0,36.0)); 
  std::vector<std::vector<std::pair<double,double>> > c_geometry = std::vector<std::vector<std::pair<double,double>> > (row_size, help);
  calorimeter_geometry(){
    this->init();
  }

  int get_R (int area, int r) {
      if (area == 0){
        return r/6;
      }
      else if (area == 1){
        return r/3 - 32;
      }
      else if (area == 2){
        return r/2 - 64;
      }
      return -10;
  }

  int get_r (int area, int R) {
      if (area == 0){
        return R*6;
      }
      else if (area == 1){
        return (R + 32)*3;
      }
      else if (area == 2){
        return (R + 64)*2;
      }
      return -10;
  }

private:
  bool init() {
    //define 1-area
    for (int i = 96; i < 288; i++){
      for (int j = 132; j< 252; j++){
        c_geometry[i][j] = std::make_pair(1, 9.0);
      }
    }
    //define 2-area
    for (int i = 144; i < 240; i++){
      for (int j = 156; j< 228; j++){
        c_geometry[i][j] = std::make_pair(2, 4.0);
      }
    }
    return true;
  }

};

}

// Declaration of the Tool Factory
DECLARE_COMPONENT( GammaPi0XGBoostTool )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
static calorimeter_geometry m_cgeom;

//=============================================================================
// Initialization
//=============================================================================
StatusCode GammaPi0XGBoostTool::initialize() {

  StatusCode sc = base_class::initialize(); // must be executed first
  if ( sc.isFailure() ) return sc;  // error printed already by GaudiAlgorithm

  if( UNLIKELY( msgLevel(MSG::DEBUG) ) ) debug() << "==> Initialize" << endmsg;

  /// Retrieve geometry of detector
  m_ecal = getDet<DeCalorimeter>( DeCalorimeterLocation::Ecal );


  const std::string paramEnv = "PARAMFILESROOT";   
  std::string paramRoot = std::string(getenv(paramEnv.c_str()));    

  m_xgb0 = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_simpl0.model");
  m_xgb1 = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_simpl1.model");
  m_xgb2 = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_simpl2.model");
  m_xgb03_b = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_bound03.model");
  m_xgb06_b = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_bound06.model");
  m_xgb14_b = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_bound14.model");
  m_xgb17_b = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_bound17.model");
  m_xgb25_b = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_bound25.model");
  m_xgb28_b = std::make_unique<XGBClassifierPhPi0>(paramRoot + "/data/GammaPi0XgbTool_bound28.model");

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
  m_xgb03_b.reset();
  m_xgb06_b.reset();
  m_xgb14_b.reset();
  m_xgb17_b.reset();
  m_xgb25_b.reset();
  m_xgb28_b.reset();

  return base_class::finalize(); // must be executed last
}

//=============================================================================
// Main execution
//=============================================================================


double GammaPi0XGBoostTool::isPhoton(const LHCb::CaloHypo* hypo){

  if ( !m_ecal ) return m_def;
  if ( LHCb::CaloMomentum(hypo).pt() < m_minPt) return m_def;
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );
  if (!cluster) return m_def;
  int area = cluster->seed().area();

  std::vector<double> rawEnergyVector(25, 0.0);

  int cluster_type = getClusterType(cluster->seed());
  bool rawEnergy = GetRawEnergy(hypo, cluster_type, rawEnergyVector);

  if(!rawEnergy) return m_def;
  
  double prediction_xgb = XGBDiscriminant(cluster_type, rawEnergyVector);
  return prediction_xgb;
}


double GammaPi0XGBoostTool::XGBDiscriminant(int cluster_type, std::vector<double>& row_energies)
{

        double value = -1e10;
        switch (cluster_type) {
            case 0: value = m_xgb0->getClassifierValues(row_energies);
                    break;
            case 1: value = m_xgb1->getClassifierValues(row_energies);
                    break;
            case 2: value = m_xgb2->getClassifierValues(row_energies);
                    break;
            case 3: value = m_xgb03_b->getClassifierValues(row_energies);
                    break;
            case 4: value = m_xgb14_b->getClassifierValues(row_energies);
                    break;
            case 5: value = m_xgb25_b->getClassifierValues(row_energies);
                    break;
            case 6: value = m_xgb06_b->getClassifierValues(row_energies);
                    break;
            case 7: value = m_xgb17_b->getClassifierValues(row_energies);
                    break;
            case 8: value = m_xgb28_b->getClassifierValues(row_energies);
                    break;
            default: if( UNLIKELY( msgLevel(MSG::WARNING) ) ) WARNING() << "GammaPi0XGBoostTool: Unsupperted cluster type" << endmsg;

        }
        return value;
}

bool GammaPi0XGBoostTool::GetRawEnergy(const LHCb::CaloHypo* hypo, int& cluster_type, std::vector<double>& rowEnergy) const{
  if( nullptr == hypo)return false;
  LHCb::CaloDigits * digits_full = getIfExists<LHCb::CaloDigits>(LHCb::CaloDigitLocation::Ecal);
  const LHCb::CaloCluster* cluster = LHCb::CaloAlgUtils::ClusterFromHypo( hypo );   // OD 2014/05 - change to Split Or Main  cluster
  
  if( nullptr == digits_full || nullptr == cluster) return false;
  
  LHCb::CaloCellID centerID = cluster->seed();

  std::vector<std::vector<double>> vector_cells (5, std::vector<double>(5, 0.0));
  
  std::vector<int> col_numbers = {(int)centerID.col() - 2, (int)centerID.col() - 1, (int)centerID.col(), (int)centerID.col() + 1, (int)centerID.col() + 2};
  std::vector<int> row_numbers = {(int)centerID.row() - 2, (int)centerID.row() - 1, (int)centerID.row(), (int)centerID.row() + 1, (int)centerID.row() + 2};

  if (cluster_type > 2)
  {
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

  else 
  { 
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
    }
  }
  return true;
}

std::vector<std::vector<double>> GammaPi0XGBoostTool::GetCluster(const LHCb::CaloCellID& centerID, const LHCb::CaloDigits * digits_full) const{
      int start_x = m_cgeom.get_r(centerID.area(), centerID.col());
      int start_y = m_cgeom.get_r(centerID.area(), centerID.row());
      int expected_area = centerID.area();
      int shift = m_cgeom.cell_size[expected_area];
      std::unordered_map<LHCb::CaloCellID, double, cellID_hash> hash_of_energy;
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
              int local_col = m_cgeom.get_R(local_area, k);           
              int local_row = m_cgeom.get_R(local_area, l);
              const auto id_ = LHCb::CaloCellID(centerID.calo(), local_area, local_row, local_col);
              std::unordered_map<LHCb::CaloCellID, double, cellID_hash>::const_iterator is_exist = hash_of_energy.find (id_);
              if (is_exist != hash_of_energy.end()){
                vector_cells[i+2][j+2] +=  is_exist -> second /m_cgeom.c_geometry[k][l].second;
                continue;
              }
              auto * test = digits_full->object(id_);
              if (test){
                vector_cells[i+2][j+2] +=  double(test->e())/m_cgeom.c_geometry[k][l].second;
                hash_of_energy.insert(std::make_pair(id_, double(test->e())));
              }
            }
          }
        }
      }
      return vector_cells;
}