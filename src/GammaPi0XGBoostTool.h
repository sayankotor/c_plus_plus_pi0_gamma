#ifndef GammaPi0XGBoostTool_H
#define GammaPi0XGBoostTool_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "CaloDet/DeCalorimeter.h"
#include "CaloInterfaces/IGammaPi0SeparationTool.h"

// Math
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "LHCbMath/Line.h"
#include "LHCbMath/GeomFun.h"

//using xgb in c++

#include "XGBClassifierPhPi0.h"


/** @class GammaPi0XGBoostTool GammaPi0XGBoostTool.h
 *
 *
 *  @author @sayankotor
 *  @date   2018-03-24
 */

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

  int get_R (int area, int r) {
      //std::cout<<std::endl;
      //std::cout<<area<<" "<<r<<std::endl;
      if (area == 0){
        //std::cout<<0<<std::endl;
        return r/6;
      }
      else if (area == 1){
        //std::cout<<1<<std::endl;
        return r/3 - 32;
      }
      else if (area == 2){
        //std::cout<<2<<std::endl;
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

};


class GammaPi0XGBoostTool : public extends<GaudiTool, IGammaPi0SeparationTool>{
public:
  functor_cell comparer;
  /// Standard constructor
  GammaPi0XGBoostTool( const std::string& type,
                          const std::string& name,
                          const IInterface* parent);

  StatusCode initialize() override;
  StatusCode finalize() override;

  //double isPhoton(const LHCb::Particle* gamma);
  double isPhoton(const LHCb::CaloHypo* hypo) override;

  bool GetRawEnergy(const LHCb::CaloHypo* hypo, bool isBorder, std::vector<double>& rowEnergy);
  std::vector<std::vector<double>> GetCluster(LHCb::CaloCellID centerID, LHCb::CaloDigits * digits_full);

  double inputData(std::string data) override { //@TODO: const-ify
    // try ecal data
    auto it = m_data.find(data);
    if( it != m_data.end() )return it->second;
    // else try prs data :
    auto itp = m_prsdata.find(data);
    if( itp != m_prsdata.end() )return itp->second;
    // yapa
    return 0.;
  }
  std::map<std::string,double> inputDataMap() override {return m_data;}
  std::map<std::string,double> inputPrsDataMap() override {return m_prsdata;}
  bool ClusterVariables(const LHCb::CaloHypo* hypo, double& fr2, double& fasym, double& fkappa, double& fr2r4, double& etot,
  double& Eseed, double& E2, int& area);
  bool PrsVariables(const LHCb::CaloHypo* hypo, double& r2PS, double& asymPS, double& kappaPS, double& r2r4PS,                                           double& eSumPS, double& ePrs, double& eMaxPS, double& e2ndPS, double& ecornerPS, 
    int& multiPS, int& multiPS15, int& multiPS30, int& multiPS45){return true;};
  double isPhoton(const double* v){return 0.0;};


private:

  double m_minPt = 2000.;

  std::unique_ptr<XGBClassifierPhPi0> m_xgb0;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb1;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb2;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb0_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb1_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb2_b;

  const DeCalorimeter* m_ecal = nullptr; 


  double XGBDiscriminant(int area, bool isBorder, std::vector<double>& row_energies);

  std::map<std::string,double> m_data;
  std::map<std::string,double> m_prsdata;
  calorimeter_geometry m_cgeom;
  double m_def = -1.e+06;
};
#endif // GammaPi0XGBoostTool_H
