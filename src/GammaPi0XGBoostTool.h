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
      if ((int)cell_a.cellID().col() == (int)cell_b.cellID().col()) {
        return (int)cell_a.cellID().row() < (int)cell_b.cellID().row();
      }
      return ((int)cell_a.cellID().col() < (int)cell_b.cellID().col());
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

  bool GetRawEnergy(const LHCb::CaloHypo* hypo, std::vector<double>& rowEnergy);

  bool ClusterVariables(const LHCb::CaloHypo* hypo,
                        double& fr2, double& fasym, double& fkappa, double& fr2r4, double& etot,
                        double& Eseed, double& E2, int& area) override;

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
  bool PrsVariables(const LHCb::CaloHypo* hypo, double& r2PS, double& asymPS, double& kappaPS, double& r2r4PS,                                           double& eSumPS, double& ePrs, double& eMaxPS, double& e2ndPS, double& ecornerPS, 
    int& multiPS, int& multiPS15, int& multiPS30, int& multiPS45){return true;};
  double isPhoton(const double* v){return 0.0;};


private:

  double m_minPt = 2000.;

  std::unique_ptr<XGBClassifierPhPi0> m_xgb0;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb1;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb2;

  const DeCalorimeter* m_ecal = nullptr; 

  double XGBDiscriminant(int area, std::vector<double>& row_energies);

  std::map<std::string,double> m_data;
  std::map<std::string,double> m_prsdata;
  double m_def = -1.e+06;
};
#endif // GammaPi0XGBoostTool_H
