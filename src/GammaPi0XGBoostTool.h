#ifndef GammaPi0XGBoostTool_H
#define GammaPi0XGBoostTool_H

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


class GammaPi0XGBoostTool : public extends<GaudiTool, IGammaPi0SeparationTool>{
public:
  /// Standard constructor
  GammaPi0XGBoostTool( const std::string& type,
                          const std::string& name,
                          const IInterface* parent);

  StatusCode initialize() override;
  StatusCode finalize() override;

  //double isPhoton(const LHCb::Particle* gamma);
  double isPhoton(const LHCb::CaloHypo* hypo) override;

  bool GetRawEnergy(const LHCb::CaloHypo* hypo, int& cluster_type, std::vector<double>& rowEnergy) const;
  std::vector<std::vector<double>> GetCluster(const LHCb::CaloCellID & centerID, const LHCb::CaloDigits * digits_full) const;

  double inputData(std::string) override { //@TODO: const-ify
    /* // try ecal data */
    /* auto it = m_data.find(data); */
    /* if( it != m_data.end() )return it->second; */
    /* // else try prs data : */
    /* auto itp = m_prsdata.find(data); */
    /* if( itp != m_prsdata.end() )return itp->second; */
    /* // yapa */
    return 0.;
  }
  
  std::map<std::string,double> inputDataMap() override {return std::map<std::string,double>();}
  std::map<std::string,double> inputPrsDataMap() override {return std::map<std::string,double>();}
  bool ClusterVariables(const LHCb::CaloHypo*, double&, double&, double&, double&, double&, 
      double&, double&, int&) override {return true;} 
  bool PrsVariables(const LHCb::CaloHypo*, double&, double&, double&, double&,
        double&, double&, double&, double&, double&, 
        int&, int&, int&, int&) override {return true;}
  double isPhoton(const double*) override {return 0.0;};

private:

  double m_minPt = 2000.;

  std::unique_ptr<XGBClassifierPhPi0> m_xgb0;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb1;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb2;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb03_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb14_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb25_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb06_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb17_b;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb28_b;

  const DeCalorimeter* m_ecal = nullptr;


  double XGBDiscriminant(int cluster_type, std::vector<double>& row_energies);

  const double m_def = -1.e+06;
};
#endif // GammaPi0XGBoostTool_H