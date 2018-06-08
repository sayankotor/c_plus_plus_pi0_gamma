#ifndef GAMMAPI0SEPARATIONTOOL_H
#define GAMMAPI0SEPARATIONTOOL_H

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

//using namespace LHCb;

#include "TMV_MLP_inner.C"
#include "TMV_MLP_middle.C"
#include "TMV_MLP_outer.C"


/** @class GammaPi0SeparationTool GammaPi0SeparationTool.h
 *
 *
 *  @author Miriam Calvo Gomez
 *  @date   2010-03-29
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

class GammaPi0SeparationTool1 : public extends<GaudiTool, IGammaPi0SeparationTool>{
public:
  functor_cell comparer;
  /// Standard constructor
  GammaPi0SeparationTool1( const std::string& type,
                          const std::string& name,
                          const IInterface* parent);

  StatusCode initialize() override;
  StatusCode finalize() override;

  //double isPhoton(const LHCb::Particle* gamma);
  double isPhoton(const LHCb::CaloHypo* hypo) override;

  bool ClusterVariables(const LHCb::CaloHypo* hypo,
                        double& fr2, double& fasym, double& fkappa, double& fr2r4, double& etot,
                        double& Eseed, double& E2, int& area) override;

  bool PrsVariables(const LHCb::CaloHypo* hypo,
                    double& r2PS, double& asymPS, double& kappaPS, double& r2r4PS,
                    double& eSumPS, double& ePrs, double& eMaxPS, double& e2ndPS, double& ecornerPS,
                    int& multiPS, int& multiPS15, int& multiPS30, int& multiPS45) override;

  bool GetRawEnergy(const LHCb::CaloHypo* hypo, std::vector<double>& rowEnergy);

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
  double isPhoton(const double* v) override {
    return photonDiscriminant(int(v[0]),
                              v[1],v[2],v[3],
                              v[4],v[5],v[6],
                              v[7],v[8],v[9],v[10],
                              int(v[11]),int(v[12]),int(v[13]),int(v[14]));
  }

private:

  double m_minPt = 2000.;

  std::unique_ptr<IClassifierReader> m_reader0;
  std::unique_ptr<IClassifierReader> m_reader1;
  std::unique_ptr<IClassifierReader> m_reader2;

  std::unique_ptr<XGBClassifierPhPi0> m_xgb0;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb1;
  std::unique_ptr<XGBClassifierPhPi0> m_xgb2;

  const DeCalorimeter* m_ecal = nullptr;
  //const LHCb::CaloDigits* digits_full = nullptr;

  double photonDiscriminant(int area,
                            double r2, double r2r4, double asym,
                            double kappa, double Eseed, double E2,
                            double r2PS, double asymPS, double eMaxPS, double e2ndPS,
                            int multiPS, int multiPS15, int multiPS30, int multiPS45);

  double XGBDiscriminant(int area, std::vector<double>& row_energies);

  std::map<std::string,double> m_data;
  std::map<std::string,double> m_prsdata;
  double m_def = -1.e+06;
};
#endif // GAMMAPI0SEPARATIONTOOL1_H
