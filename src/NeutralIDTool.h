// $Id: $
#ifndef NEUTRALIDTOOL_H
#define NEUTRALIDTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "CaloDet/DeCalorimeter.h"
#include "CaloInterfaces/INeutralIDTool.h"
#include "CaloInterfaces/ICaloHypoEstimator.h"
// Math
#include "GaudiKernel/Point3DTypes.h"
#include "GaudiKernel/Vector3DTypes.h"
#include "LHCbMath/Line.h"
#include "LHCbMath/GeomFun.h"

//using namespace LHCb;

#include "TMV_MLP_H.C"
#include "TMV_MLP_E.C"

/** @class neutralIDTool neutralIDTool.h
 *
 *
 *  @author Mostafa HOBALLAH
 *  @date   2013-07-25
 */
class NeutralIDTool : public GaudiTool, virtual public INeutralIDTool{
public:

  // Return the interface ID
  //  static const InterfaceID& interfaceID() { return IID_NeutralIDTool; }

  /// Standard constructor
  NeutralIDTool( const std::string& type,
                          const std::string& name,
                          const IInterface* parent);

  StatusCode initialize() override;
  StatusCode finalize() override;

  //double isPhoton(const LHCb::Particle* gamma);
  double isNotE(const LHCb::CaloHypo* hypo,ICaloHypoEstimator* e=NULL) override;
  double isNotH(const LHCb::CaloHypo* hypo,ICaloHypoEstimator* e=NULL) override;
  void Variables(const LHCb::CaloHypo* hypo, double& clmatch,
                         double& prse, double& e19, double& hclecl, double& prse19,double& prse49,
                         double& sprd, double& prse4mx, double& prsm, double& spdm) override;


  double isNotE(const double* v) override {
    return photonDiscriminantE(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9]);
  }
  double isNotH(const double* v) override {
    return photonDiscriminantH(v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9]);
  }

private:

  Gaudi::Property<float> m_minPt {this, "MinPt", 75.};

  std::unique_ptr<IClassifierReader> m_reader0;
  std::unique_ptr<IClassifierReader> m_reader1;



  //Gaudi::XYZPoint m_vertex;
  //Gaudi::Plane3D m_planePrs;

  double photonDiscriminantE(double clmatch, double prse, double e19, double hclecl, double prse19,
                             double prse49, double sprd, double prse4mx, double prsm, double spdm);

  double photonDiscriminantH(double clmatch, double prse, double e19, double hclecl, double prse19,
                              double prse49, double sprd, double prse4mx, double prsm, double spdm);


  ICaloHypoEstimator* m_estimator = nullptr;

};
#endif // NEUTRALIDTOOL_H
