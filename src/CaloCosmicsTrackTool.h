#ifndef CALOCOSMICSTRACKTOOL_H
#define CALOCOSMICSTRACKTOOL_H 1

// Include files
// from Gaudi
#include "GaudiKernel/IEventTimeDecoder.h"
#include "GaudiAlg/GaudiTupleTool.h"

// from LHCb
#include "CaloInterfaces/ICaloCosmicsTrackTool.h"            // Interface

/** @class CaloCosmicsTrackTool CaloCosmicsTrackTool.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2008-05-17
 */
class CaloCosmicsTrackTool : public extends<GaudiTupleTool, ICaloCosmicsTrackTool> {
public:
  /// Standard constructor
  CaloCosmicsTrackTool( const std::string& type,
                        const std::string& name,
                        const IInterface* parent);

  StatusCode initialize() override;
  //
  StatusCode processing() override;
  StatusCode tupling(unsigned int unit) override;
  double phi() override {return m_phi;}
  double phiVariance() override {return m_sPhi;}
  double theta() override {return m_theta;}
  double thetaVariance() override {return m_sTheta;}
  const Gaudi::XYZVector slopes() override {return m_slopes;}
  const Gaudi::XYZPoint  referencePoint() override { return m_refPt;}
  const Gaudi::SymMatrix3x3 slopesCovariance() override {return m_slopesCov;}
  const Gaudi::SymMatrix3x3 referencePointCovariance() override { return m_refPtCov;}
  bool   forward() override {return (m_dir == 1) ? true : false;}
  bool   tracked() override {return m_tracked;}
  double time() override {return m_time;}
  double timeVariance() override {return m_stime;}
  bool   timed() override {return  m_timed;}
  ICaloCosmicsTool* ecal() override {return m_eCosmics;}
  ICaloCosmicsTool* hcal() override {return m_hCosmics;}
  const LHCb::Track& track() override {return m_track;}
  StatusCode propagate(Gaudi::Plane3D plane) override;
  StatusCode propagate(double z ) override;

private:
  StatusCode matching();
  StatusCode fit3D();
  double dist2D(Gaudi::XYZPoint,Gaudi::XYZPoint);
  StatusCode buildTrack();

  //
  ICaloCosmicsTool* m_eCosmics = nullptr;
  ICaloCosmicsTool* m_hCosmics = nullptr;
  IEventTimeDecoder* m_odin = nullptr;
  //
  long m_run;
  long m_evt;
  long m_bx;
  double m_tx;
  double m_ty;
  Gaudi::XYZVector m_slopes;
  Gaudi::SymMatrix3x3 m_slopesCov;
  double m_stx;
  double m_sty;
  double m_theta;
  double m_sTheta;
  double m_phi;
  double m_sPhi;
  double m_chi2;
  Gaudi::Property<bool> m_intern {this, "UseInternalPlanes", false};
  Gaudi::Property<bool> m_extern {this, "UseExternalPlanes", false};
  double m_time;
  double m_stime;
  Gaudi::XYZPoint m_refPoint[2];
  Gaudi::SymMatrix3x3 m_refPointCov[2];
  Gaudi::XYZPoint m_refPt;
  Gaudi::SymMatrix3x3 m_refPtCov;
  int m_ref;

  int m_dir;
  bool m_tracked;
  bool m_timed;
  LHCb::Track m_track;
  // properties
  Gaudi::Property<std::string> m_cosmics {this, "CosmicsTool", "CaloCosmicsTool"};
  Gaudi::Property<float> m_chi2max {this, "MaxChi2", 15};
  Gaudi::Property<float> m_chi2min {this, "MinChi2", 0};
  Gaudi::Property<float> m_fac {this, "Factor", 1.5};
  Gaudi::Property<bool> m_tuple {this, "Ntupling", false, "produce ntuple"};
  Gaudi::Property<std::string> m_timer {this, "Timer", "EcalElseHcal"};
};
#endif // CALOCOSMICSTRACKTOOL_H
