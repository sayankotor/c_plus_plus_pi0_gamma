#ifndef CALOCOSMICSTOOL_H
#define CALOCOSMICSTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTupleTool.h"
// from LHCb
#include "CaloInterfaces/ICaloCosmicsTool.h"            // Interface
#include "CaloDAQ/ICaloDataProvider.h"
#include "CaloKernel/CaloVector.h"


/** @class CaloCosmicsTool CaloCosmicsTool.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2008-04-07
 */
class CaloCosmicsTool : public extends<GaudiTupleTool, ICaloCosmicsTool> {
public:
  /// Standard constructor
  CaloCosmicsTool( const std::string& type,
               const std::string& name,
               const IInterface* parent);


  StatusCode initialize() override;

  Gaudi::XYZPoint referencePoint() override {return m_refPoint;}
  Gaudi::XYZPoint referencePointVariance() override {return m_eRefPoint;}
  const std::pair<Gaudi::XYZPoint,Gaudi::XYZPoint>& extrema() override {return m_bound;}
  double deposit() override {return m_adcSum;}
  double phi() override {return m_phi;}
  double phiVariance() override {return m_sPhi;}
  double asymmetry() override {return (m_maxLoc != "none") ? m_Rmean[m_maxLoc] : -99.;}
  double asymmetryVariance() override {return (m_maxLoc !="none") ? m_Rmean2[m_maxLoc] : -99.;}
  double slot() override {return m_offset;}
  double time() override {return m_time;}
  double timeVariance() override {return m_stime;}
  double timeDispersion() override {return (m_maxLoc !="none") ? m_td[m_maxLoc] : -99.;;}
  DeCalorimeter* det() override {return m_calo;}
  StatusCode processing() override;
  StatusCode tupling(unsigned int unit) override;
  bool tracked() override {return m_tracked;}
  bool timed() override {return m_timed;}
  double kernel() override {return m_kernel;};


protected:

private:
  // private methods
  StatusCode getBanks();
  StatusCode zSup();
  StatusCode fit2D();
  StatusCode timing();
  //
  DeCalorimeter* m_calo = nullptr;
  Gaudi::Plane3D m_plane;
  double m_delta;
  std::map<std::string, ICaloDataProvider*> m_daqs;          // decoding tool per time-slot
  LHCb::CaloAdc m_max;                        // highest deposit (sum over time-slots)
  std::vector<LHCb::CaloAdc> m_zsupADCs;      // vector of ADC sum over time-slots (after Zsup)
  std::vector<LHCb::CaloAdc> m_cosmicsADCs;   // vector of cellID associated to cosmics 2D track

  std::map<std::string,double> m_Rmean; // mean asymmetry
  std::map<std::string,double> m_Rmean2;// dispersion of the asymmetry
  std::map<std::string,double> m_R;     // global asymetry
  std::map<std::string,double> m_slotSum;
  std::map<std::string,double> m_td;


  double m_time;  // m_Rmean converted to time
  double m_stime;
  Gaudi::XYZPoint  m_refPoint;    // reference point on the cosmics 2D track
  Gaudi::XYZPoint  m_eRefPoint;   // error on reference point
  std::pair<Gaudi::XYZPoint,Gaudi::XYZPoint> m_bound; // extrema of the 2D track segment
  double m_phi;
  double m_sPhi;
  bool m_timed;
  bool m_tracked;
  long m_adcSum;
  std::vector<std::string> m_slots;        // full list of requested time slots
  std::string m_maxLoc;
  double m_offset;

  // properties
  std::string m_det;                                // detector name
  std::string m_readoutTool;                       // Name of the readout tool
  std::vector<std::string> m_seq;                  // sequence of time-slots to be used for trajectory computation
  std::vector<std::string> m_kern;                  // sequence of time-slots to be removed from kernel
  std::map<std::string, std::vector<std::string> > m_asy;    // pairs of time-slots to be used for asymmetry computation
  //
  long m_zSup   ;
  long m_zInf   ;
  long m_minD   ;
  long m_maxD   ;
  long m_minM   ;
  long m_maxM   ;
  double m_tol  ;
  double m_minR ;
  double m_maxR ;
  std::vector<double> m_par;
  double m_tRes;
  bool m_tuple;
  int m_maxAdc;
  bool m_full;
  double m_kernel;
};
#endif // CALOCOSMICSTOOL_H
