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


  double m_time = -999;  // m_Rmean converted to time
  double m_stime = -999;
  Gaudi::XYZPoint m_refPoint;    // reference point on the cosmics 2D track
  Gaudi::XYZPoint m_eRefPoint;   // error on reference point
  std::pair<Gaudi::XYZPoint,Gaudi::XYZPoint> m_bound; // extrema of the 2D track segment
  double m_phi = -999;
  double m_sPhi = 999;
  bool m_timed = false;
  bool m_tracked = false;
  long m_adcSum;
  std::vector<std::string> m_slots;        // full list of requested time slots
  std::string m_maxLoc;
  double m_offset = -99999;
  double m_kernel = -1.;

  // properties
  Gaudi::Property<std::string> m_det {this, "Detector"};
  
  Gaudi::Property<std::string> m_readoutTool 
    {this, "ReadoutTool", "CaloDataProvider", "Name of the readout tool"};
  
  Gaudi::Property<std::vector<std::string>> m_seq
    {this, "TrajectorySlots", {"Prev1", "T0", "Next1"},
    "sequence of time-slots to be used for trajectory computation"};

  Gaudi::Property<std::vector<std::string>> m_kern
    {this, "RemoveSlotForKernel", {"T0"}, 
    "sequence of time-slots to be removed from kernel"};

  Gaudi::Property<std::map<std::string, std::vector<std::string> >> m_asy
    {this, "AsymmetrySlots",
    { 
      {"-25.", {"Prev1", "T0"}}, 
      {"0.", {"T0", "Next1"}} 
    }, 
    "pairs of time-slots to be used for asymmetry computation"};

  //
  Gaudi::Property<long> m_zSup    {this, "ZeroSuppression", 0, "Zero suppression threshold for hits ADC (sum over BX)"};
  Gaudi::Property<long> m_zInf    {this, "MaxSuppression", 99999, "Remove largest ADC"};
  Gaudi::Property<long> m_minD    {this, "MinCosmicsDeposit", 0, "minimal ADC sum over hits in cosmics track"};
  Gaudi::Property<long> m_maxD    {this, "MaxCosmicsDeposit", 99999, "maximal ADC sum over hits in cosmics track"};
  Gaudi::Property<long> m_minM    {this, "MinCosmicsMult", 0, "minimal multiplicity of hits in cosmics track"};
  Gaudi::Property<long> m_maxM    {this, "MaxCosmicsMult", 6156, "minimal multiplicity of hits in cosmics track"};
  Gaudi::Property<float> m_tol   {this, "MaxDistance2Line", 2.0, "maximal distance between hit and cosmics track (cellSize unit)"};

  // Timing
  Gaudi::Property<float> m_minR  {this, "MinR", 0., "minimal asymmetry range to compute time (absolute value)"};
  Gaudi::Property<float> m_maxR  {this, "MaxR", 0.8, "maximal asymmetry range to compute time (absolute value)"};
  Gaudi::Property<std::vector<float>> m_par {this, "RtoTime", {1.4, -0.7, 25.0, 0.19}, "parameters to convert R to time "};
  Gaudi::Property<float> m_tRes {this, "TRes", 0, "time resolution parameter"};

  // Tupling setup
  Gaudi::Property<bool> m_tuple  {this, "Ntupling"    , false, "produce ntuple"};
  Gaudi::Property<int>  m_maxAdc {this, "MaxArraySize", 500, "ntuple max array (# cosmics ADCs)"};
  Gaudi::Property<bool> m_full   {this, "AllDigits"   , false, "fill digit vector with all 0-sup ADC"};
};
#endif // CALOCOSMICSTOOL_H
