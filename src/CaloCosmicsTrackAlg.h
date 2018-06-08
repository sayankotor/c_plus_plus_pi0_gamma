#ifndef CALOCOSMICSTRACKALG_H 
#define CALOCOSMICSTRACKALG_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTupleAlg.h"
// from LHCb
#include "CaloInterfaces/ICaloCosmicsTrackTool.h"       


/** @class CaloCosmicsTrackAlg CaloCosmicsTrackAlg.h
 *  
 *
 *  @author Olivier Deschamps
 *  @date   2008-05-17
 */
class CaloCosmicsTrackAlg : public GaudiTupleAlg {
public: 
  /// Standard constructor
  using GaudiTupleAlg::GaudiTupleAlg;

  StatusCode initialize() override;    ///< Algorithm initialization
  StatusCode execute() override;    ///< Algorithm execution

private:
  ICaloCosmicsTrackTool* m_caloTrack = nullptr;

  Gaudi::Property<std::string> m_trackToolType 
    {this, "TrackTool", "CaloCosmicsTrackTool"};
  
  Gaudi::Property<std::string> m_forward 
    {this, "ForwardTrackContainer", LHCb::TrackLocation::CaloCosmicsForward};

  Gaudi::Property<std::string> m_backward 
    {this, "BackwardTrackContainer", LHCb::TrackLocation::CaloCosmicsBackward};

  Gaudi::Property<bool> m_monitor {this, "Monitor", false};

};
#endif // CALOCOSMICSTRACKALG_H
