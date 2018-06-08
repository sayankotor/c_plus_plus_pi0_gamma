#ifndef CHECKCALOHYPOREF_H
#define CHECKCALOHYPOREF_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiAlgorithm.h"
#include "CaloInterfaces/ICounterLevel.h"

/** @class CheckCaloHypoRef CheckCaloHypoRef.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2012-05-14
 */
class CheckCaloHypoRef : public GaudiAlgorithm {
public:
  /// Standard constructor
  CheckCaloHypoRef( const std::string& name, ISvcLocator* pSvcLocator );

  StatusCode execute() override;    ///< Algorithm execution

private:
  Gaudi::Property<std::vector<std::string>> m_inputs {this, "CaloHypos"};
  ICounterLevel* counterStat = nullptr;
};
#endif // CHECKCALOHYPOREF_H
