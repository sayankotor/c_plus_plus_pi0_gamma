// ============================================================================
#ifndef CALOGETTERTOOL_H
#define CALOGETTERTOOL_H 1
// ============================================================================
// Include files
// ============================================================================
// STD & STL
// ============================================================================
#include <algorithm>
// ============================================================================
// GaudiKernel
// ============================================================================
#include "GaudiKernel/IIncidentListener.h"
// ============================================================================
// GaudiAlg
// ============================================================================
#include "GaudiAlg/GaudiTool.h"
// ============================================================================
// Event
// ============================================================================
#include "Event/CaloDigit.h"
#include "Event/CaloCluster.h"
#include "Event/CaloHypo.h"
// ============================================================================
// Calo
// ============================================================================
#include "CaloDAQ/ICaloDataProvider.h"
#include "CaloInterfaces/ICaloGetterTool.h"            // Interface
#include "CaloInterfaces/ICounterLevel.h"
// ============================================================================
/** @class CaloGetterTool CaloGetterTool.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2009-04-17
 */
class CaloGetterTool : public extends<GaudiTool , ICaloGetterTool , IIncidentListener >
{
public:
  /// Standard constructor
  CaloGetterTool( const std::string& type,
                  const std::string& name,
                  const IInterface* parent);

  // ==========================================================================
  StatusCode initialize() override;
  StatusCode finalize() override;
  // ==========================================================================

  // update
  void update() override;
  int detectorsMask() override {return m_detMask;}



  // setters
  void addToDigits ( const std::string&loc,bool clear ) override
  {
    if ( clear )
    {
      m_digiLoc.clear () ;
      m_digits.clear  () ;
    }
    //
    if ( m_digiLoc.end() == std::find ( m_digiLoc.begin() ,
                                        m_digiLoc.end  () , loc ) )
    { m_digiLoc.value().push_back( loc ) ; }
    //
    if ( m_digits.end() == m_digits.find( loc ) )
    { m_digits[loc] = 0 ; }
    //
  }
  void addToClusters(const std::string& loc,bool clear ) override
  {
    if ( clear )
    {
      m_clusLoc.clear  () ;
      m_clusters.clear () ;
    }
    //
    if ( m_clusLoc.end() == std::find ( m_clusLoc.begin () ,
                                        m_clusLoc.end   () ,  loc ) )
    { m_clusLoc.value().push_back( loc ) ; }
    //
    if ( m_clusters.end() == m_clusters.find( loc ) )
    { m_clusters[loc] = 0 ; }
    //
  }
  void addToHypos  (const std::string& loc,bool clear ) override
  {
    if ( clear )
    {
      m_hypoLoc.clear () ;
      m_hypos.clear   () ;
    }
    //
    if ( m_hypoLoc.end() == std::find ( m_hypoLoc.begin () ,
                                        m_hypoLoc.end   () , loc ) )
    { m_hypoLoc.value().push_back( loc ) ; }
    //
    if ( m_hypos.end() == m_hypos.find( loc ) )
    { m_hypos[loc] = 0 ; }
    //
  }

  //getters
  LHCb::CaloDigits*   digits   ( const std::string& loc ) override;
  LHCb::CaloClusters* clusters ( const std::string& loc ) override;
  LHCb::CaloHypos*    hypos    ( const std::string& loc ) override;

  bool hasData(const std::string& det) override
  {
    std::map<std::string,bool>::iterator it = m_prov.find(det);
    return (it == m_prov.end()) ? false : (*it).second;
  }
  ICaloDataProvider* provider( const std::string& det) override
  {
    std::map<std::string,ICaloDataProvider*>::iterator it = m_provider.find(det);
    return (it == m_provider.end()) ? NULL : (*it).second;
  }

public:
  // =========================================================================
  /// Inform that a new incident has occurred
  void handle(const Incident& /* inc */ ) override { nullify() ; }
  // =========================================================================
protected:
  // =========================================================================
  /// nullify the pointers
  void nullify() ;
  // =========================================================================
private:
  ICounterLevel* counterStat = nullptr;
  // =========================================================================
  std::map<std::string,LHCb::CaloDigits*>   m_digits;
  std::map<std::string,LHCb::CaloClusters*> m_clusters;
  std::map<std::string,LHCb::CaloHypos*>    m_hypos;
  std::map<std::string,ICaloDataProvider*>  m_provider;
  std::map<std::string,bool> m_prov;
  Gaudi::Property<bool> m_digiUpd {this, "GetDigits", true};
  Gaudi::Property<bool> m_clusUpd {this, "GetClusters", false};
  Gaudi::Property<bool> m_hypoUpd {this, "GetHypos", false};
  //  bool m_provUpd;
  Gaudi::Property<std::vector<std::string>> m_digiLoc {this, "DigitLocations"};
  Gaudi::Property<std::vector<std::string>> m_clusLoc {this, "ClusterLocations"};
  Gaudi::Property<std::vector<std::string>> m_hypoLoc {this, "HypoLocations"};
  Gaudi::Property<int> m_detMask {this, "DetectorMask", 0xF};
};
// ============================================================================
// The END
// ============================================================================
#endif // CALOGETTERTOOL_H
// ============================================================================
