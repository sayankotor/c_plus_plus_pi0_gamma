// Include files

// from Gaudi
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"


// local
#include "CaloRelationsGetter.h"

//-----------------------------------------------------------------------------
// Implementation file for class : CaloRelationsGetter
//
// 2013-10-04 : Olivier Deschamps
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_COMPONENT( CaloRelationsGetter )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
CaloRelationsGetter::CaloRelationsGetter( const std::string& type,
                                          const std::string& name,
                                          const IInterface* parent )
: GaudiTool ( type, name , parent )
{
  declareInterface<ICaloRelationsGetter>(this);
  declareInterface<IIncidentListener>(this);
}

//=============================================================================


StatusCode CaloRelationsGetter::initialize(){
  StatusCode sc = GaudiTool::initialize();
  if( UNLIKELY( msgLevel(MSG::DEBUG) ) )debug() << "Initialize CaloRelationsGetter tool " << endmsg;
  // subscribe to the incidents
  IIncidentSvc* inc = incSvc() ;
  if ( inc )inc -> addListener  ( this , IncidentType::BeginEvent ) ;
  return sc;
}

StatusCode CaloRelationsGetter::finalize(){
  IIncidentSvc* inc = incSvc() ;
  if ( inc ) { inc -> removeListener  ( this ) ; }
  return GaudiTool::finalize();
}

LHCb::Calo2Track::ITrHypoTable2D* CaloRelationsGetter::getTrHypoTable2D(std::string location){
  auto it = m_trHypo2D.find(location);
  if( it == m_trHypo2D.end() ) m_trHypo2D[location] = getIfExists<LHCb::Calo2Track::ITrHypoTable2D> (location);
  return m_trHypo2D[location];
}

LHCb::Calo2Track::IHypoEvalTable* CaloRelationsGetter::getHypoEvalTable(std::string location){
  auto it = m_hypoEval.find(location);
  if( it == m_hypoEval.end()) m_hypoEval[location] = getIfExists<LHCb::Calo2Track::IHypoEvalTable> (location);
  return m_hypoEval[location];
}

LHCb::Calo2Track::IClusTrTable* CaloRelationsGetter::getClusTrTable(std::string location){
  auto it = m_clusTr.find(location);
  if( it == m_clusTr.end()) m_clusTr[location] = getIfExists<LHCb::Calo2Track::IClusTrTable> (location);
  return m_clusTr[location];
}

void CaloRelationsGetter::clean(){
  m_clusTr.clear();
  m_hypoEval.clear();
  m_trHypo2D.clear();

}

