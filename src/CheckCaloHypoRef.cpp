// Include files

// from  LHCb
#include "CaloUtils/CaloAlgUtils.h"
#include "Event/CaloHypo.h"
// local
#include "CheckCaloHypoRef.h"

//-----------------------------------------------------------------------------
// Implementation file for class : CheckCaloHypoRef
//
// 2012-05-14 : Olivier Deschamps
//-----------------------------------------------------------------------------

// Declaration of the Algorithm Factory
DECLARE_COMPONENT( CheckCaloHypoRef )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
CheckCaloHypoRef::CheckCaloHypoRef( const std::string& name,
                                    ISvcLocator* pSvcLocator)
: GaudiAlgorithm ( name , pSvcLocator )
{
  using namespace LHCb::CaloAlgUtils;
  m_inputs.value() = {
    CaloHypoLocation("Photons"   , context()),
    CaloHypoLocation("Electrons" , context()),
    CaloHypoLocation("MergedPi0s", context()) 
  };
}

//=============================================================================
// Main execution
//=============================================================================
StatusCode CheckCaloHypoRef::execute() {

 if ( msgLevel(MSG::DEBUG) ) debug() << "==> Execute" << endmsg;

 for (const auto& loc :  m_inputs ) {
    const LHCb::CaloHypos* hypos = getIfExists<LHCb::CaloHypos> (loc);
    if ( !hypos ) continue;
    if(counterStat->isQuiet()) counter("#Hypos in " + loc) += hypos->size();
    int bLink=0;
    for (const auto& h : *hypos ) {
      for( const auto&  cluster : h->clusters() ) {
        if( !cluster ) bLink++;
        else if(counterStat->isVerbose())counter("Cluster energy " +loc)+=cluster->e();
      }
    }
    if(counterStat->isQuiet())counter("Broken SmarRef " +loc) += bLink;
    if(bLink != 0)Warning("CaloHypo -> CaloCluster* SmartReference is broken for "+loc,StatusCode::SUCCESS).ignore();
 }
  return StatusCode::SUCCESS;
}

//=============================================================================
