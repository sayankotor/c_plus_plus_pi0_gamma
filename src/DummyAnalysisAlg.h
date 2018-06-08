// $Id: DummyAnalysisAlg.h,v 1.1 2009-10-02 07:24:53 cattanem Exp $
#ifndef DUMMYANALYSISALG_H
#define DUMMYANALYSISALG_H 1

// Include files
// from DaVinci, this is a specialized GaudiAlgorithm
#include "Kernel/DaVinciTupleAlgorithm.h"
#include "LoKi/IMCDecay.h"
#include "CaloInterfaces/IGammaPi0SeparationTool.h"

/** @class DummyAnalysisAlg DummyAnalysisAlg.h
 *
 *  Dummy algorithm to allow creation of empty library in Tutorial example
 *
 *  @author Marco Cattaneo
 *  @date   2009-10-02
 */
class DummyAnalysisAlg : public DaVinciTupleAlgorithm
{

public:
  
  /// Standard constructor
  DummyAnalysisAlg( const std::string& name, ISvcLocator* pSvcLocator );

  virtual ~DummyAnalysisAlg( ); ///< Destructor

  StatusCode initialize() override;    ///< Algorithm initialization
  StatusCode execute   () override;    ///< Algorithm execution
  StatusCode finalize  () override;    ///< Algorithm finalization

protected:

private:
  Decays::IMCDecay::Finder m_finder_BKPiGamma;
  Decays::IMCDecay::Finder m_finder_BKPiPi0;
  IGammaPi0SeparationTool* m_GammaPi0;
  IGammaPi0SeparationTool* m_GammaPi0_XGB;

};
#endif // DUMMYANALYSISALG_H
