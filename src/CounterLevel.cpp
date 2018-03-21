// Include files 

// from Gaudi
#include "GaudiKernel/ToolFactory.h" 

// local
#include "CounterLevel.h"

//-----------------------------------------------------------------------------
// Implementation file for class : CounterLevel
//
// 2016-08-13 : Olivier Deschamps
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( CounterLevel )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
CounterLevel::CounterLevel( const std::string& type,
                            const std::string& name,
                            const IInterface* parent )
: base_class ( type, name , parent )
{
  declareInterface<ICounterLevel>(this);
  auto p = declareProperty("SetLevel", m_clevel);
  p->declareUpdateHandler( 
    [=](const Property&) {
        this->m_isQuiet   = ( this->m_clevel > 0 );
        this->m_isVerbose = ( this->m_clevel > 1 );
  });
  p->useUpdateHandler(); // sync m_isQuiet and m_isVerbose with m_clevel
}

//=============================================================================
