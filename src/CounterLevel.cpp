#include "CounterLevel.h"

//-----------------------------------------------------------------------------
// Implementation file for class : CounterLevel
//
// 2016-08-13 : Olivier Deschamps
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_COMPONENT( CounterLevel )


//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
CounterLevel::CounterLevel( const std::string& type,
                            const std::string& name,
                            const IInterface* parent )
: base_class ( type, name , parent )
{
  declareInterface<ICounterLevel>(this);

  // sync m_isQuiet and m_isVerbose with m_clevel
  m_clevel.declareUpdateHandler(
    [=](const Property&) {
        this->m_isQuiet   = ( this->m_clevel > 0 );
        this->m_isVerbose = ( this->m_clevel > 1 );
  });
  m_clevel.useUpdateHandler();
}

//=============================================================================
