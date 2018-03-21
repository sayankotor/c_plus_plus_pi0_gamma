#ifndef COUNTERLEVEL_H 
#define COUNTERLEVEL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "CaloInterfaces/ICounterLevel.h"            // Interface


/** @class CounterLevel CounterLevel.h
 *  
 *
 *  @author Olivier Deschamps
 *  @date   2016-08-13
 */
class CounterLevel final : public extends<GaudiTool, ICounterLevel> {
public: 
  /// Standard constructor
  CounterLevel( const std::string& type, 
                const std::string& name,
                const IInterface* parent);

  bool isQuiet()      const override {return m_isQuiet;     };
  bool isVerbose()    const override {return m_isVerbose;   };
  bool isLevel(int l) const override {return m_clevel >= l; };
  int  level()        const override {return m_clevel     ; };
    
private:

  int m_clevel = 1; // quiet mode is the default
  bool m_isQuiet = true;
  bool m_isVerbose = false;
};
#endif // COUNTERLEVEL_H
