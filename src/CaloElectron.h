#ifndef CALOELECTRON_H
#define CALOELECTRON_H 1

// Include files
#include "Part2Calo.h"

//from LHCb
#include "CaloUtils/ICaloElectron.h"
#include "CaloUtils/CaloMomentum.h"

// Forward declarations
namespace LHCb
{
  class ProtoParticle;
}


/** @class CaloElectron CaloElectron.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2006-11-30
 */
class CaloElectron : public extends<Part2Calo, ICaloElectron> {
public:
  /// Standard constructor
  CaloElectron( const std::string& type,
              const std::string& name,
              const IInterface* parent);

  bool  set(const  LHCb::Particle* particle,
            std::string det = DeCalorimeterLocation::Ecal,
            CaloPlane::Plane plane = CaloPlane::ShowerMax,
            double delta =0 ) override;
  bool  set(const  LHCb::ProtoParticle* proto,
            std::string det = DeCalorimeterLocation::Ecal,
            CaloPlane::Plane plane = CaloPlane::ShowerMax,
            double delta =0 ) override;

  LHCb::CaloHypo*    electron() override;
  LHCb::CaloHypo*    bremstrahlung() override;
  LHCb::CaloMomentum bremCaloMomentum() override;
  double ecalE() override;
  double eOverP() override;
  using ICaloElectron::closestState;
  LHCb::State closestState(std::string toWhat = "hypo") override;
  double caloTrajectoryZ(CaloPlane::Plane refPlane = CaloPlane::ShowerMax ,std::string toWhat = "hypo") override;
  double caloTrajectoryL(CaloPlane::Plane refPlane = CaloPlane::ShowerMax ,std::string toWhat = "hypo") override;



protected:
  bool caloSetting ();
private:
  LHCb::CaloHypo*            m_electron = nullptr;
  LHCb::CaloHypo*            m_bremstrahlung = nullptr;
  const LHCb::CaloPosition*  m_calopos = nullptr;
  double m_zOffset = 0;
};
#endif // CALOELECTRON_H
