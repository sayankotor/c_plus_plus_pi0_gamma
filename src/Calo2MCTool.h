#ifndef CALO2MCTOOL_H
#define CALO2MCTOOL_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "CaloInterfaces/ICalo2MCTool.h"            // Interface
#include "CaloInterfaces/ICounterLevel.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"
// from LHCb : the calo->particle stream
#include "Event/CaloDigit.h"
#include "Event/CaloCluster.h"
#include "Event/CaloHypo.h"
#include "Event/ProtoParticle.h"
#include "Event/Particle.h"
#include "Event/MCParticle.h"
#include "Kernel/IParticlePropertySvc.h"
#include "Kernel/ParticleProperty.h"
#include "CaloUtils/CaloMomentum.h"
// from LHCb : some  typedef and utilities
#include "CaloUtils/Calo2MC.h"
#include "CaloUtils/CaloMCTools.h"



/** @class Calo2MCTool Calo2MCTool.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2009-07-27
 */
class Calo2MCTool : public extends<GaudiTool , ICalo2MCTool , IIncidentListener > {

public:
  // category enum

  enum MCCategory  {
    /* Single cluster/hypo categories */
    UnMatched           = 0,
    Photon              = 1,
    BremStrahlung       = 2,
    Electron            = 3,
    ConvertedPhoton     = 4,
    MergedPi0           = 5,
    ChargedHadron       = 6,
    NeutralHadron       = 7,
    Spillover           = 8
  };


  /// Standard constructor
  Calo2MCTool( const std::string& type,
               const std::string& name,
               const IInterface* parent);

  StatusCode initialize( ) override;

  //
  ICalo2MCTool* from(const LHCb::CaloDigit*     digit    ) override;
  ICalo2MCTool* from(const LHCb::CaloCluster*   cluster  ) override;
  ICalo2MCTool* from(const LHCb::CaloHypo*      hypo     ) override;
  ICalo2MCTool* from(const LHCb::ProtoParticle* proto    ) override;
  ICalo2MCTool* from(const LHCb::Particle*      particle ) override;
  const LHCb::MCParticle* bestMC() const override;
  const LHCb::MCParticle* maxMC() const override;
  double weight(const LHCb::MCParticle*) const override;
  const LHCb::MCParticle* findMC(LHCb::ParticleID id, double threshold = 0 ) const override;
  const LHCb::MCParticle* findMCOrBest(LHCb::ParticleID id, double threshold = 0 ) const override;
  const LHCb::MCParticle* findMC(std::string name, double threshold = 0 ) const override;
  const LHCb::MCParticle* findMCOrBest(std::string name, double threshold = 0 ) const override;
  double quality(const LHCb::MCParticle*) const override;
  std::string descriptor() const override;
  bool isCalo(LHCb::Particle* particle) const override {
    return particle && LHCb::CaloParticle( particle ).isCalo();
  }
  bool isPureNeutralCalo(const LHCb::Particle* particle) const override {  // SHOULD BE IN CALOPARTICLE
    return particle && LHCb::CaloParticle( const_cast<LHCb::Particle*>(particle) ).isPureNeutralCalo();// non pure calorimetric object
  }



  /// Inform that a new incident has occurred
  void handle(const Incident& /* inc */ ) override { clear();}


  // TO BE INTERFACED :
  int MCCategory(){return m_category;};
  ICalo2MCTool* fragment(unsigned int i);
  unsigned int numberOfFragments(){return m_nFrag;}
  //ostream << category
  // clusters()
  // hypos()
  // protos()

private:
  ICounterLevel* counterStat = nullptr;
  StatusCode process();
  void addDigit   (const LHCb::CaloDigit* digit);
  void addCluster (const LHCb::CaloCluster* cluster);
  void addHypo    (const LHCb::CaloHypo* hypo);
  void addProto   (const LHCb::ProtoParticle* proto, const LHCb::Particle* parent = NULL);
  void clear();
  void mcDigest();
  void mcTree(const LHCb::MCParticle* part, std::vector<const LHCb::MCParticle*>& tree , std::string& sTree);


  std::vector<const LHCb::CaloDigit*>     m_digits;
  std::vector<const LHCb::CaloCluster*>   m_clusters;
  std::vector<const LHCb::CaloHypo*>      m_hypos;
  std::vector<const LHCb::ProtoParticle*> m_protos;
  std::vector<const LHCb::Particle*>      m_parts;
  LHCb::CaloDigit* m_digit = nullptr;
  LHCb::CaloCluster* m_cluster = nullptr;
  LHCb::CaloHypo* m_hypo = nullptr;
  LHCb::ProtoParticle* m_proto = nullptr;
  LHCb::Particle* m_part = nullptr;
  CaloMCTools::CaloMCMap m_mcMap;
  //
  Gaudi::Property<std::string> m_cluster2MCLoc 
    {this, "Cluster2MCTable", "", "Cluster->MC relation table location"};

  Gaudi::Property<std::string> m_digit2MCLoc 
    {this, "Digit2MCTable", "Relations/" + LHCb::CaloDigitLocation::Default, 
    "Digit->MC relation table location"};

  LHCb::Calo2MC::IClusterTable* m_cluster2MC  = nullptr;
  LHCb::Calo2MC::IDigitTable*   m_digit2MC  = nullptr;
  double m_sum = 0.;
  const LHCb::MCParticle* m_maxMC = nullptr;
  const LHCb::MCParticle* m_bestMC = nullptr;
  std::map<std::string,std::vector<const LHCb::MCParticle*> > m_treeMap;
  SmartIF<LHCb::IParticlePropertySvc> m_ppsvc;

  //
  Gaudi::Property<bool> m_hypo2Cluster 
    {this, "Hypo2Cluster", false, "(part->protoP)->hypo->cluster cascade ( or use  hypo->MC linker tables)"};

  Gaudi::Property<bool> m_cluster2Digit 
    {this, "Cluster2Digits", false, "(part->protoP)->hypo->cluster->digit cascade ( or use cluster->MC relation tables)"};

  Gaudi::Property<bool> m_merged2Split 
    {this, "Merged2Split", false, "expert usage (merged->split->cluster - NOT IMPLEMENTED so far)"};

  Gaudi::Property<int> m_sFilter 
    {this, "DigitStatusFilter", LHCb::CaloDigitStatus::UseForEnergy,
    "digit filter in case of ..->digit->MC table is used"};

  int m_category = 0;
  int m_depth = -1;
  unsigned int m_nFrag = 0;
  ICalo2MCTool* m_tool = nullptr;
};
#endif // CALO2MCTOOL_H
