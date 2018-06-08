#ifndef CALOHYPOESTIMATOR_H
#define CALOHYPOESTIMATOR_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "GaudiKernel/IIncidentListener.h"
#include "GaudiKernel/IIncidentSvc.h"
#include "GaudiKernel/Incident.h"
#include "Event/CaloDataFunctor.h"
#include "CaloUtils/ClusterFunctors.h"
#include "CaloUtils/CaloMomentum.h"
#include "CaloUtils/CaloAlgUtils.h"
#include "CaloUtils/Calo2Track.h"
#include "Relations/Relation2D.h"
#include "Relations/IRelationWeighted.h"
#include "Relations/IRelationWeighted2D.h"
#include "Event/Track.h"
#include "CaloInterfaces/ICaloHypoEstimator.h"            // Interface
#include "CaloInterfaces/ICounterLevel.h"
#include "CaloInterfaces/IGammaPi0SeparationTool.h"
#include "CaloInterfaces/INeutralIDTool.h"
#include "CaloUtils/ICaloElectron.h"
#include "CaloInterfaces/ICaloRelationsGetter.h"
#include "CaloDet/DeCalorimeter.h"

/** @class CaloHypoEstimator CaloHypoEstimator.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2010-08-18
 */


class CaloHypoEstimator : public GaudiTool, virtual public ICaloHypoEstimator, virtual public IIncidentListener {
public:
  /// Standard constructor
  CaloHypoEstimator( const std::string& type,
                     const std::string& name,
                     const IInterface* parent);

  StatusCode initialize() override;
  StatusCode finalize() override;

  double data(const LHCb::CaloCluster* cluster ,CaloDataType::DataType type, double def = CaloDataType::Default) override;
  double data(const LHCb::CaloHypo* hypo ,CaloDataType::DataType type, double def = CaloDataType::Default) override;

  void handle(const Incident&  ) override {
    if( UNLIKELY( msgLevel(MSG::DEBUG) ) )
      debug() << "IIncident Svc reset" << endmsg;
    clean();
  }
  ICaloHypo2Calo* hypo2Calo() override {return m_toCalo;};
  const LHCb::Track* toTrack(CaloMatchType::MatchType match) override {
    caloMatchType::iterator it = m_track.find( match );
    if( it == m_track.end() )return NULL;
    return (*it).second;
  }
  const LHCb::CaloCluster* toCluster(CaloClusterType::ClusterType clus=CaloClusterType::SplitOrMain) override {
    caloClusterType::iterator it = m_clusters.find( clus );
    if( it == m_clusters.end() )return NULL;
    return (*it).second;
  }


  StatusCode  _setProperty(const std::string& p,const std::string& v) override {return  setProperty(p,v);};
  bool status() override {return m_status;}


private:
  ICounterLevel* counterStat = nullptr;
  bool estimator(const LHCb::CaloCluster* cluster,const LHCb::CaloHypo* fromHypo=NULL);
  bool estimator(const LHCb::CaloHypo* hypo);
  void clean();
  
  Gaudi::Property<bool> m_extrapol {this, "Extrapolation", true};
  Gaudi::Property<bool> m_seed {this, "AddSeed", false};
  Gaudi::Property<bool> m_neig {this, "AddNeighbors", false};
  
  ICaloHypo2Calo* m_toCalo = nullptr;
  caloDataType m_data;
  caloMatchType m_track;
  caloClusterType m_clusters;
  LHCb::CaloHypo* m_hypo = nullptr;
  LHCb::CaloCluster* m_cluster = nullptr;
  
  Gaudi::Property<std::map<std::string,std::string>> m_pidLoc {this, "NeutralIDLocations"};
  Gaudi::Property<std::string> m_cmLoc {this, "ClusterMatchLocation"};
  Gaudi::Property<std::string> m_emLoc {this, "ElectronMatchLocation"};
  Gaudi::Property<std::string> m_bmLoc {this, "BremMatchLocation"};
  Gaudi::Property<bool> m_skipC {this, "SkipChargedID", false};
  Gaudi::Property<bool> m_skipN {this, "SkipNeutralID", false};
  Gaudi::Property<bool> m_skipCl {this, "SkipClusterMatch", false};

  std::map<std::string,LHCb::Calo2Track::IHypoEvalTable*> m_idTable;
  bool m_status = true;
  ICaloElectron * m_electron = nullptr;
  IGammaPi0SeparationTool* m_GammaPi0 = nullptr;
  IGammaPi0SeparationTool* m_GammaPi0XGB = nullptr;
  INeutralIDTool* m_neutralID = nullptr;
  ICaloRelationsGetter*    m_tables = nullptr;
  DeCalorimeter* m_ecal = nullptr;
};
#endif // CALOHYPOESTIMATOR_H
