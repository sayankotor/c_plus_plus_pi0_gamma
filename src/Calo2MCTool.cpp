// Include files
#include "boost/algorithm/string.hpp"

// from LHCb
#include "Relations/RelationWeighted1D.h"
#include "CaloUtils/CaloParticle.h"
#include "CaloUtils/Calo2MC.h"
#include "GaudiKernel/IRegistry.h"
#include "Linker/LinkedTo.h"
#include "Linker/LinkedFrom.h"

// local
#include "Calo2MCTool.h"

//-----------------------------------------------------------------------------
// Implementation file for class : Calo2MCTool
//
// 2009-07-27 : Olivier Deschamps
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( Calo2MCTool )

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
Calo2MCTool::Calo2MCTool( const std::string& type,
                          const std::string& name,
                          const IInterface* parent )
  : base_class ( type, name , parent )
{
  declareInterface<ICalo2MCTool>(this);

  declareProperty ( "Cluster2MCTable" , m_cluster2MCLoc ) ; // Cluster->MC relation table location
  declareProperty ( "Digit2MCTable"   , m_digit2MCLoc = "Relations/" + LHCb::CaloDigitLocation::Default ) ;  // Digit->MC relation table location
  declareProperty ( "Hypo2Cluster"    , m_hypo2Cluster  = false   ); // (part->protoP)->hypo->cluster cascade ( or use  hypo->MC linker tables)
  declareProperty ( "Cluster2Digits"  , m_cluster2Digit = false  );  // (part->protoP)->hypo->cluster->digit cascade ( or use cluster->MC relation tables)
  declareProperty ( "Merged2Split"    , m_merged2Split  = false  );  // expert usage (merged->split->cluster - NOT IMPLEMENTED so far)
  declareProperty ( "DigitStatusFilter"  , m_sFilter    = LHCb::CaloDigitStatus::UseForEnergy       ) ; // digit filter in case of ..->digit->MC table is used

  m_cluster2MCLoc = "Relations/" + ( boost::iequals( context(), "HLT")
                                     ? LHCb::CaloClusterLocation::DefaultHlt
                                     : LHCb::CaloClusterLocation::Default );
}

//=============================================================================


StatusCode Calo2MCTool::initialize(){
  StatusCode sc = base_class::initialize();
  if (sc.isFailure()) return Error("Failed to initialize", sc);
  m_ppsvc = service( "LHCb::ParticlePropertySvc", true);
  if( !m_ppsvc )return StatusCode::FAILURE;

  // incidentSvc
  IIncidentSvc* inc = incSvc() ;
  if ( 0 != inc )inc -> addListener  ( this , IncidentType::BeginEvent ) ;


  // init
  counterStat = tool<ICounterLevel>("CounterLevel");
  m_digit   = nullptr;
  m_cluster = nullptr;
  m_hypo    = nullptr;
  m_proto   = nullptr;
  m_part    = nullptr;
  m_digit2MC= nullptr;
  m_cluster2MC=nullptr;
  m_category = 0 ;
  m_depth    = -1 ;
  //m_hypo2MC = nullptr;
  // get fragment tool and propagate properties
  m_tool = tool<ICalo2MCTool>("Calo2MCTool/CaloFragment2MCTool"); // CAUTION : THE TOOL CANNOT BE PRIVATE !!
  //@TODO: dynamic_cast to IProperty... remove IProperty from ICalo2MCTool...
  m_tool->setProperty ( "Cluster2MCTable" , m_cluster2MCLoc ) ;
  m_tool->setProperty ( "Digit2MCTable"   , m_digit2MCLoc   ) ;
  m_tool->setProperty ( "Cluster2Digits"  , Gaudi::Utils::toString( m_cluster2Digit ) );
  m_tool->setProperty ( "Hypo2Cluster"    , Gaudi::Utils::toString( m_hypo2Cluster  ) );
  m_tool->setProperty ( "Merged2Split"    , Gaudi::Utils::toString( m_merged2Split  ) );
  m_tool->setProperty ( "DigitStatusFilter"  , Gaudi::Utils::toString( m_sFilter    ) ) ;

  m_sum = 0.;
  m_maxMC  = nullptr ;
  m_bestMC = nullptr ;

  //
  if( m_hypo2Cluster) {
    Warning(" ... Hypo->Cluster reference is used in Calo2MC  ", StatusCode::SUCCESS).ignore();
    Warning(" ... CaloCluster  will be re-processed (require full DST)", StatusCode::SUCCESS).ignore();
    Warning(" ... assume an identical reconstruction version  ", StatusCode::SUCCESS).ignore();
  }



  //
  clear();
  return StatusCode::SUCCESS;
}




/*-------------------------- from Part to MC  ------------------------*/
// associate single particle
ICalo2MCTool* Calo2MCTool::from(const LHCb::Particle*   part  ){
  if( part == m_part)return this; // process only if needed

  clear();
  m_depth = 4; // particle level
  if( !part )return this;

  // check the particle is full calorimetric (using caloParticle functionalities)
  LHCb::CaloParticle cPart = LHCb::CaloParticle( (LHCb::Particle*) part );
  if( !cPart.isCalo() ){
    Warning("Cannot associate non-pure calorimetric particle to MC").ignore();
    return this;
  }


  // register final state particles
  m_parts = cPart.caloEndTree();
  m_part = const_cast<LHCb::Particle*>(part);
  if( m_parts.empty())m_parts.push_back( part );



  // particle->protoparticle cascade (mandatory)
  m_depth=3; // protoparticle level
  for(const auto& fs : m_parts ) {
    const LHCb::ProtoParticle* proto = fs->proto();
    if( !proto){
      Warning("ProtoParticle point to NULL (should not)").ignore();
      continue;
    }
    addProto( proto , fs );
  }
 StatusCode sc = process();
  if ( sc.isFailure() )Warning("Processing Calo2MCTool from Particle failed").ignore();
  return this;
}


/*-------------------------- from Proto to MC  ------------------------*/
// Accumulate protos->hypo()s [->clusters]
void Calo2MCTool::addProto(const LHCb::ProtoParticle* proto, const LHCb::Particle* parent ){
  if( !proto)return;

  // register proto if not yet done
  bool ok = std::none_of( m_protos.begin(), m_protos.end(),
                          [&](const LHCb::ProtoParticle* p) { return p == proto; } );
  if(!ok){
    if( UNLIKELY( msgLevel(MSG::DEBUG) ) )
      debug() << "ProtoParticle appears twice in the same decay chain" << endmsg;
    if(counterStat->isQuiet())counter("ProtoParticle appears twice in the same decay chain" ) += 1.;
    return;
  }
  m_protos.push_back( proto ) ;


  // proto->hypo cascade (mandatory)
  m_depth = 2; // hypo level
  const auto& hypos = proto->calo();
  if( hypos.empty() )return;

  // special treatment for Bremstrahlung electrons
  bool charged = ( proto->track() );/* charged part->proto */
  if (!charged) {
    for( const auto& hypo : hypos ) addHypo( hypo );
  } else {
    bool brem = ( parent && ( std::abs(parent->momentum().P() - proto->track()->firstState().p())> 1E-4 ) );/* bremstrahlung corrected particle */
    for( const auto& hypo : hypos ) {
       if ( hypo->hypothesis() == LHCb::CaloHypo::EmCharged  ||
          ( brem && hypo->hypothesis() == LHCb::CaloHypo::Photon ) ) addHypo(hypo);
    }
  }
}

// associate single protoParticle
ICalo2MCTool* Calo2MCTool::from(const LHCb::ProtoParticle*   proto  ){
  if( proto == m_proto)return this; // process only if needed
  clear();
  m_depth = 3; // proto level
  if( !proto  )return this;
  m_proto = const_cast<LHCb::ProtoParticle*>(proto);
  addProto( proto );
  StatusCode sc = process();
  if ( sc.isFailure() )Warning("Processing Calo2MCTool from ProtoParticle failed").ignore();
  return this;
}


/*-------------------------- from Hypo to MC  ------------------------*/
// Accumulate hypos [->clusters]
void Calo2MCTool::addHypo(const LHCb::CaloHypo* hypo){
  if( !hypo  )return;
  // register hypo if not yet done
  bool ok = std::none_of( m_hypos.begin(), m_hypos.end(),
                          [&](const auto& h) { return h==hypo; } );
  if(!ok){
    if( UNLIKELY( msgLevel(MSG::DEBUG) ) )
      debug() << "CaloHypo appears twice in the same decay chain" << endmsg;
    if(counterStat->isQuiet())counter("CaloHypo appears twice in the same decay chain" )+=1;
    return;
  }
  m_hypos.push_back( hypo ) ;
  m_depth=2;

  // ----- hypo->MC association :
  if( ! m_hypo2Cluster ){
    // 2 - get the relevant linker
    std::string loc = ( hypo->parent() && hypo->parent()->registry() ) ? hypo->parent()->registry()->identifier() : "";
    LHCb::Calo2MC::HypoLinkTo linker( evtSvc() , msgSvc() , loc ) ;
    if ( linker.notFound() ){
      Warning( "No Hypo2MC link at '" + loc + "' ",StatusCode::SUCCESS,1 ).ignore() ;
      Warning(" ... try using Hypo->Cluster reference  ", StatusCode::SUCCESS,1).ignore();
      Warning(" ... CaloCluster  will be re-processed (require full DST)", StatusCode::SUCCESS,1).ignore();
      Warning(" ... assume an identical reconstruction version  ", StatusCode::SUCCESS,1).ignore();
      if(counterStat->isQuiet())counter("!! Hypo->Cluster")+=1;
      m_hypo2Cluster=true;
    } // ===> force hypo->cluster cascade
    else{
      //debug() << "Linker found at " << loc << endmsg;
      // - built (particle,weight) map
      for ( const LHCb::MCParticle* particle = linker.first( hypo ) ; particle; particle = linker.next() ) {
        m_mcMap[particle] += linker.weight();
      }
      m_sum += hypo->e();
      if(counterStat->isQuiet())counter("CaloHypo->MC matching") += 1;
    }
  }

  // hypo->cluster->MC cascade if requested/needed
  if( m_hypo2Cluster ){
    m_depth = 1; // cluster level
    const SmartRefVector<LHCb::CaloCluster> clusters = hypo->clusters();
    if( clusters.empty() )return;
    const LHCb::CaloCluster* cluster =
      ( clusters.size() > 1 && hypo->hypothesis() == LHCb::CaloHypo::PhotonFromMergedPi0 ) ?
      *(clusters.begin() + 1 ) : *(clusters.begin()); //@ToDO use cluster::Type (when defined)
    if( !cluster )return ;
    if( counterStat->isQuiet() &&
        clusters.size() !=1 && hypo->hypothesis() != LHCb::CaloHypo::PhotonFromMergedPi0 )counter("ill-defined CaloHypo") += 1;
    addCluster( cluster );
  }
  return;
}
// associate single hypo
ICalo2MCTool* Calo2MCTool::from(const LHCb::CaloHypo*   hypo  ){
  if( hypo == m_hypo)return this; // process only if needed
  clear();
  m_depth = 2; // hypo level
  if( !hypo )return this;
  m_hypo = const_cast<LHCb::CaloHypo*>(hypo);
  // special case for MergedPi0
  if( hypo->hypothesis() == LHCb::CaloHypo::Pi0Merged && m_merged2Split){
    const auto& hyps = hypo->hypos();
    if( hyps.empty() )return this;
    addHypo( hyps.front() ); // splitPhoton1
    addHypo( *std::next(hyps.begin()) ); // splitPhoton2
  } else{
    addHypo( hypo ); // mother CaloHypo
  }
  StatusCode sc = process();
  if ( sc.isFailure() )Warning("Processing Calo2MCTool from CaloHypo failed").ignore();
  return this;
}

/*-------------------------- from Cluster to MC  ------------------------*/
// Accumulate clusters [->digits]
void Calo2MCTool::addCluster(const LHCb::CaloCluster*   cluster  ){
  if( !cluster )return;

  // register cluster (if not yet done)
  bool ok = std::none_of( m_clusters.begin(), m_clusters.end(),
                          [&](const auto& c) { return c == cluster; } );
  if(!ok){
    if( UNLIKELY( msgLevel(MSG::DEBUG) ) )debug() << "Warning : CaloCluster appears twice in the same decay chain" << endmsg;
    if(counterStat->isQuiet())counter( "Warning : CaloCluster appears twice in the same decay chain" )+=1;
    return;
  }
  m_clusters.push_back( cluster ) ;

  m_depth=1;
  // -- Cluster->MC mapping
  if( ! m_cluster2Digit ){
    if(  exist<LHCb::Calo2MC::IClusterTable>(m_cluster2MCLoc ) )m_cluster2MC  = get<LHCb::Calo2MC::IClusterTable>( m_cluster2MCLoc );
    if( !m_cluster2MC ){
      Warning("No Cluster2MC table at " + m_cluster2MCLoc , StatusCode::SUCCESS,1);
      if(counterStat->isQuiet())counter("!! Cluster->Digit")+=1;
      m_cluster2Digit = true;
    } // ===> force clsuter->digit cascade
    else{
      // - built (particle,weight) map
      for( const auto& ir : m_cluster2MC->relations( cluster )) {
          m_mcMap[ir.to()] += ir.weight();
      }
      m_sum += cluster->e();
      if(counterStat->isQuiet())counter("CaloCluster->MC matching") += 1;
    }
  }

  // cluster->digit cascade if requested
  if(m_cluster2Digit){
    m_depth = 0;
    for( const auto& entry : cluster->entries() ) {
      if ( m_sFilter >= 0 && ( entry.status() & m_sFilter) == 0 ){ continue ; }
      const LHCb::CaloDigit* digit = entry.digit();
      if( !digit )continue;
      // check digit is not yet registered
      ok = std::none_of( m_digits.begin(), m_digits.end(),
                         [&](const auto& d) { return d == digit; } );
      if(ok)addDigit(digit);
    }
  }
  return;
}
// associate single cluster
ICalo2MCTool* Calo2MCTool::from(const LHCb::CaloCluster*   cluster  ){
  if( cluster == m_cluster)return this; // process only if needed
  clear();
  m_depth = 1;// cluster level
  if( !cluster )return this;
  m_cluster = (LHCb::CaloCluster*) cluster;
  //
  addCluster( cluster );
  //
  StatusCode sc = process();
  if ( sc.isFailure() )Warning("Processing Calo2MCTool from CaloCluster failed").ignore();
  return this;
}

/*-------------------------- from Digit to MC  ------------------------*/
void Calo2MCTool::addDigit(const LHCb::CaloDigit*   digit  ){
  if( !digit )return;
  m_digits.push_back( digit );
  // -- CaloDigit->MC association
  m_digit2MC  = getIfExists<LHCb::Calo2MC::DigitTable>( m_digit2MCLoc );
  if( !m_digit2MC ){
    Warning(" Digit <-> MC relation table not found at " + m_digit2MCLoc , StatusCode::FAILURE,10).ignore();
    return;
  }
  // - built (particle,weight) map
  for( const auto& ir : m_digit2MC->relations( digit ) ){
    m_mcMap[ir.to()] += ir.weight();
  }
  m_sum += digit->e();
}


// associate single digit
ICalo2MCTool* Calo2MCTool::from(const LHCb::CaloDigit*   digit  ){
  if( digit == m_digit)return this; // process only if needed
  clear();
  m_depth = 0; // digit level
  if( !digit )return this;
  m_digit = const_cast<LHCb::CaloDigit*>(digit);
  m_digits.push_back( digit ) ;
  StatusCode sc = process();
  if ( sc.isFailure() )Warning("Processing Calo2MCTool from CaloDigit failed").ignore();
  return this;
}

/*-------------------------- Generic processing ------------------------*/
StatusCode Calo2MCTool::process(){
  mcDigest();
  verbose() << " Processing Calo2MCTool " << std::endl << descriptor() << endmsg;
  return StatusCode::SUCCESS;
}

void Calo2MCTool::clear(){
  m_mcMap.clear();
  m_treeMap.clear();
  m_digits.clear();
  m_clusters.clear();
  m_hypos.clear();
  m_protos.clear();
  m_parts.clear();
  m_nFrag = 0;
  m_sum = 0.;
  m_maxMC  = nullptr ;
  m_bestMC = nullptr ;
  m_digit   = nullptr;
  m_cluster = nullptr;
  m_hypo    = nullptr;
  m_proto   = nullptr;
  m_part    = nullptr;
  m_digit2MC= nullptr;
  m_cluster2MC=nullptr;
  m_category = 0 ;
  m_depth    = -1 ;
}


void Calo2MCTool::mcDigest(){

  double mcMax = 0.;
  double mcBest = 0.;
  m_maxMC  = nullptr ;
  m_bestMC = nullptr ;

  if( m_sum <= 0 ) return;
  // loop over contributing particle :
  for( const auto& imap : m_mcMap ) {
    const LHCb::MCParticle* mcPart = imap.first;
    double w  = weight( mcPart );
    double q  = quality( mcPart );
    double m = mcPart->momentum().M();
    // the most contributing MCParticle (with smallest mass when several MCPart with same weight)
    if( w >= mcMax ){
      bool ok = true;
      if( m_maxMC  &&  w == mcMax && m > m_maxMC->momentum().M() )ok= false;
      if(ok){
        mcMax = w;
        m_maxMC = mcPart ;
      }
    }
    // the best matching MCParticle
    if( q >= mcBest ){
      mcBest = q;
      m_bestMC = mcPart;
    }

  } // end loop over MCParticles
  // build MC tree
  // 1- get related MC particle (seed) without any descendant listed in the related mcParticles
  std::vector<const LHCb::MCParticle*> seeds;
  for( const auto& imap : m_mcMap ) {
    const LHCb::MCParticle* mcPart = imap.first;
    int _pID =  mcPart->particleID().abspid() ;

    auto hasProd = [&]( const LHCb::MCParticle* mcp ) {
      const auto& vertices = mcp->endVertices();
      return std::any_of( vertices.begin(), vertices.end(),
                                [&](const auto& v) {
        return std::any_of( v->products().begin(), v->products().end(),
                            [&](const auto& p)
                            { return m_mcMap.find(p)!=m_mcMap.end(); } );
      });
    };

    // include electron with brems and converted photons as seeds
    if ( _pID == 11 || _pID == 22 || !hasProd(mcPart) ) seeds.push_back( mcPart );
  }

  // 2- build the seed upstream tree
  for( const LHCb::MCParticle* seed : seeds ) {
    std::vector<const LHCb::MCParticle*> tree;
    std::string sTree ;
    mcTree( seed , tree, sTree);
    std::stringstream ss;
    ss << format(" %6.1f %% from : ", weight( seed ) *100. )
       << sTree
       << " ( " << format(" %6.1f %% of the MC particle energy contributing",
                          (weight( seed ) == 0)?0: quality(seed) /weight( seed )*100.) << " )";
    m_treeMap[ ss.str() ] = std::move(tree);
  }
}

void Calo2MCTool::mcTree(const LHCb::MCParticle* part, std::vector<const LHCb::MCParticle*>& tree , std::string& sTree){
  if( !part ) return;
  tree.push_back( part );
  const LHCb::ParticleProperty* prop = m_ppsvc->find( part->particleID() );
  sTree = ( prop ? prop->name() : "??" ) + sTree;
  if( part->mother() ){
    sTree = " -> " + sTree;
    mcTree( part->mother() , tree, sTree );
  }
}

double Calo2MCTool::weight(const LHCb::MCParticle* part) const {
  return ( part && m_sum>0 ) ? m_mcMap[part]/m_sum : 0.;
}

double Calo2MCTool::quality(const LHCb::MCParticle* part) const {
  return ( part && part->momentum().E() != 0) ? weight(part) * m_mcMap[part]/part->momentum().E() : 0.;
}


std::string Calo2MCTool::descriptor() const {
  std::stringstream ss;
  ss  << "\n     ---------- Calo MC contribution " ;
  if( m_part  ) ss << "to particle (pid = " << m_part->particleID().pid() << ")\n" ;
  if( m_parts.size() > 1 ) ss << " -> to " << m_parts.size() << " particle(s) -------- \n" ;
  if( !m_protos.empty() ) ss << "to " << m_protos.size() << " protoParticle(s) -------- \n" ;
  if( !m_hypos.empty() ) ss << "to " << m_hypos.size() << " hypo(s) -------- \n" ;
  if( !m_digits.empty() ) ss << "to " << m_digits.size() << " digit(s) -------- \n" ;
  if( !m_clusters.empty() ) ss << "to " << m_clusters.size() << " cluster(s) ------- \n" ;
  ss  << "     ---- Total calo energy deposit : " << m_sum << " MeV \n" ;
  for( const auto& im : m_treeMap ) ss  << "        -- " << im.first << '\n';

  if( bestMC() ){
    const LHCb::ParticleProperty* prop = m_ppsvc->find( bestMC()->particleID() );
    std::string p = ( prop ? prop->name() : "??" );
    ss << "      --> Best matching MCParticle : [" << p << "] ==  (Quality/Weight : "
       << quality( bestMC() ) <<" / " << weight( bestMC() )<< ")\n" ;
  }

  if( maxMC() ){
    const LHCb::ParticleProperty* prop = m_ppsvc->find( maxMC()->particleID() );
    std::string p = ( prop ? prop->name() : "??" );
    ss << "      --> Maximum weight MCParticle : [" << p << "] == (Quality/Weight : "
       << quality( maxMC() ) <<" / " << weight( maxMC() )<<")\n" ;
  }

  ss << "      -------------------------------- ";
  return ss.str();
}


const LHCb::MCParticle* Calo2MCTool::findMCOrBest(LHCb::ParticleID id, double threshold ) const {
  const LHCb::MCParticle* found = findMC(id,threshold);
  return found ? found : bestMC();
}
const LHCb::MCParticle* Calo2MCTool::findMCOrBest(std::string name, double threshold ) const {
  const LHCb::MCParticle* found = findMC(name,threshold);
  return found ? found : bestMC();
}
const LHCb::MCParticle* Calo2MCTool::findMC(std::string name, double threshold ) const{
  const LHCb::ParticleProperty* prop = m_ppsvc->find( name );
  return prop ? findMC( prop->particleID() , threshold) : nullptr;
}
const LHCb::MCParticle* Calo2MCTool::findMC(LHCb::ParticleID id, double threshold ) const {
  double t = threshold;
  const LHCb::MCParticle* best  = nullptr;
  for( const auto& imap : m_mcMap ) {
    const LHCb::MCParticle* mcPart = imap.first;
    if( mcPart->particleID().abspid() != id.abspid() )continue;
    double q = quality( mcPart );
    if( q < t )continue;
    t = q;
    best =  mcPart;
  }
  return best;
}
const LHCb::MCParticle* Calo2MCTool::bestMC() const {
  return m_bestMC;
}
const LHCb::MCParticle* Calo2MCTool::maxMC() const {
  return m_maxMC;
}


ICalo2MCTool* Calo2MCTool::fragment(unsigned int i){
  m_nFrag = 0;
  /* CaloHypo level */
  if( m_depth == 2 ){
    m_nFrag = m_hypos.size();
    if( i >= m_nFrag)return nullptr;
    if(m_nFrag == 1)return this;
    return m_tool->from( m_hypos[i] );
  }
  /* CaloCluster level */
  if( m_depth == 1 ){
    m_nFrag = m_clusters.size();
    if( i >= m_nFrag)return nullptr;
    if(m_nFrag == 1)return this;
    return m_tool->from( m_clusters[i] );
  }
  /* CaloDigit level */
  if( m_depth == 0 ){
    m_nFrag = m_digits.size();
    if( i >= m_nFrag)return nullptr;
    if(m_nFrag == 1)return this;
    return m_tool->from( m_digits[i] );
  }
  return nullptr;
}
