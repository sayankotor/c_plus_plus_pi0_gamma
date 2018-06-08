#ifndef L0CALO2CALOTOOL_H
#define L0CALO2CALOTOOL_H 1
// ============================================================================
#include <string>
#include <iostream>
#include <string>
#include <iostream>
#include "GaudiAlg/GaudiTool.h"
#include "CaloKernel/CaloVector.h"
#include "CaloDet/DeCalorimeter.h"

#include "CaloInterfaces/IL0Calo2Calo.h"
#include "CaloInterfaces/ICaloClusterization.h"

#include "CaloDAQ/ICaloDataProvider.h"
#include "CaloDAQ/ICaloEnergyFromRaw.h"

// forward declarations
class DeCalorimeter;

/** @class L0Calo2CaloTool L0Calo2CaloTool.h
 * Tool to get a list of CaloClusters (owned by TES) in the vicinity of the input L0CaloCandidate(s),
 * if necessary invoking decoding and clusterization.
 *
 * Remarks:
 * - the returned clusters are owned by TES and should not be deleted by the user;
 * - in case of multiple calls of the tool on intersecting zones of Ecal
 *   - there might appear duplicated clusters on TES,
 *   - unless the input/output   std::vector<LHCb::CaloCluster*> clusters
 *     is clear()'ed by the user before calling L0Calo2CaloTool::clusterize(),
 *     this vector also might contain duplicates
 *
 *  @author Dmitry Golubkov
 *  @date   2009-07-29
 */
class L0Calo2CaloTool : public extends<GaudiTool,IL0Calo2Calo>
{

public:
  /// Standard constructor
  L0Calo2CaloTool( const std::string& type,
                   const std::string& name,
                   const IInterface* parent);
  // ==========================================================================
  /** obtain CaloClusters corresponding to L0CaloCandidates
   *
   * Get a list of CaloClusters in the vicinity of the L0CaloCandidate,
   * invoke decoding and clusterization if necessary.
   * The output clusters are stored in TES location CaloDigitLocation, the
   * created CaloDigits are stored on TES in CaloDigitLocation.
   *
   * @param clusters  (OUTPUT) vector of pointers of Calo clusters (output)
   * @param candidate (INPUT) pointer to L0CaloCandidate
   * @param level     (INPUT) number of neigbour levels around the candidate->id() cell
   *                  for the CaloClusterization tool
   * @return FAILURE when anything goes wrong, otherwise the return StatusCode of the ICaloClusterization tool used
   */
  StatusCode clusterize
  ( std::vector<LHCb::CaloCluster*>&      clusters ,
    const LHCb::L0CaloCandidate*          candidate,
    const unsigned int                    level     ) const override;
  // ==========================================================================
  /** obtain CaloClusters corresponding to L0CaloCandidate
   *
   * Get a list of CaloClusters in the vicinity of an L0CaloCandidate,
   * if necessary invoke decoding and clusterization.
   *
   * @param clusters (OUTPUT) vector of pointers of Calo clusters
   * @param candidate (INPUT) pointer to L0CaloCandidate
   */
  StatusCode clusterize
  ( std::vector<LHCb::CaloCluster*>&      clusters,
    const LHCb::L0CaloCandidate*          candidate ) const override;
  // ==========================================================================
  /** obtain CaloClusters corresponding to L0CaloCandidates
   *
   * Get a list of CaloClusters in the vicinity of an ObjectVector of L0CaloCandidates,
   * invoke decoding and clusterization if necessary.
   * @param clusters (OUTPUT) vector of pointers of Calo clusters
   * @param candidates (INPUT) pointer to L0CaloCandidates
   * @param level (INPUT) number of neigbour levels around cell for the CaloClusterization tool
   */
  StatusCode clusterize
  ( std::vector<LHCb::CaloCluster*>&      clusters,
    const LHCb::L0CaloCandidates*         candidates,
    const unsigned int                    level       ) const override
  {
    (void) clusters, (void) candidates, (void) level; // avoid compiler warning
    fatal() << "L0Calo2CaloTool::clusterize(..., const LHCb::L0CaloCandidates*, ...) NOT IMPLEMENTED" << endmsg;
    return StatusCode::FAILURE;
  }
  // ==========================================================================
  /** obtain CaloClusters corresponding to L0CaloCandidates
   *
   * Get a list of CaloClusters in the vicinity of an ObjectVector of L0CaloCandidates,
   * invoke decoding and clusterization if necessary.
   * @param clusters (OUTPUT) vector of pointers of Calo clusters
   * @param candidates (INPUT) pointer to L0CaloCandidates
   */
  StatusCode clusterize
  ( std::vector<LHCb::CaloCluster*>&      clusters,
    const LHCb::L0CaloCandidates*         candidates  ) const override { return clusterize(clusters, candidates, m_neighbourLevel); }
  // ==========================================================================
  /** obtain CaloClusters around a CaloCellID
   *
   * Get a list of CaloClusters in the vicinity of the CaloCellID,
   * if necessary invoke decoding and clusterization.
   *
   * @param clusters (OUTPUT) vector of pointers of Calo clusters
   * @param cellID   (INPUT)  pointer to CaloCellID
   * @param level    (INPUT)  number of neigbour levels around the cell for the ICaloClusterization tool
   */
  StatusCode clusterize
  ( std::vector<LHCb::CaloCluster*>&      clusters ,
    const LHCb::CaloCellID&               cellID,
    const unsigned int                    level     ) const override;
  // ==========================================================================
  /** obtain CaloClusters around a CaloCellID
   *
   * Get a list of CaloClusters in the vicinity of the CaloCellID,
   * if necessary invoke decoding and clusterization.
   *
   * @param clusters (OUTPUT) vector of pointers of Calo clusters
   * @param cellID   (INPUT)  pointer to CaloCellID
   */
  StatusCode clusterize
  ( std::vector<LHCb::CaloCluster*>&      clusters,
    const LHCb::CaloCellID&               cellID   ) const override;
  // ==========================================================================
  /** Interface to ICaloClusterizationTool::iterations() */
  virtual unsigned int iterations() const { return m_clusterizationTool ? m_clusterizationTool->iterations() : 0; }
  // ==========================================================================
  StatusCode initialize() override;
  // ==========================================================================
protected:
  // ==========================================================================
  /** obtain CaloClusters corresponding to the digits using cellID as a seed
   *
   * Get a list of CaloClusters in the vicinity of cellID, see makeDigits()
   * for explanation of the region used for cluster search.
   *
   * @param clusters (OUTPUT) vector of pointers to Calo clusters
   * @param digits (INPUT) pointer to the CaloDigits
   * @param cellID (INPUT) of the L0CaloCandidate seed cell
   * @param level (INPUT) the radius of the area to be clusterized
   */
  StatusCode makeClusters
  ( std::vector<LHCb::CaloCluster*>&      clusters,
    const LHCb::CaloDigits*               digits,
    const LHCb::CaloCellID&               cellID,
    const unsigned int                    level   ) const ;
  // ==========================================================================
  /** put the vector CaloClutsters on TES
   *
   * @param clusters (INPUT/OUTPUT) vector of pointers to Calo clusters
   *
   * NB: in case of multiple calls of the L0Calo2CaloTool on intersecting
   * zones of Ecal, there might appear duplicates!
   */
  StatusCode putClustersOnTES
  ( std::vector<LHCb::CaloCluster*>&      clusters ) const ;
  // ==========================================================================
  /** prepare CaloDigits around the given cellID and store them on TES
   *
   * @param cellID of the seed cell
   * @param level defines the size of the digitized area (see below)
   * @return container of CaloDigits corresponding to the TES
   *
   *  The decoding is done for the Tell1's which correspond to the cells
   *  in the region centered around the seed cellID and having the area of the
   *  the size ~ 2*(1+Level)x2*(1+Level) as steered by the level parameter.
   *  The same region of cells is used by the CaloClusterization tool
   *  invoked by makeClusters(). The area is constructed by recurrent addition of
   *  DeCalorimeter::neighborCell(...) to the the area starting from the L0CaloCandidate.
   *  Typical value which should be OK for the NeighbourLevel is ~2.
   *
   *  Below is an illustration for NeighbourLevel = 3:
   * -+-+-+-+-+-+-+-+-+-+-+-+-
   * .|.|.|.|.|.|.|.|.|.|.|.|.
   * -+-+-+-+-+-+-+-+-+-+-+-+-
   * .|.|.|*|*|*|*|*|*|*|*|.|.<--
   * -+-+-+-+-+-+-+-+-+-+-+-+- |
   * .|.|.|*|*|*|*|*|*|*|*|.|. d
   * -+-+-+-+-+-+-+-+-+-+-+-+- |
   * .|.|.|*|*|*|*|*|*|*|*|.|.<--   d     = level (defaults to the NighbourLevel property)
   * -+-+-+-+-+-+-+-+-+-+-+-+-      .     = Ecal cells
   * .|.|.|*|*|*|#|#|*|*|*|.|.      *     = cells included in the digitization
   * -+-+-+-+-+-+-+-+-+-+-+-+-      {@,#} = 2x2 cell region of the L0CaloCandidate
   * .|.|.|*|*|*|@|#|*|*|*|.|.      @     = seed cell
   * -+-+-+-+-+-+-+-+-+-+-+-+-
   * .|.|.|*|*|*|*|*|*|*|*|.|. <--
   * -+-+-+-+-+-+-+-+-+-+-+-+-  |
   * .|.|.|*|*|*|*|*|*|*|*|.|.  d
   * -+-+-+-+-+-+-+-+-+-+-+-+-  |
   * .|.|.|*|*|*|*|*|*|*|*|.|. <--
   * -+-+-+-+-+-+-+-+-+-+-+-+-
   * .|.|.|.|.|.|.|.|.|.|.|.|.
   *      ^     ^   ^     ^
   *      |<-d->|   |<-d->|
   */
  const LHCb::CaloDigits* makeDigits( const LHCb::CaloCellID &cellID, const unsigned int level ) const ;
  // ==========================================================================
  /** check if cell ID is usable
   *
   * A cell is accepted if it corresponds to Ecal
   * and if DeCalorimeter::valid( cellID ) is true
   * @param cellID of the cell in question
   * @return true if the cell is accepted / false otherwise
   */
  bool isUsable( const LHCb::CaloCellID &cellID ) const ;
  // ==========================================================================
  /** check if the L0CaloCandidate is usable
   *
   * An L0CaloCandidate is accepted if its type() is Electron or Photon (see L0DUBase::CaloType )
   * @param  candidate pointer to the tested L0CaloCandidate
   * @return true if the candidate is accepted / false otherwise
   */
  bool isUsable( const LHCb::L0CaloCandidate *candidate ) const ;
  // ==========================================================================

  /** collect the Tell1 numbers corresponding to the given L0CaloCandidate seed cellID
   *
   * @param tell1s (OUTPUT) set of tell1 IDs
   * @param cellID (INPUT) the seed cell of a 2x2 cell L0CaloCandidate
   * @param level (INPUT) the neighbour level which defines the size of the area, see makeDigits(...)
   */
  void collectTell1s
  ( std::set<int>&          tell1s,
    const LHCb::CaloCellID& cellID,
    const unsigned int      level  ) const ;
  // ==========================================================================
private:
  Gaudi::Property<std::string> m_clusterizationToolName {this, "CaloClusterizationTool", "CaloClusterizationTool"};
  Gaudi::Property<std::string> m_dataProviderToolName   {this, "CaloDataProviderTool"  , "CaloDataProvider"};

  ICaloClusterization* m_clusterizationTool = nullptr;
  ICaloDataProvider*   m_dataProviderTool = nullptr;
  DeCalorimeter*       m_calo = nullptr;

  unsigned int         m_ecalCaloNum;

  Gaudi::Property<std::string>  m_digitLocation  {this, "CaloDigitLocation"  , LHCb::CaloDigitLocation::Hlt1Ecal};
  Gaudi::Property<std::string>  m_clusterLocation{this, "CaloClusterLocation", LHCb::CaloClusterLocation::EcalHlt1};
  
  Gaudi::Property<unsigned int> m_neighbourLevel
    {this, "NeighbourLevel", 2, "Level parameter for the CaloClusterizationTool, search clusters in (1+2*Level)x(1+2*Level) region around the seed cell"};

  Gaudi::Property<bool> m_sort
    {this, "Sort", false, "sort the clusters due to energy"};

  Gaudi::Property<bool> m_sortET
    {this, "SortET", false, "if Sort: sort the clusters due to transverse energy"};

  Gaudi::Property<bool> m_decodeFullEcal
    {this, "DecodeFullEcal", false, "false = decode only the Tell1s around the L0CaloCandidate cellID"};

  Gaudi::Property<bool> m_clusOnTES
    {this, "ClusterOnTES", false};

  mutable std::set<int>  m_decodedSources;
};
#endif // L0CALO2CALOTOOL_H
