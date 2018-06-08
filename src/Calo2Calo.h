#ifndef CALO2CALO_H
#define CALO2CALO_H 1

// Include files
// from Gaudi
#include "GaudiAlg/GaudiTool.h"
#include "CaloInterfaces/ICalo2Calo.h"            // Interface
#include "CaloDet/DeCalorimeter.h"
#include  "CaloInterfaces/ICaloGetterTool.h"
/** @class Calo2Calo Calo2Calo.h
 *
 *
 *  @author Olivier Deschamps
 *  @date   2007-05-29
 */

class Calo2Calo : public extends<GaudiTool, ICalo2Calo> {
public:
  /// Standard constructor
  Calo2Calo( const std::string& type,
             const std::string& name,
             const IInterface* parent);

  StatusCode initialize() override;
  // setting
  void setCalos ( const std::string& fromCalo ,
                  const std::string& toCalo   ) override;
  // CaloCellIDs
  const std::vector<LHCb::CaloCellID>&
  cellIDs ( const LHCb::CaloCluster& fromCluster ,
            const std::string&       toCalo      ) override;
  const std::vector<LHCb::CaloCellID>&
  cellIDs ( const LHCb::CaloCellID&  fromId     ,
            const std::string&       toCalo, bool init=true) override;
  const std::vector<LHCb::CaloCellID>& cellIDs() override {return m_cells;};
  // Digits
  const std::vector<LHCb::CaloDigit*>& digits
  ( const LHCb::CaloCellID& fromId     ,
    const std::string&      toCalo     ) override;
  const std::vector<LHCb::CaloDigit*>& digits
  ( const LHCb::CaloCluster& fromCluster ,
    const std::string&       toCalo    ) override;
  const std::vector<LHCb::CaloDigit*>& digits() override {return m_digits;};
  // Energy
  double energy ( const LHCb::CaloCellID&  fromId ,
                  const std::string&       toCalo ) override;
  double energy ( const LHCb::CaloCluster& fromCluster ,
                  const std::string&       toCalo ) override;
  double energy () override {return m_energy;};
  // multiplicity
  int multiplicity( const LHCb::CaloCellID&  fromId ,
                    const std::string&       toCalo ) override;
  int multiplicity( const LHCb::CaloCluster& fromCluster ,
                    const std::string&       toCalo ) override;
  int multiplicity() override {return m_count;};
  // Additional
  bool isLocalMax ( const LHCb::CaloDigit& digit) override;

protected:
  //
  void reset();
  const std::vector<LHCb::CaloCellID>& addCell
  ( const LHCb::CaloCellID& id, const std::string& toCalo);
  // Calo Maps
  std::map<std::string,DeCalorimeter*> m_det;
  std::map<std::string,std::string> m_loc;
  std::map<std::string,double>m_refSize;
  std::map<std::string,Gaudi::Plane3D> m_plane;
  //
  Gaudi::Property<std::string> m_fromCalo {this, "FromCalo", "??"};
  Gaudi::Property<std::string> m_toCalo {this, "ToCalo", "??"};
  std::vector<LHCb::CaloCellID> m_cells;
  std::vector<LHCb::CaloDigit*>  m_digits;
  double m_energy = 0.;
  int m_count = 0;
  DeCalorimeter* m_fromDet = nullptr;
  DeCalorimeter* m_toDet = nullptr;
  double m_fromSize = 0.;
  double m_toSize = 0.;
  std::string m_toLoc;
  Gaudi::Plane3D m_toPlane;
  LHCb::CaloDigits* m_digs = nullptr;
  ICaloGetterTool* m_getter = nullptr;
  bool m_ok = true;
private:
  Gaudi::Property<bool> m_geo {this, "IdealGeometry", true};
  Gaudi::Property<std::string> m_getterName {this, "GetterName", "CaloGetter"};

};
#endif // CALO2CALO_H
