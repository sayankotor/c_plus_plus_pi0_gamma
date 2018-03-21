#ifndef XGBCLASSIFIER_H
#define XGBCLASSIFIER_H

#include <string>
#include <vector>
#include "include/xgboost/c_api.h"
//#include "/afs/cern.ch/user/v/vchekali/public/xgboost/include/xgboost/c_api.h"

class IClassifier
{
public:

    /**    * @brief      Main classification method    *
           * Takes a vector of values of features and returns the corresponding MVA
           * output.
           *    * @param[in]  featureValues  A vector of feature values
           *    * @return     MVA classifier value    */
    virtual double getClassifierValues(const std::vector<double>& featureValues) = 0;
};

class XGBClassifierPhPi0 : public IClassifier {
  public:
    XGBClassifierPhPi0();
    XGBClassifierPhPi0(const std::string& path);
    double getClassifierValues(const std::vector<double>& featureValues);
    void setPath(const std::string& path);

  private:
    std::string xgb_path ="def_path";
    std::vector<float> m_predictionsCache = {0, 0};
    DMatrixHandle m_cache_matrix,
                  m_feature_matrix;
    BoosterHandle m_booster;

};


#endif // XGBCLASSIFIER_H
