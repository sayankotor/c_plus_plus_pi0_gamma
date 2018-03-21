#include "XGBClassifierPhPi0.h"
#include <iostream>
#include <exception>

XGBClassifierPhPi0::XGBClassifierPhPi0()
{
  XGDMatrixCreateFromMat(
      &m_predictionsCache.at(0),
      1,
      2,
      0.5,
      &m_cache_matrix);
}

XGBClassifierPhPi0::XGBClassifierPhPi0(const std::string& path)
{
  xgb_path = path;
  //std::cout<<"PATH str: "<< xgb_path << std::endl;
  XGDMatrixCreateFromMat(
      &m_predictionsCache.at(0),
      1,
      2,
      0.5,
      &m_cache_matrix);
  XGBoosterCreate(&m_cache_matrix, 1, &m_booster);
  XGBoosterLoadModel(m_booster, xgb_path.c_str());
}

void XGBClassifierPhPi0::setPath(const std::string& path)
{
  XGBoosterCreate(&m_cache_matrix, 1, &m_booster);
  try{
    XGBoosterLoadModel(m_booster, path.c_str());
  }
  catch(const std::exception& e){
    std::cout << "Standard exception: " << e.what() << std::endl;
  }
}

double XGBClassifierPhPi0::getClassifierValues(const std::vector<double>& featureValues)
{
  // currently XGBoost only supports float
  std::vector<float> features_f(featureValues.begin(), featureValues.end());

  // fill the feature vector into a XGBoost DMatrix
  XGDMatrixCreateFromMat(&features_f.at(0),
                         1,
                         features_f.size(),
                         0,
                         &m_feature_matrix);

  // XGBoost returns the predictions into a arrayish object and return
  // its size
  unsigned long predictions_length;
  const float *predictions;
  try{
  XGBoosterPredict(m_booster,
                   m_feature_matrix,
                   0,
                   0,
                   &predictions_length,
                   &predictions);
  }
  //catch(const std::exception& e){
  catch(...){
    std::cout << "Some exception: " << std::endl;
  }
  std::vector<double> vec_predictions;
  for (int i =0; i<predictions_length; i++){
      vec_predictions.push_back(predictions[i]);
  }
    return vec_predictions[0];
}