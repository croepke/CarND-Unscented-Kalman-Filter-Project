#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */

   VectorXd rmse = VectorXd(4);
   rmse << 0,0,0,0;

   if (estimations.size() != ground_truth.size()||estimations.size()==0) {
     return rmse;
   }

   for (unsigned int i=0; i<estimations.size(); ++i) {
     VectorXd residuals;
     residuals = estimations[i]-ground_truth[i];
     residuals = residuals.array()*residuals.array();
     rmse+=residuals;
   }

   rmse /= estimations.size();

   return rmse.array().sqrt();

}
