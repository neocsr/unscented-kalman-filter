#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    throw std::invalid_argument("Estimations or ground truth data size error");
  }

  // Calculate squared residuals
  for (int i = 0; i < estimations.size(); i++) {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // Calculate the mean and squared root
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}
