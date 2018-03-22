#include <iostream>
#include "tools.h"
#include <string>
using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}
//Function to handle the debug message
void Tools::Debug(string msg1){
  bool show_debug = false;
  if(show_debug)
  {
    cout<<msg1<<endl;
  }
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  //DEFINE THE RMSE VECTOR
  VectorXd rmse(4);
  rmse<<0,0,0,0;
  
  /*if((estimations.size() != ground_truth.size()) || estimations.size() == 0){
    printf("Error in the estimation and ground truth size.;");
    //RETURN 0 
    return rmse;
  }*/
  
  for(int i=0;i<estimations.size();i++){
    VectorXd temp = estimations[i] - ground_truth[i];
    temp = temp.array() * temp.array();
    rmse += temp;
  }
  //CALCULATING MEAN BY DIVING THE SQUARED DIFFERENCE BY NUMBER OF ELEMNETS
  rmse = rmse / estimations.size();
  //CALCULATING THE SQUARE ROOT OF THE ELEMENT
  rmse = rmse.array().sqrt();
  
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  float px, py, vx, vy;
  px = x_state[0];
  py = x_state[1];
  vx = x_state[2];
  vy = x_state[3];
  
  //Calculating repeated terms
  float t1 = px*px + py*py;
  float t2 = sqrt(t1);
  float t3 = t1 * t2;
  //Handling the case where the px or py might be zero.
 /* if(t1 == 0 || t1<0.000001)
    t1 = 0.000001;
  if(t2 == 0 || t2<0.000001)
    t2 = 0.000001;
  if(t3 == 0 || t3<0.000001)
    t3 = 0.000001;*/
  //Calculate Jacobian Matrix
  Hj<<px/t2, py/t2, 0, 0,
  -py/(t1), px/t1, 0, 0,
  py*(vx*py - vy*px)/t3, px*(vy*px - vx*py)/t3, px/t2, py/t2;
  
  return Hj;
}
