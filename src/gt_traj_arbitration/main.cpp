#include "ros/ros.h"
#include <gt_traj_arbitration/gt_traj_arbitration.h>



bool solveRiccati(const Eigen::MatrixXd &A,
                  const Eigen::MatrixXd &B,
                  const Eigen::MatrixXd &Q,
                  const Eigen::MatrixXd &R, Eigen::MatrixXd &P)
{

  const uint dim_x = A.rows();
  const uint dim_u = B.cols();

  Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
  Ham << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();

  ROS_INFO_STREAM("ham\n"<<Ham);
  
  Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);

  Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
  int j = 0;
  for (int i = 0; i < 2 * dim_x; ++i) {
    if (Eigs.eigenvalues()[i].real() < 0.) {
      eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
      ++j;
    }
    else
      ROS_ERROR_STREAM("eigenvalue with positive real value ! "<<Eigs.eigenvalues()[i].real() );
  }

  Eigen::MatrixXcd Vs_1, Vs_2;
  Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
  Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
  P = (Vs_2 * Vs_1.inverse()).real();

  return true;
}

int main(int argc, char **argv)
{
  ros::init(argc, argv, "main");
  ros::NodeHandle n;
  
  double m,c,k;
  m = 10;
  c = 25;
  k = 0;
  
  Eigen::MatrixXd A;
  Eigen::MatrixXd B;
  Eigen::MatrixXd C;
  
  Eigen::MatrixXd Qh;
  Eigen::MatrixXd Qr;
  Eigen::MatrixXd Rh;
  Eigen::MatrixXd Rr;
  
  Eigen::MatrixXd Q_alpha;
  Eigen::MatrixXd R_alpha;
  
  A.resize(3,3);
  B.resize(3,2);
  C.resize(3,3);
  
  Qh.resize(3,3);
  Qr.resize(3,3);
  Rh.resize(1,1);
  Rr.resize(1,1);
  
  Q_alpha.resize(3,3);
  R_alpha.resize(2,2);
  
  A .setZero();
  B .setZero();
  C .setZero();
  Qh.setZero();
  Qr.setZero();
  Rh.setZero();
  Rr.setZero();
  Q_alpha.setZero();
  R_alpha.setZero();
  
  
  A <<  0,    0,    1,
        0,    0,    1,
        0,  -k/m, -c/m;
  
  B <<   0,   0,
         0,   0,
        1/m, 1/m;
  
  C = Eigen::MatrixXd::Identity(3,3);
  
  Qh.diagonal() << 0, 1, 0.1;
  Qr.diagonal() << 1, 0, 0.1;
  Rh << 0.0005;
  Rr << 0.0005;
  
  double alpha = 0.5;
  
  Q_alpha = alpha*Qh+(1-alpha)*Qr;
  
  R_alpha.diagonal() << alpha*Rh , (1-alpha)*Rr;
  
  ROS_INFO_STREAM("A\n"<<A);
  ROS_INFO_STREAM("B\n"<<B);
  ROS_INFO_STREAM("C\n"<<C);
  
  ROS_INFO_STREAM("Qh\n"<<Qh);
  ROS_INFO_STREAM("Qr\n"<<Qr);
  ROS_INFO_STREAM("Rh\n"<<Rh);
  ROS_INFO_STREAM("Rr\n"<<Rr);
  
  ROS_INFO_STREAM("Q_alpha\n"<<Q_alpha);
  ROS_INFO_STREAM("R_alpha\n"<<R_alpha);
  
  std::vector<double> mult = {1,10,100,1000,10000,100000,1000000,10000000,100000000,1000000000};
  
  for (int i=0;i<mult.size();i++)
  {
    Eigen::MatrixXd P ;
    solveRiccati(A,B,Q_alpha*mult[i],R_alpha*mult[i],P);
    
    ROS_INFO_STREAM("mult: "<<i<<", P\n"<<P);
    
    
    Eigen::MatrixXd K = -R_alpha.inverse()*B.transpose()*P; 
    
//     ROS_INFO_STREAM("K\n"<<K);
  }
  return 0;
}
