#pragma once

#include <cmath>
#include <Eigen/Core>
#include <ros/time.h>
#include <geometry_msgs/WrenchStamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <sensor_msgs/JointState.h>

#include <state_space_filters/filtered_values.h>
#include <cnr_controller_interface/cnr_joint_command_controller_interface.h>
#include <cnr_hardware_interface/posveleff_command_interface.h>
#include <cnr_hardware_interface/veleff_command_interface.h>
#include <control_msgs/GripperCommandAction.h>
#include <actionlib/client/simple_action_client.h>
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/neldermeadsolver.h>

#include <std_msgs_stamped/Float32Stamped.h>

namespace ect = eigen_control_toolbox;

namespace cnr
{
namespace control
{


/**
 * @brief The GTTrajArbitration class
 */
class GTTrajArbitration: public cnr::control::JointCommandController<
        hardware_interface::JointHandle, hardware_interface::VelocityJointInterface>
{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

  GTTrajArbitration();
  ~GTTrajArbitration();
  bool doInit();
  bool doUpdate  (const ros::Time& time, const ros::Duration& period);
  bool doStarting(const ros::Time& time);
  bool doStopping(const ros::Time& time);
  
  enum Control
  {
    CGT   = 0
    ,LQR  = 1
    ,NCGT = 2
    ,ICGT = 3
  };
  const std::map<std::string, GTTrajArbitration::Control> control_map_ = 
  {
    {"cgt"             ,GTTrajArbitration::Control::CGT}
    ,{"cooperative"    ,GTTrajArbitration::Control::CGT}
    ,{"lqr"            ,GTTrajArbitration::Control::LQR}
    ,{"ncgt"           ,GTTrajArbitration::Control::NCGT}
    ,{"non_cooperative",GTTrajArbitration::Control::NCGT}
    ,{"inverse_cgt"    ,GTTrajArbitration::Control::ICGT}
  }; 
  const std::map<int, std::string> control_map_reverse_ = 
  {
    { 0 , "cgt"            }
    ,{0 , "cooperative"    }
    ,{1 , "lqr"            }
    ,{2 , "ncgt"           }
    ,{2 , "non_cooperative"}
    ,{3 , "inverse_cgt"    }
  }; 

protected:

  std::mutex m_mtx;
  std::mutex gains_mtx_;
  
  double alpha_;
  double alpha_max_;
  double alpha_min_;
  double alpha_switch_;
  
  std::thread* gain_thread_;
  
  int ctr_switch_;
  bool new_gain_available_;
  
  
  ect::FilteredVectorXd wrench_fitler_;

  Eigen::VectorXd dq_sp_;
  Eigen::VectorXd q_sp_;
  Eigen::VectorXd ddq_;
  Eigen::VectorXd dq_;
  Eigen::VectorXd q_;

  geometry_msgs::PoseStamped robot_pose_sp_;
  geometry_msgs::PoseStamped human_pose_sp_;

  Eigen::Affine3d T_robot_base_targetpose_;
  Eigen::Affine3d T_human_base_targetpose_;
  
  Eigen::MatrixXd A_;
  Eigen::MatrixXd B_;
  Eigen::MatrixXd B_single_;
  
  Eigen::MatrixXd Qhh_;
  Eigen::MatrixXd Qhr_;
  Eigen::MatrixXd Qrh_;
  Eigen::MatrixXd Qrr_;
  
  Eigen::MatrixXd Qh_;
  Eigen::MatrixXd Qr_;
  
  Eigen::MatrixXd Rh_;
  Eigen::MatrixXd Rr_;
  
  
  Eigen::MatrixXd Q_gt_;
  Eigen::MatrixXd R_gt_;

  Eigen::MatrixXd P_;
  Eigen::MatrixXd CGT_gain_;
  
  Eigen::MatrixXd Kh_lqr_;
  Eigen::MatrixXd Kr_lqr_;
  Eigen::MatrixXd Kh_nc_ ;
  Eigen::MatrixXd Kr_nc_;
  
  Eigen::MatrixXd K_cgt;
  Eigen::MatrixXd Kh_lqr;
  Eigen::MatrixXd Kr_lqr;
  Eigen::MatrixXd Kh_nc;
  Eigen::MatrixXd Kr_nc;
  
  std::string control_type_;
  std::string base_control_type_;
  
  
  bool use_cartesian_reference_;
  bool robot_active_;

  bool first_cycle_;
  bool new_sp_available_;
  
  int n_dofs_;
  
  bool w_b_init_;
  bool use_filtered_wrench_;
  
  Eigen::Vector6d w_b_filt_;
  Eigen::Vector6d w_b_;
  Eigen::Vector6d w_b_0_;
  Eigen::Vector6d wrench_deadband_;

  rosdyn::ChainPtr chain_bs_;
  rosdyn::ChainPtr chain_bt_;

  size_t filtered_wrench_base_pub_;
  size_t wrench_base_pub_;
  size_t wrench_tool_pub_;
  size_t robot_wrench_pub_;
  size_t nominal_h_wrench_pub_;
  size_t current_pose_pub_;
  size_t current_vel_pub_;
  size_t delta_pub_;
  size_t reference_pose_pub_;
  size_t robot_ref_pose_pub_;
  size_t delta_W_pub_;
  size_t human_ref_pose_pub_;
  size_t human_wrench_pub_ ;  
  
  ros::Subscriber sub_;

  Eigen::Vector6d M_;
  Eigen::Vector6d M_inv_;
  Eigen::Vector6d D_;
  Eigen::Vector6d K_;
  
  Eigen::Vector6d mask_;
  
  bool getImpedanceParams(Eigen::Vector6d& M, Eigen::Vector6d& C, Eigen::Vector6d& K );
  Eigen::Vector6d getMask();
  bool getWeightMatrix(const std::string & param, const int & size, Eigen::MatrixXd& W);
  bool getSSMatrix(const int dofs, const Eigen::Vector6d& M_inv, const Eigen::Vector6d& C, const Eigen::Vector6d& K, Eigen::MatrixXd& A, Eigen::MatrixXd& B);
  
  void computeGains();
  void updateGTMatrices(const double& alpha );
  
  void wrenchCallback             (const geometry_msgs::WrenchStampedConstPtr& msg );
  void setRobotTargetPoseCallback (const geometry_msgs::PoseStampedConstPtr&   msg );
  void setHumanTargetPoseCallback (const geometry_msgs::PoseStampedConstPtr&   msg );
  void setTargetJointsCallback    (const sensor_msgs::JointStateConstPtr&      msg );
  void setAlpha                   (const std_msgs_stamped::Float32StampedConstPtr&  msg );
  
  bool eigVecToWrenchMsg(const Eigen::Vector6d& vec, geometry_msgs::Wrench&       msg);
  bool eigToTwistMsgs   (const Eigen::Vector6d& ev , geometry_msgs::TwistStamped& msg);
  
  Eigen::MatrixXd solveRiccati( const Eigen::MatrixXd &A,
                                const Eigen::MatrixXd &B,
                                const Eigen::MatrixXd &Q,
                                const Eigen::MatrixXd &R,
                                      Eigen::MatrixXd &P) ;
  
  void solveNashEquilibrium( const Eigen::MatrixXd &A,
            const Eigen::MatrixXd &B1,
            const Eigen::MatrixXd &B2,
            const Eigen::MatrixXd &Q1,
            const Eigen::MatrixXd &Q2,
            const Eigen::MatrixXd &R1,
            const Eigen::MatrixXd &R2, 
            const Eigen::MatrixXd &R12,
            const Eigen::MatrixXd &R21, 
            Eigen::MatrixXd &P1,Eigen::MatrixXd &P2);
  
  
  
};


}
}
