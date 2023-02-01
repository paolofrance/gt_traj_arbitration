#include <gt_traj_arbitration/gt_traj_arbitration.h>
#include <gt_traj_arbitration/utils.h>
#include <state_space_filters/filtered_values.h>
#include <eigen_matrix_utils/overloads.h>
#include <geometry_msgs/PoseStamped.h>
#include <pluginlib/class_list_macros.h>
#include <Eigen/Dense>
#include <eigen_conversions/eigen_msg.h>
#include <std_msgs/Float32.h>
#include <sensor_msgs/JointState.h>
#include <rosdyn_core/primitives.h>
#include <name_sorting/name_sorting.h>
#include <tf2_eigen/tf2_eigen.h>
#include <tf/transform_broadcaster.h>
#include <tf_conversions/tf_eigen.h>

PLUGINLIB_EXPORT_CLASS(cnr::control::GTTrajArbitration  , controller_interface::ControllerBase)



namespace cnr
{
namespace control
{


/**
 * @brief GTTrajArbitration::GTTrajArbitration
 */
GTTrajArbitration::GTTrajArbitration()
{
}
GTTrajArbitration::~GTTrajArbitration()
{
  gain_thread_->join();
}


bool GTTrajArbitration::getImpedanceParams(Eigen::Vector6d& M, Eigen::Vector6d& C, Eigen::Vector6d& K )
{
  std::vector<double> M_r(6,0), D_r(6,0), K_r(6,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "M_r", M_r, 6 , "<=" );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "K_r", K_r, 6 , "<"  );
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "D_r", D_r, 6 , "<"  );
  bool is_damping_ratio;
  GET_AND_RETURN( this->getControllerNh(), "damping_is_ratio", is_damping_ratio);
  
  if (is_damping_ratio)
  {
    Eigen::Vector6d d_tmp;
      for (int i=0; i<D_r.size();i++)
        d_tmp(i) = D_r.data()[i] * 2 * std::sqrt( M_r.data()[i] * K_r.data()[i] );
      C = d_tmp;
  }
  else
      C = Eigen::Vector6d( D_r.data() );
  
  M = Eigen::Vector6d( M_r.data() );
  K = Eigen::Vector6d( K_r.data() );

  return true;
}

Eigen::Vector6d GTTrajArbitration::getMask()
{
  std::vector<double> mask(6,0);
  if (!this->getControllerNh().getParam("mask", mask))
  {
    CNR_INFO(this->logger(),"mask not found ! default all active (mask = [1,1,1,1,1,1]) ");
    mask = {1,1,1,1,1,1};
  }
  
  for (size_t i=0; i<mask.size();i++)
  {
    if (mask.at(i) !=0.0 && mask.at(i) != 1.0)
    {
      CNR_WARN(this->logger(),"mask at "<<i<<" is not zero nor one. setting it to one");
      mask.at(i)=1;
    }
  }
  
  return Eigen::Vector6d( mask.data() );
}


bool GTTrajArbitration::getWeightMatrix(const std::string & param, const int & size, Eigen::MatrixXd & W)
{  
  std::vector<double> w(size,0);
  GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), param, w, size, "" );
  Eigen::VectorXd W_vec = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(w.data(), w.size());
  W = W_vec.asDiagonal();
  return true;
}

bool GTTrajArbitration::getSSMatrix(const int dofs, const Eigen::Vector6d& M_inv, const Eigen::Vector6d& C, const Eigen::Vector6d& K, Eigen::MatrixXd& A, Eigen::MatrixXd& B)
{
  
  Eigen::VectorXd km = - K.cwiseProduct(M_inv).segment(0,dofs);
  Eigen::VectorXd cm = - C.cwiseProduct(M_inv).segment(0,dofs);
  
  A.topLeftCorner    (dofs, dofs) = Eigen::MatrixXd::Zero(dofs, dofs);
  A.topRightCorner   (dofs,dofs)  = Eigen::MatrixXd::Identity(dofs, dofs);
  A.bottomLeftCorner (dofs,dofs)  = km.asDiagonal();
  A.bottomRightCorner(dofs,dofs)  = cm.asDiagonal();
  
//   B.topLeftCorner    (2*dofs, 2*dofs) = Eigen::MatrixXd::Zero(2*dofs, 2*dofs);
//   B.bottomLeftCorner (dofs,dofs)  = M_inv.segment(0,dofs).asDiagonal();
//   B.bottomRightCorner(dofs,dofs)  = M_inv.segment(0,dofs).asDiagonal();
  
  B.resize(2*n_dofs_,n_dofs_); B_single_.setZero();
  B.topLeftCorner   (n_dofs_, n_dofs_) = Eigen::MatrixXd::Zero(n_dofs_, n_dofs_);
  B.bottomLeftCorner(n_dofs_, n_dofs_) = M_inv_.segment(0,n_dofs_).asDiagonal();
  
  return true;
}

// void GTTrajArbitration::updateGTMatrices(const double& alpha )
// {
//   Qh_ = alpha*Qhh_ + (1-alpha)*Qrh_;
//   Qr_ = alpha*Qhr_ + (1-alpha)*Qrr_;
//   
//   Q_gt_ = Qh_ + Qr_;
//   
//   R_gt_.topLeftCorner    (Rh_.cols(),Rh_.cols()) = alpha*Rh_;
//   R_gt_.bottomRightCorner(Rh_.cols(),Rh_.cols()) = (1-alpha)*Rr_;
// }

void GTTrajArbitration::firstCycleInit()
{
  q_sp_ = this->getPosition();
  dq_sp_ = this->getVelocity();
  T_human_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  T_robot_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  robot_pose_sp_.pose = tf2::toMsg (T_robot_base_targetpose_);
  human_pose_sp_.pose = tf2::toMsg (T_robot_base_targetpose_);
  first_cycle_ = false;
  count_=0;
  
  Eigen::Affine3d T_b_t_init = chain_bt_->getTransformation(q_); 
  Eigen::Vector3d cp = T_b_t_init.translation();
  initial_pose_ = cp(0);
}

bool GTTrajArbitration::eigVecToWrenchMsg(const Eigen::Vector6d& vec, geometry_msgs::Wrench& msg)
{
  msg.force.x = vec(0);
  msg.force.y = vec(1);
  msg.force.z = vec(2);
  msg.torque.x = vec(3);
  msg.torque.y = vec(4);
  msg.torque.z = vec(5);
  return true;
}

bool GTTrajArbitration::eigToTwistMsgs(const Eigen::Vector6d& ev, geometry_msgs::TwistStamped& msg)
{
  msg.twist.linear.x  = ev[0];
  msg.twist.linear.y  = ev[1];
  msg.twist.linear.z  = ev[2];
  msg.twist.angular.x = ev[3];
  msg.twist.angular.y = ev[4];
  msg.twist.angular.z = ev[5];
  return true;
}


Eigen::VectorXd GTTrajArbitration::getReference(const Eigen::Affine3d& targetpose, const Eigen::Affine3d& T_bt)
{
  Eigen::VectorXd ret;  ret.resize(2*n_dofs_); 
  Eigen::Matrix<double,6,1> cartesian_error;
  rosdyn::getFrameDistance(targetpose , T_bt, cartesian_error);
    
  ret.segment(0,n_dofs_) = cartesian_error.segment(0,n_dofs_);
  ret.segment(n_dofs_,n_dofs_) = Eigen::VectorXd::Zero(n_dofs_);
  
  return ret;
  
}

// void GTTrajArbitration::computeGains()
// {
//   while(ros::ok())
//   {
//     gains_mtx_.lock();
//     switch(ctr_switch_)
//     {
//       case GTTrajArbitration::Control::CGT :
//       {
//         CGT_gain_ = solveRiccati(A_,B_,Q_gt_,R_gt_,P_);
//         break;
//       }
//       case GTTrajArbitration::Control::LQR :
//       {
//         Eigen::MatrixXd Ph_lq,Pr_lq;    
//         Kh_lqr_ = solveRiccati(A_,B_single_,Qh_,Rh_,Ph_lq);
//         Kr_lqr_ = solveRiccati(A_,B_single_,Qr_,Rr_,Pr_lq);
//         break;
//       }
//       case GTTrajArbitration::Control::NCGT :
//       {
//         Eigen::MatrixXd Ph_nc,Pr_nc;
//         solveNashEquilibrium(A_,B_single_,B_single_,Qh_,Qr_,Rh_,Rr_,Eigen::MatrixXd::Zero(n_dofs_, n_dofs_),Eigen::MatrixXd::Zero(n_dofs_, n_dofs_),Ph_nc,Pr_nc);
//         Kh_nc_ = Rh_.inverse()*B_single_.transpose()*Ph_nc;
//         Kr_nc_ = Rr_.inverse()*B_single_.transpose()*Pr_nc;
//         break;
//       }
//       case GTTrajArbitration::Control::ICGT :
//       {
//         CGT_gain_ = solveRiccati(A_,B_,Q_gt_,R_gt_,P_);
//         break;
//       }
//       default:
//       {
//         CGT_gain_ = solveRiccati(A_,B_,Q_gt_,R_gt_,P_);
//         break;
//       }
//     }
//     gains_mtx_.unlock();
//     
//     new_gain_available_ = true;
//     
//     ros::Duration(0.01).sleep();
//   }
//   return;
// }

bool GTTrajArbitration::doInit()
{
  
try
{
  
  ROS_INFO_STREAM(cnr_logger::BLUE()<< "qui"<<this->jointNames().size());
  
  for(auto n:this->jointNames())
    ROS_INFO_STREAM(cnr_logger::BLUE()<< "joint name: "<<n);
  
  
  
  q_sp_ .resize(this->jointNames().size());  q_sp_ .setZero();
  dq_sp_.resize(this->jointNames().size());  dq_sp_.setZero();
  q_    .resize(this->jointNames().size());  q_    .setZero();
  dq_   .resize(this->jointNames().size());  dq_   .setZero();
  ddq_  .resize(this->jointNames().size());  ddq_  .setZero();

  ROS_INFO_STREAM(cnr_logger::BLUE()<< "q_sp_ .: "<<q_sp_ .transpose());  
  ROS_INFO_STREAM(cnr_logger::BLUE()<< "dq_sp_.: "<<dq_sp_.transpose());  
  ROS_INFO_STREAM(cnr_logger::BLUE()<< "q_    .: "<<q_    .transpose());  
  ROS_INFO_STREAM(cnr_logger::BLUE()<< "dq_   .: "<<dq_   .transpose());  
  ROS_INFO_STREAM(cnr_logger::BLUE()<< "ddq_  .: "<<ddq_  .transpose());  
  
  
  
  
  
  
  
  
  
  GET_AND_RETURN(this->getControllerNh(),"n_dofs",n_dofs_);
  
  A_   .resize(2*n_dofs_,2*n_dofs_);   A_   .setZero();
  B_   .resize(2*n_dofs_,n_dofs_  );   B_   .setZero();
  Qhh_ .resize(2*n_dofs_,2*n_dofs_);   Qhh_ .setZero();
  Qhr_ .resize(2*n_dofs_,2*n_dofs_);   Qhr_ .setZero();
  Qrh_ .resize(2*n_dofs_,2*n_dofs_);   Qrh_ .setZero();
  Qrr_ .resize(2*n_dofs_,2*n_dofs_);   Qrr_ .setZero();
  Qh_  .resize(2*n_dofs_,2*n_dofs_);   Qh_  .setZero();
  Qr_  .resize(2*n_dofs_,2*n_dofs_);   Qr_  .setZero();  
  Q_gt_.resize(2*n_dofs_,2*n_dofs_);   Q_gt_.setZero();
  R_gt_.resize(2*n_dofs_,2*n_dofs_);   R_gt_.setZero();
  
  CGT_gain_.resize(2*n_dofs_,2*n_dofs_); CGT_gain_.setZero();
//   P_       .resize(2*n_dofs_,2*n_dofs_); P_       .setZero();
  
  
  //INIT PUB/SUB
  {
    std::string external_wrench_topic ;
    GET_AND_RETURN( this->getControllerNh(), "external_wrench_topic"  , external_wrench_topic );
    this->template add_subscriber<geometry_msgs::WrenchStamped>(
          external_wrench_topic,5,boost::bind(&GTTrajArbitration::wrenchCallback,this,_1), false);
    
    this->template add_subscriber<std_msgs::Float32>("/alpha",5,boost::bind(&GTTrajArbitration::setAlpha,this,_1), false);
    this->template add_subscriber<geometry_msgs::PoseStamped>("/human_target",5,boost::bind(&GTTrajArbitration::setHumanTargetPoseCallback,this,_1), false);
    
    GET_AND_DEFAULT( this->getControllerNh(), "robot_active", robot_active_, true);
    GET_AND_DEFAULT( this->getControllerNh(), "use_cartesian_reference", use_cartesian_reference_, false);
    if(use_cartesian_reference_)
    {
      std::string pose_target;
      GET_AND_RETURN( this->getControllerNh(), "pose_target"  , pose_target);
      this->template add_subscriber<geometry_msgs::PoseStamped>(pose_target,5,boost::bind(&GTTrajArbitration::setRobotTargetPoseCallback,this,_1), false);
    }
    else
    {
      std::string joint_target;
      GET_AND_RETURN( this->getControllerNh(), "joint_target_topic"  , joint_target);
      this->template add_subscriber<sensor_msgs::JointState>(joint_target,5,boost::bind(&GTTrajArbitration::setTargetJointsCallback,this,_1), false);
    }
  }
  
  service_ = this->getControllerNh().advertiseService("/update_gt", &GTTrajArbitration::updateGTSrv, this);

  {  // URDF parsing and chain creation
    urdf::Model urdf_model;
    if ( !urdf_model.initParam ( "/robot_description" ) ) 
    {
        ROS_ERROR ( "Urdf robot_description '%s' does not exist", (  this->getControllerNamespace()+"/robot_description" ).c_str() );
        return false;
    }
    
    Eigen::Vector3d gravity;  gravity << 0, 0, -9.806;
    std::string robot_base_frame; std::string  robot_tip_frame; std::string  force_sensor_frame;
    GET_AND_RETURN( this->getControllerNh(), "robot_tip_frame"   , robot_tip_frame);
    GET_AND_RETURN( this->getControllerNh(), "robot_base_frame"  , robot_base_frame);
    GET_AND_RETURN( this->getControllerNh(), "force_sensor_frame", force_sensor_frame);
    
    chain_bs_ = rosdyn::createChain ( urdf_model,robot_base_frame, force_sensor_frame, gravity );
    chain_bt_ = rosdyn::createChain ( urdf_model,robot_base_frame, robot_tip_frame   , gravity );
  }

  {  // wrench filter params
    wrench_deadband_.setZero();
    w_b_            .setZero();
    
    GET_AND_DEFAULT(this->getControllerNh(),"use_filtered_wrench",use_filtered_wrench_,false);
    std::vector<double> wrench_deadband(6,0);
    GET_PARAM_VECTOR_AND_RETURN ( this->getControllerNh(), "wrench_deadband", wrench_deadband, 6, "<=" );
    
    double omega;
    GET_AND_DEFAULT(this->getControllerNh(),"omega_wrench",omega,10.0);
    ect::FilteredVectorXd::Value dead_band;
    ect::FilteredVectorXd::Value saturation;
    ect::FilteredVectorXd::Value init_value;
    
    dead_band  = Eigen::Vector6d( wrench_deadband.data() );
    saturation = 1000.0 * dead_band;
    init_value = 0.0 * dead_band;
    
    if(!wrench_fitler_.activateFilter ( dead_band, saturation, (omega / (2 * M_PI)), this->m_sampling_period, init_value ))
    {
      CNR_RETURN_FALSE(this->logger());
    }
    
    wrench_deadband_   = Eigen::Vector6d( wrench_deadband.data() );
    w_b_filt_ = wrench_fitler_.getUpdatedValue();
  }
  
  
  {  // IMpedance PArameters and mass inverse
    getImpedanceParams(M_,D_,K_);
    for (unsigned int iAx=0;iAx<6;iAx++)
      M_inv_(iAx)=1.0/M_(iAx);  
  }
  
  mask_ = getMask();
  
  std::vector<double> wrench_gains(6,0);
  if (!this->getControllerNh().getParam("wrench_gains", wrench_gains))
  {
    wrench_gains = {1,1,1,1,1,1};
    CNR_INFO(this->logger(),"wrench_gains not found ! default (wrench_gains = 1,1,1,1,1,1) ");
  }
  wrench_gains_ = Eigen::Vector6d (wrench_gains.data());
  
  getSSMatrix(n_dofs_,M_inv_,D_,K_,A_,B_);
  
  getWeightMatrix("Qhh",2*n_dofs_,Qhh_);
  getWeightMatrix("Qhr",2*n_dofs_,Qhr_);
  getWeightMatrix("Qrh",2*n_dofs_,Qrh_);
  getWeightMatrix("Qrr",2*n_dofs_,Qrr_);

  getWeightMatrix("Rh",n_dofs_,Rh_);
  getWeightMatrix("Rr",n_dofs_,Rr_);
  
  GET_AND_DEFAULT( this->getControllerNh(), "alpha", alpha_,0.5);
  GET_AND_DEFAULT( this->getControllerNh(), "alpha_max", alpha_max_,0.95);
  GET_AND_DEFAULT( this->getControllerNh(), "alpha_min", alpha_min_,0.05);
  GET_AND_DEFAULT( this->getControllerNh(), "alpha_switch", alpha_switch_,0.5);
  
  GET_AND_DEFAULT( this->getControllerNh(), "arbitrate", arbitrate_,true);
  
  GET_AND_DEFAULT( this->getControllerNh(), "use_same_reference", use_same_reference_,false);
  
  
//   B_single_.resize(2*n_dofs_,n_dofs_); B_single_.setZero();
//   B_single_.topLeftCorner   (n_dofs_, n_dofs_) = Eigen::MatrixXd::Zero(n_dofs_, n_dofs_);
//   B_single_.bottomLeftCorner(n_dofs_, n_dofs_) = M_inv_.segment(0,n_dofs_).asDiagonal();
  
  
  {  // INIT CGT
    cgt_= new CoopGT(n_dofs_, this->m_sampling_period);
    cgt_->setSysParams(A_,B_);
    cgt_->setCostsParams(Qhh_,Qhr_,Qrh_,Qrr_,Rh_,Rr_);
    cgt_->computeCooperativeGains(alpha_);

    CGT_gain_ = cgt_->getCooperativeGains();    
    ROS_INFO_STREAM("CGT_gain:\n "<<CGT_gain_);
  }
  {  // INIT NCGT
    if(!cgt_->getCostMatrices(Qh_,Qr_,Rh_,Rr_))
    {
      CNR_ERROR(this->logger(),"cost parameters not set!");
      return false;
    }
    ncgt_= new NonCoopGT(n_dofs_, this->m_sampling_period);
    ncgt_->setSysParams(A_,B_);
    ncgt_->setCostsParams(Qh_,Qr_,Rh_,Rr_);
    ncgt_->computeNonCooperativeGains();
    ncgt_->getNonCooperativeGains(Kh_nc_,Kr_nc_);    
    ROS_INFO_STREAM("NCGT_gain:\n "<<Kh_nc_);
    ROS_INFO_STREAM("NCGT_gain:\n "<<Kr_nc_);
  }

  
  

  CNR_INFO(this->logger(),CYAN<<"A\n"   << A_);
  CNR_INFO(this->logger(),CYAN<<"B\n"   << B_);
  CNR_INFO(this->logger(),CYAN<<"Q_gt\n"<< Q_gt_);
  CNR_INFO(this->logger(),CYAN<<"Q_h\n" << Qh_);
  CNR_INFO(this->logger(),CYAN<<"Q_r\n" << Qr_);
  CNR_INFO(this->logger(),CYAN<<"R_Gt\n"<< R_gt_);
  CNR_INFO(this->logger(),CYAN<<"K_gt\n"<< CGT_gain_);
  
  
  w_b_init_ = false;
  first_cycle_ = true;  
  
  GET_AND_DEFAULT(this->getControllerNh(),"control_switch",switch_control_type_,"ncgt");
  GET_AND_DEFAULT(this->getControllerNh(),"base_control",base_control_type_,"cgt");
  
  
  
//   Eigen::MatrixXd Ph_lq,Pr_lq;    
//   Kh_lqr_ = solveRiccati(A_,B_single_,Qh_,Rh_,Ph_lq);
//   Kr_lqr_ = solveRiccati(A_,B_single_,Qr_,Rr_,Pr_lq);
  
//   Eigen::MatrixXd Ph_nc,Pr_nc;
//   solveNashEquilibrium(A_,B_single_,B_single_,Qh_,Qr_,Rh_,Rr_,Eigen::MatrixXd::Zero(n_dofs_, n_dofs_),Eigen::MatrixXd::Zero(n_dofs_, n_dofs_),Ph_nc,Pr_nc);
//   Kh_nc_ = Rh_.inverse()*B_single_.transpose()*Ph_nc;
//   Kr_nc_ = Rr_.inverse()*B_single_.transpose()*Pr_nc;
    
  filtered_wrench_base_pub_ = this->template add_publisher<geometry_msgs::WrenchStamped>("/filtered_wrench_base",5);
  wrench_base_pub_          = this->template add_publisher<geometry_msgs::WrenchStamped>("/wrench_base",5);
  wrench_tool_pub_          = this->template add_publisher<geometry_msgs::WrenchStamped>("/wrench_tool",5);
  robot_wrench_pub_         = this->template add_publisher<geometry_msgs::WrenchStamped>("/robot_wrench",5);
  nominal_h_wrench_pub_     = this->template add_publisher<geometry_msgs::WrenchStamped>("/nominal_h_wrench",5);
  current_pose_pub_         = this->template add_publisher<geometry_msgs::PoseStamped>  ("/current_pose",5);
  current_vel_pub_          = this->template add_publisher<geometry_msgs::TwistStamped> ("/current_velocity",5);
  delta_pub_                = this->template add_publisher<geometry_msgs::TwistStamped> ("/delta_error",5);
  reference_pose_pub_       = this->template add_publisher<geometry_msgs::PoseStamped>  ("/gt_reference",5);
  robot_ref_pose_pub_       = this->template add_publisher<geometry_msgs::PoseStamped>  ("/robot_ref_pos",5);
  human_ref_pose_pub_       = this->template add_publisher<geometry_msgs::PoseStamped>  ("/human_ref_pos",5);
  human_wrench_pub_         = this->template add_publisher<geometry_msgs::WrenchStamped>("/human_wrench",5);
  delta_W_pub_              = this->template add_publisher<geometry_msgs::WrenchStamped>("/delta_force",5);
  kp_pub_              = this->template add_publisher<std_msgs::Float32>("/Kp",5);
  kv_pub_              = this->template add_publisher<std_msgs::Float32>("/Kv",5);
  alpha_pub_              = this->template add_publisher<std_msgs::Float32>("/alpha_gt",5);
  
//   gain_thread_= new std::thread(&GTTrajArbitration::computeGains,this);
  
  CNR_INFO(this->logger(),"intialized !!");
  CNR_RETURN_TRUE(this->logger());
  
  
}
catch(std::exception& e)
{
  CNR_ERROR( this->logger(), cnr_logger::RED() << "got exception !" );
  
  CNR_RETURN_FALSE(this->logger());

  
}

}



bool GTTrajArbitration::doStarting(const ros::Time& time)
{
  CNR_TRACE_START(this->logger(),"Starting Controller");

  q_sp_  = this->getPosition();
  dq_sp_ = this->getVelocity();
  q_  = q_sp_;
  dq_ = dq_sp_;
  this->setCommandPosition(q_);
  this->setCommandVelocity(dq_);
  T_robot_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  robot_pose_sp_.pose = tf2::toMsg (T_robot_base_targetpose_);
  
  CNR_WARN(this->logger(),robot_pose_sp_);

  dq_sp_ = 0 * this->getVelocity();

  CNR_RETURN_TRUE(this->logger());
}


bool GTTrajArbitration::doStopping(const ros::Time& time)
{
  CNR_TRACE_START(this->logger(),"Stopping Controller");
  CNR_RETURN_TRUE(this->logger());
}


bool GTTrajArbitration::doUpdate(const ros::Time& time, const ros::Duration& period)
{
  auto start = std::chrono::steady_clock::now();

  CNR_TRACE_START_THROTTLE_DEFAULT(this->logger());
  std::stringstream report;
  std::lock_guard<std::mutex> lock(m_mtx);

  if (first_cycle_)
    firstCycleInit();
  
  Eigen::Vector6d wrench = (use_filtered_wrench_) ? w_b_filt_ : w_b_;
  
  if(use_cartesian_reference_)
  {
    tf2::fromMsg (robot_pose_sp_.pose, T_robot_base_targetpose_);
  }
  else
  {
    T_robot_base_targetpose_ = chain_bt_->getTransformation(q_sp_);
  }
  
  tf2::fromMsg (human_pose_sp_.pose, T_human_base_targetpose_);
  
  Eigen::Affine3d T_b_t = chain_bt_->getTransformation(q_); 
  
  { // for showing TFs on Rviz
    tf::Transform t;
    tf::transformEigenToTF(T_robot_base_targetpose_, t);
    
    tf::Transform tbt;
    tf::transformEigenToTF(T_b_t, tbt);
    
    static tf::TransformBroadcaster br;
    br.sendTransform(tf::StampedTransform(t, ros::Time::now(), "base_link", "robot_target_pose"));
    br.sendTransform(tf::StampedTransform(tbt, ros::Time::now(), "base_link", "TBT")); 
  }

  ctr_switch_ = (alpha_<alpha_switch_) ? control_map_.at(switch_control_type_) : control_map_.at(base_control_type_);
  CNR_INFO_THROTTLE(this->logger(),5.0,"control type active: " << control_map_reverse_.at(ctr_switch_)<<" . alpha : " << alpha_);
  
  Eigen::VectorXd control(2*n_dofs_); control.setZero();
  double Kp,Kv;
  
  Eigen::Matrix6Xd J_of_t_in_b  = chain_bt_->getJacobian(q_);
  Eigen::Vector6d cart_vel_of_t_in_b  = J_of_t_in_b*dq_;
  Eigen::Vector6d cart_acc_nl_of_t_in_b  = chain_bt_->getDTwistNonLinearPartTool(q_,dq_); // DJ*Dq
  Eigen::Vector6d cart_acc_of_t_in_b; cart_acc_of_t_in_b.setZero();
  
  Eigen::Vector3d cart_pos = T_b_t.translation();
  Eigen::VectorXd X; X.resize(2*n_dofs_);
  X << cart_pos.segment(0,n_dofs_), cart_vel_of_t_in_b.segment(0,n_dofs_); 
  
  Eigen::VectorXd ref_h = T_human_base_targetpose_.translation().segment(0,n_dofs_);
  Eigen::VectorXd ref_r = T_robot_base_targetpose_.translation().segment(0,n_dofs_);

  if (use_same_reference_)
    Eigen::VectorXd ref_r = ref_h;
  
  switch(ctr_switch_)
  {
    case GTTrajArbitration::Control::CGT :
    {
      cgt_->computeCooperativeGains(alpha_);
      CGT_gain_ = cgt_->getCooperativeGains();
      cgt_->setPosReference(ref_h,ref_r);
      cgt_->setCurrentState(X);
      control = cgt_->computeControlInputs();
      Kp = CGT_gain_(n_dofs_,0);
      Kv = CGT_gain_(n_dofs_,n_dofs_);
      break;
    }
    case GTTrajArbitration::Control::NCGT :
    {      
      cgt_->updateGTMatrices(alpha_);
      if(!cgt_->getCostMatrices(Qh_,Qr_,Rh_,Rr_))
      {
        CNR_ERROR(this->logger(),"cost parameters not set!");
        return false;
      }

      ncgt_->setCostsParams(Qh_,Qr_,Rh_,Rr_);
      ncgt_->computeNonCooperativeGains();
      ncgt_->getNonCooperativeGains(Kh_nc_,Kr_nc_);
      ncgt_->setPosReference(ref_h,ref_r);
      ncgt_->setCurrentState(X);
      
      Kp = Kr_nc_(0,0);
      Kv = Kr_nc_(n_dofs_,n_dofs_);
      
      control = ncgt_->computeControlInputs();
      
      break;
    }
  }
    

    
  Eigen::Vector6d human_wrench_ic = mask_.cwiseProduct(wrench);
  
  Eigen::Vector6d nominal_human_wrench_ic; nominal_human_wrench_ic.setZero(); 
  nominal_human_wrench_ic.segment(0,n_dofs_) = control.segment(0,n_dofs_);
    
  Eigen::Vector6d robot_wrench_ic; robot_wrench_ic.setZero(); 
  
  
  Eigen::Matrix<double,6,1> robot_cartesian_error_actual_target_in_b;
  rosdyn::getFrameDistance(T_robot_base_targetpose_ , T_b_t, robot_cartesian_error_actual_target_in_b);
  Eigen::Vector6d global_reference;
  if(robot_active_)
  {
    global_reference.segment(0,n_dofs_) = cgt_->getReference().segment(0,n_dofs_);
    global_reference.segment(n_dofs_,6-n_dofs_) = robot_cartesian_error_actual_target_in_b.segment(n_dofs_,6-n_dofs_);
    robot_wrench_ic.segment(0,n_dofs_) = control.segment(n_dofs_,n_dofs_);  
}
  else
  {
    global_reference = robot_cartesian_error_actual_target_in_b;
  }
  
  
  CNR_DEBUG_THROTTLE(this->logger(),2.0,"human_wrench_ic \n" <<human_wrench_ic .transpose());
  CNR_DEBUG_THROTTLE(this->logger(),2.0,"robot_wrench_ic \n" <<robot_wrench_ic .transpose());
  
  cart_acc_of_t_in_b = (M_inv_).cwiseProduct(
                        K_.cwiseProduct(global_reference) +
                        D_.cwiseProduct(-cart_vel_of_t_in_b) +
                        human_wrench_ic + robot_wrench_ic);
  
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(J_of_t_in_b, Eigen::ComputeThinU | Eigen::ComputeThinV);
  if (svd.singularValues()(svd.cols()-1)==0)
    ROS_WARN_THROTTLE(1,"SINGULARITY POINT");
  else if (svd.singularValues()(0)/svd.singularValues()(svd.cols()-1) > 1e2)
    ROS_WARN_THROTTLE(1,"SINGULARITY POINT");
  
  
  // INVERSE KINEMATICS
  
  ddq_ = svd.solve(cart_acc_of_t_in_b-cart_acc_nl_of_t_in_b);  
  q_  += dq_  * period.toSec() + ddq_*std::pow(period.toSec(),2.0)*0.5;
  dq_ += ddq_ * period.toSec();


  this->setCommandPosition( q_ );
  this->setCommandVelocity( dq_);
  
  
  { // publishing section
    geometry_msgs::Pose cp = tf2::toMsg (T_b_t);
    
    ros::Time stamp = ros::Time::now();
    
    geometry_msgs::PoseStamped ps;
    ps.header.stamp = stamp;
    ps.pose = cp;  
    this->publish(current_pose_pub_,ps);
    
    {
      geometry_msgs::TwistStamped cv;
      cv.header.stamp = stamp;
      eigToTwistMsgs(cart_vel_of_t_in_b,cv);
      this->publish(current_vel_pub_,cv);
    }
    {
      geometry_msgs::PoseStamped p;
      p.pose = tf2::toMsg (T_robot_base_targetpose_);
      p.header.stamp = stamp;
      this->publish(robot_ref_pose_pub_,p);
    }
    {
      geometry_msgs::PoseStamped p;
      p.pose = tf2::toMsg (T_human_base_targetpose_);
      p.header.stamp = stamp;
      this->publish(human_ref_pose_pub_,p);
    }
    {
      geometry_msgs::WrenchStamped w;
      eigVecToWrenchMsg(wrench,w.wrench);
      w.header.stamp = stamp;
      this->publish(human_wrench_pub_,w);
    }
    
    {
      geometry_msgs::PoseStamped ref;
      ref.header.stamp = stamp;
      
      Eigen::Vector3d world_reference = Q_gt_.inverse() * ( Qh_*T_human_base_targetpose_.translation() + Qr_*T_robot_base_targetpose_.translation() );
      ref.pose.position.x = world_reference(0);
      ref.pose.position.y = world_reference(1);
      ref.pose.position.z = world_reference(2);
      
      ref.pose.orientation = ps.pose.orientation;
      
      this->publish(reference_pose_pub_,ref);
    }
    
    {
      std_msgs::Float32 m;
      m.data = Kp;
      this->publish(kp_pub_,m);
    }
    {
      std_msgs::Float32 m;
      m.data = Kv;
      this->publish(kv_pub_,m);
    }
    {
      std_msgs::Float32 m;
      m.data = alpha_;
      this->publish(alpha_pub_,m);
    }
    
    
    geometry_msgs::WrenchStamped human_w,robot_w;
    eigVecToWrenchMsg(nominal_human_wrench_ic,human_w.wrench);
    eigVecToWrenchMsg(robot_wrench_ic,robot_w.wrench);
    
    human_w.header.stamp = stamp;
    robot_w.header.stamp = stamp;
    
    this->publish(robot_wrench_pub_,robot_w);
    this->publish(nominal_h_wrench_pub_,human_w);
  }
  
  
  auto mid = std::chrono::steady_clock::now();
  CNR_INFO_COND(this->logger(),std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count()>=8000
                 ,RED<<"too much time to command: "<<std::chrono::duration_cast<std::chrono::microseconds>(mid - start).count());

  CNR_RETURN_TRUE_THROTTLE_DEFAULT(this->logger());

  }


  
  void GTTrajArbitration::wrenchCallback(const geometry_msgs::WrenchStampedConstPtr& msg )
  {
    if(!w_b_init_)
    {
      w_b_0_ ( 0 ) = msg->wrench.force.x;
      w_b_0_ ( 1 ) = msg->wrench.force.y;
      w_b_0_ ( 2 ) = msg->wrench.force.z;
      w_b_0_ ( 3 ) = msg->wrench.torque.x;
      w_b_0_ ( 4 ) = msg->wrench.torque.y;
      w_b_0_ ( 5 ) = msg->wrench.torque.z;

      w_b_init_ = true;
    }
    
    Eigen::Vector6d wrench_s;
    wrench_s( 0 ) = msg->wrench.force.x  - w_b_0_ ( 0 );
    wrench_s( 1 ) = msg->wrench.force.y  - w_b_0_ ( 1 );
    wrench_s( 2 ) = msg->wrench.force.z  - w_b_0_ ( 2 );
    wrench_s( 3 ) = msg->wrench.torque.x - w_b_0_ ( 3 );
    wrench_s( 4 ) = msg->wrench.torque.y - w_b_0_ ( 4 );
    wrench_s( 5 ) = msg->wrench.torque.z - w_b_0_ ( 5 );

    Eigen::Affine3d T_bs = chain_bs_->getTransformation ( this->getPosition() );
    Eigen::Affine3d T_bt = chain_bt_->getTransformation ( this->getPosition() );
    Eigen::Affine3d T_ts = T_bt.inverse() * T_bs;
    Eigen::Vector6d w_t = rosdyn::spatialDualTranformation ( wrench_s , T_ts );
    Eigen::Vector6d wrench;
    wrench = rosdyn::spatialRotation ( w_t, T_bt.linear() );

    for ( unsigned int idx=0; idx<6; idx++ )
    {
      if ( ( wrench ( idx ) >wrench_deadband_ ( idx ) ) )
      {
          w_b_ ( idx ) = wrench ( idx )-wrench_deadband_ ( idx );
      }
      else if ( ( wrench ( idx ) <-wrench_deadband_ ( idx ) ) )
      {
          w_b_ ( idx ) = wrench ( idx )+wrench_deadband_ ( idx );
      }
      else
      {
          w_b_ ( idx ) =0;
      }
    }

    geometry_msgs::WrenchStamped tool_w;

    tool_w.header.frame_id = "robotiq_ft_frame_id";
    tool_w.header.stamp = ros::Time::now();
    tool_w.wrench.force.x  = wrench_s( 0 );
    tool_w.wrench.force.y  = wrench_s( 1 );
    tool_w.wrench.force.z  = wrench_s( 2 );
    tool_w.wrench.torque.x = wrench_s( 3 );
    tool_w.wrench.torque.y = wrench_s( 4 );
    tool_w.wrench.torque.z = wrench_s( 5 );

    geometry_msgs::WrenchStamped base_w;

//     double grav_force = M_l_*g_;
    
    base_w.header.frame_id = "ur5_base_link";
    base_w.header.stamp = ros::Time::now();
    base_w.wrench.force.x  = w_b_( 0 );
    base_w.wrench.force.y  = w_b_( 1 );
    base_w.wrench.force.z  = w_b_( 2 );// - grav_force;
    base_w.wrench.torque.x = w_b_( 3 );
    base_w.wrench.torque.y = w_b_( 4 );
    base_w.wrench.torque.z = w_b_( 5 );

    
    Eigen::Vector6d old = w_b_filt_ ;
    
    wrench_fitler_.update(wrench);
    w_b_filt_ = wrench_gains_.cwiseProduct( wrench_fitler_.getUpdatedValue() );
    
//     for(int i=0;i<6;i++)
//       w_b_filt_(i) = w_b_filt_(i)*wrench_gains_(i);
    
    CNR_INFO_THROTTLE(this->logger(),5.0, cnr_logger::BLUE()<<"subscribing wrench : "<< w_b_filt_.transpose());
    
    Eigen::Vector6d dW = (w_b_filt_ - old)/this->m_sampling_period;
    
    geometry_msgs::WrenchStamped filter_base_w, deltaW;    
    filter_base_w.header.frame_id = "ur5_base_link";
    filter_base_w.header.stamp = ros::Time::now();
    deltaW.header.stamp = ros::Time::now();
    eigVecToWrenchMsg(w_b_filt_, filter_base_w.wrench);
    eigVecToWrenchMsg(dW, deltaW.wrench);
    
    this->publish(wrench_base_pub_,base_w);
    this->publish(filtered_wrench_base_pub_ ,filter_base_w);
    this->publish(wrench_tool_pub_,tool_w);
    this->publish(delta_W_pub_,deltaW);
  }

//   Eigen::MatrixXd GTTrajArbitration::solveRiccati(const Eigen::MatrixXd &A,
//                                                   const Eigen::MatrixXd &B,
//                                                   const Eigen::MatrixXd &Q,
//                                                   const Eigen::MatrixXd &R, Eigen::MatrixXd &P)
//   {
//     const uint dim_x = A.rows();
//     const uint dim_u = B.cols();
// 
//     Eigen::MatrixXd Ham = Eigen::MatrixXd::Zero(2 * dim_x, 2 * dim_x);
//     Ham << A, -B * R.inverse() * B.transpose(), -Q, -A.transpose();
// 
//     Eigen::EigenSolver<Eigen::MatrixXd> Eigs(Ham);
// 
//     Eigen::MatrixXcd eigvec = Eigen::MatrixXcd::Zero(2 * dim_x, dim_x);
//     int j = 0;
//     for (int i = 0; i < 2 * dim_x; ++i) {
//       if (Eigs.eigenvalues()[i].real() < 0.) {
//         eigvec.col(j) = Eigs.eigenvectors().block(0, i, 2 * dim_x, 1);
//         ++j;
//       }
//     }
// 
//     Eigen::MatrixXcd Vs_1, Vs_2;
//     Vs_1 = eigvec.block(0, 0, dim_x, dim_x);
//     Vs_2 = eigvec.block(dim_x, 0, dim_x, dim_x);
//     P = (Vs_2 * Vs_1.inverse()).real();
//     
//     return R.inverse()*B.transpose()*P;
//   }
  
//   void GTTrajArbitration::solveNashEquilibrium( const Eigen::MatrixXd &A,
//                                                 const Eigen::MatrixXd &B1,
//                                                 const Eigen::MatrixXd &B2,
//                                                 const Eigen::MatrixXd &Q1,
//                                                 const Eigen::MatrixXd &Q2,
//                                                 const Eigen::MatrixXd &R1,
//                                                 const Eigen::MatrixXd &R2, 
//                                                 const Eigen::MatrixXd &R12,
//                                                 const Eigen::MatrixXd &R21, 
//                                                 Eigen::MatrixXd &P1,Eigen::MatrixXd &P2)
//   {
//     Eigen::MatrixXd S1  = B1 * R1.inverse() * B1.transpose();
//     Eigen::MatrixXd S2  = B2 * R2.inverse() * B2.transpose();
//     Eigen::MatrixXd S12 = B1 * R1.inverse() * R21 * R1.inverse() * B1.transpose();
//     Eigen::MatrixXd S21 = B2 * R2.inverse() * R12 * R2.inverse()* B2.transpose();
// 
//     solveRiccati(A,B1,Q1,R1,P1);
//     solveRiccati(A,B2,Q2,R2,P2);
//     
//     Eigen::MatrixXd P1_prev = P1;
//     Eigen::MatrixXd P2_prev = P2;
//     double err_1 = 1;
//     double err_2 = 1;
//     double toll = 0.00001;
// 
//     while (err_1>toll && err_2>toll)
//     {    
//       Eigen::MatrixXd A1 = A - S2*P2;
//       Eigen::MatrixXd A2 = A - S1*P1;
//       
//       Eigen::MatrixXd Q_1 = Q1 + P1*S21*P1;
//       solveRiccati(A1,B1,Q_1,R1,P1);
//       Eigen::MatrixXd Q_2 = Q2 + P2*S12*P2;
//       solveRiccati(A2,B2,Q_2,R2,P2);
//     
//       err_1 = (P1-P1_prev).norm();
//       err_2 = (P2-P2_prev).norm();
//       
//       P1_prev = P1;
//       P2_prev = P2;
//     }
//     
//     return;
//   }
  
  void GTTrajArbitration::setRobotTargetPoseCallback(const geometry_msgs::PoseStampedConstPtr& msg)
  {
    try
    {
      robot_pose_sp_ = *msg;
      new_sp_available_ = true;
    }
    catch(...)
    {
      ROS_ERROR("Something wrong in target callback");
    }
  }
  void GTTrajArbitration::setHumanTargetPoseCallback(const geometry_msgs::PoseStampedConstPtr& msg)
  {
    try
    {
      human_pose_sp_ = *msg;
      
    }
    catch(...)
    {
      ROS_ERROR("Something wrong in target callback");
    }
  }

  void GTTrajArbitration::setTargetJointsCallback(const sensor_msgs::JointStateConstPtr& msg)
  {
    try
    {
      sensor_msgs::JointState tmp_msg=*msg;
      if (!name_sorting::permutationName(this->jointNames(),tmp_msg.name,tmp_msg.position,tmp_msg.velocity,tmp_msg.effort))
      {
        CNR_ERROR(this->logger(),"joints not found");
        return;
      }
      for (unsigned int iAx=0;iAx<q_sp_.rows();iAx++)
      {
        q_sp_(iAx)=tmp_msg.position.at(iAx);
        dq_sp_(iAx)=tmp_msg.velocity.at(iAx);
      }

    }
    catch(...)
    {
      CNR_ERROR(this->logger(),"Something wrong in target callback");
    }
  }
  
  void GTTrajArbitration::setAlpha(const std_msgs::Float32ConstPtr& msg )
  {
    alpha_ = msg->data;
    
    if ( ( alpha_ > alpha_max_ ) )
    {
      alpha_ = alpha_max_ ;
      CNR_INFO_THROTTLE(this->logger(),2.0,"saturating alpha to max: "<<alpha_);
    }
    else if ( alpha_ < alpha_min_  )
    {
      alpha_ = alpha_min_;
      CNR_INFO_THROTTLE(this->logger(),2.0,"saturating alpha to min: "<<alpha_);
    }    
  }
  
  
  bool GTTrajArbitration::updateGTSrv(pbo_service::updateGT::Request  &req,
                                      pbo_service::updateGT::Response &res)
  {
    
    ROS_INFO_STREAM("New values available for the Rr = "<< req.Rr);
    ROS_INFO_STREAM("New values available for the alpha = "<< req.alpha);
    
    Rr_ = Eigen::MatrixXd::Identity(n_dofs_, n_dofs_) * req.Rr;
    ROS_INFO_STREAM("New Rr\n "<< Rr_);
    alpha_ = req.alpha;
    cgt_->setCostsParams(Qhh_,Qhr_,Qrh_,Qrr_,Rh_,Rr_);

    res.res = true;
    
    return true;
  }
  

}
}
