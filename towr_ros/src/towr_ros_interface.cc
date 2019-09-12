/******************************************************************************
Copyright (c) 2018, Alexander W. Winkler. All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <towr_ros/towr_ros_interface.h>

#include <std_msgs/Int32.h>

#include <xpp_states/convert.h>
#include <xpp_msgs/topic_names.h>
#include <xpp_msgs/TerrainInfo.h>

#include <towr/terrain/height_map.h>
#include <towr/variables/euler_converter.h>
#include <towr_ros/topic_names.h>
#include <towr_ros/towr_xpp_ee_map.h>
#include <ctime>



namespace towr {


TowrRosInterface::TowrRosInterface ()
{
  ::ros::NodeHandle n;

  user_command_sub_ = n.subscribe(towr_msgs::user_command, 1,
                                  &TowrRosInterface::UserCommandCallback, this);

  controller_sub_ = n.subscribe(towr_msgs::controller_command, 1,
                                  &TowrRosInterface::ReplanningCallback, this);

  initial_state_pub_  = n.advertise<xpp_msgs::RobotStateCartesian>
                                          (xpp_msgs::robot_state_desired, 1);

  robot_parameters_pub_  = n.advertise<xpp_msgs::RobotParameters>
                                    (xpp_msgs::robot_parameters, 1);

  trajectory_ = n.advertise<xpp_msgs::RobotStateCartesianTrajectory>("/xpp/trajectory_des",1);

  dwltrajectory_=n.advertise<dwl_msgs::WholeBodyTrajectory>("hyq/plan",1);

  recompute_sub = n.subscribe("/hyq/recompute", 1, &TowrRosInterface::RecomputePlan, this);

  solver_ = std::make_shared<ifopt::IpoptSolver>();

  visualization_dt_ = 0.01;
  offsetBF<<0.0229786, 5.2e-5, -0.0397;
}

BaseState
TowrRosInterface::GetGoalState(const TowrCommandMsg& msg) const
{
  BaseState goal;
  goal.lin.at(kPos) = xpp::Convert::ToXpp(msg.goal_lin.pos);
  goal.lin.at(kVel) = xpp::Convert::ToXpp(msg.goal_lin.vel);
  goal.ang.at(kPos) = xpp::Convert::ToXpp(msg.goal_ang.pos);
  goal.ang.at(kVel) = xpp::Convert::ToXpp(msg.goal_ang.vel);

  return goal;
}

BaseState TowrRosInterface::GetInitialState()
{
  return initialBaseState;
}
BaseState TowrRosInterface::GetInitialStateCoM()
{
  return initialCoMState;
}

void TowrRosInterface::ReplanningCallback(const dwl_msgs::WholeBodyController & msg){
  
  initialBaseState.ang.at(kPos)(0) = msg.actual.base[0].position;
  initialBaseState.ang.at(kPos)(1) = msg.actual.base[1].position;
  initialBaseState.ang.at(kPos)(2) = msg.actual.base[2].position;
  
  initialBaseState.ang.at(kVel)(0) = msg.actual.base[0].velocity;
  initialBaseState.ang.at(kVel)(1) = msg.actual.base[1].velocity;
  initialBaseState.ang.at(kVel)(2) = msg.actual.base[2].velocity;

  initialBaseState.lin.at(kPos)(0) = msg.actual.base[3].position;
  initialBaseState.lin.at(kPos)(1) = msg.actual.base[4].position;
  initialBaseState.lin.at(kPos)(2) = msg.actual.base[5].position;

  initialBaseState.lin.at(kVel)(0) = msg.actual.base[3].velocity;
  initialBaseState.lin.at(kVel)(1) = msg.actual.base[4].velocity;
  initialBaseState.lin.at(kVel)(2) = msg.actual.base[5].velocity;
  //initialBaseState.ang.at(kPos) = msg.actual.base;
  //initialBaseState.ang.at(kVel) = msg.actual.base;


  initial_foot_lf_B(0) = msg.actual.contacts[0].position.x;
  initial_foot_lf_B(1) = msg.actual.contacts[0].position.y;
  initial_foot_lf_B(2) = msg.actual.contacts[0].position.z;
 
  initial_foot_rf_B(0) = msg.actual.contacts[1].position.x;
  initial_foot_rf_B(1) = msg.actual.contacts[1].position.y;
  initial_foot_rf_B(2) = msg.actual.contacts[1].position.z;

  initial_foot_lh_B(0) = msg.actual.contacts[2].position.x;
  initial_foot_lh_B(1) = msg.actual.contacts[2].position.y;
  initial_foot_lh_B(2) = msg.actual.contacts[2].position.z;

  initial_foot_rh_B(0) = msg.actual.contacts[3].position.x;
  initial_foot_rh_B(1) = msg.actual.contacts[3].position.y;
  initial_foot_rh_B(2) = msg.actual.contacts[3].position.z;

  /* this is needed to decouple the planning from the base Z coordinate given by 
  the DLS framework's state estimator */
  double average_foot_height = (initial_foot_lf_B(2) + initial_foot_rf_B(2) + initial_foot_lh_B(2) + initial_foot_rh_B(2))/4.0;
  double robot_height = - average_foot_height;
  initialBaseState.lin.at(kPos)(2) = robot_height;//+foot_radius_;

  //std::cout<<"Average foot height is : "<<average_foot_height<<std::endl;
  //std::cout<<"Average robot's height is: "<<robot_height<<std::endl;

  //std::cout<<"initial foot pos B CALLBACK"<<initial_foot_rh_B.transpose()<<std::endl;


  //const EulerAngles euler_angles = initialBaseState.ang.at(kPos);
  auto base_angular=EulerConverter(solution.base_angular_);

  Eigen::Matrix3d w_R_b = base_angular.GetRotationMatrixBaseToWorld(initialBaseState.ang.at(kPos));
  initialCoMState.lin.at(kPos) = initialBaseState.lin.at(kPos)+ w_R_b*offsetBF;
  

  initial_foot_lf_W = w_R_b*initial_foot_lf_B + initialBaseState.lin.at(kPos);
  initial_foot_rf_W = w_R_b*initial_foot_rf_B + initialBaseState.lin.at(kPos);
  initial_foot_lh_W = w_R_b*initial_foot_lh_B + initialBaseState.lin.at(kPos);
  initial_foot_rh_W = w_R_b*initial_foot_rh_B + initialBaseState.lin.at(kPos);
  



  //std::cout<<"initial foot pos W CALLBACK"<<initial_foot_rh_W.transpose()<<std::endl;

}

void
TowrRosInterface::RecomputePlan(const geometry_msgs::Vector3& msg)
{
  // robot model
  formulation_.model_ = RobotModel(RobotModel::Hyq);
  auto robot_params_msg = BuildRobotParametersMsg(formulation_.model_);
  robot_parameters_pub_.publish(robot_params_msg);

  // terrain
  //auto terrain_id = static_cast<HeightMap::TerrainID>(msg.terrain);
  //formulation_.terrain_ = HeightMap::MakeTerrain(terrain_id);

 int n_ee = formulation_.model_.kinematic_model_->GetNumberOfEndeffectors();
 formulation_.params_ = GetTowrParametersReplanningCallback(n_ee, time_);
 

  

  
  std::vector<Eigen::Vector3d> initial_feet_pos;
  initial_feet_pos.push_back(initial_foot_lf_W);
  initial_feet_pos.push_back(initial_foot_rf_W);
  initial_feet_pos.push_back(initial_foot_lh_W);
  initial_feet_pos.push_back(initial_foot_rh_W);

  std::cout<<"initial foot pos LF WF: "<<initial_foot_lf_W.transpose()<<std::endl;
  std::cout<<"initial foot pos RF WF: "<<initial_foot_rf_W.transpose()<<std::endl;
  std::cout<<"initial foot pos LH WF: "<<initial_foot_lh_W.transpose()<<std::endl;
  std::cout<<"initial foot pos RH WF: "<<initial_foot_rh_W.transpose()<<std::endl;

  SetTowrInitialState(initial_feet_pos); 
  std::cout<<"initial foot pos LF BF: "<<initial_foot_lf_B.transpose()<<std::endl;
  std::cout<<"initial foot pos RF BF: "<<initial_foot_rf_B.transpose()<<std::endl;
  std::cout<<"initial foot pos LH BF: "<<initial_foot_lh_B.transpose()<<std::endl;
  std::cout<<"initial foot pos RH BF: "<<initial_foot_rh_B.transpose()<<std::endl;
  
  formulation_.initial_base_.lin.at(kPos)=initialCoMState.lin.at(kPos);
  formulation_.initial_base_.lin.at(kVel)=initialBaseState.lin.at(kVel);
  formulation_.initial_base_.ang.at(kPos)=initialBaseState.ang.at(kPos);
  formulation_.initial_base_.ang.at(kVel)=initialBaseState.ang.at(kVel);

  //formulation_.initial_ee_W_.at(0)=initial_foot_lf_W;
  //formulation_.initial_ee_W_.at(1)=initial_foot_rf_W;
  //formulation_.initial_ee_W_.at(2)=initial_foot_lh_W;
  //formulation_.initial_ee_W_.at(3)=initial_foot_rh_W

  // visualization
  PublishInitialState();

  formulation_.final_base_.lin.at(kPos).x()=msg.x;
  formulation_.final_base_.lin.at(kPos).y()=msg.y;
  formulation_.final_base_.lin.at(kPos).z()=msg.z;
  std::cout<<"goal is: "<<std::endl;
  std::cout<<formulation_.final_base_.lin.at(kPos)<<std::endl;

  // Defaults to /home/user/.ros/
  solver_ = std::make_shared<ifopt::IpoptSolver>();
  // solver parameters
  //SetIpoptParameters(msg);
  
  time_t ct;
  ct = time(NULL);
  struct tm *localTime = localtime(&ct);
  char buffer [80];
  strftime(buffer, 80, "%d-%m-%Y-%H-%M-%S", localTime);
  std::string bag_file = "towr_trajectory-" + std::string(buffer)+ ".bag";
  //std::string bag_file = "towr_trajectory.bag";
  //if (msg.optimize || msg.play_initialization) {
    nlp_ = ifopt::Problem();
    for (auto c : formulation_.GetVariableSets(solution))
      nlp_.AddVariableSet(c);
    for (auto c : formulation_.GetConstraints(solution,nlp_))
      nlp_.AddConstraintSet(c);
    for (auto c : formulation_.GetCosts())
      nlp_.AddCostSet(c);

    solver_->Solve(nlp_);
   SaveOptimizationAsRosbag(bag_file, robot_params_msg, false);
   int success = system(("rosbag play --topics "
       + xpp_msgs::robot_state_desired + " "
       + xpp_msgs::terrain_info
       + " -r " + std::to_string(1.0)
       + " --quiet " + bag_file).c_str());
    
  //}

  // playback using terminal commands
  //if (msg.replay_trajectory || msg.play_initialization || msg.optimize) {


  //if (msg.plot_trajectory) {
  //  int success = system(("killall rqt_bag; rqt_bag " + bag_file + "&").c_str());
  //}

  // to publish entire trajectory (e.g. to send to controller)
  xpp_msgs::RobotStateCartesianTrajectory xpp_msg = xpp::Convert::ToRos(GetTrajectory());

  //dwl_msgs::WholeBodyTrajectory wbtraj = ToRos();

  trajectory_.publish(xpp_msg);
  
  char a;
  std::cout<<"do you want to publish the trajectory? y/n"<<std::endl;
  std::cin>>a;
  if(a=='y'){
    std::cout<<"New trajectory being published."<<std::endl;
    ToRosAndPublish();
    //dwltrajectory_.publish(wbtraj);
  }else{
    if(a=='n'){  
      std::cout<<"New trajectory will not be published."<<std::endl;
    }else{
      std::cout<<"Please re-optimize and then enter a valid digit (either y or n)."<<std::endl;
    }
  }

  int success1 = system (("mv " + bag_file + " ~/misc_ws/src/bag_files").c_str()); //mio computer
  // int success1 = system (("mv"+ bag_file + " ~/catkin_ws/bag_files")c_str()); //hyq_furious 
  
}

void
TowrRosInterface::UserCommandCallback(const TowrCommandMsg& msg)
{
  // robot model
  formulation_.model_ = RobotModel(static_cast<RobotModel::Robot>(msg.robot));
  auto robot_params_msg = BuildRobotParametersMsg(formulation_.model_);
  robot_parameters_pub_.publish(robot_params_msg);

  // terrain
  auto terrain_id = static_cast<HeightMap::TerrainID>(msg.terrain);
  formulation_.terrain_ = HeightMap::MakeTerrain(terrain_id);

  int n_ee = formulation_.model_.kinematic_model_->GetNumberOfEndeffectors();
  formulation_.params_ = GetTowrParameters(n_ee, msg);
  formulation_.final_base_ = GetGoalState(msg);
  time_=msg.total_duration;
  SetTowrDefaultState();


  // visualization
  PublishInitialState();

  // Defaults to /home/user/.ros/
  solver_ = std::make_shared<ifopt::IpoptSolver>();
  // solver parameters
  SetIpoptParameters(msg);
  
  time_t ct;
  ct = time(NULL);
  struct tm *localTime = localtime(&ct);
  char buffer [80];
  strftime(buffer, 80, "%d-%m-%Y-%H-%M-%S", localTime);
  std::string bag_file = "towr_trajectory-" + std::string(buffer)+ ".bag";
  //std::string bag_file = "towr_trajectory.bag";
  if (msg.optimize || msg.play_initialization) {
    nlp_ = ifopt::Problem();
    for (auto c : formulation_.GetVariableSets(solution))
      nlp_.AddVariableSet(c);
    for (auto c : formulation_.GetConstraints(solution,nlp_))
      nlp_.AddConstraintSet(c);
    for (auto c : formulation_.GetCosts())
      nlp_.AddCostSet(c);

    solver_->Solve(nlp_);
    SaveOptimizationAsRosbag(bag_file, robot_params_msg, msg, false);
        // int success1 = system (("mv"+ bag_file + " ~/catkin_ws/bag_files")c_str()); //hyq_furious

    ToRos();


  }

  // playback using terminal commands
  if (msg.replay_trajectory || msg.play_initialization || msg.optimize) {
    int success = system(("rosbag play --topics "
        + xpp_msgs::robot_state_desired + " "
        + xpp_msgs::terrain_info
        + " -r " + std::to_string(msg.replay_speed)
        + " --quiet " + bag_file).c_str());
  }

  if (msg.plot_trajectory) {
    int success = system(("killall rqt_bag; rqt_bag " + bag_file + "&").c_str());
  }

  // to publish entire trajectory (e.g. to send to controller)
  xpp_msgs::RobotStateCartesianTrajectory xpp_msg = xpp::Convert::ToRos(GetTrajectory());

  dwl_msgs::WholeBodyTrajectory wbtraj;

  trajectory_.publish(xpp_msg);
  if (msg.optimize){
    int success1 = system (("mv " + bag_file + " ../misc_ws/src/bag_files").c_str()); //mio computer
    char a;
    std::cout<<"do you want to publish the trajectory? y/n"<<std::endl;
    std::cin>>a;
    if(a=='y'){
      std::cout<<"New trajectory being published."<<std::endl;
      dwltrajectory_.publish(wbtraj);
    }else{
      if(a=='n'){  
        std::cout<<"New trajectory will not be published."<<std::endl;
      }else{
        std::cout<<"Please re-optimize and then enter a valid digit (either y or n)."<<std::endl;
      }
    }
  }
}

void
TowrRosInterface::PublishInitialState()
{
  int n_ee = formulation_.initial_ee_W_.size();
  xpp::RobotStateCartesian xpp(n_ee);
  xpp.base_.lin.p_ = formulation_.initial_base_.lin.p();
  xpp.base_.ang.q  = EulerConverter::GetQuaternionBaseToWorld(formulation_.initial_base_.ang.p());

  for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
    int ee_xpp = ToXppEndeffector(n_ee, ee_towr).first;
    xpp.ee_contact_.at(ee_xpp)   = true;
    xpp.ee_motion_.at(ee_xpp).p_ = formulation_.initial_ee_W_.at(ee_towr);
    xpp.ee_forces_.at(ee_xpp).setZero(); // zero for visualization
  }

  initial_state_pub_.publish(xpp::Convert::ToRos(xpp));
}

std::vector<TowrRosInterface::XppVec>
TowrRosInterface::GetIntermediateSolutions ()
{
  std::vector<XppVec> trajectories;

  for (int iter=0; iter<nlp_.GetIterationCount(); ++iter) {
    nlp_.SetOptVariables(iter);
    trajectories.push_back(GetTrajectory());
  }

  return trajectories;
}

TowrRosInterface::XppVec
TowrRosInterface::GetTrajectory () const
{
  XppVec trajectory;
  double t = 0.0;
  double T = solution.base_linear_->GetTotalTime();

  EulerConverter base_angular(solution.base_angular_);

  while (t<=T+1e-5) {
    int n_ee = solution.ee_motion_.size();
    xpp::RobotStateCartesian state(n_ee);

    state.base_.lin = ToXpp(solution.base_linear_->GetPoint(t));

    state.base_.ang.q  = base_angular.GetQuaternionBaseToWorld(t);
    state.base_.ang.w  = base_angular.GetAngularVelocityInWorld(t);
    state.base_.ang.wd = base_angular.GetAngularAccelerationInWorld(t);

    for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
      int ee_xpp = ToXppEndeffector(n_ee, ee_towr).first;

      state.ee_contact_.at(ee_xpp) = solution.phase_durations_.at(ee_towr)->IsContactPhase(t);
      state.ee_motion_.at(ee_xpp)  = ToXpp(solution.ee_motion_.at(ee_towr)->GetPoint(t));
      state.ee_forces_ .at(ee_xpp) = solution.ee_force_.at(ee_towr)->GetPoint(t).p();
    }

    state.t_global_ = t;
    trajectory.push_back(state);
    t += visualization_dt_;
  }

  return trajectory;
}

xpp_msgs::RobotParameters
TowrRosInterface::BuildRobotParametersMsg(const RobotModel& model) const
{
  xpp_msgs::RobotParameters params_msg;
  auto max_dev_xyz = model.kinematic_model_->GetMaximumDeviationFromNominal();
  params_msg.ee_max_dev = xpp::Convert::ToRos<geometry_msgs::Vector3>(max_dev_xyz);

  auto nominal_B = model.kinematic_model_->GetNominalStanceInBase();
  int n_ee = nominal_B.size();
  for (int ee_towr=0; ee_towr<n_ee; ++ee_towr) {
    Vector3d pos = nominal_B.at(ee_towr);
    params_msg.nominal_ee_pos.push_back(xpp::Convert::ToRos<geometry_msgs::Point>(pos));
    params_msg.ee_names.push_back(ToXppEndeffector(n_ee, ee_towr).second);
  }

  params_msg.base_mass = model.dynamic_model_->m();

  return params_msg;
}

void
TowrRosInterface::SaveOptimizationAsRosbag (const std::string& bag_name,
                                   const xpp_msgs::RobotParameters& robot_params,
                                   const TowrCommandMsg user_command_msg,
                                   bool include_iterations)
{
  rosbag::Bag bag;
  bag.open(bag_name, rosbag::bagmode::Write);
  ::ros::Time t0(1e-6); // t=0.0 throws ROS exception

  // save the a-priori fixed optimization variables
  bag.write(xpp_msgs::robot_parameters, t0, robot_params);
  bag.write(towr_msgs::user_command+"_saved", t0, user_command_msg);

  // save the trajectory of each iteration
  if (include_iterations) {
    auto trajectories = GetIntermediateSolutions();
    int n_iterations = trajectories.size();
    for (int i=0; i<n_iterations; ++i)
      SaveTrajectoryInRosbag(bag, trajectories.at(i), towr_msgs::nlp_iterations_name + std::to_string(i));

    // save number of iterations the optimizer took
    std_msgs::Int32 m;
    m.data = n_iterations;
    bag.write(towr_msgs::nlp_iterations_count, t0, m);
  }

  // save the final trajectory
  auto final_trajectory = GetTrajectory();
  SaveTrajectoryInRosbag(bag, final_trajectory, xpp_msgs::robot_state_desired);

  bag.close();
}

void
TowrRosInterface::SaveTrajectoryInRosbag (rosbag::Bag& bag,
                                 const XppVec& traj,
                                 const std::string& topic) const
{
  for (const auto state : traj) {
    auto timestamp = ::ros::Time(state.t_global_ + 1e-6); // t=0.0 throws ROS exception

    xpp_msgs::RobotStateCartesian msg;
    msg = xpp::Convert::ToRos(state);
    bag.write(topic, timestamp, msg);

    xpp_msgs::TerrainInfo terrain_msg;
    for (auto ee : state.ee_motion_.ToImpl()) {
      Vector3d n = formulation_.terrain_->GetNormalizedBasis(HeightMap::Normal, ee.p_.x(), ee.p_.y());
      terrain_msg.surface_normals.push_back(xpp::Convert::ToRos<geometry_msgs::Vector3>(n));
      terrain_msg.friction_coeff = formulation_.terrain_->GetFrictionCoeff();
    }

    bag.write(xpp_msgs::terrain_info, timestamp, terrain_msg);
  }
}

void TowrRosInterface::ToRosAndPublish()
{

  auto base_angular=EulerConverter(solution.base_angular_);
  double speed_factor=0.5;
  double sampling_time=0.004*speed_factor;

  Eigen::MatrixXd oldFootPosDesWF(3,4);
  Eigen::MatrixXd oldFootVelDesWF(3,4);
  oldFootVelDesWF.setZero();
  ros::Rate loop_rate_(250.0);
  for(int ee=0; ee<solution.ee_motion_.size(); ee++){
  oldFootPosDesWF.col(ee) = formulation_.initial_ee_W_.at(ee);
  }
  
  for(int i=0; i<solution.base_linear_->GetTotalTime()/sampling_time; i++)
  {
    dwl_msgs::WholeBodyTrajectory planned_wt;
    double t=(double)i*sampling_time;
    dwl_msgs::WholeBodyState planned_wbs_msg;
    for(int ee=0; ee<solution.ee_motion_.size(); ee++)
    {
      dwl_msgs::ContactState contact;
      Eigen::Matrix3d w_R_b = base_angular.GetRotationMatrixBaseToWorld(t);
      //Eigen::Vector3d footPosDesCoM = w_R_b.transpose()*(solution.ee_motion_.at(ee)->GetPoint(t).p()-solution.base_linear_->GetPoint(t).p());
      Eigen::Vector3d footPosDesWF = solution.ee_motion_.at(ee)->GetPoint(t).p();
      contact.position.x = footPosDesWF(0) - solution.ee_motion_.at(ee)->GetPoint(0.0).p().x();
      contact.position.y = footPosDesWF(1) - solution.ee_motion_.at(ee)->GetPoint(0.0).p().y();
      contact.position.z = footPosDesWF(2) - solution.ee_motion_.at(ee)->GetPoint(0.0).p().z(); 
      //Eigen::Vector3d footVelDesCoM = w_R_b.transpose()*(solution.ee_motion_.at(ee)->GetPoint(t).v() - solution.base_linear_->GetPoint(t).v());
      Eigen::Vector3d footVelDesWF = solution.ee_motion_.at(ee)->GetPoint(t).v()*speed_factor;
      //Eigen::Vector3d footVelDesWF = (footPosDesWF - oldFootPosDesWF.col(ee))/sampling_time;
      contact.velocity.x = footVelDesWF(0);
      contact.velocity.y = footVelDesWF(1);
      contact.velocity.z = footVelDesWF(2);
      oldFootPosDesWF.col(ee) = footPosDesWF;

      //Eigen::Vector3d footAccDesCoM = w_R_b.transpose()*(solution.ee_motion_.at(ee)->GetPoint(t).a() - solution.base_linear_->GetPoint(t).a());
      Eigen::Vector3d footAccDesWF = solution.ee_motion_.at(ee)->GetPoint(t).a();
      //Eigen::Vector3d footAccDesWF = (footVelDesWF - oldFootVelDesWF.col(ee))/sampling_time;
      contact.acceleration.x = footAccDesWF(0)*pow(speed_factor,2);
      contact.acceleration.y = footAccDesWF(1)*pow(speed_factor,2);
      contact.acceleration.z = footAccDesWF(2)*pow(speed_factor,2);
      oldFootVelDesWF.col(ee) = footVelDesWF;
      switch(ee)
      {
       case 0: contact.name = "01_lf_foot"; break;
       case 1: contact.name = "02_rf_foot"; break;
       case 2: contact.name = "03_lh_foot"; break;
       case 3: contact.name = "04_rh_foot"; break;
      }
      if(solution.phase_durations_.at(ee)->IsContactPhase(t))
      {
       contact.wrench.force.x=100;
       contact.wrench.force.y=100;
       contact.wrench.force.z=100;
      }
      else
      {
        contact.wrench.force.x=0;
        contact.wrench.force.y=0;
        contact.wrench.force.z=0;
      }
      //if (i=0)
      //{std::cout<<footPosDesCoM.transpose()<<std::endl;}
    planned_wbs_msg.contacts.push_back(contact);
    }

    //unsigned int contact_counter = 1;

    //std::cout<<"i am here1"<<std::endl;
    auto pos=solution.base_angular_->GetPoint(t).p();
    auto vel=base_angular.GetAngularVelocityInWorld(t);
    auto accel=base_angular.GetAngularAccelerationInWorld(t);
    for(int base=0; base<3; base++)
    {
      dwl_msgs::BaseState base_state;
      base_state.name="floating_base";
      base_state.position =pos(base);
      base_state.velocity =vel(base)*speed_factor;
      base_state.acceleration=accel(base)*pow(speed_factor,2);
      switch(base)
      {
       case 0: base_state.id = base_state.AX; break;
       case 1: base_state.id = base_state.AY; break;
       case 2: base_state.id = base_state.AZ; break;
      }
      planned_wbs_msg.base.push_back(base_state);
    }
    Eigen::Vector3d pos_lin=solution.base_linear_->GetPoint(t).p()-solution.base_linear_->GetPoint(0.0).p();
    //std::cout << "t=" << t << "\n";
    //std::cout<<solution.base_linear_->GetPoint(t).p()<<std::endl;
    //std::cout<<" "<<std::endl;
    //std::cout<<solution.base_linear_->GetPoint(0).p()<<std::endl;
    //std::cout<<" "<<std::endl;
    //std::cout<<"sottrazione "<<std::endl;
    //std::cout<<solution.base_linear_->GetPoint(t).p()-solution.base_linear_->GetPoint(0).p()<<std::endl;
    //std::cout<<"pos_lin"<<std::endl;
    //std::cout<<pos_lin<<std::endl;

    auto vel_lin=solution.base_linear_->GetPoint(t).v();;
    auto acc_lin=solution.base_linear_->GetPoint(t).a();
    for(int base=0; base<3; base++)
    {
      dwl_msgs::BaseState base_state;
      base_state.name="floating_base";
      base_state.position =pos_lin(base);
      base_state.velocity =vel_lin(base)*speed_factor;
      base_state.acceleration=acc_lin(base)*pow(speed_factor,2);
      switch(base)
      {
       case 0: base_state.id = base_state.LX; break;
       case 1: base_state.id = base_state.LY; break;
       case 2: {
                base_state.id = base_state.LZ;
                //base_state.position += foot_radius_;
                break;}
      }
      planned_wbs_msg.base.push_back(base_state);

    }
    planned_wt.trajectory.push_back(planned_wbs_msg);

    dwltrajectory_.publish(planned_wt);
    loop_rate_.sleep();

  }
  
}

void TowrRosInterface::ToRos()
{

  time_t ct;
  ct = time(NULL);
  struct tm *localTime = localtime(&ct);
  char buffer [80];
  strftime(buffer, 80, "%d-%m-%Y-%H-%M-%S", localTime);
  std::string bag_file_name = "whole-body-trajectory-" + std::string(buffer)+ ".bag";
  rosbag::Bag bag;
  bag.open(bag_file_name, rosbag::bagmode::Write);
  const std::string wb_topic = "hyq/plan";

  //planned_wt.resize(solution.base_linear_->GetTotalTime()/0.04);
  auto base_angular=EulerConverter(solution.base_angular_);
  double speed_factor=0.5;
  double sampling_time=0.004*speed_factor;

  Eigen::MatrixXd oldFootPosDesWF(3,4);
  Eigen::MatrixXd oldFootVelDesWF(3,4);
  oldFootVelDesWF.setZero();
  for(int ee=0; ee<solution.ee_motion_.size(); ee++){
  oldFootPosDesWF.col(ee) = formulation_.initial_ee_W_.at(ee);
  }
  
  for(int i=0; i<solution.base_linear_->GetTotalTime()/sampling_time; i++)
  {

    dwl_msgs::WholeBodyTrajectory planned_wt;

    double t=(double)i*sampling_time;
    dwl_msgs::WholeBodyState planned_wbs_msg;
    for(int ee=0; ee<solution.ee_motion_.size(); ee++)
    {
      dwl_msgs::ContactState contact;
      Eigen::Matrix3d w_R_b = base_angular.GetRotationMatrixBaseToWorld(t);
      //Eigen::Vector3d footPosDesCoM = w_R_b.transpose()*(solution.ee_motion_.at(ee)->GetPoint(t).p()-solution.base_linear_->GetPoint(t).p());
      Eigen::Vector3d footPosDesWF = solution.ee_motion_.at(ee)->GetPoint(t).p();
      contact.position.x = footPosDesWF(0) - solution.ee_motion_.at(ee)->GetPoint(0.0).p().x();
      contact.position.y = footPosDesWF(1) - solution.ee_motion_.at(ee)->GetPoint(0.0).p().y();
      contact.position.z = footPosDesWF(2) - solution.ee_motion_.at(ee)->GetPoint(0.0).p().z(); 
      //Eigen::Vector3d footVelDesCoM = w_R_b.transpose()*(solution.ee_motion_.at(ee)->GetPoint(t).v() - solution.base_linear_->GetPoint(t).v());
      Eigen::Vector3d footVelDesWF = solution.ee_motion_.at(ee)->GetPoint(t).v()*speed_factor;
      //Eigen::Vector3d footVelDesWF = (footPosDesWF - oldFootPosDesWF.col(ee))/sampling_time;
      contact.velocity.x = footVelDesWF(0);
      contact.velocity.y = footVelDesWF(1);
      contact.velocity.z = footVelDesWF(2);
      oldFootPosDesWF.col(ee) = footPosDesWF;

      //Eigen::Vector3d footAccDesCoM = w_R_b.transpose()*(solution.ee_motion_.at(ee)->GetPoint(t).a() - solution.base_linear_->GetPoint(t).a());
      Eigen::Vector3d footAccDesWF = solution.ee_motion_.at(ee)->GetPoint(t).a();
      //Eigen::Vector3d footAccDesWF = (footVelDesWF - oldFootVelDesWF.col(ee))/sampling_time;
      contact.acceleration.x = footAccDesWF(0)*pow(speed_factor,2);
      contact.acceleration.y = footAccDesWF(1)*pow(speed_factor,2);
      contact.acceleration.z = footAccDesWF(2)*pow(speed_factor,2);
      oldFootVelDesWF.col(ee) = footVelDesWF;
      switch(ee)
      {
       case 0: contact.name = "01_lf_foot"; break;
       case 1: contact.name = "02_rf_foot"; break;
       case 2: contact.name = "03_lh_foot"; break;
       case 3: contact.name = "04_rh_foot"; break;
      }
      if(solution.phase_durations_.at(ee)->IsContactPhase(t))
      {
       contact.wrench.force.x=100;
       contact.wrench.force.y=100;
       contact.wrench.force.z=100;
      }
      else
      {
        contact.wrench.force.x=0;
        contact.wrench.force.y=0;
        contact.wrench.force.z=0;
      }
      //if (i=0)
      //{std::cout<<footPosDesCoM.transpose()<<std::endl;}
    planned_wbs_msg.contacts.push_back(contact);
    }

    //unsigned int contact_counter = 1;

    //std::cout<<"i am here1"<<std::endl;
    auto pos=solution.base_angular_->GetPoint(t).p();
    auto vel=base_angular.GetAngularVelocityInWorld(t);
    auto accel=base_angular.GetAngularAccelerationInWorld(t);
    for(int base=0; base<3; base++)
    {
      dwl_msgs::BaseState base_state;
      base_state.name="floating_base";
      base_state.position =pos(base);
      base_state.velocity =vel(base)*speed_factor;
      base_state.acceleration=accel(base)*pow(speed_factor,2);
      switch(base)
      {
       case 0: base_state.id = base_state.AX; break;
       case 1: base_state.id = base_state.AY; break;
       case 2: base_state.id = base_state.AZ; break;
      }
      planned_wbs_msg.base.push_back(base_state);
    }
    Eigen::Vector3d pos_lin=solution.base_linear_->GetPoint(t).p()-solution.base_linear_->GetPoint(0.0).p();
    //std::cout << "t=" << t << "\n";
    //std::cout<<solution.base_linear_->GetPoint(t).p()<<std::endl;
    //std::cout<<" "<<std::endl;
    //std::cout<<solution.base_linear_->GetPoint(0).p()<<std::endl;
    //std::cout<<" "<<std::endl;
    //std::cout<<"sottrazione "<<std::endl;
    //std::cout<<solution.base_linear_->GetPoint(t).p()-solution.base_linear_->GetPoint(0).p()<<std::endl;
    //std::cout<<"pos_lin"<<std::endl;
    //std::cout<<pos_lin<<std::endl;

    auto vel_lin=solution.base_linear_->GetPoint(t).v();;
    auto acc_lin=solution.base_linear_->GetPoint(t).a();
    for(int base=0; base<3; base++)
    {
      dwl_msgs::BaseState base_state;
      base_state.name="floating_base";
      base_state.position =pos_lin(base);
      base_state.velocity =vel_lin(base)*speed_factor;
      base_state.acceleration=acc_lin(base)*pow(speed_factor,2);
      switch(base)
      {
       case 0: base_state.id = base_state.LX; break;
       case 1: base_state.id = base_state.LY; break;
       case 2: {
                base_state.id = base_state.LZ;
                //base_state.position += foot_radius_;
                break;}
      }
      planned_wbs_msg.base.push_back(base_state);

    }
    planned_wt.trajectory.push_back(planned_wbs_msg);
    
    //timestamp = t;
    auto timestamp = ::ros::Time(t + 1e-6);
    bag.write(wb_topic, timestamp, planned_wt);
    }
    //double t = 0.0;
    //        while (t<=solution.base_linear_->GetTotalTime() + 1e-5)
    //         {
    //          std::cout << "t=" << t << "\n";
    //          std::cout<<solution.base_linear_->GetPoint(t).p()<<std::endl;
    //          std::cout<<" "<<std::endl;
    //          std::cout<<solution.base_linear_->GetPoint(0).p()<<std::endl;//
    //          std::cout << "LF ee_motion position x,y,z:   \t";
    //          //std::cout << solution.ee_motion_.at(0)->GetPoint(t).p().transpose() << "\t[m]" << std::endl;
    //          //Eigen::Matrix3d w_R_b = base_angular.GetRotationMatrixBaseToWorld(t);
    //          //Eigen::Vector3d basePosWF = solution.base_linear_->GetPoint(t).p()-w_R_b*offsetBF;
    //          //Eigen::Vector3d footPosDesBase = w_R_b.transpose()*(solution.ee_motion_.at(0)->GetPoint(t).p() - basePosWF);
    //          //std::cout << "LF ee_motion position x,y,z in base frame:   \t";
    //          //std::cout << footPosDesBase << "\t[m]" << std::endl;
    //          //std::cout <<"foot pos LF wf "<< solution.ee_motion_.at(0)->GetPoint(t).p() << std::endl;
    //          //std::cout << solution.ee_motion_.at(0)->GetPoint(t).p().transpose() << "\t[m]" << std::endl;
    //          t += 0.1;
    //}


    ////std::cout<<"i am here2"<<std::endl;
    bag.close();

    //return planned_wt;
  
}

  void
TowrRosInterface::SaveOptimizationAsRosbag (const std::string& bag_name,
                                   const xpp_msgs::RobotParameters& robot_params,
                                   //const TowrCommandMsg user_command_msg,
                                   bool include_iterations)
{
  rosbag::Bag bag;
  bag.open(bag_name, rosbag::bagmode::Write);
  ::ros::Time t0(1e-6); // t=0.0 throws ROS exception

  // save the a-priori fixed optimization variables
  bag.write(xpp_msgs::robot_parameters, t0, robot_params);
  //bag.write(towr_msgs::user_command+"_saved", t0, user_command_msg);

  // save the trajectory of each iteration
  if (include_iterations) {
    auto trajectories = GetIntermediateSolutions();
    int n_iterations = trajectories.size();
    for (int i=0; i<n_iterations; ++i)
      SaveTrajectoryInRosbag(bag, trajectories.at(i), towr_msgs::nlp_iterations_name + std::to_string(i));

    // save number of iterations the optimizer took
    std_msgs::Int32 m;
    m.data = n_iterations;
    bag.write(towr_msgs::nlp_iterations_count, t0, m);
  }

  // save the final trajectory
  auto final_trajectory = GetTrajectory();
  SaveTrajectoryInRosbag(bag, final_trajectory, xpp_msgs::robot_state_desired);

  bag.close();
}
} /* namespace towr */

