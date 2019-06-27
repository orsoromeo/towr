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



namespace towr {


TowrRosInterface::TowrRosInterface ()
{
  ::ros::NodeHandle n;

  user_command_sub_ = n.subscribe(towr_msgs::user_command, 1,
                                  &TowrRosInterface::UserCommandCallback, this);

  initial_state_pub_  = n.advertise<xpp_msgs::RobotStateCartesian>
                                          (xpp_msgs::robot_state_desired, 1);

  robot_parameters_pub_  = n.advertise<xpp_msgs::RobotParameters>
                                    (xpp_msgs::robot_parameters, 1);

  trajectory_ = n.advertise<xpp_msgs::RobotStateCartesianTrajectory>("/xpp/trajectory_des",1);

  dwltrajectory_=n.advertise<dwl_msgs::WholeBodyTrajectory>("hyq/plan",1);

  solver_ = std::make_shared<ifopt::IpoptSolver>();

  visualization_dt_ = 0.01;
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

  SetTowrInitialState();

  // solver parameters
  SetIpoptParameters(msg);

  // visualization
  PublishInitialState();

  // Defaults to /home/user/.ros/
  std::string bag_file = "towr_trajectory.bag";
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

  dwl_msgs::WholeBodyTrajectory wbtraj = ToRos();

  trajectory_.publish(xpp_msg);
  
  dwltrajectory_.publish(wbtraj);
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

dwl_msgs::WholeBodyTrajectory TowrRosInterface::ToRos()
{
  dwl_msgs::WholeBodyTrajectory planned_wt;
  //planned_wt.resize(solution.base_linear_->GetTotalTime()/0.04);
  auto base_angular=EulerConverter(solution.base_angular_);
  for(int i=0; i<solution.base_linear_->GetTotalTime()/0.04; i++)
  {

    double t=i*0.04;

    dwl_msgs::WholeBodyState planned_wbs_msg;
    for(int ee=0; ee<solution.ee_motion_.size(); ee++)
    {
      dwl_msgs::ContactState contact;
      contact.position.x = solution.ee_motion_.at(ee)->GetPoint(t).p().x();
      contact.position.y = solution.ee_motion_.at(ee)->GetPoint(t).p().y();
      contact.position.z = solution.ee_motion_.at(ee)->GetPoint(t).p().z();
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
      base_state.velocity =vel(base);
      base_state.acceleration=accel(base);
      switch(base)
      {
       case 0: base_state.id = base_state.AX; break;
       case 1: base_state.id = base_state.AY; break;
       case 2: base_state.id = base_state.AZ; break;
      }
      planned_wbs_msg.base.push_back(base_state);
    }
    auto pos_lin=solution.base_linear_->GetPoint(t).p();
    auto vel_lin=solution.base_linear_->GetPoint(t).v();;
    auto acc_lin=solution.base_linear_->GetPoint(t).a();
    for(int base=0; base<3; base++)
    {
      dwl_msgs::BaseState base_state;
      base_state.name="floating_base";
      base_state.position =pos_lin(base);
      base_state.velocity =vel_lin(base);
      base_state.acceleration=acc_lin(base);
      switch(base)
      {
       case 0: base_state.id = base_state.LX; break;
       case 1: base_state.id = base_state.LY; break;
       case 2: base_state.id = base_state.LZ; break;

      }
      planned_wbs_msg.base.push_back(base_state);
    }
    planned_wt.trajectory.push_back(planned_wbs_msg);}
    std::cout<<"i am here2"<<std::endl;
    return planned_wt;
    }
    //dwl_msgs::BaseState base_state_0, base_state_1, base_state_2, base_state_3, base_state_4, base_state_5;
    //
    //base_state_0.name="floating_base";
    //base_state_0.id = base_state_0.AX;
    //base_state_0.position = solution.base_angular_->GetPoint(t).p().x();
    //base_state_0.velocity = base_angular.GetAngularVelocityInWorld(t).x();
    //base_state_0.acceleration = base_angular.GetAngularAccelerationInWorld(t).x();
    //
    //base_state_1.name="floating_base";
    //base_state_1.id = base_state_1.AY;
    //base_state_1.position = solution.base_angular_->GetPoint(t).p().y();
    //base_state_1.velocity = base_angular.GetAngularVelocityInWorld(t).y();
    //base_state_1.acceleration = base_angular.GetAngularAccelerationInWorld(t).y();
    //
    //base_state_2.name="floating_base";
    //base_state_2.id = base_state_2.AZ;
    //base_state_2.position = solution.base_angular_->GetPoint(t).p().z();
    //base_state_2.velocity = base_angular.GetAngularVelocityInWorld(t).z();
    //base_state_2.acceleration = base_angular.GetAngularAccelerationInWorld(t).z();
    //
    //base_state_3.name="floating_base";
    //base_state_3.id = base_state_3.LX;
    //base_state_3.position = solution.base_linear_->GetPoint(t).p().x();
    //base_state_3.velocity = solution.base_linear_->GetPoint(t).v().x();
    //base_state_3.acceleration = solution.base_linear_->GetPoint(t).a().x();
    //
    //base_state_4.name="floating_base";
    //base_state_4.id = base_state_4.LY;
    //base_state_4.position = solution.base_linear_->GetPoint(t).p().y();
    //base_state_4.velocity = solution.base_linear_->GetPoint(t).v().y();
    //base_state_4.acceleration = solution.base_linear_->GetPoint(t).a().y();
    //
    //base_state_5.name="floating_base";
    //base_state_5.id = base_state_5.LZ;
    //base_state_5.position = solution.base_linear_->GetPoint(t).p().z();
    //base_state_5.velocity = solution.base_linear_->GetPoint(t).v().z();
    //base_state_5.acceleration = solution.base_linear_->GetPoint(t).a().z();
    //
    //planned_wbs_msg.base.push_back(base_state_0);
    //planned_wbs_msg.base.push_back(base_state_1);
    //planned_wbs_msg.base.push_back(base_state_2);
    //planned_wbs_msg.base.push_back(base_state_3);
    //planned_wbs_msg.base.push_back(base_state_4);
    //planned_wbs_msg.base.push_back(base_state_5);
    //for (unsigned int i = 0; i < fbs_->getJointDoF(); i++) {
    //for (unsigned int i = 0; i < 12; i++) {
    //        // Getting the joint id
    //        unsigned int joint_id = getDWLJointId(JointIdentifiers(i));
    //        // Converting the actual whole-body states
    //        planned_ws_.setJointPosition(des_q_(i), joint_id);
    //        planned_ws_.setJointVelocity(des_qd_(i), joint_id);
    //        planned_ws_.setJointAcceleration(des_qdd_(i), joint_id);
    //        //no torques are sent
    //        planned_ws_.setJointEffort(0., joint_id);
    //    }
    //
    //  planned_ws_.setBaseRPY(solution.base_angular_->GetPoint(t).p());
    //  planned_ws_.setBaseRPYVelocity_W(base_angular.GetAngularAccelerationInWorld(t));
    //  planned_ws_.setBaseRPYAcceleration_W(base_angular.GetAngularAccelerationInWorld(t));
    //
    //
    //  planned_ws_.setBasePosition(solution.base_linear_->GetPoint(t).p());
    //  planned_ws_.setBaseVelocity_W(solution.base_linear_->GetPoint(t).v());
    //  planned_ws_.setBaseAcceleration_W(solution.base_linear_->GetPoint(t).a());//TODO accel
    //
    //  //set the desired foot positions in the planned trajectory
    //  planned_ws_.setContactPosition_B("01_lf_foot", solution.ee_motion_.at(0)->GetPoint(t).p());
    //  planned_ws_.setContactPosition_B("02_rf_foot", solution.ee_motion_.at(1)->GetPoint(t).p());
    //  planned_ws_.setContactPosition_B("03_lh_foot", solution.ee_motion_.at(2)->GetPoint(t).p());
    //  planned_ws_.setContactPosition_B("04_rh_foot", solution.ee_motion_.at(3)->GetPoint(t).p());
    //
    //  //send the desired foot velocities in the planned trajectory
    //  planned_ws_.setContactVelocity_B("01_lf_foot", solution.ee_motion_.at(ee_towr)->GetPoint(t).v());
    //  planned_ws_.setContactVelocity_B("02_rf_foot", solution.ee_motion_.at(ee_towr)->GetPoint(t).v());
    //  planned_ws_.setContactVelocity_B("03_lh_foot", solution.ee_motion_.at(ee_towr)->GetPoint(t).v());
    //  planned_ws_.setContactVelocity_B("04_rh_foot", solution.ee_motion_.at(ee_towr)->GetPoint(t).v());
    //
    //  //send the desired stance legs in the planned trajectory
    //  planned_ws_.setContactCondition("01_lf_foot",gl.stance_legs[LF]);
    //  planned_ws_.setContactCondition("02_rf_foot",gl.stance_legs[RF]);
    //  planned_ws_.setContactCondition("03_lh_foot",gl.stance_legs[LH]);
    //  planned_ws_.setContactCondition("04_rh_foot",gl.stance_legs[RH]);
    //
    //planned_wt.resize(1);
    //planned_wt.trajectory.push_back(planned_wbs_msg);}
  //std::cout<<"i am here2"<<std::endl;
  //return planned_wt;
//}

} /* namespace towr */

