
#include <gtest/gtest.h>
#include <iostream>
#include <cmath>

#include <towr/nlp_formulation.h>
#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/constraints/spline_acc_constraint.h>
#include <towr/terrain/examples/height_map_examples.h>
#include <towr/nlp_formulation.h>
#include <towr/variables/variable_names.h>
#include <towr/constraints/force_constraint.h>
#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/base_acc_constraint_range_lin.h>
#include <towr/constraints/base_acc_constraint_range_ang.h>
#include <towr/constraints/base_motion_constraint.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/constraints/force_constraint.h>
#include <towr/constraints/range_of_motion_constraint.h>
#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/terrain_constraint.h>
#include <towr/constraints/total_duration_constraint.h>
#include <ifopt/ipopt_solver.h>
#include <ifopt/composite.h>


using namespace towr;
using namespace std;
using namespace ifopt;
using VariablePtrVec   = std::vector<ifopt::VariableSet::Ptr>;
using ContraintPtrVec  = std::vector<ifopt::ConstraintSet::Ptr>;
using CostPtrVec       = std::vector<ifopt::CostTerm::Ptr>;
using EEPos            = std::vector<Eigen::Vector3d>;
using Vector3d         = Eigen::Vector3d;
using Ptr              = std::shared_ptr<IpoptSolver>;

/* The four following tests verify that the values of motion phase nodes and force phase nodes
* are actually splined from the inital desired value to the final desired value using the InitializeNodesTowardsGoal method */
TEST(ComputeWrenchTest, WrenchOftheRobot)
{
  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  //const std::vector<double> initial_durations = {0.3, 0.4, 0.4,0.4,0.4, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    ContraintPtrVec constraints;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.5;
    auto com=solution.base_linear_->GetPoint(t);
    //auto rp=com.p();
    //auto ra=com.a();
    //auto acc=solution.base_angular_->GetPoint(t).a();
    //auto w_R_b=EulerConverter(solution.base_angular_).GetRotationMatrixBaseToWorld(t);
    //auto alfa=EulerConverter(solution.base_angular_).GetAngularAccelerationInWorld(t);
    //Eigen::Matrix3d I_b;
    //I_b<<1.2, 0, 0.2, 0, 5.5, 0.01, 0.2, 0.01, 6;
    //auto I_w=w_R_b*I_b*w_R_b.transpose();


    Eigen::VectorXd expected(6);
    expected<<5.80092e-12,0,196.133,1.332264e-14,-1979.15849,2.1633414e-15;

    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);

    Derivative der(formulation.terrain_, solution, 1, formulation.params_.ee_phase_durations_);

    BaseAccConstraintRangeLin base (formulation.model_.dynamic_model_,
                                    formulation.params_.GetTotalTime(),
                                    formulation.params_.dt_constraint_base_acc_,
                                    solution.base_linear_, id::base_lin_nodes,
                                    formulation.terrain_,
                                    solution,
                                    formulation.params_.GetEECount(),
                                    formulation.params_.ee_phase_durations_,
                                    geom,
                                    der);
   auto calculated=base.ComputeWrench(com,t);
   auto b=expected-calculated;
   for (int i=0; i <6; i++)
   {EXPECT_LE(b(i), 0.0018);
    EXPECT_GE(b(i), -0.0018);} //the fifth value is -1979,15849, so an error of 0.0018 is acceptable!
}

TEST(ComputeWrenchTest, EdgesWithFlatGround)
{
  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.2; double ee=0;
    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);
    auto normal=geom.ReturnNormalTerrain(t);
    normal.resize(3,1);
    Eigen::Vector3d axis1;
    axis1<< normal(1), -normal(0), 0.0;
    Eigen::VectorXd axis=1/normal.norm()*axis1;
    double angle=geom.ComputeRotationAngle(normal);
    Eigen::MatrixXd expected(4,3);
    expected<< 0.5, -0.5, sqrt(2)/2,
           0.5,  0.5, sqrt(2)/2,
          -0.5,  0.5, sqrt(2)/2,
          -0.5, -0.5, sqrt(2)/2;

    auto calculated=geom.ComputeLinearPartOfTheCone(axis,angle);
    for (int i=0; i <4; i++)
    for (int j=0; j<3; j++)
      {
        auto l=calculated(i,j)-expected(i,j);
        EXPECT_LE(l, 0.0018);
        EXPECT_GE(l, -0.0018);}
}
TEST(ComputeWrenchTest, EdgesWithNoFlatGround)
{
  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);
    Eigen::Vector3d normal;
    normal<<0.5, 0.83, sqrt(2)/6;
    Eigen::Vector3d axis1;
    axis1<< normal(1), -normal(0), 0.0;
    Eigen::VectorXd axis=1/axis1.norm()*axis1;
    double angle=geom.ComputeRotationAngle(normal);
    Eigen::MatrixXd expected(4,3);
    expected<<  0.2140, -0.9771,  -3.7728e-5,
               -0.1233, -0.5390,   0.8335,
               -0.9213, -0.2017,   0.334,
               -0.5839, -0.6398,  -0.5002;

    auto calculated=geom.ComputeLinearPartOfTheCone(axis,angle);
    for (int i=0; i<4; i++)
    for (int j=0; j<3; j++)
      {
        auto l=calculated(i,j)-expected(i,j);
        EXPECT_LE(l, 0.0018);
        EXPECT_GE(l, -0.0018);
      }

}

TEST(ComputeWrenchTest, EdgesInSwingPhase)
{
  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.5; double ee=0;
    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);
    auto calculated=geom.ComputeCone(t);
    Eigen::MatrixXd expected(4,3);
    expected.setZero();
    for (int i=0; i <4; i++)
    for (int j=0; j<3; j++)
      {
        auto l=calculated(i,j)-expected(i,j);
        EXPECT_LE(l, 0.0018);
        EXPECT_GE(l, -0.0018);}
}

TEST (ComputeWrenchTest, TotalEdges)
{

  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.2;

    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);
    Eigen::MatrixXd calculated(6,4);
    calculated=geom.ComputeCone(t);
    Eigen::MatrixXd expected;
    expected.resize(4,6);
    expected.setZero();
    expected<< 0.5, -0.5, sqrt(2)/2,  0, -7.097230763, -5.0185,
               0.5,  0.5, sqrt(2)/2,  0, -7.097230763,  5.0185,
              -0.5,  0.5, sqrt(2)/2,  0, -7.097230763,  5.0185,
              -0.5, -0.5, sqrt(2)/2,  0, -7.097230763, -5.0185;
    Eigen::MatrixXd expected1=expected.transpose();
    for (int i=0; i<6; i++)
    for (int j=0; j<4; j++)
      {
        auto l=calculated(i,j)-expected1(i,j);
        EXPECT_LE(l, 0.0018);
        EXPECT_GE(l, -0.0018);}
}
TEST (ComputeWrenchTest, TotalEdgesInSwingPhase)
{

  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.5;

    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);
    Eigen::MatrixXd calculated(6,4);
    calculated=geom.ComputeCone(t);
    Eigen::MatrixXd expected;
    expected.resize(6,4);
    expected.setZero();
    for (int i=0; i<6; i++)
    for (int j=0; j<4; j++)
      {
        auto l=expected(i,j)-calculated(i,j);
        EXPECT_LE(l, 0.0018);
        EXPECT_GE(l, -0.0018);}
}

TEST (ComputeWrenchTest, Constraint)
{

  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.2;
    auto com=solution.base_linear_->GetPoint(t);
    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);

    Derivative der(formulation.terrain_, solution, 1, formulation.params_.ee_phase_durations_);
    BaseAccConstraintRangeLin base (formulation.model_.dynamic_model_,
                                    formulation.params_.GetTotalTime(),
                                    formulation.params_.dt_constraint_base_acc_,
                                    solution.base_linear_, id::base_lin_nodes,
                                    formulation.terrain_,
                                    solution,
                                    formulation.params_.GetEECount(),
                                    formulation.params_.ee_phase_durations_,
                                    geom,
                                    der);

    auto calculated=base.FillConstraint(com,t);
    Eigen::VectorXd expected(6);
    expected<<-1.55154e-11, 0, 190.476145751, 0, -1911.682153896, 0;
    for (int i=0; i<6; i++)

      {
        auto l=calculated(i)-expected(i);
        EXPECT_LE(l, 0.002);
        EXPECT_GE(l, -0.002);}
}

TEST (ComputeWrenchTest, ConstraintInSwingPhase)
{

  NlpFormulation formulation;
  formulation.terrain_ = std::make_shared<FlatGround>(0.0);
  formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
  formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
  formulation.final_base_.ang.at(towr::kPos) << 0.5, 0, 0;
  formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

  // Kinematic limits and dynamic parameters of the hopper
  formulation.model_ = RobotModel(RobotModel::Monoped);

  const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
  bool is_first_phase_contact = true;
  double min_phase_duration = 0.1;
  double max_phase_duration = 2.0;

  formulation.params_.force_limit_in_normal_direction_=2000;

   //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
   formulation.params_.ee_phase_durations_.push_back(initial_durations);
   formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
   PhaseDurations phase_durations(0,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

    ifopt::Problem nlp;
    SplineHolder solution;
    for (auto c : formulation.GetVariableSets(solution))
    nlp.AddVariableSet(c);
    double t=0.5;
    auto com=solution.base_linear_->GetPoint(t);
    Geometry geom(formulation.model_.dynamic_model_,
                   formulation.terrain_,
                   solution,
                   1,
                 formulation.params_.ee_phase_durations_);

    Derivative der(formulation.terrain_, solution, 1, formulation.params_.ee_phase_durations_);
    BaseAccConstraintRangeLin base (formulation.model_.dynamic_model_,
                                    formulation.params_.GetTotalTime(),
                                    formulation.params_.dt_constraint_base_acc_,
                                    solution.base_linear_, id::base_lin_nodes,
                                    formulation.terrain_,
                                    solution,
                                    formulation.params_.GetEECount(),
                                    formulation.params_.ee_phase_durations_,
                                    geom,
                                    der);

    auto calculated=base.FillConstraint(com,t);
    Eigen::VectorXd expected(6);
    expected<<5.80092e-12,0,196.133,1.332264e-14,-1979.15849,2.1633414e-15;//equal to the value of the first test
    for (int i=0; i<6; i++)

      {
        auto l=calculated(i)-expected(i);
        EXPECT_LE(l, 0.002);
        EXPECT_GE(l, -0.002);}
}
