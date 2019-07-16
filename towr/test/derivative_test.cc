
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

TEST (DerivativeTest, JacobianWrtEE_motion)
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
    geom.ComputeCone(t);
    NodeSpline::Jacobian calculated(3,calculated.cols());
    calculated.setZero();
    int ee=0;
    if (geom.IsInTouch(t,ee))
    calculated=base.FillJacobianLinWrenchWrtEENodes(t,ee);
    NodeSpline::Jacobian expected(3,calculated.cols());
    for (int  i=0; i<3; i++)
    for (int j=0; j<expected.cols(); j++)

      {
        auto l=calculated.coeffRef(i,j)-expected.coeffRef(i,j);
        EXPECT_LE(l, 0.002);
        EXPECT_GE(l, -0.002);}
}

TEST (DerivativeTest, JacobianWrtAngEE)
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
    geom.ComputeCone(t);
    NodeSpline::Jacobian calculated(3,calculated.cols());
    calculated.setZero();
    int ee=0;
    if (geom.IsInTouch(t,ee))
    calculated=base.FillJacobianAngWrenchWrtEENodes(t,ee);
    NodeSpline::Jacobian expected(3,calculated.cols());
    expected.coeffRef(0,1)=-5.656854249;
    expected.coeffRef(1,0)=5.656854249;
    for (int  i=0; i<3; i++)
    for (int j=0; j<expected.cols(); j++)

      {
        auto l=calculated.coeffRef(i,j)-expected.coeffRef(i,j);
        EXPECT_LE(l, 0.002);
        EXPECT_GE(l, -0.002);}
}
TEST (DerivativeTest, JacobianWrtAngEESwingPhase)
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
    der.GetNumberOfFeetInTouch(0.2); //serve per riempire ee_motion_in_touch se no è vuoto!
    geom.ComputeCone(t);
    NodeSpline::Jacobian calculated(3,calculated.cols());
    calculated.setZero();
    int ee=0;
    if (geom.IsInTouch(t,ee))
    calculated=base.FillJacobianAngWrenchWrtEENodes(t,ee);

    NodeSpline::Jacobian expected(3,calculated.cols());
    for (int  i=0; i<3; i++)
    for (int j=0; j<expected.cols(); j++)

      {
        auto l=calculated.coeffRef(i,j)-expected.coeffRef(i,j);
        EXPECT_LE(l, 0.002);
        EXPECT_GE(l, -0.002);}
}
TEST (DerivativeTest, JacobianWrtLambda)
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
    geom.ComputeCone(t);
    auto calculated1=base.FillJacobianEdgesWrtLambda(t,2);
    Eigen::MatrixXd expected1(6,calculated1.cols());
    expected1.setZero();
    Eigen::MatrixXd l(6,4);
    l<<            1,            1,             -1,           -1,
                  -1,            1,              1,           -1,
             sqrt(2),       sqrt(2),        sqrt(2),      sqrt(2),
                   0,             0,            -0,            -0,
         -14.1944615,   -14.1944615,   -14.1944615,    -14.1944615,
             -10.037,        10.037,        10.037,       -10.037;
    expected1.block(0,16,6,4)=l/2.0;
    Eigen::MatrixXd expected=Eigen::MatrixXd(expected1);
    Eigen::MatrixXd calculated=Eigen::MatrixXd(calculated1);
    
    //std::cout<<calculated<<std::endl;
    //std::cout<<"uu "<<solution.lambda_->GetPoint(t).p()<<std::endl;

    EXPECT_EQ(calculated.cols(), 96);
    
    for (int  i=0; i<6; i++)  
    {
      for (int  j=0; j<calculated.cols(); j++)
      { //std::cout<<i<<" "<<j<<std::endl;
        double p=calculated(i,j)-expected(i,j);
        //std::cout<<calculated(i,j)<<" "<<expected(i,j)<<" "<<p<<std::endl;
        EXPECT_LE(p, 0.002);
        EXPECT_GE(p, -0.002);
      }
    }

}
TEST (DerivativeTest, UpdateJacobianAtInstance)
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
    geom.ComputeCone(t);
    ifopt::Component::Jacobian calculated1;
    int variables_num = 96;
    int constraints_num = 12*4;
    //ConstraintSet::GetRows();
    calculated1.resize(constraints_num,variables_num);
    calculated1.setZero();
    std::string vars_set = id::lambda_;
    std::cout<<vars_set<<std::endl;
    base.UpdateJacobianAtInstance(0.2, 2, vars_set, calculated1);
    Eigen::MatrixXd expected1(6,calculated1.cols());
    expected1.setZero();
    Eigen::MatrixXd l(6,4);
    l<<            1,            1,             -1,           -1,
                  -1,            1,              1,           -1,
             sqrt(2),       sqrt(2),        sqrt(2),      sqrt(2),
                   0,             0,            -0,            -0,
         -14.1944615,   -14.1944615,   -14.1944615,    -14.1944615,
             -10.037,        10.037,        10.037,       -10.037;
    expected1.block(0,16,6,4)=l/2.0;
    Eigen::MatrixXd expected=Eigen::MatrixXd(expected1);
    Eigen::MatrixXd calculated=Eigen::MatrixXd(calculated1);
    
    std::cout<<calculated<<std::endl;
    //std::cout<<"uu "<<solution.lambda_->GetPoint(t).p()<<std::endl;

    EXPECT_EQ(calculated.cols(), variables_num);
    
    for (int  i=0; i<6; i++)  
    {
      for (int  j=0; j<calculated.cols(); j++)
      { //std::cout<<i<<" "<<j<<std::endl;
        double p=calculated(i,j)-expected(i,j);
        //std::cout<<calculated(i,j)<<" "<<expected(i,j)<<" "<<p<<std::endl;
        EXPECT_LE(p, 0.002);
        EXPECT_GE(p, -0.002);
      }
    }

}
