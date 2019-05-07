/**
 * 07/11/2018
 *
 * Author: Romeo Orsolino
 *
 * email: rorsolino@ihmc.us
 *
 */

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

TEST(TOWR, optimizeTrajectory){
        NlpFormulation formulation;
	uint eeID = 0;

        //test to apply only the spline acc constraint (not costs) and swing constraint to a monoped robot which has to do only one jump

        formulation.terrain_ = std::make_shared<FlatGround>(0.0);
        formulation.initial_base_.lin.at(kPos) << 10.0, 0, 0.5;
        formulation.final_base_.lin.at(towr::kPos) << 10.2, 0, 0.5;
        formulation.final_base_.ang.at(towr::kPos) << 0, 0, 0;
        formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0, 0));

        // Kinematic limits and dynamic parameters of the hopper
        formulation.model_ = RobotModel(RobotModel::Monoped);
        
        const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
        //const std::vector<double> initial_durations = {0.3, 0.4, 0.4,0.4,0.4,0.4,0.4};
	bool is_first_phase_contact = true;
	double min_phase_duration = 0.1;
	double max_phase_duration = 2.0;

        formulation.params_.force_limit_in_normal_direction_=2000;

        //formulation.params_.dt_constraint_base_acc_=0.1; i have added it in the parameter class
        formulation.params_.ee_phase_durations_.push_back(initial_durations);
        formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
        PhaseDurations phase_durations(eeID,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

          ifopt::Problem nlp;
          SplineHolder solution;
          for (auto c : formulation.GetVariableSets(solution))
            nlp.AddVariableSet(c);
          ContraintPtrVec constraints;

          auto ee_motion_name =  "ee-motion_" + std::to_string(formulation.params_.GetEECount()-1);
          int ee_count = formulation.params_.GetEECount()-1;

          double tot_time = formulation.params_.GetTotalTime();
          formulation.params_.OptimizePhaseDurations();
          //constraints.push_back(std::make_shared<TotalDurationConstraint>(tot_time, ee_count));
          constraints.push_back(std::make_shared<TerrainConstraint>(formulation.terrain_, ee_motion_name));

          //swing_constraint
          constraints.push_back(std::make_shared<SwingConstraint>(id::EEMotionNodes(formulation.params_.GetEECount()-1)));

          //force_constraint

          constraints.push_back(std::make_shared<ForceConstraint>(formulation.terrain_,
                                                                  formulation.params_.force_limit_in_normal_direction_,
                                                                  formulation.params_.GetEECount()-1));

          //spline_acc_constraint

          constraints.push_back(std::make_shared<SplineAccConstraint>
                                (solution.base_linear_, id::base_lin_nodes)) ;

          constraints.push_back(std::make_shared<SplineAccConstraint>
                                (solution.base_angular_, id::base_ang_nodes)) ;

          constraints.push_back(std::make_shared<DynamicConstraint>(formulation.model_.dynamic_model_,
                                                                  formulation.params_.GetTotalTime(),
                                                                  formulation.params_.dt_constraint_dynamic_,
                                                                  solution));
          // base_acc_range_constraint


          constraints.push_back(std::make_shared<BaseAccConstraintRangeLin>(formulation.model_.dynamic_model_,
                                                                            formulation.params_.GetTotalTime(),
                                                                            formulation.params_.dt_constraint_base_acc_,
                                                                            solution.base_linear_, id::base_lin_nodes)) ;

          //constraints.push_back(std::make_shared<BaseAccConstraintRangeAng>(formulation.model_.dynamic_model_,
          //                                                                  formulation.params_.GetTotalTime(),
          //                                                                  formulation.params_.dt_constraint_base_acc_,
          //                                                                  solution.base_angular_, solution.base_linear_,
          //                                                                  id::base_ang_nodes)) ;
          //
          //constraints.push_back(std::make_shared<BaseMotionConstraint>(formulation.params_.GetTotalTime(),
          //                                                              formulation.params_.GetTotalTime(),
          //                                                            formulation.params_.dt_constraint_base_motion_,
          //                                                            solution));
          //BASE MOTION SEMPRE COMMENTATO!

          constraints.push_back(std::make_shared<RangeOfMotionConstraint>(formulation.model_.kinematic_model_,
                                                                  formulation.params_.GetTotalTime(),
                                                                  formulation.params_.dt_constraint_range_of_motion_,
                                                                  ee_count,
                                                                  solution));




          //devo lanciarlo per tutti gli elementi di constraints!

          for (auto l:constraints)
            nlp.AddConstraintSet(l);
          for (auto c : formulation.GetCosts())
          nlp.AddCostSet(c);

          auto solver = std::make_shared<ifopt::IpoptSolver>();

          solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
          solver->SetOption("max_cpu_time", 20.0);


          solver->Solve(nlp);



         //I copied this part from hopper_example

            double t = 0.0;
            while (t<=solution.base_linear_->GetTotalTime() + 1e-5) {
            cout << "t=" << t << "\n";
            cout << "Base linear position x,y,z:   \t";
            cout << solution.base_linear_->GetPoint(t).p().transpose() << "\t[m]" << endl;
            cout << "Base linear vel x,y,z:   \t";
            cout << solution.base_linear_->GetPoint(t).v().transpose() << "\t[m/s]" << endl;
            cout << "Base linear acc x,y,z:   \t";
            cout << solution.base_linear_->GetPoint(t).a().transpose() << "\t[m/s^2]" << endl;

            cout << "Base Euler roll, pitch, yaw:  \t";
            Eigen::Vector3d rad = solution.base_angular_->GetPoint(t).p();
            cout << (rad/M_PI*180).transpose() << "\t[deg]" << endl;
            cout << "Base Euler roll, pitch, yaw vel:  \t";
            Eigen::Vector3d ang_vel = solution.base_angular_->GetPoint(t).v();
            cout<< (ang_vel/M_PI*180).transpose() << "\t[deg/s]" << endl;

            cout << "Foot position x,y,z:          \t";
            cout << solution.ee_motion_.at(0)->GetPoint(t).p().transpose() << "\t[m]" << endl;
            cout << "Foot velocity x,y,z:          \t";
            cout << solution.ee_motion_.at(0)->GetPoint(t).v().transpose() << "\t[m/s]" << endl;

            cout << "Contact force x,y,z:          \t";
            cout << solution.ee_force_.at(0)->GetPoint(t).p().transpose() << "\t[N]" << endl;

            cout << "Contact force dertivative x,y,z:          \t";
            cout << solution.ee_force_.at(0)->GetPoint(t).v().transpose() << "\t[N/s]" << endl;


            bool contact = solution.phase_durations_.at(0)->IsContactPhase(t);
            std::string foot_in_contact = contact? "yes" : "no";
            cout << "Foot in contact:              \t" + foot_in_contact << endl;

            t += 0.1;
            }
}
