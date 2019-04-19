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

        //test to apply only the sline acc constraint (not costs) to a monoped robot which has to do only one jump

        formulation.terrain_ = std::make_shared<FlatGround>(0.0);
        formulation.initial_base_.lin.at(kPos) << 10.0, 0.0, 0.5;
        formulation.final_base_.lin.at(towr::kPos) << 11.0, 0.0, 0.5;
        formulation.initial_ee_W_.push_back(Eigen::Vector3d(10.0, 0.0, 0.0));

        // Kinematic limits and dynamic parameters of the hopper
        formulation.model_ = RobotModel(RobotModel::Monoped);

        const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
	bool is_first_phase_contact = true;
	double min_phase_duration = 0.1;
	double max_phase_duration = 2.0;
	ifopt::Composite composite("composite", false);

        formulation.params_.OptimizePhaseDurations();
        formulation.params_.ee_phase_durations_.push_back(initial_durations);
        formulation.params_.ee_in_contact_at_start_.push_back(is_first_phase_contact);
        PhaseDurations phase_durations(eeID,
                                                                initial_durations,
                                                                is_first_phase_contact,
                                                                min_phase_duration,
                                                                max_phase_duration);

          //VariablePtrVec vars;
          //SplineHolder spline_holder_;
          // initial conditions
          //auto base_motion = MakeBaseVariables();
          //NlpFormulation base_motion;
          //base_motion.MakeBaseVariables();
          //vars.insert(vars.end(), base_motion.begin(), base_motion.end());
          //
          //auto ee_motion = MakeEndeffectorVariables();
          //vars.insert(vars.end(), ee_motion.begin(), ee_motion.end());
          //
          //auto ee_force = MakeForceVariables();
          //vars.insert(vars.end(), ee_force.begin(), ee_force.end());
          //
          //auto contact_schedule = MakeContactScheduleVariables();
          //if (params_.OptimizeTimings()) {
          //  vars.insert(vars.end(), contact_schedule.begin(), contact_schedule.end());
          //}
          // spline_holder_ = SplineHolder(base_motion.at(0), // linear
          //                               base_motion.at(1), // angular
          //                               params_.GetBasePolyDurations(),
          //                               ee_motion,
          //                               ee_force,
          //                               contact_schedule,
          //                             params_.OptimizeTimings());
          ifopt::Problem nlp;
          SplineHolder solution;
          for (auto c : formulation.GetVariableSets(solution))
            nlp.AddVariableSet(c);
          ContraintPtrVec constraints;

          constraints.push_back(std::make_shared<SplineAccConstraint>
                                (solution.base_linear_, id::base_lin_nodes));

          constraints.push_back(std::make_shared<SplineAccConstraint>
                                (solution.base_angular_, id::base_ang_nodes));

          //devo lanciarlo per tutti gli elementi di constraints!
          for (auto l:constraints)
          nlp.AddConstraintSet(l);
          //SplineAccConstraint AccConstraint_;
          //AccConstraint_.FillJacobianBlock("splineacc",jac);
          //AccConstraint_.GetBounds ();
          auto solver = std::make_shared<ifopt::IpoptSolver>();
          //ifopt::IpoptSolver::Ptr solver;
          solver->SetOption("jacobian_approximation", "exact"); // "finite difference-values"
          solver->SetOption("max_cpu_time", 20.0);
          //nlp.AddVariableSet(vars);
          //nlp.AddConstraintSet(AccConstraint_);
          solver->Solve(nlp);

         //I copied this part from hopper_example

            double t = 0.0;
            while (t<=solution.base_linear_->GetTotalTime() + 1e-5) {
            cout << "t=" << t << "\n";
            cout << "Base linear position x,y,z:   \t";
            cout << solution.base_linear_->GetPoint(t).p().transpose() << "\t[m]" << endl;
            cout << "Base linear vel x,y,z:   \t";
            cout << solution.base_linear_->GetPoint(t).v().transpose() << "\t[m/s]" << endl;
            
            cout << "Base Euler roll, pitch, yaw:  \t";
            Eigen::Vector3d rad = solution.base_angular_->GetPoint(t).p();
            cout << (rad/M_PI*180).transpose() << "\t[deg]" << endl;
            cout << "Base Euler roll, pitch, yaw vel:  \t";
            Eigen::Vector3d ang_vel = solution.base_angular_->GetPoint(t).v();
            cout<< (ang_vel/M_PI*180).transpose() << "\t[deg/s]" << endl;

            cout << "Foot position x,y,z:          \t";
            cout << solution.ee_motion_.at(0)->GetPoint(t).p().transpose() << "\t[m]" << endl;

            cout << "Contact force x,y,z:          \t";
            cout << solution.ee_force_.at(0)->GetPoint(t).p().transpose() << "\t[N]" << endl;

            bool contact = solution.phase_durations_.at(0)->IsContactPhase(t);
            std::string foot_in_contact = contact? "yes" : "no";
            cout << "Foot in contact:              \t" + foot_in_contact << endl;

            //cout << "Phase durations:              \t" << endl;
            //for (int j = 0; j<3; j++){
            //cout<<solution.phase_durations_[j];
            //cout << endl;
            //} 
            t += 0.05;
            }
}
