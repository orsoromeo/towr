/**
 * 07/11/2018
 *
 * Author: Romeo Orsolino
 *
 * email: rorsolino@ihmc.us
 *
 */

#include <gtest/gtest.h>

#include <cmath>
#include <towr/nlp_formulation.h>
#include <towr/models/centroidal_model.h>
#include <ifopt/ipopt.h>
#include <ifopt/composite.h>
using namespace towr;
using VariablePtrVec   = std::vector<ifopt::VariableSet::Ptr>;
using ContraintPtrVec  = std::vector<ifopt::ConstraintSet::Ptr>;
using CostPtrVec       = std::vector<ifopt::CostTerm::Ptr>;
using EEPos            = std::vector<Eigen::Vector3d>;
using Vector3d         = Eigen::Vector3d;

TEST(TOWR, optimizeTrajectory){

	uint eeID = 0;
	const std::vector<double> initial_durations = {0.3, 0.4, 0.4};
	bool is_first_phase_contact = true;
	double min_phase_duration = 0.1;
	double max_phase_duration = 2.0;
	ifopt::Composite composite("composite", false);
	//PhaseDurations phase_durations(eeID,
	//							initial_durations,
	//							is_first_phase_contact,
	//							min_phase_duration,
	//							max_phase_duration);
    //
	//  VariablePtrVec vars;
	//  SplineHolder spline_holder_;
	//  auto base_motion = MakeBaseVariables();
	//  vars.insert(vars.end(), base_motion.begin(), base_motion.end());
    //
	//  auto ee_motion = MakeEndeffectorVariables();
	//  vars.insert(vars.end(), ee_motion.begin(), ee_motion.end());
    //
	//  auto ee_force = MakeForceVariables();
	//  vars.insert(vars.end(), ee_force.begin(), ee_force.end());
    //
	//  auto contact_schedule = MakeContactScheduleVariables();
	//  if (params_.OptimizeTimings()) {
	//    vars.insert(vars.end(), contact_schedule.begin(), contact_schedule.end());
	//  }
	//  spline_holder_ = SplineHolder(base_motion.at(0), // linear
	//                                base_motion.at(1), // angular
	//                                params_.GetBasePolyDurations(),
	//                                ee_motion,
	//                                ee_force,
	//                                contact_schedule,
	//                                params_.OptimizeTimings());
}

