/**
 * 17/11/2018
 *
 * Author: Romeo Orsolino
 *
 * email: rorsolino@ihmc.us
 *
 */
#include <iostream>
#include <gtest/gtest.h>

#include <towr/nlp_factory.h>

#include <towr/variables/variable_names.h>
#include <towr/variables/base_nodes.h>
#include <towr/variables/phase_durations.h>

#include <towr/constraints/base_motion_constraint.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/constraints/force_constraint.h>
#include <towr/constraints/range_of_motion_constraint.h>
#include <towr/constraints/swing_constraint.h>
#include <towr/constraints/terrain_constraint.h>
#include <towr/constraints/total_duration_constraint.h>
#include <towr/constraints/spline_acc_constraint.h>

#include <towr/models/examples/monoped_model.h>
#include <towr/models/centroidal_model.h>

#include <towr/costs/node_cost.h>

namespace towr {

TEST(PhaseNodesTest, motionPhaseNodesInitializationTowardsGoal){

	NlpFactory nlp;
	std::vector<PhaseNodes::Ptr> vars;

	double T = 1.5;
	int ee = 0;
	int phase_count = 2;
	int nodes_number = phase_count+1;
	int number_of_variables = nodes_number*6;
	/*
	 * TODO: Understand why the motion phase nodes passes this test only if the sequence starts in swing phase.
	 * Therefore the following test should be false for the test to pass.
	 */
	bool in_contact_start = false;
	const std::string var_name  ="vars";
	int n_polys_in_changing_phase = 2;

	auto motionNodes = std::make_shared<PhaseNodes>(phase_count, in_contact_start, var_name,
													n_polys_in_changing_phase, PhaseNodes::Motion);
	Eigen::VectorXd motion_values = motionNodes->GetValues();
	std::cout<<"default motion nodes values: "<<motion_values.transpose()<<std::endl;

	for (int idx = 0; idx<number_of_variables; idx++)
		EXPECT_EQ(0.0, motion_values(idx));

    Eigen::Vector3d initial_ee_pos_W = Eigen::Vector3d(1.0, 0.0, 0.42);
    Eigen::Vector3d final_ee_pos_W = Eigen::Vector3d(3.0, 0.0, 0.42);

    motionNodes->InitializeNodesTowardsGoal(initial_ee_pos_W, final_ee_pos_W, T);
    Eigen::VectorXd initialized_motion_values = motionNodes->GetValues();
    std::cout<<"initialized motion nodes values: "<<initialized_motion_values.transpose()<<std::endl;
    EXPECT_EQ(initial_ee_pos_W(0), initialized_motion_values(0));
    EXPECT_EQ(initial_ee_pos_W(1), initialized_motion_values(1));
    EXPECT_EQ(initial_ee_pos_W(2), initialized_motion_values(2));

}

TEST(PhaseNodesTest, forcePhaseNodesInitializationTowardsGoal){

	NlpFactory nlp;
	std::vector<PhaseNodes::Ptr> vars;

	double T = 1.5;
	int ee = 0;
	int phase_count = 2;
	int nodes_number = phase_count+1;
	int number_of_variables = nodes_number*6;
	/*
	 * TODO: Understand why the force phase nodes passes this test only if the sequence starts in stance phase.
	 * Therefore the following test should be true for the test to pass.
	 */
	bool in_contact_start = true;
	const std::string var_name  ="vars";
	int n_polys_in_changing_phase = 2;

	auto forceNodes = std::make_shared<PhaseNodes>(phase_count, in_contact_start, var_name,
													n_polys_in_changing_phase, PhaseNodes::Force);
	Eigen::VectorXd force_values = forceNodes->GetValues();
	std::cout<<"default motion nodes values: "<<force_values.transpose()<<std::endl;

	for (int idx = 0; idx<number_of_variables; idx++)
		EXPECT_EQ(0.0, force_values(idx));

    Eigen::Vector3d initial_ee_pos_W = Eigen::Vector3d(1.0, 0.0, 0.42);
    Eigen::Vector3d final_ee_pos_W = Eigen::Vector3d(3.0, 0.0, 0.42);

    forceNodes->InitializeNodesTowardsGoal(initial_ee_pos_W, final_ee_pos_W, T);
    Eigen::VectorXd initialized_force_values = forceNodes->GetValues();
    std::cout<<"initialized force nodes values: "<<initialized_force_values.transpose()<<std::endl;
    EXPECT_EQ(initial_ee_pos_W(0), initialized_force_values(0));
    EXPECT_EQ(initial_ee_pos_W(1), initialized_force_values(1));
    EXPECT_EQ(initial_ee_pos_W(2), initialized_force_values(2));

}

TEST(PhaseNodesTest, phaseNodesInitializationTowardsGoal){

	// Kinematic limits and dynamic parameters
	RobotModel model;
	model.dynamic_model_   = std::make_shared<MonopedDynamicModel>();
	model.kinematic_model_ = std::make_shared<MonopedKinematicModel>();
	// Parameters that define the motion. See c'tor for default values or
	// other values that can be modified.
	Parameters params;
	params.t_total_ = 1.6; // [s] time to reach goal state
	// here we define the initial phase durations, that can however be changed
	// by the optimizer. The number of swing and stance phases however is fixed.
	// alternating stance and swing:     ____-----_____-----_____-----_____
	params.ee_phase_durations_.push_back({0.4, 0.2, 0.4});
	params.min_phase_duration_ = 0.1;
	params.ee_in_contact_at_start_.push_back(true);

	NlpFactory nlp;
	std::vector<PhaseNodes::Ptr> vars;

	//double T = 1.5;
	int ee = 0;
	int phase_count = 3;
	int nodes_number = phase_count+1;
	//int number_of_variables = nodes_number*6;
	/*
	 * TODO: Understand why the force phase nodes passes this test only if the sequence starts in stance phase.
	 * Therefore the following test should be true for the test to pass.
	 */
	bool in_contact_start = false;
	//const std::string var_name  ="vars";
	int n_polys_in_changing_phase = 2;

	NlpFactory::VariablePtrVec variables;

	/* Create base linear nodes */
	std::vector<Nodes::Ptr> baseLinNodesList;
	auto baseLinNodes = std::make_shared<BaseNodes>(nodes_number, id::base_lin_nodes);
	baseLinNodes->InitializeNodesTowardsGoal(nlp.initial_base_.lin.p(), nlp.final_base_.lin.p(), params.t_total_);
	baseLinNodes->AddStartBound(kPos, {X,Y,Z}, nlp.initial_base_.lin.p());
	baseLinNodes->AddStartBound(kVel, {X,Y,Z}, nlp.initial_base_.lin.v());
	baseLinNodes->AddFinalBound(kPos, {X,Y},   nlp.final_base_.lin.p());
	baseLinNodes->AddFinalBound(kVel, {X,Y,Z}, nlp.final_base_.lin.v());
	baseLinNodesList.push_back(baseLinNodes);
	variables.insert(variables.end(), baseLinNodesList.begin(), baseLinNodesList.end());

	/* Create base angular nodes */
	std::vector<Nodes::Ptr> baseAngNodesList;
	auto baseAngNodes = std::make_shared<BaseNodes>(nodes_number, id::base_ang_nodes);
	baseAngNodes->InitializeNodesTowardsGoal(nlp.initial_base_.ang.p(), nlp.final_base_.ang.p(), params.t_total_);
	baseAngNodes->AddStartBound(kPos, {X,Y,Z}, nlp.initial_base_.ang.p());
	baseAngNodes->AddStartBound(kVel, {X,Y,Z}, nlp.initial_base_.ang.v());
	baseAngNodes->AddFinalBound(kPos, {X,Y},   nlp.final_base_.ang.p());
	baseAngNodes->AddFinalBound(kVel, {X,Y,Z}, nlp.final_base_.ang.v());
	baseAngNodesList.push_back(baseAngNodes);
	variables.insert(variables.end(), baseAngNodesList.begin(), baseAngNodesList.end());

	/* Create ee force nodes */
	std::vector<PhaseNodes::Ptr> forceNodesList;
	auto forceNodes = std::make_shared<PhaseNodes>(phase_count, in_contact_start, id::EEForceNodes(ee),
													n_polys_in_changing_phase, PhaseNodes::Force);
    // initialize with mass of robot distributed equally on all legs
    double m = model.dynamic_model_->m();
    double g = model.dynamic_model_->g();
    Eigen::Vector3d f_stance(0.0, 0.0, m*g);
    forceNodes->InitializeNodesTowardsGoal(f_stance, f_stance, params.t_total_);
    forceNodesList.push_back(forceNodes);
	variables.insert(variables.end(), forceNodesList.begin(), forceNodesList.end());

	/* Create ee motion nodes */
	std::vector<PhaseNodes::Ptr> motionNodesList;
	auto motionNodes = std::make_shared<PhaseNodes>(phase_count, in_contact_start, id::EEMotionNodes(ee),
													n_polys_in_changing_phase, PhaseNodes::Motion);
	/* initialize the starting and target points of the trajectory */
    Eigen::Vector3d initial_ee_pos_W = Eigen::Vector3d(1.0, 0.0, 0.42);
    Eigen::Vector3d final_ee_pos_W = Eigen::Vector3d(3.0, 0.0, 0.42);
    motionNodes->InitializeNodesTowardsGoal(initial_ee_pos_W, final_ee_pos_W, params.t_total_);
    motionNodes->AddStartBound(kPos, {X,Y,Z}, initial_ee_pos_W);
	motionNodesList.push_back(motionNodes);
	variables.insert(variables.end(), motionNodesList.begin(), motionNodesList.end());

	/* Create phase duration nodes */
	std::vector<PhaseDurations::Ptr> phaseDurationsList;
	auto phaseDurations = std::make_shared<PhaseDurations>(ee, params.ee_phase_durations_.at(ee),
	                                                params.ee_in_contact_at_start_.at(ee),
	                                                params.min_phase_duration_,
	                                                params.max_phase_duration_);
	phaseDurationsList.push_back(phaseDurations);
	variables.insert(variables.end(), phaseDurationsList.begin(), phaseDurationsList.end());

    Eigen::VectorXd initialized_force_values = forceNodes->GetValues();
    std::cout<<"initialized force nodes values: "<<initialized_force_values.transpose()<<std::endl;

	// set the initial position
	BaseState initial_base;
	initial_base.lin.at(kPos) << 1.0, 0.0, 1.0;
	// define the desired goal state
	BaseState goal;
	goal.lin.at(towr::kPos) << 3.0, 0.0, 1.0;

	Eigen::Vector3d initial_foot_pos_W = Eigen::Vector3d::Zero();


    SplineHolder spline_holder;

    //Nodes::Ptr base_lin_nodes = baseLinNodesList.at(ee);
    //Nodes::Ptr base_ang_nodes = baseAngNodesList.at(ee);
    const std::vector<double> base_poly_durations = params.GetBasePolyDurations();
    //std::vector<PhaseNodes::Ptr> ee_motion_nodes = motionNodesList;
    //std::vector<PhaseNodes::Ptr> ee_force_nodes = forceNodesList;
    //std::vector<PhaseDurations::Ptr> phase_durations = phaseDurationsList;

    std::cout<<"Define the Spline Holder "<<std::endl;

    spline_holder = SplineHolder(baseLinNodes, // linear
    									baseAngNodes, // angular
										base_poly_durations,
										motionNodesList,
										forceNodesList,
										phaseDurationsList,
    									false);

}
}
