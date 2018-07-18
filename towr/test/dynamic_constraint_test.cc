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

#include <towr/towr.h>
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

#include <towr/terrain/examples/height_map_examples.h>

#include <towr/models/examples/monoped_model.h>
#include <towr/models/centroidal_model.h>

#include <towr/costs/node_cost.h>

#include <ifopt/ipopt.h>

namespace towr {

TEST(DynamicConstraintTest, testDynamicConstraintValues){

	ifopt::Problem nlp;
	NlpFactory factory;

	// Kinematic limits and dynamic parameters
	RobotModel model;
	model.dynamic_model_   = std::make_shared<MonopedDynamicModel>();
	model.kinematic_model_ = std::make_shared<MonopedKinematicModel>();


	// set the initial position
	BaseState initial_base;
	initial_base.lin.at(kPos).z() = 0.5;

	Eigen::Vector3d initial_foot_pos_W = Eigen::Vector3d::Zero();

	// define the desired goal state
	BaseState goal;
	goal.lin.at(towr::kPos) << 1.0, 0.0, 0.5;


	// Parameters that define the motion. See c'tor for default values or
	// other values that can be modified.
	Parameters params;
	params.t_total_ = 1.6; // [s] time to reach goal state
	// here we define the initial phase durations, that can however be changed
	// by the optimizer. The number of swing and stance phases however is fixed.
	// alternating stance and swing:     ____-----_____-----_____-----_____
	params.ee_phase_durations_.push_back({0.4, 0.2, 0.4});
	params.ee_in_contact_at_start_.push_back(true);
	DynamicModel::Ptr monopedDynamicModel = std::make_shared<MonopedDynamicModel>();

    //for (auto v : factory.GetVariableSets()){
    //	Eigen::VectorXd variablesValues = v->GetValues();
    //	std::cout<<"Variables values: "<<variablesValues.transpose()<<std::endl;
    //}

	std::cout<<"I am here! "<<std::endl;

	factory.initial_base_ = initial_base;
	factory.initial_ee_W_ = {initial_foot_pos_W};
	factory.final_base_ = goal;
	factory.params_ = params;
	factory.model_ = model;
	factory.terrain_ = std::make_shared<FlatGround>();

    for (auto v : factory.GetVariableSets()){
    	std::cout<<"Variables name: "<<v->GetName()<<std::endl;
    	Eigen::VectorXd variablesValues = v->GetValues();
    	std::cout<<"Variables values: "<<variablesValues.transpose()<<std::endl;
    }

    auto dynamicConstraint = std::make_shared<DynamicConstraint>(monopedDynamicModel,
                                                          params,
														  factory.spline_holder_);
    std::cout<<"add dynamic constraint"<<std::endl;
    nlp.AddConstraintSet(dynamicConstraint);

    //Eigen::VectorXd constraintValues = dynamicConstraint->GetValues();
    //for (auto c : factory.GetConstraints()){
    	//nlp.AddConstraintSet(c);
    	//Eigen::VectorXd constraintValues = c->GetValues();
    //}

    //for (auto c : factory.GetCosts())
    	//nlp.AddCostSet(c);

    //nlp_problem.AddConstraintSet(dynamicConstraint);
    //Eigen::VectorXd constraintValues = dynamicConstraint->GetValues();
    //std::cout<<"Dynamic constraint values: "<<constraintValues.transpose()<<std::endl;
}

TEST(DynamicConstraintTest, testConstraintSetVariables){

	NlpFactory nlp_factory;
	ifopt::Problem nlp_problem;
	NlpFactory::VariablePtrVec vars = nlp_factory.GetVariableSets();

	double T = 1.5;
	int ee = 0;
	int phase_count = 2;
	int nodes_number = phase_count+1;
	int number_of_variables = nodes_number*6;

    DynamicModel::Ptr monopedDynamicModel = std::make_shared<MonopedDynamicModel>();
    Parameters params;
    std::cout<<"created dynamic model"<<std::endl;

    auto dynamicConstraint = std::make_shared<DynamicConstraint>(monopedDynamicModel,
                                                          params,
														  nlp_factory.spline_holder_);
    for (auto v : nlp_factory.GetVariableSets()){
    	Eigen::VectorXd variablesValues = v->GetValues();
    	std::cout<<"Variables values: "<<variablesValues.transpose()<<std::endl;
    }

    nlp_problem.AddConstraintSet(dynamicConstraint);
    Eigen::VectorXd constraintValues = dynamicConstraint->GetValues();
    std::cout<<"Dynamic constraint values: "<<constraintValues.transpose()<<std::endl;
}

}
