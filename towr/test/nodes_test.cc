/**
 * 07/11/2018
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

#include <towr/costs/node_cost.h>

namespace towr {

TEST(NodesTest, nodesInitializationTowardsGoal){
	std::cout<<"test initialization"<<std::endl;
//	VariablePtrVars vars;
	NlpFactory nlp;
	std::vector<PhaseNodes::Ptr> vars;

	// Endeffector Motions
	double T = 1.5;
	int ee = 0;
	int phase_count = 3;
	bool in_contact_start = true;
	const std::string var_name  ="vars";
	int n_polys_in_changing_phase = 2;

	auto motionNodes = std::make_shared<PhaseNodes>(phase_count, in_contact_start, var_name,
													n_polys_in_changing_phase, PhaseNodes::Motion);
	Eigen::VectorXd motion_values = motionNodes->GetValues();
	std::cout<<"default motion nodes values: "<<motion_values.transpose()<<std::endl;

	EXPECT_EQ(0.0, motion_values(0));

	double yaw = 0.0;
	//Eigen::Vector3d euler(0.0, 0.0, yaw);
    //Eigen::Matrix3d w_R_b = EulerConverter::GetRotationMatrixBaseToWorld(euler);

    Eigen::Vector3d initial_ee_pos_W = Eigen::Vector3d(1.0, 0.0, 0.42);
    Eigen::Vector3d final_ee_pos_W = Eigen::Vector3d(3.0, 0.0, 0.42);
    ////std::cout<<"test initialization"<<motionNodes->GetValues()<<std::endl;
    motionNodes->InitializeNodesTowardsGoal(initial_ee_pos_W, final_ee_pos_W, T);
    Eigen::VectorXd initialized_motion_values = motionNodes->GetValues();
    std::cout<<"initialized motion nodes values: "<<initialized_motion_values.transpose()<<std::endl;
    EXPECT_EQ(initial_ee_pos_W(0), initialized_motion_values(0));
    EXPECT_EQ(initial_ee_pos_W(1), initialized_motion_values(1));
    EXPECT_EQ(initial_ee_pos_W(2), initialized_motion_values(2));

    //
    //motionNodes->AddStartBound(kPos, {X,Y,Z}, nlp.initial_ee_W_.at(ee));
    //vars.push_back(motionNodes);

}
}
