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

#ifndef TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_MONOPED_MODEL_H_
#define TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_MONOPED_MODEL_H_

#include <towr/models/kinematic_model.h>
#include <towr/models/single_rigid_body_dynamics.h>

namespace towr {

/**
 * @brief The Kinematics of a one-legged hopper with HyQ leg.
 */
class MonopedKinematicModel : public KinematicModel {
public:
  MonopedKinematicModel () : KinematicModel(1)
  {
    nominal_stance_.at(0) = Eigen::Vector3d( 0.0, 0.0, -0.58);
    max_dev_from_nominal_ << 0.25, 0.15, 0.2;


    // values for debugging:
    //ThetaL_<<0.0,  0.0,  0.0,  0.0,  1.57,   1.57,   1.57,   1.57,   3.1415, 3.1415, 3.1415, 3.1415, 4.71, 4.71, 4.71, 4.71;
    //ThetaN_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;
    //ThetaR_<<0.0,  0.0,  0.0,  0.0,  1.57,   1.57,   1.57,   1.57,   3.1415, 3.1415, 3.1415, 3.1415, 4.71, 4.71, 4.71, 4.71;

    //ThetaL_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;
    //ThetaN_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;
    //ThetaR_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;

    ThetaL_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;
    ThetaN_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;
    ThetaR_<<1.57, 1.57, 1.57, 1.57, 3.1415, 3.1415, 3.1415, 3.1415, 4.71,   4.71,   4.71,   4.71,    0.0,  0.0,  0.0,  0.0;

    // values for HyQ:
    //ThetaL_<<0.19, 0.19, 2.94, 2.94, 3.33, 3.33, -0.19, -0.19, -0.49, -0.49, 3.64, 3.64, 2.64, 2.64, 0.49, 0.49;
    //ThetaN_<<0.82, 0.82, 2.31, 2.31, 3.96, 3.96, -0.82, -0.82, -0.0129, -0.0129, 3.154, 3.154, 3.12, 3.12, 0.0129, 0.0129;
    //ThetaR_<<1.18, 1.18, 1.95, 1.95, 4.32, 4.32, -1.18, -1.18, 3.62, 3.62, -0.47, -0.47, 0.47, 0.47, 2.66, 2.66;
    
    //ThetaR_<<0.82, 0.82, 2.31, 2.31, 3.96, 3.96, -0.82, -0.82, -0.0129, -0.0129, 3.154, 3.154, 3.12, 3.12, 0.0129, 0.0129;
    //Theta_cost_<<1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    
    //coeffDL_<<1000, 1000, 1000, 1000;
    //coeffDN_<<1000, 1000, 1000, 1000;
    //coeffDR_<<1000, 1000, 1000, 1000;

    coeffDL_<<1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000, 1000;
    coeffDN_<<800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800, 800;
    coeffDR_<<200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200;

    //coeffDL_<<300, 300, 300, 300, 300, 300, 300, 300, 212, 212, 212, 212, 212, 212, 212, 212;
    //coeffDN_<<300, 300, 300, 300, 300, 300, 300, 300, 212, 212, 212, 212, 212, 212, 212, 212;
    //coeffDR_<<300, 300, 300, 300, 300, 300, 300, 300, 212, 212, 212, 212, 212, 212, 212, 212;
  }
};

/**
 * @brief The Dynamics of a one-legged hopper with HyQ leg.
 */
class MonopedDynamicModel : public SingleRigidBodyDynamics {
public:
  MonopedDynamicModel()
  : SingleRigidBodyDynamics(20,                      // mass of the robot
                    1.2, 5.5, 6.0, 0.0, -0.2, -0.01, // base inertia
                    1) {}                            // number of endeffectors
};

} /* namespace towr */

#endif /* TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_MONOPED_MODEL_H_ */
