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

#ifndef TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HYQ_MODEL_H_
#define TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HYQ_MODEL_H_

#include <towr/models/kinematic_model.h>
#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/models/endeffector_mappings.h>

namespace towr {

/**
 * @brief The Kinematics of the quadruped robot HyQ.
 */
class HyqKinematicModel : public KinematicModel {
public:
  HyqKinematicModel () : KinematicModel(4)
  {
    const double x_nominal_b = 0.36743;
    const double y_nominal_b = 0.3272;
    const double z_nominal_b = -0.55;

    nominal_stance_.at(LF) <<  x_nominal_b,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RF) <<  x_nominal_b,  -y_nominal_b, z_nominal_b;
    nominal_stance_.at(LH) << -x_nominal_b,   y_nominal_b, z_nominal_b;
    nominal_stance_.at(RH) << -x_nominal_b,  -y_nominal_b, z_nominal_b;

    max_dev_from_nominal_ << 0.25, 0.10, 0.06;
    //ThetaL_<<1.07, 1.07, 2.31, 2.31, 3.96, 3.96, -0.82, -0.82, -0.0129, -0.0129, 3.154, 3.154, 3.12, 3.12, 0.0129, 0.0129;

    ThetaL_<<0.19, 0.19, 2.94, 2.94, 3.33, 3.33, -0.19, -0.19, -0.49, -0.49, 3.64, 3.64, 2.64, 2.64, 0.49, 0.49;
    ThetaN_<<0.82, 0.82, 2.31, 2.31, 3.96, 3.96, -0.82, -0.82, -0.0129, -0.0129, 3.154, 3.154, 3.12, 3.12, 0.0129, 0.0129;
    ThetaR_<<1.18, 1.18, 1.95, 1.95, 4.32, 4.32, -1.18, -1.18, 3.62, 3.62, -0.47, -0.47, 0.47, 0.47, 2.66, 2.66;
    //ThetaR_<<0.82, 0.82, 2.31, 2.31, 3.96, 3.96, -0.82, -0.82, -0.0129, -0.0129, 3.154, 3.154, 3.12, 3.12, 0.0129, 0.0129;
    //Theta_cost_<<1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0;
    
    coeffDL_<<300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300;
    coeffDN_<<300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300;
    coeffDR_<<300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300, 300;

    //coeffDL_<<300, 300, 300, 300, 300, 300, 300, 300, 212, 212, 212, 212, 212, 212, 212, 212;
    //coeffDN_<<300, 300, 300, 300, 300, 300, 300, 300, 212, 212, 212, 212, 212, 212, 212, 212;
    //coeffDR_<<300, 300, 300, 300, 300, 300, 300, 300, 212, 212, 212, 212, 212, 212, 212, 212;
 }
};

/**
 * @brief The Dynamics of the quadruped robot HyQ.
 */
class HyqDynamicModel : public SingleRigidBodyDynamics {
public:
  HyqDynamicModel() : SingleRigidBodyDynamics(86.8659,
                      3.9369, 11.1741, 12.5216, 0.1459, -0.3039, -0.0135,
                      4) {}
};

} /* namespace towr */

#endif /* TOWR_TOWR_ROS_INCLUDE_TOWR_ROS_HYQ_MODEL_H_ */
