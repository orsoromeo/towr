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

#ifndef TOWR_MODELS_KINEMATIC_MODEL_H_
#define TOWR_MODELS_KINEMATIC_MODEL_H_

#include <memory>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

namespace towr {

/**
 * @brief  Contains all the robot specific kinematic parameters.
 *
 * This class is mainly used to formulate the @ref RangeOfMotionConstraint,
 * restricting each endeffector to stay inside it's kinematic range.
 *
 * @ingroup Robots
 */
class KinematicModel {
public:
  using Ptr      = std::shared_ptr<KinematicModel>;
  using EEPos    = std::vector<Eigen::Vector3d>;
  using Vector3d = Eigen::Vector3d;

  /**
   * @brief Constructs a kinematic model of a robot with zero range of motion.
   * @param n_ee  The number of endeffectors of the robot.
   */
  KinematicModel (int n_ee)
  {
    nominal_stance_.resize(n_ee);
    max_dev_from_nominal_.setZero();
    ThetaL_.resize(4,4);
    ThetaN_.resize(4,4);
    ThetaR_.resize(4,4);
    Theta_cost_.resize(2,4);
    coeffDL_.resize(4,4);
    coeffDN_.resize(4,4);
    coeffDR_.resize(4,4);
    

  }

  virtual ~KinematicModel () = default;

  /**
   * @brief  The xyz-position [m] of each foot in default stance.
   * @returns The vector from base to each foot expressed in the base frame.
   */
  virtual EEPos GetNominalStanceInBase() const
  {
    return nominal_stance_;
  }

  /**
   * @brief How far each foot can deviate from its nominal position.
   * @return The deviation [m] expresed in the base frame.
   */
  virtual Vector3d GetMaximumDeviationFromNominal() const
  {
    return max_dev_from_nominal_;
  }

  /**
   * @returns returns the number of endeffectors of this robot.
   */
  int GetNumberOfEndeffectors() const
  {
    return nominal_stance_.size();
  }
  Eigen::MatrixXd GetThetaL () const
  { 
    return ThetaL_;
  }
  Eigen::MatrixXd GetThetaN () const
  {
    return ThetaN_;
  }
  Eigen::MatrixXd GetThetaR () const
  {
    return ThetaR_;
  }
  Eigen::MatrixXd GetThetaCost () const
  {
    return Theta_cost_;
  }
  Eigen::MatrixXd GetDL () const
  {
    return coeffDL_;
  }
  Eigen::MatrixXd GetDN () const
  {
    return coeffDN_;
  }
  Eigen::MatrixXd GetDR () const
  {
    return coeffDR_;
  }
  Eigen::VectorXd GetDCost () const
  {
    return coeffDcost_;
  }
protected:
  EEPos nominal_stance_;
  Vector3d max_dev_from_nominal_;
  public:
  Eigen::MatrixXd ThetaL_;
  Eigen::MatrixXd ThetaN_;
  Eigen::MatrixXd ThetaR_;
  Eigen::MatrixXd Theta_cost_;
  Eigen::MatrixXd coeffDL_;
  Eigen::MatrixXd coeffDN_;
  Eigen::MatrixXd coeffDR_;
  Eigen::Vector2d coeffDcost_;


};

} /* namespace towr */


#endif /* TOWR_MODELS_KINEMATIC_MODEL_H_ */
