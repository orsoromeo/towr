#ifndef TOWR_CONSTRAINTS_BASE_ACC_CONSTRAINT_RANGE_LIN_H_
#define TOWR_CONSTRAINTS_BASE_ACC_CONSTRAINT_RANGE_LIN_H_

#include <towr/variables/spline_holder.h>
#include <towr/variables/spline.h>
#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/terrain/height_map.h>
#include <towr/models/dynamic_model.h>

#include <iostream>

#include "time_discretization_constraint.h"

namespace towr {

/**
 * @brief Keeps the base acc in a specified range.
 *
 * 
 *
 * @ingroup Constraints
 */
class Geometry {
public:
  /**
   * @brief Links the base variables and sets hardcoded bounds on the state.
   * @param T  The total time of the optimization horizon.
   * @param dt The discretization interval of the constraints.
   * @param spline_holder  Holds pointers to the base variables.
   */
  Geometry (DynamicModel::Ptr model, HeightMap::Ptr terrain, const SplineHolder& spline_holder);
  virtual ~Geometry () = default;


private:
  
  Eigen::MatrixXd  LinearEdges_;
  Eigen::MatrixXd  AngularEdges_;
  Eigen::Vector3d  base_;
  HeightMap::Ptr terrain_;
  std::vector<NodeSpline::Ptr> ee_motion_;
  std::vector<NodeSpline::Ptr> ee_motion_in_touch_;
  mutable DynamicModel::Ptr model_;
  int GetNumberOfFeetInTouch (double t);
  Eigen::MatrixXd ReturnNormalTerrain (double t);
  Eigen::MatrixXd ReturnTangentTerrain(double t);
  Eigen::Matrix3d RotationMatrix (Eigen::Vector3d t, double angle) const;
  Eigen::Vector3d ComputeNextEdge (Eigen::Matrix3d M, Eigen::Vector3d n) const;
  void ComputeCone (double t);
  Eigen::MatrixXd ComputeLinearPartOfTheCone (Eigen::Vector3d axis, double angle);
  Eigen::MatrixXd ComputeAngularPartOfTheCone (double t, Eigen::Vector3d axis, double anglem, double ee);
  Eigen::Matrix3d RotationX (double angle) const;
  Eigen::Matrix3d RotationZ (double angle) const;
  double ComputeRotationAngle (Eigen::Vector3d axis) const;
  
};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_BASE_MOTION_CONSTRAINT_H_ */
