#ifndef TOWR_CONSTRAINTS_GEOMETRY_H_
#define TOWR_CONSTRAINTS_GEOMETRY_H_

#include <ifopt/problem.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/spline.h>
#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/terrain/height_map.h>
#include <towr/models/dynamic_model.h>
#include <ifopt/constraint_set.h>
#include <ifopt/composite.h>
#include "time_discretization_constraint.h"
#include <towr/constraints/total_duration_constraint.h>
#include <iostream>
#include <towr/variables/phase_durations.h>
#include <towr/constraints/derivative.h>
namespace towr {


class Geometry {
public:

  Geometry (DynamicModel::Ptr model,
            const HeightMap::Ptr &terrain,
            const SplineHolder& spline_holder, int numberofleg,
            const std::vector<std::vector<double>> & phase_durations);
  virtual ~Geometry () = default;
  Eigen::MatrixXd ComputeCone(double t) const;
  mutable Eigen::MatrixXd  LinearEdges_;
  mutable Eigen::MatrixXd  AngularEdges_;
  mutable Eigen::MatrixXd  EhatLin_;
  mutable Eigen::MatrixXd  EhatAng_;
  mutable Eigen::VectorXd  angle_;
  mutable Eigen::MatrixXd  normal_;
  mutable Eigen::MatrixXd  Edges_;
  mutable std::vector<NodeSpline::Ptr> ee_motion_in_touch_;
  Eigen::MatrixXd ReturnNormalTerrain (double t) const ;
  double ComputeRotationAngle (Eigen::Vector3d normal) const;
  int GetNumberOfFeetInTouch (double t) const;
  Eigen::MatrixXd ComputeLinearPartOfTheCone(const Eigen::Vector3d &axis, const double &angle) const;
  int numberoflegs_;
  bool IsInTouch (double t,int ee) const;
private:
  int number_of_leg_in_touch_;
  Eigen::MatrixXd  ComputeLinear_;
  Eigen::Vector3d  base_;
  HeightMap::Ptr terrain_;
  std::vector<std::vector<double>> phase_durations_;
  mutable std::vector<NodeSpline::Ptr> ee_motion_;

  mutable DynamicModel::Ptr model_;
  void SetZeroGlobalVariables() const;
  Eigen::Matrix3d RotationMatrix (Eigen::Vector3d t, double angle) const;
  Eigen::Vector3d ComputeNextEdge (Eigen::Matrix3d M, Eigen::Vector3d n) const;
  Eigen::MatrixXd ComputeAngularPartOfTheCone (double t, double ee) const;
  Eigen::Matrix3d RotationX (const double &angle) const;
  Eigen::Matrix3d RotationZ (const double &angle) const;


};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_GEOMETRY_H_ */
