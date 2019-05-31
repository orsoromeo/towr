#ifndef TOWR_CONSTRAINTS_DERIVATIVE_H_
#define TOWR_CONSTRAINTS_DERIVATIVE_H_

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
namespace towr {
class Derivative {
public:

  Derivative (const HeightMap::Ptr & terrain,
              const SplineHolder& spline_holder, int numberofleg, const std::vector<std::vector<double> > &phase_durations);
  virtual ~Derivative () = default;
  NodeSpline::Jacobian GetDerivativeofLinearEdgeWrtNodes(int ee, double t,Eigen::Vector3d normal, double angle, Eigen::Vector3d edge) const;
  NodeSpline::Jacobian GetDerivativeofAngularEdgeWrtNodes(int ee, double t, Eigen::Vector3d normal, double angle, Eigen::Vector3d edge, Eigen::Vector3d edges) const;
  mutable std::vector<NodeSpline::Ptr> ee_motion_in_touch_;
  void GetNumberOfFeetInTouch (double t) const;
private:
Eigen::Vector3d  base_;
HeightMap::Ptr   terrain_;
mutable std::vector<NodeSpline::Ptr> ee_motion_;
std::vector<std::vector<double>> phase_durations_;
int numberoflegs_;
NodeSpline::Jacobian GetDerivativeOfNormalWrtNodes (int ee, double t) const;
Eigen::Vector3d      GetDerivativeofAngleWrtNormal(Eigen::Vector3d normal, double angle) const;
NodeSpline::Jacobian GetDerivativeofAngleWrtNodes(int ee, double t, Eigen::Vector3d normal, double angle) const;
NodeSpline::Jacobian GetDerivativeofLinearEdgeWrtNormal (Eigen::Vector3d normal, double angle, Eigen::Vector3d edge) const;
Eigen::MatrixXd      GetDerivativeofRotationAxisWrtNormal(Eigen::Vector3d vector) const;
};

} /* namespace towr */

#endif /* TOWR_CONSTRAINTS_DERIVATIVE_H_ */
