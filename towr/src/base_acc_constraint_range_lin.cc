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

#include <towr/constraints/base_acc_constraint_range_lin.h>
#include <ifopt/constraint_set.h>
#include <ifopt/problem.h>
#include <politopixAPI.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/terrain/height_map.h>
#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <ifopt/composite.h>
#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/models/dynamic_model.h>
#include <iostream>
#include <cmath>
#define PI 3.14159265

namespace towr {


BaseAccConstraintRangeLin::BaseAccConstraintRangeLin (const DynamicModel::Ptr& model,
                                                      double T,
                                                      double dt,
                                                      const NodeSpline::Ptr& spline,
                                                      std::string node_variable_name,
                                                      const HeightMap::Ptr& terrain,
                                                      const SplineHolder& spline_holder,
                                                      int numberofleg)
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueLin-"), geom_(model_,terrain, spline_holder, numberofleg)

{
 numberofleg_=numberofleg;
 model_=model;
 node_variables_id_=node_variable_name;
 spline_=spline;
 node_bounds_.resize(6);
 NumberNodes_=(spline_->GetPolynomialCount()+1);
 SetRows(NumberNodes_*6);
 ee_motion_=spline_holder.ee_motion_;
 base_ << 0.0, 0.0, 1.0;
 lambda_=spline_holder.lambda_;
 //Derivative der(terrain, spline_holder);
 m_=model_->m();
 g_= model_->g();
 base_angular_ = EulerConverter(spline_holder.base_angular_);
 I_b=model_->ReturnInertia();

}

void
BaseAccConstraintRangeLin::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{

  if (k<NumberNodes_)
  {

    auto com=spline_->GetPoint(t);
    g.middleRows(GetRow(k, 0), 6) =FillConstraint(com,t);

  }
}

void
BaseAccConstraintRangeLin::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  if (k<NumberNodes_)
  {
    for (int i=0; i<6; i++)
    bounds.at(GetRow(k,i)) = ifopt::BoundZero;


  }
} //

void
BaseAccConstraintRangeLin::UpdateJacobianAtInstance (double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{
//
//   if (k<NumberNodes_)
//   {
//
//   //if (var_set==id::base_lin_nodes)
//   if (var_set ==node_variables_id_)
//    {
//
//      jac.middleRows(6*k,6)=FillJacobian(spline_,t);
//
//      //spline_->GetJacobianWrtNodes(k,0.1, kAcc);
//      //std::cout<<"lin "<<std::endl;
//      //std::cout<<jac<<std::endl;
//
//   }
// }
}
int
BaseAccConstraintRangeLin::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}
Eigen::VectorXd
BaseAccConstraintRangeLin::FillConstraint (State com, double t) const
 {
  VectorXd g;
  g.resize(6);
  g=ComputeWrench(com,t);
  auto lambda=lambda_->GetPoint(t).p();
  Eigen::MatrixXd edges=geom_.ComputeCone(t);
  auto lambdaedge=edges*lambda;
  return g-lambdaedge;
}
NodeSpline::Jacobian
BaseAccConstraintRangeLin::FillJacobian(NodeSpline::Ptr spline_,double t) const
{

  auto jac1=spline_->GetJacobianWrtNodes(t, kAcc);

  double mu=10.0;
  Jacobian g=Jacobian(6, jac1.cols());
  g.row(0)=jac1.row(0)+mu*jac1.row(2);
  g.row(1)=jac1.row(1)+mu*jac1.row(2);
  g.row(2)=jac1.row(0)-mu*jac1.row(2);
  g.row(3)=jac1.row(1)-mu*jac1.row(2);
  
  return g;

}

Eigen::Vector3d BaseAccConstraintRangeLin::ComputeLinearWrench (State com) const
{
  Eigen::Vector3d lin_wrench=m_*(com.a()+Eigen::Vector3d(0.0, 0.0, g_));
  return lin_wrench;
}

Eigen::Vector3d
BaseAccConstraintRangeLin::ComputeAngularWrench (Eigen::VectorXd acc, Jac I_w, State com) const
{

  auto Iwacc=I_w*acc;
  Eigen::Vector3d rp=com.p();
  Eigen::Vector3d ra=com.a();
  auto rmg=rp.cross(m_*(Eigen::Vector3d (0.0, 0.0, g_)));
  Eigen::Vector3d rma=rp.cross(m_*ra);
  Eigen::Vector3d ang_wrench=Iwacc+rmg+rma;
  return ang_wrench;
 }

Eigen::VectorXd
BaseAccConstraintRangeLin::ComputeWrench (State com, double t) const
{
  VectorXd g;
  g.resize(6);
  g.middleRows(0, 3)=ComputeLinearWrench(com);
  Eigen::Vector3d acc= base_angular_.GetAngularAccelerationInWorld(t);
  Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  Jac I_w = w_R_b.sparseView() * I_b* w_R_b.transpose().sparseView();
  //Jac I_w=I_w1.sparseView();
  g.middleRows(3, 3)=ComputeAngularWrench(acc,I_w,com);
  return g;
}


} /* namespace towr */

