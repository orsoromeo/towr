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
                                                      int numberofleg,
                                                      const std::vector<std::vector<double>>& phase_durations)
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueLin-"), geom_(model_,terrain, spline_holder, numberofleg, phase_durations), der_(terrain, spline_holder)

{
 numberofleg_=numberofleg;
 model_=model;
 node_variables_id_=node_variable_name;
 base_linear_=spline;
 node_bounds_.resize(6);
 NumberNodes_=(base_linear_->GetPolynomialCount()+1);
 SetRows(NumberNodes_*6);
 ee_motion_=spline_holder.ee_motion_;
 base_ << 0.0, 0.0, 1.0;
 lambda_=spline_holder.lambda_;
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

    auto com=base_linear_->GetPoint(t);
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

   if (k<NumberNodes_)
   {

    if (var_set==id::base_lin_nodes)
     {
      jac.middleRows(6*k,3)=FillJacobianLinWrenchWrtLin(t);
      jac.middleRows(6*k+3,3)=FillJacobianAngWrenchWrtLin(t,k);
     }
    if (var_set ==id::base_ang_nodes)
     {
      jac.middleRows(6*k+3,3)=FillJacobianAngWrenchWrtAng(t);
     }



      //spline_->GetJacobianWrtNodes(k,0.1, kAcc);
      //std::cout<<"lin "<<std::endl;
      //std::cout<<jac<<std::endl;

   }

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
NodeSpline::Jacobian
BaseAccConstraintRangeLin::FillJacobianLinWrenchWrtLin(double t) const
{
   Jacobian da=base_linear_->GetJacobianWrtNodes(t,kAcc);
   return m_*da;
}

NodeSpline::Jacobian
BaseAccConstraintRangeLin::FillJacobianAngWrenchWrtLin(double t, int k ) const
{
  Jac jac1=DerivativeOfrxma(t);

  Jac jac2= Jac(3, jac1.cols());
  jac2.coeffRef(0,6*k+1)=m_*g_;
  jac2.coeffRef(1,6*k)=-m_*g_;
  Jac total_jac=jac1+jac2;


  return total_jac;
}

BaseAccConstraintRangeLin::Jac
BaseAccConstraintRangeLin::DerivativeOfrxma(double t) const
{
  Jacobian dr=base_linear_->GetJacobianWrtNodes(t,kPos);
  Jacobian da=base_linear_->GetJacobianWrtNodes(t,kAcc);
  State s=base_linear_->GetPoint(t);
  Eigen::Vector3d p=s.p();
  Eigen::Vector3d a=s.a();
  Jac jac3(3,dr.cols());
  jac3.row(0)=p.y()*da.row(2)+dr.row(1)*a.z()-p.z()*da.row(1)-dr.row(2)*a.y();
  jac3.row(1)=p.z()*da.row(0)+dr.row(2)*a.x()-p.x()*da.row(2)-dr.row(0)*a.z();
  jac3.row(2)=p.x()*da.row(1)+dr.row(0)*a.y()-p.y()*da.row(0)-dr.row(1)*a.x();
  return jac3;
}

NodeSpline::Jacobian
BaseAccConstraintRangeLin::FillJacobianAngWrenchWrtAng(double t) const
{
  Eigen::Vector3d acc= base_angular_.GetAngularAccelerationInWorld(t);
  Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
  Jac I_w = w_R_b.sparseView() * I_b * w_R_b.transpose().sparseView();

  // 1st term of product rule (derivative of R)
  Eigen::Vector3d v11 = I_b*w_R_b.transpose()*acc;

  Jac jac11 = base_angular_.DerivOfRotVecMult(t, v11, false);

  // 2nd term of product rule (derivative of R^T)
  Jac jac12 = w_R_b.sparseView()*I_b*base_angular_.DerivOfRotVecMult(t, acc, true);

  // 3rd term of product rule (derivative of wd)
  Jac jac_ang_acc = base_angular_.GetDerivOfAngAccWrtEulerNodes(t);
  Jac jac13 = I_w * jac_ang_acc;


  return jac11 + jac12 + jac13;
}


} /* namespace towr */

