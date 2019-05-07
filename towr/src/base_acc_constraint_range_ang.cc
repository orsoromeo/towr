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
#include <towr/constraints/base_acc_constraint_range_ang.h>
#include <towr/constraints/dynamic_constraint.h>
#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/models/dynamic_model.h>
#include <towr/variables/euler_converter.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>
#include <iostream>



namespace towr {


BaseAccConstraintRangeAng::BaseAccConstraintRangeAng (const DynamicModel::Ptr& model, double T, double dt,
                                                      const NodeSpline::Ptr& angular, const NodeSpline::Ptr& linear,
                                                      std::string node_variable_name)
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueAng-")

{
 model_=model;
 spline_=angular;
 base_linear_=linear;
 node_variables_id_=node_variable_name;
 base_angular_ = EulerConverter(angular);
 node_bounds_.resize(4);
 SetRows((spline_->GetPolynomialCount()+1)*4);
 I_b=model_->ReturnInertia();
 m_=model_->m();
 g_= model_->g();
 //I_b << 1.2, 0.0, 0.2,
 //       0.0, 5.5, 0.0,
 //       0.2, 0.0, 6.0;
 //
}

void
BaseAccConstraintRangeAng::UpdateConstraintAtInstance (double t,int k,
                                                  VectorXd& g) const
{

  if (k<spline_->GetPolynomialCount()+1)
  {
   Eigen::Vector3d acc= base_angular_.GetAngularAccelerationInWorld(t);
   Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
   Jac I_w = w_R_b.sparseView() * I_b* w_R_b.transpose().sparseView();
   //Jac I_w=I_w1.sparseView();
   State r=base_linear_->GetPoint(t);
   g.middleRows(GetRow(k,0), 4) = FillConstraint(acc, I_w, r);


  }
}
void
BaseAccConstraintRangeAng::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  if (k<spline_->GetPolynomialCount()+1)
  {
    bounds.at(GetRow(k,0)) = ifopt::BoundGreaterZero;
    bounds.at(GetRow(k,1)) = ifopt::BoundGreaterZero;
    bounds.at(GetRow(k,2)) = ifopt::BoundSmallerZero;
    bounds.at(GetRow(k,3)) = ifopt::BoundSmallerZero;
  }
} //

void
BaseAccConstraintRangeAng::UpdateJacobianAtInstance (double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{

 if (k<spline_->GetPolynomialCount()+1)
 {
   if (var_set ==node_variables_id_)
 {
   Eigen::Vector3d acc= base_angular_.GetAngularAccelerationInWorld(t);
   Eigen::Matrix3d w_R_b = base_angular_.GetRotationMatrixBaseToWorld(t);
   jac.middleRows(4*k, 4)= FillJacobian(w_R_b,acc, k,t );
   //std::cout<<"ang "<<std::endl;
   //std::cout<<jac<<std::endl;
 }
 }
}
int
BaseAccConstraintRangeAng::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}

Eigen::VectorXd
BaseAccConstraintRangeAng::FillConstraint (Eigen::VectorXd acc, Jac I_w, State r) const
{

  auto Iwacc=I_w*acc;
  //Eigen::Vector3d grav=(0.0, 0.0, -9.81);
  Eigen::Vector3d rp=r.p();
  Eigen::Vector3d ra=r.a();
  auto rmg=rp.cross(m_*(Eigen::Vector3d (0.0, 0.0, g_)));
  Eigen::Vector3d rma=ra.cross(m_*ra);
  Eigen::Vector3d com=Iwacc+rma+rmg;
  Eigen::VectorXd constraint;
  constraint.resize(4);
  //double g[4];
  double mu=4.0;
  constraint(0)=com(0)+mu*com(2);//>0    mu1_lx(p1, p2, p3, p4)*w_lx + mu1_ly(..)*w_ly + mu1_lz(..)*w_lz + mu1_ax(..)*w_ax + mu1_ay(..)*w_ay + mu1_az(..)*w_az
  constraint(1)=com(1)+mu*com(2);//>0    mu1*x + mu2*y + mu3*z
  constraint(2)=com(0)-mu*com(2);//<0    mu1*x + mu2*y + mu3*z
  constraint(3)=com(1)-mu*com(2);//<0    mu1*x + mu2*y + mu3*z
  return constraint;
}

NodeSpline::Jacobian
BaseAccConstraintRangeAng::FillJacobian(Eigen::Matrix3d w_R_b, Eigen::Vector3d acc, int k, double t ) const
{
  //Jac I_b_sparse=I_b.sparseView();
  Jac I_w = w_R_b.sparseView() * I_b * w_R_b.transpose().sparseView();

  // 1st term of product rule (derivative of R)
  Eigen::Vector3d v11 = I_b*w_R_b.transpose()*acc;
  
  Jac jac11 = base_angular_.DerivOfRotVecMult(t, v11, false);

  // 2nd term of product rule (derivative of R^T)
  Jac jac12 = w_R_b.sparseView()*I_b*base_angular_.DerivOfRotVecMult(t, acc, true);

  // 3rd term of product rule (derivative of wd)
  Jac jac_ang_acc = base_angular_.GetDerivOfAngAccWrtEulerNodes(t);
  Jac jac13 = I_w * jac_ang_acc;
  

  Jac jac1 = jac11 + jac12 + jac13;
  //Derivative of rxma
  Jac jac3=DerivativeOfrxma(t);
  //Derivative of rxmg
  Jac jac2= Jac(3, jac1.cols());
  jac2.coeffRef(0,6*k+1)=-m_*g_;
  jac2.coeffRef(1,6*k)=+m_*g_;

  Jac total_jac=jac1+jac3-jac2;


  double mu=4.0;
  NodeSpline::Jacobian jacobian=NodeSpline::Jacobian(4, total_jac.cols());
  jacobian.row(0)=total_jac.row(0)+mu*total_jac.row(2);
  jacobian.row(1)=total_jac.row(1)+mu*total_jac.row(2);
  jacobian.row(2)=total_jac.row(0)-mu*total_jac.row(2);
  jacobian.row(3)=total_jac.row(1)-mu*total_jac.row(2);
 
  return jacobian;
}

BaseAccConstraintRangeAng::Jac
BaseAccConstraintRangeAng::DerivativeOfrxma(double t) const
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
} /* namespace towr */

