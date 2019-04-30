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

#include <towr/variables/euler_converter.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>
#include <iostream>



namespace towr {


BaseAccConstraintRangeAng::BaseAccConstraintRangeAng (double T, double dt,
                                            const NodeSpline::Ptr& spline, std::string node_variable_name) //forse non mi serve spline_holder
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueAng-")
{
 spline_=spline;
 node_variables_id_=node_variable_name;
 base_angular_ = EulerConverter(spline);
 node_bounds_.resize(4);
 SetRows(spline_->GetPolynomialCount()*4);
}

void
BaseAccConstraintRangeAng::UpdateConstraintAtInstance (double t,int k,
                                                  VectorXd& g) const
{ if (k<spline_->GetPolynomialCount())
  {

   auto com = base_angular_.GetAngularAccelerationInWorld(t);
   g.middleRows(GetRow(k, AX), 4) = FillConstraint(com);
   //std::cout<<"l' acc ang è "<<com<<std::endl;

  }
}

void
BaseAccConstraintRangeAng::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  if (k<spline_->GetPolynomialCount())
  { bounds.at(GetRow(k,0)) = ifopt::BoundGreaterZero;
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
 if (k<spline_->GetPolynomialCount())
 {

   //if (var_set==id::base_lin_nodes)
   //if (var_set ==node_variables_id_)
   // {
   //  auto jac1=spline_->base_angular_.GetDerivOfAngAccWrtEulerNodes(t);
   //  jac.middleRows(3*k,3)=FillJacobian(jac1);
   //}
 }
}
int
BaseAccConstraintRangeAng::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}

Eigen::VectorXd
BaseAccConstraintRangeAng::FillConstraint (Eigen::VectorXd com) const
{
  Eigen::VectorXd g;
  g.resize(4);
  //double g[4];
  double mu=1.0;
  g(0)=com(0)+mu*com(2);//>0
  g(1)=com(1)+mu*com(2);//>0
  g(2)=com(0)-mu*com(2);//<0
  g(3)=com(1)-mu*com(2);//<0
  return g;
}

NodeSpline::Jacobian
BaseAccConstraintRangeAng::FillJacobian(Jacobian jac1) const
{
  double mu=1.0;
  NodeSpline::Jacobian g=NodeSpline::Jacobian(4, jac1.cols());
  g.row(0)=jac1.row(0)+mu*jac1.row(2);
  g.row(1)=jac1.row(1)+mu*jac1.row(2);
  g.row(2)=jac1.row(0)-mu*jac1.row(2);
  g.row(3)=jac1.row(1)-mu*jac1.row(2);

  return g;
}
} /* namespace towr */
