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
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>


namespace towr {


BaseAccConstraintRangeLin::BaseAccConstraintRangeLin (double T, double dt,
                                            const NodeSpline::Ptr& spline, std::string node_variable_name) //forse non mi serve spline_holder
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueLin-")
{
 node_variables_id_=node_variable_name;
 spline_=spline;
 node_bounds_.resize(4);
 SetRows(spline_->GetPolynomialCount()*4);
}

void
BaseAccConstraintRangeLin::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{ if (k<spline_->GetPolynomialCount())
  {

    auto com=spline_->GetPoint(t);
    g.middleRows(GetRow(k, 0), 4) =FillConstraint(com);
    std::cout<<"l' acc lin è "<<com.a()<<std::endl;
  }
}

void
BaseAccConstraintRangeLin::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  if (k<spline_->GetPolynomialCount())
  {for (int dim=0; dim<node_bounds_.size(); ++dim)
   bounds.at(GetRow(k,dim)) = Bounds(-30,30);


  }
} //

void
BaseAccConstraintRangeLin::UpdateJacobianAtInstance (double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{
 if (k<spline_->GetPolynomialCount())
 {

   //if (var_set==id::base_lin_nodes)
   if (var_set ==node_variables_id_)
    {
      jac.middleRows(4*k,4)=FillJacobian(spline_,k); //spline_->GetJacobianWrtNodes(k,0.1, kAcc);

     }
 }
}
int
BaseAccConstraintRangeLin::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}
Eigen::VectorXd
BaseAccConstraintRangeLin::FillConstraint (State com) const
{
  VectorXd g;
  int mu=1;
  g(0)=com.a().x()+mu*com.a().z();
  g(1)=com.a().y()+mu*com.a().z();
  g(2)=com.a().x()-mu*com.a().z();
  g(3)=com.a().y()-mu*com.a().z();
  return g;
}
NodeSpline::Jacobian
BaseAccConstraintRangeLin::FillJacobian(NodeSpline::Ptr spline_,int k) const
{
  auto jac1=spline_->GetJacobianWrtNodes(k,0.1, kAcc);
  int mu=1;
  Jacobian g=spline_->GetJacobianWrtNodes(k,0.1, kAcc);
  std::cout<<"doing g"<<std::endl;
  std::cout<<"g 0: "<<g.row(0)<<std::endl;
  std::cout<<"jac1 0: "<<jac1.row(0)+mu*jac1.row(2)<<std::endl;
  g.row(0)=jac1.row(0)+mu*jac1.row(2);
  g.row(1)=jac1.row(1)+mu*jac1.row(2);
  g.row(2)=jac1.row(0)-mu*jac1.row(2);
  std::cout<<" problem! "<<std::endl;
  g.row(3)=jac1.row(1)-mu*jac1.row(2);
  std::cout<<"g done"<<std::endl;
  return g;

}
} /* namespace towr */

