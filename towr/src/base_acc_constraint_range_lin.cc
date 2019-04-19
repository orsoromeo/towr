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
#include <iostream>


namespace towr {


BaseAccConstraintRangeLin::BaseAccConstraintRangeLin (double T, double dt,
                                            const NodeSpline::Ptr& spline, std::string node_variable_name) //forse non mi serve spline_holder
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueLin-")
{
 node_variables_id_=node_variable_name;
 spline_=spline;
 node_bounds_.resize(k3D);
 SetRows(spline_->GetPolynomialCount()*k3D);
}

void
BaseAccConstraintRangeLin::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{ if (k<spline_->GetPolynomialCount())
  {

    auto com=spline_->GetPoint(t);
    g.middleRows(GetRow(k, 0), k3D) = com.a();
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
      jac.middleRows(3*k,3)=spline_->GetJacobianWrtNodes(k,0.1, kAcc);

     }
 }
}
int
BaseAccConstraintRangeLin::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}

} /* namespace towr */

