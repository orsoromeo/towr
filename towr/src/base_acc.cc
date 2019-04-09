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

#include <towr/constraints/base_acc.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>
#include <towr/variables/spline_holder.h>

namespace towr {


BaseAccConstraintRange::BaseAccConstraintRange (double T, double dt,
                                            const SplineHolder& spline_holder) //forse non mi serve spline_holder
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValue")
{
  base_linear_  = spline_holder.base_linear_;
  base_angular_ = spline_holder.base_angular_;

  SetRows(GetNumberOfNodes()*k3D);
}

void
BaseAccConstraintRange::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{
  auto com = base_linear_->GetPoint(t);
  auto com1 = base_angular_->GetPoint(t);
  g.middleRows(GetRow(k, LX), k3D) = com.a();
  g.middleRows(GetRow(k, AX), k3D) = com1.a();
  //Eigen::Matrix<int,3,1> l;
  //l << 1,1,1;
  //g.middleRows(GetRow(k, AX), k3D)=l.transpose();
}

void
BaseAccConstraintRange::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<node_bounds_.size(); ++dim)
    bounds.at(GetRow(k,dim)) = (0,100);
}

void
BaseAccConstraintRange::UpdateJacobianAtInstance (double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{
  if (var_set == id::base_ang_nodes)
    jac.middleRows(GetRow(k,AX), k3D) = base_angular_->GetJacobianWrtNodes(t, kPos);

  if (var_set == id::base_lin_nodes)
    jac.middleRows(GetRow(k,LX), k3D) = base_linear_->GetJacobianWrtNodes(t, kPos);
}

int
BaseAccConstraintRange::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}

} /* namespace towr */

