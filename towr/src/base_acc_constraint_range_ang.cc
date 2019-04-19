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
 node_bounds_.resize(k3D);
 SetRows(spline_->GetPolynomialCount()*k3D);
}

void
BaseAccConstraintRangeAng::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{ if (k<spline_->GetPolynomialCount())
  {

   auto com = base_angular_.GetAngularAccelerationInWorld(t);
   g.middleRows(GetRow(k, AX), k3D) = com;
   //std::cout<<"l' acc ang è "<<com<<std::endl;

  }
}

void
BaseAccConstraintRangeAng::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  if (k<spline_->GetPolynomialCount())
  {for (int dim=0; dim<node_bounds_.size(); ++dim)
   bounds.at(GetRow(k,dim)) = Bounds(-30,30);
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
   if (var_set ==node_variables_id_)
    {
    jac.middleRows(3*k,3)=base_angular_.GetDerivOfAngAccWrtEulerNodes(t);

    //auto com=spline_->GetPoint(t);
    //auto com1=spline_->GetPoint(t+1);
    //accelerazione x
    // jac.coeffRef(3*k,k+1)=sin(com.p().y())*cos(com.p().z())*com.v().x()/0.1;
    // jac.coeffRef(3*k,k+2)=cos(com.p().z())*com.v().y()/0.1;
    // jac.coeffRef(3*k,k+3)=-cos(com.p().y())*cos(com.p().z())/0.1;
    // jac.coeffRef(3*k,k+4)=sin(com.p().z())/0.1;
    // jac.coeffRef(3*k,k+7)=-sin(com1.p().y())*cos(com1.p().z())*com1.v().x()/0.1;
    // jac.coeffRef(3*k,k+8)=-cos(com1.p().z())*com1.v().y()/0.1;
    // jac.coeffRef(3*k,k+9)=cos(com1.p().y())*cos(com1.p().z())/0.1;
    // jac.coeffRef(3*k,k+10)=-sin(com.p().z())/0.1;


     // jac.coeffRef(3*k,6*k+3)=-1/0.1;
     // jac.coeffRef(3*k,6*(k+1)+3)=1/0.1;
     // jac.coeffRef(3*k+1,6*k+4)=-1/0.1;
     // jac.coeffRef(3*k+1,6*(k+1)+4)=1/0.1;
     // jac.coeffRef(3*k+2,6*k+5)=-1/0.1;
     // jac.coeffRef(3*k+2,6*(k+1)+5)=1/0.1;

    }

  //if (var_set == id::base_ang_nodes)
  // {
  //     jac.coeffRef(3*k,6*k+3)=-1/0.1;
  //     jac.coeffRef(3*k,6*(k+1)+3)=1/0.1;
  //     jac.coeffRef(3*k+1,6*k+4)=-1/0.1;
  //     jac.coeffRef(3*k+1,6*(k+1)+4)=1/0.1;
  //     jac.coeffRef(3*k+2,6*k+5)=5;
  //     jac.coeffRef(3*k+2,6*(k+1)+5)=1/0.1;
  //
  //     //
  //  }

//in realtà per il la parte angolare sto calcolando la derivata seconda degli angoli. Devo
//quindi lanciare il componente solo una volta, passandogli tutta la spline perchè devo
//differenziare i due jacobiani. Bisognerà cambiare anche il costruttore

 }
}
int
BaseAccConstraintRangeAng::GetRow (int node, int dim) const
{
  return node*node_bounds_.size() + dim;
}

} /* namespace towr */
