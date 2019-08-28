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

#include <towr/costs/base_cost.h>

#include <iostream>
namespace towr {

BaseCost::BaseCost (const std::string& nodes_id, Dx deriv, int dim, double weight,const SplineHolder& spline_holder)

    :CostTerm(nodes_id
               +"-dx_"+std::to_string(deriv)
               +"-dim_"+std::to_string(dim))
{

  base_lin_=spline_holder.base_linear_;
  base_ang_ = EulerConverter(spline_holder.base_angular_);
  spline_holder_=spline_holder;


}

//Eigen::VectorXd
//BaseCost::GetValues () const
//{
//  
//  
//}
  void BaseCost::InitVariableDependedQuantities(const VariablesPtr& x)

{
  SetRows(1);
}

double BaseCost::GetCost () const{

  VectorXd g = VectorXd::Zero(1);

  double cost_value=0.0;
  double t=0.0;
  for (double i=0; i<base_lin_->GetTotalTime()/0.1 + 1e-6; i++)
  {
    t=0.1*i; 
    auto tmp=base_lin_->GetPoint(t).v();
    auto tmp1=base_ang_.GetAngularVelocityInWorld(t);

    //auto tmp1=base_ang_.GetAngularAccelerationInWorld(t);
    for (int p=0; p<3; p++)
    {
       cost_value = cost_value + fabs(tmp1(p));// + fabs(tmp1(p));
       ///cost_value = cost_value +pow(tmp(p),2)+pow(tmp1(p),2);// 
      
    }
  }

  return cost_value;
}

void BaseCost::FillJacobianBlock(std::string var_set, Jacobian& jac) const 
  { 
  
    Eigen::MatrixXd intermediate_sum;
    intermediate_sum.resize(jac.rows(),jac.cols());
    intermediate_sum.setZero();
//    {
//    
//     for (double i=0; i<base_lin_->GetTotalTime()/0.1 + 1e-6; i++)
//     {
//      double t=0.1*i;
//      Eigen::Vector3d lin_vel=base_lin_->GetPoint(t).v();
//      Jacobian jac_base_lin_vel=base_lin_->GetJacobianWrtNodes(t,kVel);
//      for (int p=0; p<3; p++)
//      {
//        Eigen::MatrixXd tmp;
//        tmp.resize(1,jac.cols());
//        tmp=2*lin_vel(p)*jac_base_lin_vel.row(p);
//        if (lin_vel(p)<0)
//        {
//          tmp= -tmp;
//        }
//        //tmp=2.0*lin_acc(p)*jac_base_lin_acc.row(p);
//        intermediate_sum = intermediate_sum + tmp;
//      }
//      //std::cout<<"Get Point "<<std::endl;
//      //std::cout<<base_lin_->GetPoint(t).v()<<std::endl;
//      //std::cout<<"JacobianWrtNodes "<<std::endl;
//      //std::cout<<Eigen::MatrixXd(base_lin_->GetJacobianWrtNodes(t,kVel))<<std::endl;
//      //std::cout<<"Sum "<<std::endl;
//      //std::cout<<intermediate_sum<<std::endl;
//
//      }
//     jac=intermediate_sum.sparseView();
//    } 
//    
    if (var_set == id::base_ang_nodes) 
    {

     for (double i=0; i<base_lin_->GetTotalTime()/0.1 + 1e-6; i++)
     {
      double t=0.1*i;
      Eigen::Vector3d ang_vel=base_ang_.GetAngularVelocityInWorld(t);
      Jacobian jac_base_ang_vel=base_ang_.GetDerivOfAngVelWrtEulerNodes(t);
           // Jacobian jac_base_ang_acc=base_ang_.GetDerivOfAngAccWrtEulerNodes(t);

      for (int p=0; p<3; p++)
      {
        Eigen::MatrixXd tmp;
        tmp.resize(1,jac.cols());
        tmp=jac_base_ang_vel.row(p);
        if (ang_vel(p)<0)
        {
          tmp= -tmp;
        }
        //tmp=2*ang_vel(p)*jac_base_ang_vel.row(p);
        intermediate_sum= intermediate_sum+tmp ;
      }
    }

    jac=intermediate_sum.sparseView();
  }


 }
} /* namespace towr */
