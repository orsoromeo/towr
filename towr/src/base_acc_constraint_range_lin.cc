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
#include <towr/constraints/dynamic_constraint.h>
#include <towr/terrain/height_map.h>
#include <towr/variables/nodes_variables_phase_based.h>
#include <towr/variables/spline_holder.h>
#include <towr/variables/variable_names.h>
#include <towr/variables/cartesian_dimensions.h>

#include <towr/models/single_rigid_body_dynamics.h>
#include <towr/models/dynamic_model.h>
#include <iostream>
#include <cmath>
#define PI 3.14159265

namespace towr {


BaseAccConstraintRangeLin::BaseAccConstraintRangeLin (const DynamicModel::Ptr& model, double T, double dt,
                                                      const NodeSpline::Ptr& spline, std::string node_variable_name, HeightMap::Ptr terrain, const SplineHolder& spline_holder)
    :TimeDiscretizationConstraint(T, dt, "baseAccConstraintValueLin-")

{
 model_=model;
 node_variables_id_=node_variable_name;
 spline_=spline;
 node_bounds_.resize(4);
 NumberNodes_=(spline_->GetPolynomialCount()+1);
 SetRows(NumberNodes_*4);
 terrain_=terrain;
 ee_motion_=spline_holder.ee_motion_;
 mu_= terrain->GetFrictionCoeff();
}

void
BaseAccConstraintRangeLin::UpdateConstraintAtInstance (double t, int k,
                                                  VectorXd& g) const
{

  if (k<NumberNodes_)
  {
    auto com=spline_->GetPoint(t);
    g.middleRows(GetRow(k, 0), 4) =FillConstraint(com);
    //std::cout<<"l' acc lin è "<<com.a()<<std::endl;
  }
}

void
BaseAccConstraintRangeLin::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  if (k<NumberNodes_)
  {
   bounds.at(GetRow(k,0)) = ifopt::BoundGreaterZero;
   bounds.at(GetRow(k,1)) = ifopt::BoundGreaterZero;
   bounds.at(GetRow(k,2)) = ifopt::BoundSmallerZero;
   bounds.at(GetRow(k,3)) = ifopt::BoundSmallerZero;
  }
} //

void
BaseAccConstraintRangeLin::UpdateJacobianAtInstance (double t, int k,
                                                std::string var_set,
                                                Jacobian& jac) const
{

   if (k<NumberNodes_)
   {

   //if (var_set==id::base_lin_nodes)
   if (var_set ==node_variables_id_)
    {

      jac.middleRows(4*k,4)=FillJacobian(spline_,k);

      //spline_->GetJacobianWrtNodes(k,0.1, kAcc);
      //std::cout<<"lin "<<std::endl;
      //std::cout<<jac<<std::endl;

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
  g.resize(4);
  //double g[4];
  double mu=10.0;
  g(0)=com.a().x()+mu*(com.a().z()+model_->g());//>0
  g(1)=com.a().y()+mu*(com.a().z()+model_->g());//>0
  g(2)=com.a().x()-mu*(com.a().z()+model_->g());//<0. gravità
  g(3)=com.a().y()-mu*(com.a().z()+model_->g());//<0
  return g;
}
NodeSpline::Jacobian
BaseAccConstraintRangeLin::FillJacobian(NodeSpline::Ptr spline_,int k) const
{
  double a;
  int pol;
  if (k==0)
      {
        a=0;
        pol=k;
      }
  else
       { a=0.1;
         pol=k-1;
       }
  auto jac1=spline_->GetJacobianWrtNodes(pol,a, kAcc);

  double mu=10.0;
  Jacobian g=Jacobian(4, jac1.cols());
  g.row(0)=jac1.row(0)+mu*jac1.row(2);
  g.row(1)=jac1.row(1)+mu*jac1.row(2);
  g.row(2)=jac1.row(0)-mu*jac1.row(2);
  g.row(3)=jac1.row(1)-mu*jac1.row(2);
  
  return g;

}
// questa parte viene da cone.cc
Eigen::MatrixXd BaseAccConstraintRangeLin::ReturnNormalTerrain (double t)
{

  Eigen::MatrixXd n(4,3);
  for (int ee=0; GetNumberOfFeetInTouch (t); ee++)
  {Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
   n.row(ee)= terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
  }

 return n;
}

Eigen::MatrixXd BaseAccConstraintRangeLin::ReturnTangentTerrain (double t)
{
 Eigen::MatrixXd tang (4,3);
 for (int ee=0; GetNumberOfFeetInTouch (t); ee++)
  {   Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
      tang.row(ee) = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
  }
 return tang;
}

Eigen::Matrix3d BaseAccConstraintRangeLin::RotationMatrix (Eigen::Vector3d t, double angle) const

{
  double cmu=cos(angle);
  double smu=sin(angle);
  Eigen::Matrix3d M;
  M(0,0)=std::pow(t(0),2)*(1-cmu)+cmu;          M(0,1)=t(0)*t(1)*(1-cmu)-t(2)*smu;      M(0,2)=t(0)*t(2)*(1-cmu)+t(1)*smu;
  M(1,0)=t(0)*t(1)*(1-cmu)+t(2)*smu;            M(1,1)=std::pow(t(1),2)*(1-cmu)+cmu;    M(1,2)=t(1)*t(2)*(1-cmu)-t(0)*smu;
  M(2,0)=t(0)*t(2)*(1-cmu)-t(1)*smu;            M(2,1)=t(1)*t(2)*(1-cmu)-t(0)*smu;      M(2,2)=std::pow(t(2),2)*(1-cmu)+cmu;
  return M;
}

Eigen::Vector3d BaseAccConstraintRangeLin::GetEdge (Eigen::Matrix3d M, Eigen::Vector3d n) const
{
  Eigen::Vector3d edge=M*n;
  return edge;
}
void BaseAccConstraintRangeLin::GetAllEdges (double t)
{
  Eigen::MatrixXd E(4*GetNumberOfFeetInTouch (t),3);
  for (int ee=0; GetNumberOfFeetInTouch (t); ee++)
  {
   int f=ee*4;
   Eigen::Vector3d normal=ReturnTangentTerrain(t).row(ee);
   Eigen::Vector3d tang=ReturnNormalTerrain(t).row(ee);
   Eigen::Matrix3d M=RotationMatrix(tang,mu_);
   Eigen::Matrix3d M1=RotationMatrix(normal,PI/4);

   E.row(f++)=GetEdge(M,normal);
   E.row(f++)=GetEdge(M1,E.row(f-1));
   E.row(f++)=GetEdge(M1,E.row(f-1));
   E.row(f++)=GetEdge(M1,E.row(f-1));
  }
  for (int i=0; model_->GetEECount(); i++)
  {
    //e se non sono tutti in contatto che fai di e4?
    e1.row(i)=E.row(i);
    e2.row(i)=E.row(i+4);
    e3.row(i)=E.row(i+8);
    e4.row(i)=E.row(i+12);
  }
}
int BaseAccConstraintRangeLin::GetNumberOfFeetInTouch (double t)
{

  for (int ee=0; model_->GetEECount(); ee++)
   { Eigen::Vector3d p = ee_motion_.at(ee)->GetPoint(t).p();
     if (p.z()=terrain_->GetHeight(p.x(),p.y()))
      {
       ee_motion_in_touch_.push_back(ee_motion_.at(ee));
      }
  return ee_motion_in_touch_.size();
  }
}

//void BaseAccConstraintRangeLin::GetAllEdges (double t)
//{
//  for (int ee=0; model_->GetEECount(); ee++)
//  {
//   Eigen::Vector3d normal=ReturnTangentTerrain(t).row(ee);
//   Eigen::Vector3d tang=ReturnNormalTerrain(t).row(ee);
//   Eigen::Matrix3d M=RotationMatrix(tang,mu_);
//   Eigen::Matrix3d M1=RotationMatrix(normal,PI/4);
//   e1.row(ee)=GetEdge(M,normal);
//   e2.row(ee)=GetEdge(M1,e1.row(ee));
//   e3.row(ee)=GetEdge(M1,e2.row(ee));
//   e4.row(ee)=GetEdge(M1,e3.row(ee));
//
//  }
//}
} /* namespace towr */

