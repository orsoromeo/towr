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

#include <towr/constraints/range_of_motion_constraint.h>
#include <towr/variables/variable_names.h>

namespace towr {

RangeOfMotionConstraint::RangeOfMotionConstraint (const KinematicModel::Ptr& model,
                                                  double T, double dt,
                                                  const EE& ee,
                                                  const SplineHolder& spline_holder,
                                                  const HeightMap::Ptr& terrain)
    :TimeDiscretizationConstraint(T, dt, "rangeofmotion-" + std::to_string(ee))
{

  base_linear_  = spline_holder.base_linear_;
  base_angular_ = EulerConverter(spline_holder.base_angular_);
  ee_motion_    = spline_holder.ee_motion_.at(ee);

  max_deviation_from_nominal_ = model->GetMaximumDeviationFromNominal();
  nominal_ee_pos_B_           = model->GetNominalStanceInBase().at(ee);
  ee_ = ee;
  terrain_ = terrain;
  theta_=37*M_PI/180;
  lenght_=0.3;
  HeightToCheck_= lenght_*sin(theta_);
  //SetRows(GetNumberOfNodes()*k3D);
  SetRows(GetNumberOfNodes()*6);
}

int
RangeOfMotionConstraint::GetRow (int node, int dim) const
{
  return node*6 + dim;
  //return node*k3D + dim;
}

void
RangeOfMotionConstraint::UpdateConstraintAtInstance (double t, int k, VectorXd& g) const
{
  Vector3d base_W  = base_linear_->GetPoint(t).p();
  Vector3d pos_ee_W = ee_motion_->GetPoint(t).p();
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();

  Vector3d vector_base_to_ee_W = pos_ee_W - base_W;
  Vector3d vector_base_to_ee_B = b_R_w*(vector_base_to_ee_W);

  g.middleRows(GetRow(k, X), k3D) = vector_base_to_ee_B;
  int a=GetRow(k,3);
  if (ee_<2) //quando scende deve essere cambiato
  {  
    Vector3d one=Vector3d(1.0, 1.0,1.0);
    g.middleRows(a,3)=one;
    //g.coeffRef(a,0) = pos_ee_W(2) + HeightToCheck_ - terrain_->GetHeight(pos_ee_W(0)+distance_, pos_ee_W(1));
  }
  else
  {
    g.coeffRef(a,0) = pos_ee_W(2) + HeightToCheck_ - terrain_->GetHeight(pos_ee_W(0)+lenght_*cos(theta_), pos_ee_W(1));
    g.coeffRef(a+1,0) = pos_ee_W(2) + HeightToCheck_/3 - terrain_->GetHeight(pos_ee_W(0)+(lenght_*cos(theta_)/3), pos_ee_W(1));
    g.coeffRef(a+2,0) = pos_ee_W(2) + HeightToCheck_*2/3 - terrain_->GetHeight(pos_ee_W(0)+(2*lenght_*cos(theta_)/3), pos_ee_W(1));
  }
  
  //std::cout<<"constraint "<<ee_<<" "<<t<<"  "<<g.coeffRef(a,0)<<std::endl;

}

void
RangeOfMotionConstraint::UpdateBoundsAtInstance (double t, int k, VecBound& bounds) const
{
  for (int dim=0; dim<k3D; ++dim) {
    ifopt::Bounds b;
    b += nominal_ee_pos_B_(dim);
  if (dim==2){
    b.upper_ +=0.04;
  }
  else{
    b.upper_ += max_deviation_from_nominal_(dim);
  }
  
  b.lower_ -= max_deviation_from_nominal_(dim);
     
 
    bounds.at(GetRow(k,dim)) = b;
  }
   bounds.at(GetRow(k,3)) = ifopt::Bounds (0.02, 100);
   bounds.at(GetRow(k,4)) = ifopt::Bounds (0.02, 100);
   bounds.at(GetRow(k,5)) = ifopt::Bounds (0.02, 100);

}

void
RangeOfMotionConstraint::UpdateJacobianAtInstance (double t, int k,
                                                   std::string var_set,
                                                   Jacobian& jac) const
{
  //double n_cols=jac.cols();
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(t).transpose();
  int row_start = GetRow(k,X);

  if (var_set == id::base_lin_nodes) {
    jac.middleRows(row_start, k3D) = -1*b_R_w*base_linear_->GetJacobianWrtNodes(t, kPos);
  }

  if (var_set == id::base_ang_nodes) {
    Vector3d base_W   = base_linear_->GetPoint(t).p();
    Vector3d ee_pos_W = ee_motion_->GetPoint(t).p();
    Vector3d r_W = ee_pos_W - base_W;
    jac.middleRows(row_start, k3D) = base_angular_.DerivOfRotVecMult(t,r_W, true);
  }

  if (var_set == id::EEMotionNodes(ee_)) {
    jac.middleRows(row_start, k3D) = b_R_w*ee_motion_->GetJacobianWrtNodes(t,kPos);
    Vector3d pos_ee_W = ee_motion_->GetPoint(t).p();
    if (ee_<2)
    {
     //jac.middleRows(GetRow(k,3),0) =ee_motion_->GetJacobianWrtNodes(t,kPos).row(2)  - GetDerivativeHeightWrtNodes(jac.cols(),t,pos_ee_W(0)+distance_,pos_ee_W(1));
    }
    else 
    {
     jac.middleRows(GetRow(k,3),0) = ee_motion_->GetJacobianWrtNodes(t,kPos).row(2)- GetDerivativeHeightWrtNodes(jac.cols(),t,pos_ee_W(0)+lenght_*cos(theta_),pos_ee_W(1));
     jac.middleRows(GetRow(k,4),0) = ee_motion_->GetJacobianWrtNodes(t,kPos).row(2)- GetDerivativeHeightWrtNodes(jac.cols(),t,pos_ee_W(0)+lenght_*cos(theta_)/3,pos_ee_W(1));
     jac.middleRows(GetRow(k,5),0) = ee_motion_->GetJacobianWrtNodes(t,kPos).row(2)- GetDerivativeHeightWrtNodes(jac.cols(),t,pos_ee_W(0)+lenght_*cos(theta_)*2/3,pos_ee_W(1));

    }
    //std::cout<<"Jacobian "<<jac.middleRows(row_start,4)<<std::endl;
  }

  if (var_set == id::EESchedule(ee_)) {
    jac.middleRows(row_start, k3D) = b_R_w*ee_motion_->GetJacobianOfPosWrtDurations(t);
  }
}

NodeSpline::Jacobian
RangeOfMotionConstraint::GetDerivativeHeightWrtNodes (double jac_cols, double t,double posx, double posy) const
{
  Jacobian jac1;
  jac1.resize(1, jac_cols);

  Jacobian DerPosWrtNodes=ee_motion_->GetJacobianWrtNodes(t,kPos);
  auto JacWrtNodesX= DerPosWrtNodes.row(0);
  auto JacWrtNodesY= DerPosWrtNodes.row(1);
  //
  double DerHeightWrtPosX;
  double DerHeightWrtPosY;
  
  DerHeightWrtPosX= terrain_->GetHeightDerivWrtX(posx,posy);
  DerHeightWrtPosY= terrain_->GetHeightDerivWrtY(posx,posy);
  

  jac1= DerHeightWrtPosX*JacWrtNodesX + DerHeightWrtPosY * JacWrtNodesY;
  

  return jac1;
}



} /* namespace xpp */

