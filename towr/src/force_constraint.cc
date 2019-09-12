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

#include <towr/constraints/force_constraint.h>

#include <towr/variables/variable_names.h>

namespace towr {


ForceConstraint::ForceConstraint (const KinematicModel::Ptr& robot_model,
                                  const HeightMap::Ptr& terrain,
                                  double force_limit,
                                  EE ee,
                                  const SplineHolder& spline
                                  )
    :ifopt::ConstraintSet(kSpecifyLater, "force-" + id::EEForceNodes(ee))
{
  max_deviation_from_nominal_ = robot_model->GetMaximumDeviationFromNominal();
  nominal_ee_pos_B_           = robot_model->GetNominalStanceInBase().at(ee);
  base_linear_  = spline.base_linear_;
  base_angular_ = EulerConverter(spline.base_angular_);
  terrain_ = terrain;
  fn_max_  = force_limit;
  mu_      = terrain->GetFrictionCoeff();
  ee_      = ee;
  n_constraints_per_node_ = 1 + 2*k2D; // positive normal force + 4 friction pyramid constraints
  //n_constraints_per_node_ = 1 + 2*k2D+6; // positive normal force + 4 friction pyramid constraints
  ee_force_node_=spline.ee_force_.at(ee);
  ee_motion_node_=spline.ee_motion_.at(ee);
  
}

void
ForceConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
{
  ee_force_  = x->GetComponent<NodesVariablesPhaseBased>(id::EEForceNodes(ee_));
  ee_motion_ = x->GetComponent<NodesVariablesPhaseBased>(id::EEMotionNodes(ee_));
  T_=ee_force_node_->GetPolyDurations();
  pure_stance_force_node_ids_ = ee_force_->GetIndicesOfNonConstantNodes();
  for (int f_node_id : pure_stance_force_node_ids_);
  int constraint_count = pure_stance_force_node_ids_.size()*n_constraints_per_node_;
  SetRows(constraint_count);
}

Eigen::VectorXd
ForceConstraint::GetValues () const
{
  VectorXd g(GetRows());
  bool sub=true;
  double time;
  int row=0;
  auto force_nodes = ee_force_->GetNodes();
  for (int f_node_id : pure_stance_force_node_ids_) 
  {

    int i;
    for (i=0; i<f_node_id+1; i++)
    {
      if (time<ee_force_node_->GetTotalTime())
      { 
        time += T_.at(i);
      }
      else {sub=false;}
    }
    if (sub)
    {time -= T_.at(f_node_id);}
    //std::cout<<i<<"  "<<f_node_id<<std::endl;
    //std::cout<<"time  "<<time<<std::endl;
    int phase  = ee_force_->GetPhase(f_node_id);

    Vector3d p = ee_motion_->GetValueAtStartOfPhase(phase); // doesn't change during stance phase
    Vector3d n = terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
    Vector3d f = force_nodes.at(f_node_id).p();
    

    // unilateral force
    g(row++) = f.transpose() * n; // >0 (unilateral forces)

    // frictional pyramid
    Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
    g(row++) = f.transpose() * (t1 - mu_*n); // t1 < mu*n
    g(row++) = f.transpose() * (t1 + mu_*n); // t1 > -mu*n

    Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());
    g(row++) = f.transpose() * (t2 - mu_*n); // t2 < mu*n
    g(row++) = f.transpose() * (t2 + mu_*n); // t2 > -mu*n
   
  }  
  std::cout<<f_polytope<<std::endl;
  std::cout<<"  "<<std::endl;
  std::cout<<d_polytope.transpose()<<std::endl;
  std::cout<<"  "<<std::endl;

  return g;
}

ForceConstraint::VecBound
ForceConstraint::GetBounds () const
{
  VecBound bounds;

  for (int f_node_id : pure_stance_force_node_ids_) {
    bounds.push_back(ifopt::Bounds(0.0, fn_max_)); // unilateral forces
    bounds.push_back(ifopt::BoundSmallerZero); // f_t1 <  mu*n
    bounds.push_back(ifopt::BoundGreaterZero); // f_t1 > -mu*n
    bounds.push_back(ifopt::BoundSmallerZero); // f_t2 <  mu*n
    bounds.push_back(ifopt::BoundGreaterZero); // f_t2 > -mu*n
    
  }

  return bounds;
}

void
ForceConstraint::FillJacobianBlock (std::string var_set,
                                    Jacobian& jac) const
{

  if (var_set == ee_force_->GetName())
   {
    int row = 0;
    for (int f_node_id : pure_stance_force_node_ids_) 
    {
      // unilateral force
      int phase   = ee_force_->GetPhase(f_node_id);
      Vector3d p  = ee_motion_->GetValueAtStartOfPhase(phase); // doesn't change during phase
      Vector3d n  = terrain_->GetNormalizedBasis(HeightMap::Normal,   p.x(), p.y());
      Vector3d t1 = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
      Vector3d t2 = terrain_->GetNormalizedBasis(HeightMap::Tangent2, p.x(), p.y());

      for (auto dim : {X,Y,Z}) 
      {
        int idx = ee_force_->GetOptIndex(NodesVariables::NodeValueInfo(f_node_id, kPos, dim));

        int row_reset=row;

        jac.coeffRef(row_reset++, idx) = n(dim);              // unilateral force
        jac.coeffRef(row_reset++, idx) = t1(dim)-mu_*n(dim);  // f_t1 <  mu*n
        jac.coeffRef(row_reset++, idx) = t1(dim)+mu_*n(dim);  // f_t1 > -mu*n
        jac.coeffRef(row_reset++, idx) = t2(dim)-mu_*n(dim);  // f_t2 <  mu*n
        jac.coeffRef(row_reset++, idx) = t2(dim)+mu_*n(dim);  // f_t2 > -mu*n
        
       }
      row += n_constraints_per_node_;
    }

}

  
  if (var_set == ee_motion_->GetName()) 
  {

    int row_for_polytope=0;
    int row = 0; double time=0;
    auto force_nodes = ee_force_->GetNodes();
    for (int f_node_id : pure_stance_force_node_ids_) 
  {
      int i;
    for (i=0; i<f_node_id+1; i++)
    {
      if (time<ee_force_node_->GetTotalTime())
      { //std::cout<<i<<"  "<<f_node_id<<std::endl;
        time += T_.at(i);
      }
    }
    time -= T_.at(0);
      int phase  = ee_force_->GetPhase(f_node_id);
      int ee_node_id = ee_motion_->GetNodeIDAtStartOfPhase(phase);

      Vector3d p = ee_motion_->GetValueAtStartOfPhase(phase); // doesn't change during pahse
      Vector3d f = force_nodes.at(f_node_id).p();
      int row_reset=0;
      for (auto dim : {X_,Y_}) 
      {
        Vector3d dn  = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal,   dim, p.x(), p.y());
        Vector3d dt1 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent1, dim, p.x(), p.y());
        Vector3d dt2 = terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Tangent2, dim, p.x(), p.y());

        int idx = ee_motion_->GetOptIndex(NodesVariables::NodeValueInfo(ee_node_id, kPos, dim));
        row_reset=row;

        // unilateral force
        jac.coeffRef(row_reset++, idx) = f.transpose()*dn;

        // friction force tangent 1 derivative
        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt1-mu_*dn);
        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt1+mu_*dn);

        // friction force tangent 2 derivative
        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt2-mu_*dn);
        jac.coeffRef(row_reset++, idx) = f.transpose()*(dt2+mu_*dn);
       
      }
      

    row +=n_constraints_per_node_;
    time=0;
 
 }

}
//
//
//
//
//if (var_set == id::EESchedule(ee_)) 
//{
//          std::cout<<"aaa"<<std::endl;
//
//  int row=0;
//  double time=0;
//    for (int f_node_id : pure_stance_force_node_ids_) 
//  {
//      int i;
//    for (i=0; i<f_node_id+1; i++)
//    {
//      if (time<ee_force_node_->GetTotalTime())
//      { //std::cout<<i<<"  "<<f_node_id<<std::endl;
//        time += T_.at(i);
//      }
//    }
//    time -= T_.at(0); double coeff=0;
//    EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(time).transpose();
//    Jacobian JacPosBWrtEE=b_R_w*ee_motion_node_->GetJacobianOfPosWrtDurations(time);
//     Vector3d vector_base_to_ee_B=ComputeBasetoEEB(time);
//    if (vector_base_to_ee_B(0)>nominal_ee_pos_B_(0)) 
//      {
//        coeff=ComputeCoeffForJacR(coeff2_, nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
//      }
//    else 
//      {
//        coeff=ComputeCoeffForJacL(coeff1_, nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
//      }
//
//    jac.middleRows(row*5+5, 1) = -coeff*JacPosBWrtEE.row(0);//+JacPosBWrtEE.row(1)+JacPosBWrtEE.row(2);
//    jac.middleRows(row*5+6, 1) = -coeff*JacPosBWrtEE.row(0);//+JacPosBWrtEE.row(1)+JacPosBWrtEE.row(2);
//    jac.middleRows(row*5+7, 1) = -coeff*JacPosBWrtEE.row(0);//+JacPosBWrtEE.row(1)+JacPosBWrtEE.row(2);
//    jac.middleRows(row*5+8, 1) = -coeff*JacPosBWrtEE.row(0);//+JacPosBWrtEE.row(1)+JacPosBWrtEE.row(2);
//    jac.middleRows(row*5+9, 1) = -coeff*JacPosBWrtEE.row(0);//+JacPosBWrtEE.row(1)+JacPosBWrtEE.row(2);
//    jac.middleRows(row*5+10, 1) = -coeff*JacPosBWrtEE.row(0);//+JacPosBWrtEE.row(1)+JacPosBWrtEE.row(2);
//    row++;
//  }
// }    
}

} /* namespace towr */
