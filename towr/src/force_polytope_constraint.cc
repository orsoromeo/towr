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

#include <towr/constraints/force_polytope_constraint.h>

#include <towr/variables/variable_names.h>

namespace towr {


ForcePolytopeConstraint::ForcePolytopeConstraint (const KinematicModel::Ptr& robot_model,
                                  const HeightMap::Ptr& terrain,
                                  double force_limit,
                                  EE ee,
                                  const SplineHolder& spline
                                  )
    :ifopt::ConstraintSet(kSpecifyLater, "force-" + id::EEForceNodes(ee))
{
  max_deviation_from_nominal_ = robot_model->GetMaximumDeviationFromNominal();
  nominal_ee_pos_B_           = robot_model->GetNominalStanceInBase().at(ee);
  base_to_hip_distance           = robot_model->GetDistanceBaseToHip().at(ee);
  base_linear_  = spline.base_linear_;
  base_angular_ = EulerConverter(spline.base_angular_);
  terrain_ = terrain;
  fn_max_  = force_limit;
  mu_      = terrain->GetFrictionCoeff();
  ee_      = ee;
  n_constraints_per_node_ = 6;
  ee_force_node_=spline.ee_force_.at(ee);
  ee_motion_node_=spline.ee_motion_.at(ee);
  InitializeQuantities(robot_model,ee);
}

void
ForcePolytopeConstraint::InitVariableDependedQuantities (const VariablesPtr& x)
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
ForcePolytopeConstraint::GetValues () const
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

    Vector3d vector_base_to_ee_B=ComputeBasetoEEB(time);
    double baseToHipX = vector_base_to_ee_B(0) - base_to_hip_distance(0);
    double baseToHipZ = vector_base_to_ee_B(2) - base_to_hip_distance(2);
    std::cout<<baseToHipX<<" "<<baseToHipZ<<std::endl;
    double squaredHip2FootDistance = pow(baseToHipX, 2) + pow(baseToHipZ, 2);
    double hip2FootDistance = sqrt(squaredHip2FootDistance);
    double legPitchAngle = atan2(baseToHipZ, baseToHipX);

    double right_side=nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0);
    double left_side=nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0);
    int j;
    double thetax=0;
    for (j=0; j<4; j++)
    {

       if (vector_base_to_ee_B(0)>nominal_ee_pos_B_(0)) 
     {
      thetax=ComputeBoundR(coeffN_(j),coeffR_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), right_side); 
      //f_polytope(0,j)=cos(thetax);
      //f_polytope(2,j)=sin(thetax);
      d_polytope(j)=ComputeBoundR(coeffDN_(j),coeffDR_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), right_side);
     }
      else 
    {
       thetax=ComputeBoundL(coeffL_(j), coeffN_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), left_side); 
      //f_polytope(0,j)=cos(thetax);
      //f_polytope(2,j)=sin(thetax);
       d_polytope(j)=ComputeBoundL(coeffDL_(j), coeffDN_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), left_side);
    }
   }
  
  //f_polytope(0,j)=cos(thetax);
  //f_polytope(2,j)=sin(thetax);
  f_polytope(0,j)=cos(legPitchAngle)*cos(thetax) - sin(legPitchAngle)*sin(thetax); //f_polytope(0,j)=cos(thetax + legPitchAngle);
  f_polytope(2,j)=sin(legPitchAngle)*cos(thetax) + cos(legPitchAngle)*sin(thetax); //f_polytope(2,j)=cos(thetax + legPitchAngle);

    f_polytope(1,4)=1;
    f_polytope(1,5)=-1;
    
    
    d_polytope(4)= 90;
    d_polytope(5)= 90;

    int plane; 
    for (plane=0; plane<6; plane++)
     {
      g(row++) = f.transpose() * f_polytope.col(plane)-d_polytope(plane);
     }
    
  time=0;
  }  
  std::cout<<f_polytope<<std::endl;
  std::cout<<"  "<<std::endl;
  std::cout<<d_polytope.transpose()<<std::endl;
  std::cout<<"  "<<std::endl;

  return g;
}

ForcePolytopeConstraint::VecBound
ForcePolytopeConstraint::GetBounds () const
{
  VecBound bounds;

  for (int f_node_id : pure_stance_force_node_ids_) {
    
    bounds.push_back(ifopt::BoundSmallerZero);
    bounds.push_back(ifopt::BoundSmallerZero);
    bounds.push_back(ifopt::BoundSmallerZero);
    bounds.push_back(ifopt::BoundSmallerZero);
    bounds.push_back(ifopt::BoundSmallerZero);
    bounds.push_back(ifopt::BoundSmallerZero);
  }

  return bounds;
}

void
ForcePolytopeConstraint::FillJacobianBlock (std::string var_set,
                                    Jacobian& jac) const
{

  if (var_set == id::base_lin_nodes) 
  {
    double time=0;
    int row=0;
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
    EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(time).transpose();

    Jacobian JacPosBWrtBaseLin= -1*b_R_w*base_linear_->GetJacobianWrtNodes(time, kPos);
    Vector3d vector_base_to_ee_B=ComputeBasetoEEB(time);
    Vector3d f = force_nodes.at(f_node_id).p();
    double coeff1=0; double coeff2=0; double thetax=0;
    int j;
    int row_reset=row;
    for (j=0; j<4; j++)
   {
    if (vector_base_to_ee_B(0)>nominal_ee_pos_B_(0)) 
      {

        thetax=ComputeBoundR(coeffN_(j),coeffR_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0)); 
        coeff1=ComputeCoeffForJacR(coeffN_(j),coeffR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
        coeff2=ComputeCoeffForJacR(coeffDN_(j),coeffDR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
        
      }
    else 
    {

      thetax=ComputeBoundL(coeffL_(j), coeffN_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0)); 
      coeff1=ComputeCoeffForJacL(coeffL_(j),coeffN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
      coeff2=ComputeCoeffForJacL(coeffDL_(j),coeffDN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));

    }


    jac.middleRows(row_reset++, 1) = ((-sin(thetax)*f(0)+cos(thetax)*f(2))*coeff1-coeff2)*JacPosBWrtBaseLin.row(0); //+ JacPosBWrtBaseLin.row(1)+JacPosBWrtBaseLin.row(2);
    std::cout<<"JacPosBWrtBaseLin row 0"<<JacPosBWrtBaseLin.row(0)<<std::endl;
    std::cout<<"jac.middleRows(row_reset++, 1)"<<JacPosBWrtBaseLin.row(0)<<std::endl;
  }

    row += n_constraints_per_node_;
    time=0;


 }
   
}


  if (var_set == id::base_ang_nodes) 
  {
   

    int row=0; double time=0;
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
    time -= T_.at(0); double coeff=0;
    Vector3d base_W   = base_linear_->GetPoint(time).p();
    Vector3d ee_pos_W = ee_motion_node_->GetPoint(time).p();
    Vector3d r_W = ee_pos_W - base_W;
    Jacobian JacPosBWrtBaseAng=base_angular_.DerivOfRotVecMult(time,r_W, true);
     
    Vector3d vector_base_to_ee_B=ComputeBasetoEEB(time);
    Vector3d f = force_nodes.at(f_node_id).p();

   double coeff1=0.0; double coeff2=0.0; double thetax=0.0;
    int j;
    int row_reset=row;

    for (j=0; j<4; j++)
   {
    if (vector_base_to_ee_B(0)>nominal_ee_pos_B_(0)) 
      {
        thetax=ComputeBoundR(coeffN_(j),coeffR_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0)); 
        coeff1=ComputeCoeffForJacR(coeffN_(j),coeffR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
        coeff2=ComputeCoeffForJacR(coeffDN_(j),coeffDR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
            
      }
    else 
    {
      thetax=ComputeBoundL(coeffL_(j), coeffN_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0)); 
      coeff1=ComputeCoeffForJacL(coeffL_(j),coeffN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
      coeff2=ComputeCoeffForJacL(coeffDL_(j),coeffDN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
      std::cout<<" coeff0 "<<coeffL_(j)<<"  "<<" coeff1 "<<coeffN_(j)<<std::endl;
      std::cout<<""<<std::endl;
      std::cout<<"coeff1 "<<coeff1<<std::endl;
          

    }
    std::cout<<sin(thetax);
    std::cout<<" "<<std::endl;
    std::cout<<f.transpose();
    std::cout<<" "<<std::endl;
    std::cout<<cos(thetax);
    std::cout<<" "<<std::endl;

    std::cout<<(-sin(thetax)*f(0)+cos(thetax)*f(2))*coeff1+coeff2<<std::endl;
    jac.middleRows(row_reset++, 1) = ((-sin(thetax)*f(0)+cos(thetax)*f(2))*coeff1-coeff2)*JacPosBWrtBaseAng.row(0); //+ JacPosBWrtBaseLin.row(1)+JacPosBWrtBaseLin.row(2);
    
  }
    
    row += n_constraints_per_node_;
    time=0;
   

 }
}
//}
//
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
        
        int plane=0;
        for (plane=0; plane<6; plane++)
         {
          jac.coeffRef(row_reset++,idx) = f_polytope(dim,plane);
         }
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
      
      EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(time).transpose();

      Jacobian JacPosBWrtNodes=b_R_w*ee_motion_node_->GetJacobianWrtNodes(time,kPos);
      Vector3d vector_base_to_ee_B=ComputeBasetoEEB(time);
      double coeff1=0; double coeff2=0; double thetax=0;
    int j;
    for (j=0; j<4; j++)
   {
    if (vector_base_to_ee_B(0)>nominal_ee_pos_B_(0)) 
      {
        thetax=ComputeBoundR(coeffN_(j),coeffR_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0)); 
        coeff1=ComputeCoeffForJacR(coeffN_(j),coeffR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
        coeff2=ComputeCoeffForJacR(coeffDN_(j),coeffDR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
        
      }
    else 
    {
      thetax=ComputeBoundL(coeffL_(j), coeffN_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0)); 
      coeff1=ComputeCoeffForJacL(coeffL_(j),coeffN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
      coeff2=ComputeCoeffForJacL(coeffDL_(j),coeffDN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));

    }


    jac.middleRows(row_reset++, 1) = ((-sin(thetax)*f(0)+cos(thetax)*f(2))*coeff1-coeff2)*JacPosBWrtNodes.row(0); //+ JacPosBWrtBaseLin.row(1)+JacPosBWrtBaseLin.row(2);
    
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
//
double ForcePolytopeConstraint::ComputeBoundL (double coeff0, double coeff1,  double Posx, double Pn,double ls) const
{
  double bound;
  bound=(Posx-ls)*(coeff1-coeff0)/(Pn-ls)+coeff0;
  return bound;
}
double ForcePolytopeConstraint::ComputeBoundR (double coeff0, double coeff1, double Posx, double Pn,double rs) const
{
  double bound;
  bound=(Posx-Pn)*(coeff1-coeff0)/(rs-Pn)+coeff0;
  return bound;
}
double ForcePolytopeConstraint::ComputeCoeffForJacL (double coeff0, double coeff1, double Pn,double ls) const
{
  double bound;

  bound=(coeff1-coeff0)/(Pn-ls);
  return bound;
}
double ForcePolytopeConstraint::ComputeCoeffForJacR (double coeff0, double coeff1, double Pn,double rs) const
{
  double bound;
  bound=(coeff1-coeff0)/(rs-Pn);
  return bound;
}

Eigen::Vector3d ForcePolytopeConstraint::ComputeBasetoEEB (double time) const
{
 Vector3d base_W  = base_linear_->GetPoint(time).p();
  Vector3d pos_ee_W = ee_motion_node_->GetPoint(time).p();
  EulerConverter::MatrixSXd b_R_w = base_angular_.GetRotationMatrixBaseToWorld(time).transpose();

  Vector3d vector_base_to_ee_W = pos_ee_W - base_W;
  Vector3d vector_base_to_ee_B = b_R_w*(vector_base_to_ee_W);
  return vector_base_to_ee_B;
}
void ForcePolytopeConstraint::InitializeQuantities (const KinematicModel::Ptr& robot_model, double ee) 
{
  
  coeffL_.resize(4);
  Eigen::MatrixXd pp=robot_model->GetThetaL();
  coeffL_=robot_model->GetThetaL().col(ee);
  coeffN_.resize(4);
  coeffN_=robot_model->GetThetaN().col(ee);
  coeffR_.resize(4);
  coeffR_=robot_model->GetThetaR().col(ee);
  
  coeffDL_.resize(4);
  coeffDL_=robot_model->GetDL().col(ee);
  coeffDN_.resize(4);
  coeffDN_=robot_model->GetDN().col(ee);
  coeffDR_.resize(4);
  coeffDR_=robot_model->GetDR().col(ee);


  f_polytope.resize(3,6);
  f_polytope.setZero();

  d_polytope.resize(6);
  d_polytope.setZero();

}
} /* namespace towr */
