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
#include <math.h>       /* fabs */
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
    std::cout<<"base to hip distance x y are"<<baseToHipX<<" "<<baseToHipZ<<std::endl;
    double squaredHip2FootDistance = pow(baseToHipX, 2) + pow(baseToHipZ, 2);
    double hip2FootDistance = sqrt(squaredHip2FootDistance);

    std::cout<<"base to hip radius is"<<hip2FootDistance<<std::endl;
    double legPitchAngle = 0.0;//atan2(baseToHipX, baseToHipZ);
    std::cout<<"legs pitch angle is"<<legPitchAngle<<std::endl;

    double nominalHip2FootDistance = fabs(nominal_ee_pos_B_(2)); // nominalHip2FootDistance must always be positive
    double max_extension=nominalHip2FootDistance + max_deviation_from_nominal_(2);
    double max_retraction=nominalHip2FootDistance - max_deviation_from_nominal_(2);

    std::cout<<"nominalHip2FootDistance"<<nominalHip2FootDistance<<std::endl;
    std::cout<<"max_extension"<<max_extension<<std::endl;
    std::cout<<"max_retraction"<<max_retraction<<std::endl;
    
    double thetax=0;

    int halfspacesNumber = 4;
    for (int j=0; j<halfspacesNumber; j++)
    {

      if (hip2FootDistance>=fabs(nominalHip2FootDistance)) {
        std::cout<<" hip2FootDistance is in the lower half (leg is extended) "<<std::endl;

        //double ForcePolytopeConstraint::ComputeBound(double coeff0, double coeff1, double Posx, double P0, double P1) const
        //bound=(Posx-P0)/(P1-P0)*(coeff1-coeff0)+coeff0;

       thetax=ComputeBound(coeffN_(j), coeffL_(j), hip2FootDistance, nominalHip2FootDistance, max_extension);
        //bound=(hip2FootDistance-nominalHip2FootDistance)/(max_extension-nominalHip2FootDistance)*(coeffR_-coeffN_)+coeffN_;
       d_polytope(j)=ComputeBound(coeffDN_(j), coeffDL_(j), hip2FootDistance, max_extension, nominalHip2FootDistance);
      }

     if (hip2FootDistance<fabs(nominalHip2FootDistance)){
          std::cout<<" hip2FootDistance is in the upper half (leg is retracted)"<<std::endl;
       //bound=(hip2FootDistance-max_retraction)/(nominalHip2FootDistance-max_retraction)*(coeffN_- coeffL_)+coeffL_;
        thetax=ComputeBound(coeffR_(j), coeffN_(j), hip2FootDistance, max_retraction, nominalHip2FootDistance);
        d_polytope(j)=ComputeBound(coeffDR_(j), coeffDN_(j), hip2FootDistance, max_retraction, nominalHip2FootDistance);
      
      }

//      if ((hip2FootDistance>=fabs(nominalHip2FootDistance))&&(hip2FootDistance<=fabs(max_extension))) {
//        std::cout<<" hip2FootDistance is in the lower half (leg is extended) "<<std::endl;
//       thetax=ComputeBoundL(coeffL_(j), coeffN_(j), hip2FootDistance, nominalHip2FootDistance, max_retraction);
//       //ForcePolytopeConstraint::ComputeBoundL (double coeff0, double coeff1,  double Posx, double Pn,double ls) const
//       //bound=(Posx-ls)*(coeff1-coeff0)/(Pn-ls)+coeff0;
//       d_polytope(j)=ComputeBoundL(coeffDL_(j), coeffDN_(j), hip2FootDistance, nominalHip2FootDistance, max_retraction);
//      }
//
//     if ((hip2FootDistance<fabs(nominalHip2FootDistance))&&(hip2FootDistance>=fabs(max_retraction))){
//          std::cout<<" hip2FootDistance is in the upper half (leg is retracted)"<<std::endl;
//      
//      thetax=ComputeBoundR(coeffN_(j), coeffR_(j), hip2FootDistance, nominalHip2FootDistance, max_extension); 
//      d_polytope(j)=ComputeBoundR(coeffDN_(j),coeffDR_(j), hip2FootDistance, nominalHip2FootDistance, max_extension);
//
//      }
//
//      if(hip2FootDistance<fabs(max_retraction)){
//          std::cout<<" hip2FootDistance reached the upper limit (max leg retraction)"<<std::endl;
//      thetax=ComputeBoundR(coeffR_(j), coeffR_(j), hip2FootDistance, nominalHip2FootDistance, max_extension); 
//      d_polytope(j)=ComputeBoundR(coeffDR_(j),coeffDR_(j), hip2FootDistance, nominalHip2FootDistance, max_extension);
//      }
//      if(hip2FootDistance>fabs(max_extension)){
//       thetax=ComputeBoundL(coeffL_(j), coeffL_(j), hip2FootDistance, nominalHip2FootDistance, max_retraction);
//       //ForcePolytopeConstraint::ComputeBoundL (double coeff0, double coeff1,  double Posx, double Pn,double ls) const
//       //bound=(Posx-ls)*(coeff1-coeff0)/(Pn-ls)+coeff0;
//       d_polytope(j)=ComputeBoundL(coeffDL_(j), coeffDL_(j), hip2FootDistance, nominalHip2FootDistance, max_retraction);
//      }
        
    std::cout<<"theta x is "<<thetax<<" and d is "<< d_polytope(j) <<std::endl;
    //f_polytope(0,j)=cos(thetax);
    //f_polytope(2,j)=sin(thetax);
    f_polytope(0,j)=cos(legPitchAngle)*cos(thetax) - sin(legPitchAngle)*sin(thetax); //f_polytope(0,j)=cos(thetax + legPitchAngle);
    f_polytope(2,j)=sin(legPitchAngle)*cos(thetax) + cos(legPitchAngle)*sin(thetax); //f_polytope(2,j)=sin(thetax + legPitchAngle);
   }

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
  //std::cout<<f_polytope<<std::endl;
  //std::cout<<"  "<<std::endl;
  //std::cout<<d_polytope.transpose()<<std::endl;
  //std::cout<<"  "<<std::endl;

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
    double baseToHipX = vector_base_to_ee_B(0) - base_to_hip_distance(0);
    double baseToHipZ = vector_base_to_ee_B(2) - base_to_hip_distance(2);
    double squaredHip2FootDistance = pow(baseToHipX, 2) + pow(baseToHipZ, 2);
    double hip2FootDistance = sqrt(squaredHip2FootDistance);
    double legPitchAngle = 0.0;// = atan2(baseToHipX, baseToHipZ);
    double nominalHip2FootDistance = fabs(nominal_ee_pos_B_(2)); // nominalHip2FootDistance must always be positive
    double max_extension=nominalHip2FootDistance + max_deviation_from_nominal_(2);
    double max_retraction=nominalHip2FootDistance - max_deviation_from_nominal_(2);

    Vector3d f = force_nodes.at(f_node_id).p();
    double theta_coeff=0; double d_coeff=0; double thetax=0;
    int j;
    int row_reset=row;
    for (j=0; j<4; j++)
   {

    if (hip2FootDistance>=fabs(nominalHip2FootDistance)) {
    
  //  double ForcePolytopeConstraint::ComputeCoeffForJac(double coeff0, double coeff1, double P0, double P1) const
  //  bound=(coeff1-coeff0)/(P1-P0);
      thetax=ComputeBound(coeffN_(j), coeffL_(j), hip2FootDistance, nominalHip2FootDistance, max_extension);
      theta_coeff=ComputeCoeffForJac(coeffN_(j), coeffL_(j), nominalHip2FootDistance, max_extension);
       d_coeff=ComputeCoeffForJac(coeffDN_(j), coeffDL_(j), nominalHip2FootDistance, max_extension);
      }

     if (hip2FootDistance<fabs(nominalHip2FootDistance)){
        //std::cout<<" hip2FootDistance is in the upper half (leg is retracted)"<<std::endl;
        thetax=ComputeBound(coeffR_(j), coeffN_(j), hip2FootDistance, max_retraction, nominalHip2FootDistance);
        theta_coeff=ComputeCoeffForJac(coeffR_(j), coeffN_(j), max_retraction, nominalHip2FootDistance);
        d_coeff=ComputeCoeffForJac(coeffDR_(j),coeffDN_(j), max_retraction, nominalHip2FootDistance);
     }

//    if (vector_base_to_ee_B(0)>nominal_ee_pos_B_(0)) 
//      {
//
//        thetax=ComputeBoundR(coeffN_(j),coeffR_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0)); 
//        theta_coeff=ComputeCoeffForJacR(coeffN_(j),coeffR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
//        d_coeff=ComputeCoeffForJacR(coeffDN_(j),coeffDR_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)+max_deviation_from_nominal_(0));
//        
//      }
//    else 
//    {
//
//      thetax=ComputeBoundL(coeffL_(j), coeffN_(j), vector_base_to_ee_B(0), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0)); 
//      theta_coeff=ComputeCoeffForJacL(coeffL_(j),coeffN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
//      d_coeff=ComputeCoeffForJacL(coeffDL_(j),coeffDN_(j), nominal_ee_pos_B_(0), nominal_ee_pos_B_(0)-max_deviation_from_nominal_(0));
//
//    }


    double dlegPitchAngle_dpx = - vector_base_to_ee_B(2)/squaredHip2FootDistance;
    double dlegPitchAngle_dpz = vector_base_to_ee_B(0)/squaredHip2FootDistance;
    double dr_dpx = vector_base_to_ee_B(0)/hip2FootDistance;
    double dr_dpz = vector_base_to_ee_B(2)/hip2FootDistance;
    double dthetax_dpx = theta_coeff*dr_dpx;
    double dthetax_dpz = theta_coeff*dr_dpx;

    double da_dpx = -sin(thetax + legPitchAngle)*(dlegPitchAngle_dpx + dthetax_dpx);
    double da_dpz = -sin(thetax + legPitchAngle)*(dlegPitchAngle_dpz + dthetax_dpz);
    double dc_dpx = cos(thetax + legPitchAngle)*(dlegPitchAngle_dpx + dthetax_dpx);
    double dc_dpz = cos(thetax + legPitchAngle)*(dlegPitchAngle_dpz + dthetax_dpz);
    double dd_dpx = d_coeff;
    double dd_dpz = d_coeff;
    jac.middleRows(row_reset++, 1) = (da_dpx*f(0)+dc_dpx*f(2)-dd_dpx)*JacPosBWrtBaseLin.row(0) + (da_dpz*f(0)+dc_dpz*f(2)-dd_dpz)*JacPosBWrtBaseLin.row(2);

    //double da_dpx = -sin(thetax)*theta_coeff;
    //Jacobian da_dbaseLin = da_dpx*JacPosBWrtBaseLin.row(0);
    //double dc_dpx = cos(thetax)*theta_coeff;
    //Jacobian dc_dbaseLin = dc_dbaseLin*JacPosBWrtBaseLin.row(0);
    //double dd_dpx = d_coeff;
    //Jacobian dd_dbaseLin =  dd_dpx*JacPosBWrtBaseLin.row(0);

    //jac.middleRows(row_reset++, 1) = ((-sin(thetax)*f(0)+cos(thetax)*f(2))*theta_coeff-d_coeff)*JacPosBWrtBaseLin.row(0); //+ JacPosBWrtBaseLin.row(1)+JacPosBWrtBaseLin.row(2);

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
    double baseToHipX = vector_base_to_ee_B(0) - base_to_hip_distance(0);
    double baseToHipZ = vector_base_to_ee_B(2) - base_to_hip_distance(2);
    double squaredHip2FootDistance = pow(baseToHipX, 2) + pow(baseToHipZ, 2);
    double hip2FootDistance = sqrt(squaredHip2FootDistance);
    double legPitchAngle = 0.0;// = atan2(baseToHipX, baseToHipZ);
    double nominalHip2FootDistance = fabs(nominal_ee_pos_B_(2)); // nominalHip2FootDistance must always be positive
    double max_extension=nominalHip2FootDistance + max_deviation_from_nominal_(2);
    double max_retraction=nominalHip2FootDistance - max_deviation_from_nominal_(2);

    Vector3d f = force_nodes.at(f_node_id).p();

   double theta_coeff=0.0; double d_coeff=0.0; double thetax=0.0;
    int j;
    int row_reset=row;

    for (j=0; j<4; j++)
   {
    if (hip2FootDistance>=fabs(nominalHip2FootDistance)) {
    
  //  double ForcePolytopeConstraint::ComputeCoeffForJac(double coeff0, double coeff1, double P0, double P1) const
  //  bound=(coeff1-coeff0)/(P1-P0);
      thetax=ComputeBound(coeffN_(j), coeffL_(j), hip2FootDistance, nominalHip2FootDistance, max_extension);
      theta_coeff=ComputeCoeffForJac(coeffN_(j), coeffL_(j), nominalHip2FootDistance, max_extension);
       d_coeff=ComputeCoeffForJac(coeffDN_(j), coeffDL_(j), nominalHip2FootDistance, max_extension);
      }

     if (hip2FootDistance<fabs(nominalHip2FootDistance)){
        //std::cout<<" hip2FootDistance is in the upper half (leg is retracted)"<<std::endl;
        thetax=ComputeBound(coeffR_(j), coeffN_(j), hip2FootDistance, max_retraction, nominalHip2FootDistance);
        theta_coeff=ComputeCoeffForJac(coeffR_(j), coeffN_(j), max_retraction, nominalHip2FootDistance);
        d_coeff=ComputeCoeffForJac(coeffDR_(j),coeffDN_(j), max_retraction, nominalHip2FootDistance);
     }

    double dlegPitchAngle_dpx = - vector_base_to_ee_B(2)/squaredHip2FootDistance;
    double dlegPitchAngle_dpz = vector_base_to_ee_B(0)/squaredHip2FootDistance;
    double dr_dpx = vector_base_to_ee_B(0)/hip2FootDistance;
    double dr_dpz = vector_base_to_ee_B(2)/hip2FootDistance;
    double dthetax_dpx = theta_coeff*dr_dpx;
    double dthetax_dpz = theta_coeff*dr_dpx;

    double da_dpx = -sin(thetax + legPitchAngle)*(dlegPitchAngle_dpx + dthetax_dpx);
    double da_dpz = -sin(thetax + legPitchAngle)*(dlegPitchAngle_dpz + dthetax_dpz);
    double dc_dpx = cos(thetax + legPitchAngle)*(dlegPitchAngle_dpx + dthetax_dpx);
    double dc_dpz = cos(thetax + legPitchAngle)*(dlegPitchAngle_dpz + dthetax_dpz);
    double dd_dpx = d_coeff;
    double dd_dpz = d_coeff;
    jac.middleRows(row_reset++, 1) = (da_dpx*f(0)+dc_dpx*f(2)-dd_dpx)*JacPosBWrtBaseAng.row(0) + (da_dpz*f(0)+dc_dpz*f(2)-dd_dpz)*JacPosBWrtBaseAng.row(2);
    //jac.middleRows(row_reset++, 1) = ((-sin(thetax)*f(0)+cos(thetax)*f(2))*coeff1-coeff2)*JacPosBWrtBaseAng.row(0); //+ JacPosBWrtBaseLin.row(1)+JacPosBWrtBaseLin.row(2);
    
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
      double baseToHipX = vector_base_to_ee_B(0) - base_to_hip_distance(0);
      double baseToHipZ = vector_base_to_ee_B(2) - base_to_hip_distance(2);
      double squaredHip2FootDistance = pow(baseToHipX, 2) + pow(baseToHipZ, 2);
      double hip2FootDistance = sqrt(squaredHip2FootDistance);
      double nominalHip2FootDistance = fabs(nominal_ee_pos_B_(2)); // nominalHip2FootDistance must always be positive
      double max_extension=nominalHip2FootDistance + max_deviation_from_nominal_(2);
      double max_retraction=nominalHip2FootDistance - max_deviation_from_nominal_(2);
      double legPitchAngle = 0.0;// = atan2(baseToHipX, baseToHipZ);
      double theta_coeff=0; double d_coeff=0; double thetax=0;
    int j;
    for (j=0; j<4; j++)
     {
    if (hip2FootDistance>=fabs(nominalHip2FootDistance)) {
    
  //  double ForcePolytopeConstraint::ComputeCoeffForJac(double coeff0, double coeff1, double P0, double P1) const
  //  bound=(coeff1-coeff0)/(P1-P0);
      thetax=ComputeBound(coeffN_(j), coeffL_(j), hip2FootDistance, nominalHip2FootDistance, max_extension);
      theta_coeff=ComputeCoeffForJac(coeffN_(j), coeffL_(j), nominalHip2FootDistance, max_extension);
       d_coeff=ComputeCoeffForJac(coeffDN_(j), coeffDL_(j), nominalHip2FootDistance, max_extension);
      }

     if (hip2FootDistance<fabs(nominalHip2FootDistance)){
        //std::cout<<" hip2FootDistance is in the upper half (leg is retracted)"<<std::endl;
        thetax=ComputeBound(coeffR_(j), coeffN_(j), hip2FootDistance, max_retraction, nominalHip2FootDistance);
        theta_coeff=ComputeCoeffForJac(coeffR_(j), coeffN_(j), max_retraction, nominalHip2FootDistance);
        d_coeff=ComputeCoeffForJac(coeffDR_(j),coeffDN_(j), max_retraction, nominalHip2FootDistance);
     }

    double dlegPitchAngle_dpx = - vector_base_to_ee_B(2)/squaredHip2FootDistance;
    double dlegPitchAngle_dpz = vector_base_to_ee_B(0)/squaredHip2FootDistance;
    double dr_dpx = vector_base_to_ee_B(0)/hip2FootDistance;
    double dr_dpz = vector_base_to_ee_B(2)/hip2FootDistance;
    double dthetax_dpx = theta_coeff*dr_dpx;
    double dthetax_dpz = theta_coeff*dr_dpx;

    double da_dpx = -sin(thetax + legPitchAngle)*(dlegPitchAngle_dpx + dthetax_dpx);
    double da_dpz = -sin(thetax + legPitchAngle)*(dlegPitchAngle_dpz + dthetax_dpz);
    double dc_dpx = cos(thetax + legPitchAngle)*(dlegPitchAngle_dpx + dthetax_dpx);
    double dc_dpz = cos(thetax + legPitchAngle)*(dlegPitchAngle_dpz + dthetax_dpz);
    double dd_dpx = d_coeff;
    double dd_dpz = d_coeff;
    jac.middleRows(row_reset++, 1) = (da_dpx*f(0)+dc_dpx*f(2)-dd_dpx)*JacPosBWrtNodes.row(0) + (da_dpz*f(0)+dc_dpz*f(2)-dd_dpz)*JacPosBWrtNodes.row(2);
    //jac.middleRows(row_reset++, 1) = ((-sin(thetax)*f(0)+cos(thetax)*f(2))*coeff1-coeff2)*JacPosBWrtNodes.row(0); //+ JacPosBWrtBaseLin.row(1)+JacPosBWrtBaseLin.row(2);
    
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

double ForcePolytopeConstraint::ComputeBound(double coeff0, double coeff1, double Posx, double P0, double P1) const
{
  double bound;
  bound=(Posx-P0)*(coeff1-coeff0)/(P1-P0)+coeff0;
  return bound;
}
double ForcePolytopeConstraint::ComputeCoeffForJac(double coeff0, double coeff1, double P0, double P1) const
{
  double bound;

  bound=(coeff1-coeff0)/(P1-P0);
  return bound;
}

//double ForcePolytopeConstraint::ComputeBoundL (double coeff0, double coeff1,  double Posx, double Pn,double ls) const
//{
//  double bound;
//  bound=(Posx-ls)*(coeff1-coeff0)/(Pn-ls)+coeff0;
//  return bound;
//}
//double ForcePolytopeConstraint::ComputeBoundR (double coeff0, double coeff1, double Posx, double Pn,double rs) const
//{
//  double bound;
//  bound=(Posx-Pn)*(coeff1-coeff0)/(rs-Pn)+coeff0;
//  return bound;
//}
//double ForcePolytopeConstraint::ComputeCoeffForJacL (double coeff0, double coeff1, double Pn,double ls) const
//{
//  double bound;
//
//  bound=(coeff1-coeff0)/(Pn-ls);
//  return bound;
//}
//double ForcePolytopeConstraint::ComputeCoeffForJacR (double coeff0, double coeff1, double Pn,double rs) const
//{
//  double bound;
//  bound=(coeff1-coeff0)/(rs-Pn);
//  return bound;
//}

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
  std::cout<<"leg N. "<<ee<<std::endl;
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
