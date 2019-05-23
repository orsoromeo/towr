#include <towr/constraints/geometry.h>
#include <towr/constraints/base_acc_constraint_range_lin.h>
#include <ifopt/problem.h>
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

Geometry::Geometry (DynamicModel::Ptr model,
                    const HeightMap::Ptr & terrain,
                    const SplineHolder& spline_holder,
                    int numberofleg,
                    const std::vector<std::vector<double> > &phase_durations)
{
  phase_durations_=phase_durations;
  numberoflegs_=numberofleg;
  model_=model;
  terrain_=terrain;
  ee_motion_=spline_holder.ee_motion_;
  base_ << 0.0, 0.0, 1.0;
  LinearEdges_.resize(4*numberofleg,3);
  AngularEdges_.resize(4*numberofleg,3);
  EhatLin_.resize(4*numberofleg,3);
  EhatAng_.resize(4*numberofleg,3);
}

int Geometry::GetNumberOfFeetInTouch (double t) const
{
  int m=0;
  for (int ee=0; ee<numberoflegs_; ee++)
  {
    double l=0.0;
    int z=0;
    int s=0;
    std::vector<double> vector=phase_durations_.at(ee);
    for (int i=0; i<vector.size(); i++)
     {
      double a=vector.at(i);
      l=l+a;
      if (t-l<1e-6)
       {
        z=i;
        if (t==l) {z=0;}
        break;
       }
      else {}
     }
   if ((z%2)==0) {s=1; ee_motion_in_touch_.push_back(ee_motion_.at(ee)); }
   else {s=0;}
   m=m+s;
  }
  return m;
}

Eigen::MatrixXd Geometry::ReturnNormalTerrain (double t) const
{

  Eigen::MatrixXd n(GetNumberOfFeetInTouch (t),3);
  for (int ee=0; ee<GetNumberOfFeetInTouch (t); ee++)
  {
    Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
   n.row(ee)= terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());

  }
 return n;
}


Eigen::Matrix3d Geometry::RotationMatrix (Eigen::Vector3d t, double angle) const

{
  double cmu=cos(angle);
  double smu=sin(angle);
  Eigen::Matrix3d M;
  M(0,0)=std::pow(t(0),2)*(1-cmu)+cmu;          M(0,1)=t(0)*t(1)*(1-cmu)-t(2)*smu;      M(0,2)=t(0)*t(2)*(1-cmu)+t(1)*smu;
  M(1,0)=t(0)*t(1)*(1-cmu)+t(2)*smu;            M(1,1)=std::pow(t(1),2)*(1-cmu)+cmu;    M(1,2)=t(1)*t(2)*(1-cmu)-t(0)*smu;
  M(2,0)=t(0)*t(2)*(1-cmu)-t(1)*smu;            M(2,1)=t(1)*t(2)*(1-cmu)+t(0)*smu;      M(2,2)=std::pow(t(2),2)*(1-cmu)+cmu;
  return M;
}

Eigen::Vector3d Geometry::ComputeNextEdge (Eigen::Matrix3d M, Eigen::Vector3d n) const
{
  Eigen::Vector3d edge=M*n;
  return edge;
}

Eigen::MatrixXd Geometry::ComputeCone (double t) const
{
  LinearEdges_.setZero();
  AngularEdges_.setZero();
  EhatLin_.setZero();
  EhatAng_.setZero();
  std::cout<<t<<std::endl;
  //std::cout<<ee_motion_.at(0)->GetPoint(t).p()<<std::endl;
  for (int ee=0; ee<GetNumberOfFeetInTouch (t); ee++)
   {
    Eigen::Vector3d normal=ReturnNormalTerrain(t).row(ee);
    Eigen::Vector3d axis;
    axis<< normal(1), -normal(0), 0.0;
    double angle= ComputeRotationAngle(normal);
    Eigen::MatrixXd tmp = ComputeLinearPartOfTheCone(axis,angle);
    for (int i=0; i<4; i++)
     { LinearEdges_.block(4*ee,0,4,3).row(i)= tmp.row(i);
       AngularEdges_.block(4*ee,0,4,3).row(i)=ComputeAngularPartOfTheCone(t,axis,angle,ee).row(i);
       EhatLin_.block(4*ee,0,4,3).row(i)=RotationMatrix(normal,angle).inverse()*LinearEdges_.block(4*ee,0,4,3).row(i).transpose();
       EhatAng_.block(4*ee,0,4,3).row(i)=RotationMatrix(normal,angle).inverse()*AngularEdges_.block(4*ee,0,4,3).row(i).transpose();
    }
   }

    Eigen::MatrixXd Edges(4*numberoflegs_,6);
    Edges.block(0,0,4*numberoflegs_,3)=LinearEdges_;
    Edges.block(0,3,4*numberoflegs_,3)=AngularEdges_;
    std::cout<<Edges.transpose()<<std::endl;
    return Edges.transpose();



}

Eigen::MatrixXd
Geometry::ComputeLinearPartOfTheCone (const Eigen::Vector3d& axis, const double & angle) const
{
  Eigen::MatrixXd E; E.resize(4,3); E.setZero();
  Eigen::MatrixXd LinearEdges; LinearEdges.resize(4,3); LinearEdges.setZero();
  for (int i=0; i<4; i++)
  {
    E.row(i)=ComputeNextEdge(RotationZ(PI/4+(double)i*PI/2)*RotationX(20*PI/180),base_);
    //std::cout<<E.row(i)<<std::endl;
    LinearEdges.row(i)=ComputeNextEdge(RotationMatrix(axis,angle),E.row(i));
  }

  return LinearEdges;
}

Eigen::MatrixXd Geometry::ComputeAngularPartOfTheCone (double t, Eigen::Vector3d axis, double angle,double ee) const
{

  Eigen::MatrixXd AngularEdges(4,3);
  Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
  for (int i=0; i<4; i++)
    {
      Eigen::Vector3d row=LinearEdges_.row(4*ee+i);
      AngularEdges.row(i)=p.cross(row);
    }
   return AngularEdges;
}


Eigen::Matrix3d Geometry::RotationX (const double& angle) const

{
  double cmu=cos(angle);
  double smu=sin(angle);
  Eigen::Matrix3d M; M.setZero();
  M(0,0)=1.0;
  M(1,1)=cmu;     M(1,2)=-smu;
  M(2,1)=smu;     M(2,2)=cmu;
  return M;
}

Eigen::Matrix3d Geometry::RotationZ (const double& angle) const

{
  double cmu=cos(angle);
  double smu=sin(angle);
  Eigen::Matrix3d M; M.setZero();
  M(0,0)=cmu;     M(0,1)=-smu;
  M(1,0)=smu;     M(1,1)=cmu;
  M(2,2)=1.0;
  return M;
}

double Geometry::ComputeRotationAngle (Eigen::Vector3d normal) const
{
   auto ab=normal.dot(base_);
   double angle=acos(ab/(normal.norm()*base_.norm()));
   return angle;
}


} /* namespace towr */
