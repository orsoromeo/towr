#include <towr/constraints/geometry.h>
#include <towr/constraints/base_acc_constraint_range_lin.h>
#include <politopixAPI.h>
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

Geometry::Geometry (DynamicModel::Ptr model, HeightMap::Ptr terrain, const SplineHolder& spline_holder)
{
  model_=model;
  terrain_=terrain;
  ee_motion_=spline_holder.ee_motion_;
  base_ << 0.0, 0.0, 1.0;
}

int Geometry::GetNumberOfFeetInTouch (double t)
{
  LinearEdges_(4*GetNumberOfFeetInTouch (t),3);
  AngularEdges_(4*GetNumberOfFeetInTouch (t),3);
  for (int ee=0; model_->GetEECount(); ee++)
   { Eigen::Vector3d p = ee_motion_.at(ee)->GetPoint(t).p();
     if (p.z()==terrain_->GetHeight(p.x(),p.y()))
      {
       ee_motion_in_touch_.push_back(ee_motion_.at(ee));
      }
  return ee_motion_in_touch_.size();
  }
}
// questa parte viene da cone.cc
Eigen::MatrixXd Geometry::ReturnNormalTerrain (double t)
{

  Eigen::MatrixXd n(GetNumberOfFeetInTouch (t),3);
  for (int ee=0; ee<GetNumberOfFeetInTouch (t); ee++)
  {Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
   n.row(ee)= terrain_->GetNormalizedBasis(HeightMap::Normal, p.x(), p.y());
  }

 return n;
}

Eigen::MatrixXd Geometry::ReturnTangentTerrain (double t)
{
 Eigen::MatrixXd tang (GetNumberOfFeetInTouch (t),3);
 for (int ee=0; ee<GetNumberOfFeetInTouch (t); ee++)
  {   Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
      tang.row(ee) = terrain_->GetNormalizedBasis(HeightMap::Tangent1, p.x(), p.y());
  }
 return tang;
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
void Geometry::ComputeCone (double t)
{


  for (int ee=0; ee<GetNumberOfFeetInTouch (t); ee++)
   {
    Eigen::Vector3d normal=ReturnNormalTerrain(t).row(ee);
    Eigen::Vector3d axis;
    axis<< normal(1), -normal(0), 0.0;
    double angle= ComputeRotationAngle(normal);
    LinearEdges_.middleRows(ee,4)=ComputeLinearPartOfTheCone(axis,angle);
    AngularEdges_.middleRows(ee,4)=ComputeAngularPartOfTheCone(t,axis,angle,ee);
    }

}

Eigen::MatrixXd Geometry::ComputeLinearPartOfTheCone (Eigen::Vector3d axis, double angle)
{
  Eigen::MatrixXd E;
  Eigen::MatrixXd LinearEdges;
  for (int i; i<4; i++)
  {
    E.row(i)=ComputeNextEdge(RotationZ(PI/4+i*PI/2)*RotationX(20),base_); //gli angoli vanno in radianti
    LinearEdges.row(i)=ComputeNextEdge(RotationMatrix(axis,angle),E.row(i));
  }
  return LinearEdges;
}

Eigen::MatrixXd Geometry::ComputeAngularPartOfTheCone (double t, Eigen::Vector3d axis, double angle,double ee)
{
  //Eigen::MatrixXd EA;
  Eigen::MatrixXd AngularEdges;
  Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
  for (int i=0; i<4; i++)
    {
      Eigen::Vector3d row=LinearEdges_.row(i);
      AngularEdges.row(i)=p.cross(row);
      //AngularEdges.row(i)=ComputeNextEdge(RotationMatrix(axis,angle),EA.row(i));

     }
   return AngularEdges;
}


Eigen::Matrix3d Geometry::RotationX (double angle) const

{
  double cmu=cos(angle);
  double smu=sin(angle);
  Eigen::Matrix3d M;
  M(0,0)=1;
  M(1,1)=cmu;     M(1,2)=-smu;
  M(2,1)=smu;     M(2,2)=cmu;
  return M;
}

Eigen::Matrix3d Geometry::RotationZ (double angle) const

{
  double cmu=cos(angle);
  double smu=sin(angle);
  Eigen::Matrix3d M;
  M(0,0)=cmu;     M(0,1)=-smu;
  M(1,0)=smu;     M(1,1)=cmu;
  M(2,2)=1;
  return M;
}

double Geometry::ComputeRotationAngle (Eigen::Vector3d axis) const
{
   auto ab=axis.dot(base_);
   double angle=acos(ab/(axis.norm()*base_.norm()));
   return angle;
}

}
