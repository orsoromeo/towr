#include <towr/constraints/derivative.h>

namespace towr {

Derivative::Derivative (HeightMap::Ptr &terrain,
                        const SplineHolder& spline_holder
                         )
{

  terrain_=terrain;
  ee_motion_=spline_holder.ee_motion_;
  base_ << 0.0, 0.0, 1.0;

}


NodeSpline::Jacobian Derivative::GetDerivativeOfNormalWrtNodes (int ee, double t) const
{
  //NodesVariablesPhaseBased::Ptr ee_motion_NVP = variables_->GetComponent<NodesVariablesPhaseBased>(id::EEMotionNodes(ee));
  Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
  double px=p.x();
  double py=p.y();
  HeightMap::Vector3d Dnx= terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, X_, px, py);
  HeightMap::Vector3d Dny= terrain_->GetDerivativeOfNormalizedBasisWrt(HeightMap::Normal, Y_, px, py);
  NodeSpline::Jacobian jac=ee_motion_in_touch_.at(ee)->GetJacobianWrtNodes(t,kPos);
  Eigen::MatrixXd Dn(3,3);
  Dn.col(0)=Dnx;
  Dn.col(1)=Dny;
  Dn.col(2) << 0.0, 0.0, 0.0;
  NodeSpline::Jacobian jac1=Dn.sparseView()*jac;
  return jac1;
  //int id;
  //id= Spline::GetSegmentID(t, duration_.at(ee));
  //int ee_node_id = ee_motion_NVP->GetNodeIDAtStartOfPhase(id);
  //int idx = ee_motion_NVP->GetOptIndex(NodesVariables::NodeValueInfo(ee_node_id, kPos, X));
  //Eigen::MatrixXd jac(3, ee_motion_.at(ee)->GetJacobianWrtNodes(t,kPos).cols());  //come posso determinare il numero di colonne senza chiamare quella funzione?
  //jac.col(idx)=Dnx;
  //jac.col(idx+1)=Dny;
  //NodeSpline::Jacobian jac1(3, ee_motion_.at(ee)->GetJacobianWrtNodes(t,kPos).cols());
  //NodeSpline::Jacobian jac1=jac.sparseView();
}
Eigen::Vector3d Derivative::GetDerivativeofAngleWrtNormal(Eigen::Vector3d normal, double angle) const
{//manca l'arcocoseno!
  Eigen::Vector3d jac;
  double Normx=sqrt(normal.norm());
  double norm3=std::pow(Normx,3);
  jac(0)=-normal(0)*normal(2)/(norm3);
  jac(1)=-normal(1)*normal(2)/(norm3);
  jac(2)=-std::pow(normal(2),2)/(norm3);
  return jac * (1/sqrt(1-std::pow(angle,2)));

}

NodeSpline::Jacobian Derivative::GetDerivativeofAngleWrtNodes(int ee, double t, Eigen::Vector3d normal,double angle) const
{
  NodeSpline::Jacobian jac=GetDerivativeOfNormalWrtNodes(ee,t);
  Eigen::Vector3d dang=GetDerivativeofAngleWrtNormal(normal,angle);
  NodeSpline::Jacobian jac1=dang.sparseView()*jac;
  return jac1;
}

NodeSpline::Jacobian Derivative::GetDerivativeofLinearEdgeWrtNormal (Eigen::Vector3d normal, double angle, Eigen::Vector3d edge) const
{
   NodeSpline::Jacobian jac1(3,3);
   auto nx=normal(0);
   auto ny=normal(1);
   auto ca=cos(angle);
   auto sa=sin(angle);
   auto da=GetDerivativeofAngleWrtNormal(normal,angle);
   auto dax=da(0);
   auto day=da(1);
   auto daz=da(2);
   //derivative wrt nx, 1 column
   jac1.coeffRef(0,0)=(std::pow(ny,2)*sa*dax-sa*dax)*edge(0)-(ny-ny*(ca+nx*(-sa)*dax))*edge(1)-(sa+nx*ca*dax)*edge(2);
   jac1.coeffRef(1,0)=(-ny+ny*(ca+nx*(-sa)*dax))*edge(0)+(2*nx-(2*nx*ca+std::pow(nx,2)*(-sa)*dax)-sa*dax)*edge(1)-ca*dax*edge(2);
   jac1.coeffRef(2,0)=(sa+nx*ca*dax)*edge(0)+(ny*ca*dax)*edge(1)-sa*dax*edge(2);
    //derivative wrt ny, 2 column
   jac1.coeffRef(0,1)=(2*ny-(2*ny*ca-std::pow(ny,2)*sa*day)-sa*day)*edge(0)-(nx-nx*ca+nx*sa*day)*edge(1)-(nx*ca*day)*edge(2);
   jac1.coeffRef(1,1)=(-nx+nx*ca-ny*nx*sa*day)*edge(0)+(std::pow(nx,2)*sa*day-sa*day)*edge(1)-(sa+ny*ca*day)*edge(2);
   jac1.coeffRef(2,0)=(nx*ca*day)*edge(0)+(sa+ny*ca*day)*edge(1)-sa*day*edge(2);
   //derivative wrt nz, 3 column
   jac1.coeffRef(0,2)=(std::pow(ny,2)*sa*daz-sa*daz)*edge(0)-(ny*nx*sa*daz)*edge(1)-(nx*ca*daz)*edge(2);
   jac1.coeffRef(1,2)=(-ny*nx*sa*daz)*edge(0)+(std::pow(nx,2)*sa*daz-sa*daz)*edge(1)-ny*ca*dax*edge(2);
   jac1.coeffRef(2,2)=(nx*ca*daz)*edge(0)+(ny*ca*daz)*edge(1)-sa*daz*edge(2);

   return jac1;
}

NodeSpline::Jacobian Derivative::GetDerivativeofLinearEdgeWrtNodes(int ee, double t,Eigen::Vector3d normal, double angle, Eigen::Vector3d edge) const
{
  NodeSpline::Jacobian jac1=GetDerivativeofLinearEdgeWrtNormal(normal,angle,edge);
  NodeSpline::Jacobian jac2=GetDerivativeOfNormalWrtNodes(ee,t);
  return jac1*jac2;
}

NodeSpline::Jacobian Derivative::GetDerivativeofAngularEdgeWrtNodes(int ee, double t,Eigen::Vector3d normal, double angle, Eigen::Vector3d edge) const
{
  Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
  NodeSpline::Jacobian jac_p=ee_motion_in_touch_.at(ee)->GetJacobianWrtNodes(t,kPos);
  NodeSpline::Jacobian jac_e=GetDerivativeofLinearEdgeWrtNodes(ee, t,normal, angle,edge);
  NodeSpline::Jacobian jac_ang(3,jac_p.cols());
  jac_ang.row(0)=jac_p.row(1)*edge(2)+p.y()*jac_e.row(2)-jac_p.row(2)*edge(1)-p.z()*jac_e.row(1);
  jac_ang.row(1)=jac_p.row(2)*edge(0)+p.z()*jac_e.row(0)-jac_p.row(0)*edge(2)-p.x()*jac_e.row(2);
  jac_ang.row(2)=jac_p.row(0)*edge(1)+p.x()*jac_e.row(1)-jac_p.row(1)*edge(0)-p.y()*jac_e.row(0);
  return jac_ang;
}
}/*end namespace*/
