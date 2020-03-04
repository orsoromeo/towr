#include <towr/constraints/derivative.h>

namespace towr {

Derivative::Derivative (const HeightMap::Ptr &terrain,
                        const SplineHolder& spline_holder,
                        int numberofleg,
                        const std::vector<std::vector<double> > &phase_durations

                         )
{
  phase_durations_=phase_durations;
  numberoflegs_=numberofleg;
  terrain_=terrain;
  ee_motion_=spline_holder.ee_motion_;
  base_ << 0.0, 0.0, 1.0;

}


NodeSpline::Jacobian Derivative::GetDerivativeOfNormalWrtNodes (int ee, double t) const
{

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
{
  Eigen::Vector3d jac;
  double Normx=sqrt(normal.norm());
  double norm3=std::pow(Normx,3);
  jac(0)=-normal(0)*normal(2)/(norm3);
  jac(1)=-normal(1)*normal(2)/(norm3);
  jac(2)=1/Normx-std::pow(normal(2),2)/(norm3);
  return jac * (-1/sqrt(1-std::pow(angle,2)));

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
   auto axis_norm=sqrt(std::pow(normal(0),2)+std::pow(normal(0),2));
   if (axis_norm!=0)
   {
     auto nx=normal(0);
     auto ny=normal(1);
     auto ca=cos(angle);
     auto sa=sin(angle);
     auto da=GetDerivativeofAngleWrtNormal(normal,angle);
     auto dax=da(0);
     auto day=da(1);
     auto daz=da(2);
     auto der=GetDerivativeofRotationAxisWrtNormal(normal);
     auto A=std::pow(ny,2)/(axis_norm*axis_norm);          auto DA=der.col(0);
     auto B=ny*nx/(axis_norm*axis_norm);                   auto DB=der.col(1);
     auto C=nx/axis_norm;                                  auto DC=der.col(2);
     auto D=ny/axis_norm;                                  auto DD=der.col(3);
     auto E=std::pow(nx,2)/(axis_norm*axis_norm);          auto DE=der.col(4);

     //derivative wrt nx, 1 column
     jac1.coeffRef(0,0)=(DA(0)*(1-ca)+A*sa*dax-sa*dax)*edge(0)+(-DB(0)+DB(0)*ca-B*sa*dax)*edge(1)-(DC(0)*sa+C*ca*dax)*edge(2);
     jac1.coeffRef(1,0)=(-DB(0)*(1-ca)-B*sa*dax)*edge(0)+(DE(0)-(DE(0)*ca-E*sa*dax)-sa*dax)*edge(1)-(DD(0)*sa+D*ca*dax)*edge(2);
     jac1.coeffRef(2,0)=(DC(0)*sa+C*ca*dax)*edge(0)+(DD(0)*sa+D*ca*dax)*edge(1)+DC(0)*edge(2);
     // //derivative wrt ny, 2 column
     jac1.coeffRef(0,1)=(DA(1)*(1-ca)+A*sa*day-sa*dax)*edge(0)+(-DB(1)+DB(1)*ca-B*sa*day)*edge(1)-(DC(1)*sa+C*ca*day)*edge(2);
     jac1.coeffRef(1,1)=(-DB(1)*(1-ca)-B*sa*day)*edge(0)+(DE(1)-(DE(1)*ca-E*sa*day)-sa*day)*edge(1)-(DD(1)*sa+D*ca*day)*edge(2);
     jac1.coeffRef(2,1)=(DC(1)*sa+C*ca*day)*edge(0)+(DD(1)*sa+D*ca*day)*edge(1)+DC(1)*edge(2);
     //derivative wrt nz, 3 column
     jac1.coeffRef(0,2)=(1-A)*(-sa)*daz*edge(0)-B*sa*daz*edge(1)-C*ca*daz*edge(2);
     jac1.coeffRef(1,2)=-B*sa*daz*edge(0)+(1-E)*(-sa)*daz*edge(1)-D*ca*daz*edge(2);
     jac1.coeffRef(2,2)=C*ca*daz*edge(0)+D*ca*daz*edge(1)-ca*daz*edge(2);
     ////derivative wrt nx, 1 column
     //jac1.coeffRef(0,0)=(std::pow(ny,2)*sa*dax-sa*dax)*edge(0)-(ny-ny*(ca+nx*(-sa)*dax))*edge(1)-(sa+nx*ca*dax)*edge(2);
     //jac1.coeffRef(1,0)=(-ny+ny*(ca+nx*(-sa)*dax))*edge(0)+(2*nx-(2*nx*ca+std::pow(nx,2)*(-sa)*dax)-sa*dax)*edge(1)-ca*dax*edge(2);
     //jac1.coeffRef(2,0)=(sa+nx*ca*dax)*edge(0)+(ny*ca*dax)*edge(1)-sa*dax*edge(2);
     // //derivative wrt ny, 2 column
     //jac1.coeffRef(0,1)=(2*ny-(2*ny*ca-std::pow(ny,2)*sa*day)-sa*day)*edge(0)-(nx-nx*ca+nx*sa*day)*edge(1)-(nx*ca*day)*edge(2);
     //jac1.coeffRef(1,1)=(-nx+nx*ca-ny*nx*sa*day)*edge(0)+(std::pow(nx,2)*sa*day-sa*day)*edge(1)-(sa+ny*ca*day)*edge(2);
     //jac1.coeffRef(2,0)=(nx*ca*day)*edge(0)+(sa+ny*ca*day)*edge(1)-sa*day*edge(2);
     ////derivative wrt nz, 3 column
     //jac1.coeffRef(0,2)=(std::pow(ny,2)*sa*daz-sa*daz)*edge(0)-(ny*nx*sa*daz)*edge(1)-(nx*ca*daz)*edge(2);
     //jac1.coeffRef(1,2)=(-ny*nx*sa*daz)*edge(0)+(std::pow(nx,2)*sa*daz-sa*daz)*edge(1)-ny*ca*dax*edge(2);
     //jac1.coeffRef(2,2)=(nx*ca*daz)*edge(0)+(ny*ca*daz)*edge(1)-sa*daz*edge(2);
  }

   return jac1;
}

NodeSpline::Jacobian Derivative::GetDerivativeofLinearEdgeWrtNodes(int ee, double t,Eigen::Vector3d normal, double angle, Eigen::Vector3d edge) const
{

  GetNumberOfFeetInTouch (t);
  NodeSpline::Jacobian jac1=GetDerivativeofLinearEdgeWrtNormal(normal,angle,edge);
  NodeSpline::Jacobian jac2=GetDerivativeOfNormalWrtNodes(ee,t);

  return jac1*jac2;
}

NodeSpline::Jacobian Derivative::GetDerivativeofAngularEdgeWrtNodes(int ee, double t,Eigen::Vector3d normal, double angle, Eigen::Vector3d edge, Eigen::Vector3d edges ) const
{
  GetNumberOfFeetInTouch (t);
  Eigen::Vector3d p = ee_motion_in_touch_.at(ee)->GetPoint(t).p();
  NodeSpline::Jacobian jac_p=ee_motion_in_touch_.at(ee)->GetJacobianWrtNodes(t,kPos);
  NodeSpline::Jacobian jac_e=GetDerivativeofLinearEdgeWrtNodes(ee, t,normal, angle,edge);
  NodeSpline::Jacobian jac_ang(3,jac_p.cols());
  jac_ang.row(0)=jac_p.row(1)*edges(2)+p.y()*jac_e.row(2)-jac_p.row(2)*edges(1)-p.z()*jac_e.row(1);
  jac_ang.row(1)=jac_p.row(2)*edges(0)+p.z()*jac_e.row(0)-jac_p.row(0)*edges(2)-p.x()*jac_e.row(2);
  jac_ang.row(2)=jac_p.row(0)*edges(1)+p.x()*jac_e.row(1)-jac_p.row(1)*edges(0)-p.y()*jac_e.row(0);
  return jac_ang;
}
Eigen::MatrixXd Derivative::GetDerivativeofRotationAxisWrtNormal(Eigen::Vector3d vector) const
{
  Eigen::MatrixXd der(5,2);
  double ny=vector(1);
  double nx=vector(0);
  double norm=vector.norm();
  double norm2=std::pow(norm,2);
  double nx2=std::pow(nx,2);
  double ny2=std::pow(ny,2);
  double norm4=std::pow(norm2,2);
  double norm3=std::pow(norm,3);
  der(0,0)=-2*ny2*nx/norm4;                   der(0,1)=2*ny/norm2-2*std::pow(ny,3)/norm4;
  der(1,0)=ny/norm2-2*nx2*ny/norm4;           der(1,1)=nx/norm2-2*nx*ny2/norm4;
  der(2,0)=1/norm-nx2/norm3;                  der(2,1)=nx*ny/norm3;
  der(3,0)=ny/norm3;                          der(3,1)=1/norm-ny2/norm3;
  der(4,0)=2*nx/norm2-2*std::pow(nx,3)/norm4; der(4,1)=-2*nx2*ny/norm4;
  return der;
}
void Derivative::GetNumberOfFeetInTouch (double t) const
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
        if (t-l>-1e-06) {z=0;}
        break;
       }
      else {}
     }
   if ((z%2)==0) {ee_motion_in_touch_.push_back(ee_motion_.at(ee)); }

   }
}
}/*end namespace*/
