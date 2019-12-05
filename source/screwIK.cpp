#include <math.h>
#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <Eigen/Dense>

namespace py = pybind11;

using namespace Eigen;
using namespace std;

#define l0 0.1932
#define l1 0.305
#define l2 0.1625
#define l3 0.059742
#define l4 0.037508
#define l5 0.34692
#define l6 0.32992
#define l7 0.215
#define l8 0.090
#define l9 0.092
#define l10 0.330
#define l11 0.300
#define l12 0.123005
#define l13 0.146
#define l14 0.018
#define l15 0.026
#define l16 0.0175

Matrix3d skewsim(Vector3d w){
    Matrix3d m;     
    m <<   0 ,-w(2), w(1),
         w(2),   0 ,-w(0),
        -w(1), w(0),   0 ;
    return m;
}

Matrix4d LieExp(VectorXd xi, double theta){
    
    Vector3d w;
    Vector3d v;
    Matrix3d R;
    Vector3d p;
    Matrix4d g;

    w << xi(0), xi(1), xi(2);
    v << xi(3), xi(4), xi(5);
    
    if(w.norm()==0){
        R = Matrix3d::Identity();
        p = v*theta;      
    }else{
        R = Matrix3d::Identity()+skewsim(w)*sin(theta)+(skewsim(w)*skewsim(w))*(1-cos(theta));
        Matrix3d tp;
        tp << w(0)*w(0), w(1)*w(0), w(2)*w(0),
              w(0)*w(1), w(1)*w(1), w(2)*w(1),
              w(0)*w(2), w(1)*w(2), w(2)*w(2);
        
        p = (Matrix3d::Identity()-R)*(w.cross(v))+tp*(v*theta);
    }

   g << R(0,0), R(0,1), R(0,2), p(0),
        R(1,0), R(1,1), R(1,2), p(1),
        R(2,0), R(2,1), R(2,2), p(2),
           0  ,    0  ,    0  ,   1 ; 

   return g;
}


MatrixXd LieAd(Matrix4d g){
    MatrixXd Ad_g(6,6);
    Matrix3d R = g.block(0,0,3,3) ;
    Vector3d p = g.block(0,3,3,1);
    Ad_g = Matrix<double, 6, 6>::Identity();
    Ad_g.block(0,0,3,3)=R;
    Ad_g.block(3,3,3,3)=R;
    Ad_g.block(0,3,3,3)=skewsim(p)*R;
    return Ad_g;
}

MatrixXd LieAdInv(Matrix4d g){
    MatrixXd Ad_g(6,6);
    Matrix3d R = g.block(0,0,3,3) ;
    Vector3d p = g.block(0,3,3,1);
    Ad_g = Matrix<double, 6, 6>::Identity();
    Ad_g.block(0,0,3,3)=R.transpose();
    Ad_g.block(3,3,3,3)=R.transpose();
    Ad_g.block(0,3,3,3)=-R.transpose()*skewsim(p);
    return Ad_g;
}

Vector4d cart2hom(Vector3d p){
    Vector4d ph;
    ph << p(0), p(1), p(2), 1;
    return ph;
}

double PK1(VectorXd xi, Vector3d p, Vector3d q, Vector3d r){
    Vector3d w = xi.block(0,0,3,1);
    Vector3d u = p-r;
    Vector3d v = q-r;
    Vector3d up = u - (w*w.transpose())*u;
    Vector3d vp = v - (w*w.transpose())*v;
    return atan2((w.transpose()*up.cross(vp)), (up.transpose()*vp));
}

Vector4d PK2(VectorXd xi1, VectorXd xi2, Vector3d p, Vector3d q, Vector3d r){
    Vector3d w1 = xi1.block(0,0,3,1);
    Vector3d w2 = xi2.block(0,0,3,1);
    
    Vector3d u = p-r;
    Vector3d v = q-r;
    
   
    
    
    double alpha = (
                    ((w1.transpose()*w2)(0)*(w2.transpose()*u)(0)-(w1.transpose()*v)(0))
                    /
                    ((pow((w1.transpose()*w2)(0),2)-1))
                   );
      
    double beta  = (
                    ((w1.transpose()*w2)(0)*(w1.transpose()*v)(0)-(w2.transpose()*u)(0))
                    /
                    ((pow((w1.transpose()*w2)(0), 2) - 1 ))
                   );


    double gamma2 = (
                     (pow(u.norm(),2)-pow(alpha,2)-pow(beta,2)-2*alpha*beta*(w1.transpose()*w2)(0))
                     /
                     (pow(w1.cross(w2).norm(), 2)) 
                    );


    Vector3d d = r+alpha*w1+beta*w2+sqrt(gamma2)*w1.cross(w2);
    Vector3d c = r+alpha*w1+beta*w2-sqrt(gamma2)*w1.cross(w2);

    Vector3d up = u - w2*(w2.transpose()*u);
    Vector3d vp = v - w1*(w1.transpose()*v);
    
    Vector3d m = c-r;
    Vector3d m1p = m - w1*(w1.transpose()*m);
    Vector3d m2p = m - w2*(w2.transpose()*m);
    
    Vector3d n = d-r;
    Vector3d n1p = n - w1*(w1.transpose()*n);
    Vector3d n2p = n - w2*(w2.transpose()*n);
    
    
    double theta2_1=atan2((w2.transpose()*up.cross(m2p)),(up.transpose()*m2p));
    double theta1_1=atan2((w1.transpose()*m1p.cross(vp)),(m1p.transpose()*vp));
    
    double theta2_2=atan2((w2.transpose()*up.cross(n2p)),(up.transpose()*n2p));
    double theta1_2=atan2((w1.transpose()*n1p.cross(vp)),(n1p.transpose()*vp));
    
    Vector4d solution;
    solution <<  theta1_1, theta1_2, theta2_1, theta2_2;

    return solution;

}

Vector2d PK3(VectorXd xi, Vector3d p, Vector3d q, Vector3d r, double delta){
    Vector3d w = xi.block(0,0,3,1);
    
    Vector3d u = p-r;
    Vector3d v = q-r;
    
    Vector3d up = u - w*(w.transpose()*u);
    Vector3d vp = v - w*(w.transpose()*v);
    
    double deltap2=pow(delta,2)-pow((w.transpose()*(p-q)).norm(),2);
    
    double alpha = atan2((w.transpose()*up.cross(vp)), (up.transpose()*vp));
    
    double beta = acos((pow(up.norm(),2)+pow(vp.norm(),2)-deltap2)/(2*up.norm()*vp.norm()));
    
    Vector2d solution;
    solution << alpha+beta, alpha-beta;
    
    return solution;
}

Vector3d ll_fk_com_from_foot(VectorXd angles){

    VectorXd xi1(6), xi2(6), xi3(6), xi4(6), xi5(6), xi6(6);
    //          w              v
    xi1  << 1., 0., 0.,   0., l12, 0.;
    xi2  << 0., 1., 0.,   -l12, 0., 0. ;
    xi3  << 0., 1., 0.,   -(l11 + l12), 0., 0. ;
    xi4  << 0., 1., 0.,   -(l10 + l11 + l12), 0., 0. ;
    xi5  << 1., 0., 0.,   0., (l10 + l11 + l12), (l16) ;
    xi6  << 0., 0., 1.,   -l16, 0., 0. ;

    Matrix4d Hst0;
    Hst0 << 1, 0, 0, 0,
            0, 1, 0, -l16,
            0, 0, 1, l9+l10+l11+l12,
            0, 0, 0, 1;

    Matrix4d noap;

    noap = LieExp(xi6,angles(5))*LieExp(xi5,angles(4))*LieExp(xi4,angles(3))*LieExp(xi3,angles(2))*LieExp(xi2,angles(1))*LieExp(xi1,angles(0))*Hst0;

    return noap.block(0,3,3,1);

}

// Left leg with respect to ground
VectorXd ll_com_from_foot(Vector3d com, Vector3d eulerangles){
    VectorXd solution(6);
    solution << 1,2,3,4,5,6;
    VectorXd xi1(6), xi2(6), xi3(6), xi4(6), xi5(6), xi6(6);
    //          w              v
    xi1  << 1., 0., 0.,   0., l12, 0.;
    xi2  << 0., 1., 0.,   -l12, 0., 0. ;
    xi3  << 0., 1., 0.,   -(l11 + l12), 0., 0. ;
    xi4  << 0., 1., 0.,   -(l10 + l11 + l12), 0., 0. ;
    xi5  << 1., 0., 0.,   0., (l10 + l11 + l12), (l16) ;
    xi6  << 0., 0., 1.,   -l16, 0., 0. ;
     
    Matrix3d Rx;      
    Rx << 1,       0              ,           0           ,
          0,   cos(eulerangles(0)),   -sin(eulerangles(0)),
          0,   sin(eulerangles(0)),    cos(eulerangles(0));
    
    Matrix3d Ry;  
    Ry << cos(eulerangles(1)), 0 ,  sin(eulerangles(1)),
                   0         , 1 ,            0        ,
         -sin(eulerangles(1)), 0 ,  cos(eulerangles(1));
    
    Matrix3d Rz;
    Rz << cos(eulerangles(2)), -sin(eulerangles(2)) , 0,
          sin(eulerangles(2)),  cos(eulerangles(2)) , 0,
                  0          ,            0         , 1;
    Matrix3d R;
    R = Rz*Ry*Rx;
    
    Vector3d com2hip(-0.00867, -l13+0.00003, -0.01406);
    Vector3d hip2p = R*com2hip;

    Vector3d tool_p(com(0)+hip2p(0), com(1)+(l13-l16)+hip2p(1),com(2)+hip2p(2));
    
    MatrixXd noap(4,4); 
    noap << R(0,0),R(0,1),R(0,2), tool_p(0),
            R(1,0),R(1,1),R(1,2), tool_p(1),
            R(2,0),R(2,1),R(2,2), tool_p(2),
               0  ,   0  ,   0  ,     1    ;

    Matrix4d Hst0;
    Hst0 << 1, 0, 0, 0,
            0, 1, 0, -l16,
            0, 0, 1, l9+l10+l11+l12,
            0, 0, 0, 1;
            
    Matrix4d Hst0Inv;
    Hst0Inv << 1, 0, 0, 0,
               0, 1, 0, l16,
               0, 0, 1, -(l9+l10+l11+l12),
               0, 0, 0, 1;
    
    Vector3d p(0, -l16, l9+l10+l11+l12);
    Vector3d f(0, -l16,    l10+l11+l12);
    Vector3d r(0, -l16,        l11+l12);
    Vector3d k(0,    0,            l12);
    Vector3d o(0,    0,              0);
    
    float delta = ((((noap*Hst0Inv)*cart2hom(f)).head(3))-k).norm();
    Vector2d thetas3 = PK3(xi3,f,k,r,delta);

    float theta3;
    if(thetas3(0)>-0.2872664625997 && thetas3(0)<1.4381513036){
        theta3=thetas3(0);
    }else if (thetas3(1)>-0.2872664625997 && thetas3(1)<1.4381513036){
        theta3=thetas3(1);
    }else{
        //           1    2    3    4    5    6  
        solution << NAN, NAN, NAN, NAN, NAN, NAN;
        return solution;
    }

    
    
    Vector3d kp=(noap*Hst0Inv*cart2hom(f)).head(3);
    Vector3d fp=(LieExp(xi3,theta3)*cart2hom(f)).head(3);
    Vector4d thetas1and2 = PK2(xi1,xi2,fp,kp,k);
    float theta1;
    float theta2;
    
    if(thetas1and2(0)>-0.3473205211 && thetas1and2(0)<0.7417649321){
        theta1=thetas1and2(0);
        theta2=thetas1and2(2);
    }else if(thetas1and2(1)>-0.3473205211 && thetas1and2(1)<0.7417649321){
        theta1=thetas1and2(1);
        theta2=thetas1and2(3);
    }else{
        //           1    2     3      4    5    6  
        solution << NAN, NAN, theta3, NAN, NAN, NAN;
        return solution;
    }

    
    Vector3d kpp=(LieExp(xi3,-theta3)*LieExp(xi2,-theta2)*LieExp(xi1,-theta1)*noap*Hst0Inv*cart2hom(p)).head(3);
    Vector4d thetas4and5 = PK2(xi4,xi5,p,kpp,f);

    float theta4;
    float theta5;
    if(thetas4and5(0)>-0.5515240436 && thetas4and5(0)<0.7382742736){
        theta4=thetas4and5(0);
        theta5=thetas4and5(2);
    }else if(thetas4and5(1)>-0.5515240436 && thetas4and5(1)<0.7382742736){
        theta4=thetas4and5(1);
        theta5=thetas4and5(3);
    }else{
        //            1       2       3      4    5    6  
        solution << theta1, theta2, theta3, NAN, NAN, NAN;
        return solution;
    }
    
    
    Vector3d kppp=(LieExp(xi5,-theta5)*LieExp(xi4,-theta4)*LieExp(xi3,-theta3)*LieExp(xi2,-theta2)*LieExp(xi1,-theta1)*noap*Hst0Inv*cart2hom(o)).head(3);
    float theta6=PK1(xi6,o,kppp,p);
    if(theta6>-0.5724679947 && theta6<0.4869468613){
        solution << theta1, theta2, theta3, theta4, theta5, theta6;
        return solution;
    }else{
        //            1       2       3       4       5      6  
        solution << theta1, theta2, theta3, theta4, theta5, NAN;
        return solution;
    }

}


Vector3d rl_fk_com_from_foot(VectorXd angles){

    VectorXd xi1(6), xi2(6), xi3(6), xi4(6), xi5(6), xi6(6);
    //          w              v
    xi6  << 0., 0., 1., l16, 0., 0. ;
    xi5  << 1., 0., 0., 0., (l10 + l11 + l12), -(l16)  ;
    xi4  << 0., 1., 0.,-(l10 + l11 + l12), 0., 0. ;
    xi3 << 0., 1., 0.,-(l11 + l12), 0., 0. ;
    xi2 << 0., 1., 0.,-l12, 0., 0. ;
    xi1 << 1., 0., 0., 0., l12, 0. ;


    Matrix4d Hst0;
    Hst0 << 1, 0, 0, 0,
            0, 1, 0, -l16,
            0, 0, 1, l9+l10+l11+l12,
            0, 0, 0, 1;

    Matrix4d noap;

    noap = LieExp(xi6,angles(5))*LieExp(xi5,angles(4))*LieExp(xi4,angles(3))*LieExp(xi3,angles(2))*LieExp(xi2,angles(1))*LieExp(xi1,angles(0))*Hst0;

    Matrix3d R;
    Vector3d p;

    R =  noap.block(0,0,3,3);
    p =  noap.block(0,3,3,1);


    Vector3d com2hip(-0.00867, l13+0.00003, -0.01406);
    Vector3d hip2p = R*com2hip;

    Vector3d tool_p(p(0)+hip2p(0), p(1)-(l13-l16)+hip2p(1),p(2)+hip2p(2));

    return noap.block(0,3,3,1);

}

// Right leg with respect to to ground
VectorXd rl_com_from_foot(Vector3d com, Vector3d eulerangles){
    VectorXd solution(6);
    VectorXd xi7(6), xi8(6), xi9(6), xi10(6), xi11(6), xi12(6); 
    xi7  << 0., 0., 1., l16, 0., 0. ;
    xi8  << 1., 0., 0., 0., (l10 + l11 + l12), -(l16)  ;
    xi9  << 0., 1., 0.,-(l10 + l11 + l12), 0., 0. ;
    xi10 << 0., 1., 0.,-(l11 + l12), 0., 0. ;
    xi11 << 0., 1., 0.,-l12, 0., 0. ;
    xi12 << 1., 0., 0., 0., l12, 0. ;


    float theta12;
    float theta11;
    float theta10;
    float theta9;
    float theta8;
    float theta7;
           

    Matrix3d Rx;      
    Rx << 1,       0              ,           0           ,
          0,   cos(eulerangles(0)),   -sin(eulerangles(0)),
          0,   sin(eulerangles(0)),    cos(eulerangles(0));
    
    Matrix3d Ry;  
    Ry << cos(eulerangles(1)), 0 ,  sin(eulerangles(1)),
                   0         , 1 ,            0        ,
         -sin(eulerangles(1)), 0 ,  cos(eulerangles(1));
    
    Matrix3d Rz;
    Rz << cos(eulerangles(2)), -sin(eulerangles(2)) , 0,
          sin(eulerangles(2)),  cos(eulerangles(2)) , 0,
                  0          ,            0         , 1;
    Matrix3d R;
    R = Rz*Ry*Rx;
    
    Vector3d com2hip(-0.00867, l13+0.00003, -0.01406);
    Vector3d hip2p = R*com2hip;

    Vector3d tool_p(com(0)+hip2p(0), com(1)-(l13-l16)+hip2p(1),com(2)+hip2p(2));
    
    MatrixXd noap(4,4); 
    noap << R(0,0),R(0,1),R(0,2), tool_p(0),
            R(1,0),R(1,1),R(1,2), tool_p(1),
            R(2,0),R(2,1),R(2,2), tool_p(2),
               0  ,   0  ,   0  ,     1    ;

            
    Matrix4d Hst0;
    Hst0 << 1, 0, 0, 0,
            0, 1, 0, l16,
            0, 0, 1, l9+l10+l11+l12,
            0, 0, 0, 1;
            
    Matrix4d Hst0Inv;
    Hst0Inv << 1, 0, 0, 0,
               0, 1, 0, -l16,
               0, 0, 1, -(l9+l10+l11+l12),
               0, 0, 0, 1;


    Vector3d p(0,  l16, l9+l10+l11+l12);  
    Vector3d f(0,  l16,    l12+l11+l10);
    Vector3d r(0,  l16,        l12+l11);
    Vector3d k(0,    0,            l12);  
    Vector3d o(0,    0,              0);

    float delta = ((((noap*Hst0Inv)*cart2hom(f)).head(3))-k).norm();
    Vector2d thetas10 = PK3(xi10,f,k,r,delta);
    if(thetas10(0)>-0.2872664625997 && thetas10(0)<1.4381513036){
        theta10=thetas10(0);
    }else if (thetas10(1)>-0.2872664625997 && thetas10(1)<1.4381513036){
        theta10=thetas10(1);
    }else{
        //          12   11   10    9    8    7
        solution << NAN, NAN, NAN, NAN, NAN, NAN;
        return solution;
    }

    
    
    Vector3d kp=(noap*Hst0Inv*cart2hom(f)).head(3);
    Vector3d fp=(LieExp(xi10,theta10)*cart2hom(f)).head(3);
    Vector4d theta12and11 = PK2(xi12,xi11,fp,kp,k);
    if(theta12and11(0)>-0.7417649321 && theta12and11(0)<0.3473205211){
        theta12=theta12and11(0);
        theta11=theta12and11(2);
    }else if(theta12and11(1)>-0.7417649321 && theta12and11(1)<0.3473205211){
        theta12=theta12and11(1);
        theta11=theta12and11(3);
    }else{
        //          12   11     10      9    8    7
        solution << NAN, NAN, theta10, NAN, NAN, NAN;
        return solution;
    }
   

    Vector3d kpp=(LieExp(xi10,-theta10)*LieExp(xi11,-theta11)*LieExp(xi12,-theta12)*noap*Hst0Inv*cart2hom(p)).head(3);
    Vector4d theta9and8 = PK2(xi9,xi8,p,kpp,f);
    if(theta9and8(0)>-0.5515240436 && theta9and8(0)<0.7382742736){
        theta9=theta9and8(0);
        theta8=theta9and8(2);
    }else if(theta9and8(1)>-0.5515240436 && theta9and8(1)<0.7382742736){
        theta9=theta9and8(1);
        theta8=theta9and8(3);
    }else{
        //            12        11       10     9    8    7
        solution << theta12, theta11, theta10, NAN, NAN, NAN;
        return solution;
    }

    
    Vector3d kppp=(LieExp(xi8,-theta8)*LieExp(xi9,-theta9)*LieExp(xi10,-theta10)*LieExp(xi11,-theta11)*LieExp(xi12,-theta12)*noap*Hst0Inv*cart2hom(o)).head(3);
    theta7=PK1(xi7,o,kppp,p);
    if(theta7>-0.4869468613 && theta7<0.5724679947){
        solution << theta12, theta11, theta10, theta9, theta8, theta7;
        return solution;
    }else{
        //             12      11       10       9       8      7
        solution << theta12, theta11, theta10, theta9, theta8, NAN;
        return solution;
    }
}


PYBIND11_MODULE(screwIK, m) {
    m.doc() = "Screw-based Inverse Kinematics for Teo Robot"; // optional module docstring

    m.def("skewsim", &skewsim, "Hat operator");
    m.def("LieExp", &LieExp, "Converts a unit twist and an a twist angle into a 4x4 rigid transformation matrix");
    m.def("LieAd", &LieAd, "Converts a 4x4 transformation matrix into its adjoint 6x6 form");
    m.def("LieAdInv", &LieAdInv, "Converts a 4x4 transformation matrix into the adjoint 6x6 form of its inverse");
    m.def("cart2hom", &cart2hom, "Appends 1 to a vector to obtain the hoogeneous representation of a point");
    m.def("PK1", &PK1, "Paden-Kahan subproblem 1");
    m.def("PK2", &PK2, "Paden-Kahan subproblem 2");
    m.def("PK3", &PK3, "Paden-Kahan subproblem 3");
    m.def("ll_com_from_foot", &ll_com_from_foot, "Inverse kinematics solution from the left foot to the hip");
    m.def("rl_com_from_foot", &rl_com_from_foot, "Inverse kinematics solution from the right foot to the hip");
    m.def("ll_fk_com_from_foot", &ll_fk_com_from_foot, "Forward kinematics solution from the left foot to the hip");
    m.def("rl_fk_com_from_foot", &rl_fk_com_from_foot, "Forward kinematics solution from the right foot to the hip");
}
