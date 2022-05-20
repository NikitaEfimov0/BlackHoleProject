//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_ISOHRONDERIVATIVE_H
#define FIRSTVERSION_ISOHRONDERIVATIVE_H
#include "MultMatrix/matrix.hpp"
class IsohronDerivative{
    Matrix dXdP;
    Matrix dFdX;
public:
    Matrix setInitialValues(){
        /* df = {v1, v2, v3, a1, a1, a3}
         * dx = {x, y, z, vx, vy, vz}*/
        dXdP = {{1, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,   0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,   0.0f},
                {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f}};
        return dXdP;
    }

    Matrix updateMatrix(double x, double y, double z, double m){
        dFdX = {{0,                 0,                 0,                 1, 0, 0, 0},
                {0,                 0,                 0,                 0, 1, 0, 0},
                {0,                 0,                 0,                 0, 0, 1, 0},
                {dvxdx(x, y, z, m), dvxdy(x, y, z, m), dvxdz(x, y, z, m), 0, 0, 0, dvxdm(x, y, z, m)},
                {dvydx(x, y, z, m), dvydy(x, y, z, m), dvydz(x, y, z, m), 0, 0, 0, dvydm(x, y, z, m)},
                {dvzdx(x, y, z, m), dvzdy(x, y, z, m), dvzdz(x, y, z, m), 0, 0, 0, dvzdm(x, y, z, m)},
                {0,                 0,                 0,                 0, 0, 0, 1}};
        return dFdX*dXdP;
    }



    //Производные ускорения по X

    double dvxdx(double x, double y, double z, double m){
        return -((-2*m*pow(x, 4)+m*pow(y, 4)+m*pow(z, 4) - m*pow(x, 2)*pow(y, 2)-m*pow(x, 2)*pow(z, 2)+2*m*pow(y, 2)*pow(z, 2))
        /(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 3)));
    }


    double dvxdy(double x, double y, double z, double m){
        return (3*m*x*y)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }


    double dvxdz(double x, double y, double z, double m){
        return (3*m*x*z)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }

    double dvxdm(double x, double y, double z, double m){
        return -(x)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }



    //Производные ускорения по Y


    double dvydx(double x, double y, double z, double m){
        return (3*m*x*y)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }

    double dvydy(double x, double y, double z, double m){
        return -((m*pow(x, 4)-2*m*pow(y, 4)+m*pow(z, 4) - m*pow(x, 2)*pow(y, 2)+2*m*pow(x, 2)*pow(z, 2)-m*pow(y, 2)*pow(z, 2))
                 /(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 3)));
    }

    double dvydz(double x, double y, double z, double m){
        return (3*m*z*y)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }

    double dvydm(double x, double y, double z, double m){
        return -(y)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }




    //Производные ускорения по Z


    double dvzdx(double x, double y, double z, double m){
        return (3*m*z*x)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));

    }

    double dvzdy(double x, double y, double z, double m){
        return (3*m*z*y)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));

    }

    double dvzdz(double x, double y, double z, double m){
        return -((m*pow(x, 4)+m*pow(y, 4)-2*m*pow(z, 4) +2*m*pow(x, 2)*pow(y, 2)-m*pow(x, 2)*pow(z, 2)-m*pow(y, 2)*pow(z, 2))
                 /(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 3)));
    }

    double dvzdm(double x, double y, double z, double m){
        return -(z)/(pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 2)*pow((pow(x, 2)+pow(y, 2)+pow(z, 2)), 0.5));
    }

 /*X(t) -> f(X(t))
  *
  * dxdp -> (dx/dp)' = df/dx*dxdp*/

//    double axdx(double x, double y, double z, double m){
//        return
//    }
};
#endif //FIRSTVERSION_ISOHRONDERIVATIVE_H
