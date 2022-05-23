//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_ISOHRONDERIVATIVE_H
#define FIRSTVERSION_ISOHRONDERIVATIVE_H
#include "MultMatrix/matrix.hpp"
class IsohronDerivative{
public:


    double dvxdm(double x, double y, double z, double m) {
        return (-x)/(pow(sqrt(x*x+y*y+z*z), 3));
    }

    double dvydm(double x, double y, double z, double m) {
        return (-y)/(pow(sqrt(x*x+y*y+z*z), 3));;
    }

    double dvzdm(double x, double y, double z, double m) {
        return (-z)/(pow(sqrt(x*x+y*y+z*z), 3));;
    }

    void updateMatrix(double x, double y, double z, double m){
        Matrix dFdGm = Matrix({{0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, dvxdm(x, y, z, m)},
                               {0, 0, 0, 0, 0, 0, dvydm(x, y, z, m)},
                               {0, 0, 0, 0, 0, 0, dvzdm(x, y, z, m)}});

        static Matrix dXdP = Matrix({
                       {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                       {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                       {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                       {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,   0.0f},
                       {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,   0.0f},
                       {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f}});
        std::cout<<"before:\n";
        dXdP.DebugPrint();
        std::cout<<"\n";
        Matrix dFdX = Matrix({
                {0.0f,                 0.0f,                 0.0f,        1.0f, 0.0f, 0.0f},
                {0.0f,                 0.0f,                 0.0f,        0.0f, 1.0f, 0.0f},
                {0.0f,                 0.0f,                 0.0f,        0.0f, 0.0f, 1.0f},
                {dvxdx(x, y, z, m), dvxdy(x, y, z, m), dvxdz(x, y, z, m), 1.0f, 0.0f, 0.0f},
                {dvydx(x, y, z, m), dvydy(x, y, z, m), dvydz(x, y, z, m), 0.0f, 1.0f, 0.0f},
                {dvzdx(x, y, z, m), dvzdy(x, y, z, m), dvzdz(x, y, z, m), 0.0f, 0.0f, 1.0f}});
        dFdX.DebugPrint();
        dXdP = dFdGm+dFdX*dXdP;
        std::cout<<"after:\n";
        dXdP.DebugPrint();
        std::cout<<std::endl;
    }



    //Производные ускорения по X
    double dvxdx(double x, double y, double z, double m){
        return (((-m)/(pow((sqrt(x*x+y*y+z*z)), 3)))+((x*x*m*3)/(pow((sqrt(x*x+y*y+z*z)), 5))));
    }

    double dvxdy(double x, double y, double z, double m){
        return (3*x*y*m)/(pow((sqrt(x*x+y*y+z*z)), 5));
    }

    double dvxdz(double x, double y, double z, double m){
        return (3*x*z*m)/(pow((sqrt(x*x+y*y+z*z)), 5));
    }



    //Производные ускорения по Y

    double dvydx(double x, double y, double z, double m){
        return (3*x*y*m)/(pow((sqrt(x*x+y*y+z*z)), 5));
    }

    double dvydy(double x, double y, double z, double m){
        return (((-m)/(pow((sqrt(x*x+y*y+z*z)), 3)))+((y*y*m*3)/(pow((sqrt(x*x+y*y+z*z)), 5))));
    }

    double dvydz(double x, double y, double z, double m){
        return (3*z*y*m)/(pow((sqrt(x*x+y*y+z*z)), 5));
    }


    //Производные ускорения по Z

    double dvzdx(double x, double y, double z, double m){
        return (3*x*z*m)/(pow((sqrt(x*x+y*y+z*z)), 5));
    }

    double dvzdy(double x, double y, double z, double m){
        return (3*z*y*m)/(pow((sqrt(x*x+y*y+z*z)), 5));
    }

    double dvzdz(double x, double y, double z, double m){
        return (((-m)/(pow((sqrt(x*x+y*y+z*z)), 3)))+((z*z*m*3)/(pow((sqrt(x*x+y*y+z*z)), 5))));
    }







 /*X(t) -> f(X(t))
  *
  * dxdp -> (dx/dp)' = df/dx*dxdp*/

//    double axdx(double x, double y, double z, double m){
//        return
//    }
};
#endif //FIRSTVERSION_ISOHRONDERIVATIVE_H
