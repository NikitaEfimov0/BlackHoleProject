//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_ISOHRONDERIVATIVE_H
#define FIRSTVERSION_ISOHRONDERIVATIVE_H
#include "MultMatrix/matrix.hpp"
#include "TmpStarObject.h"
class IsohronDerivative{
    std::vector<std::pair<int, Matrix>>allDeriv;
public:
    Matrix dXdPRes = Matrix(12, 7);
    double dvxdm(double x, double y, double z, double m) {
        return (-x*m*m)/(pow(sqrt(x*x+y*y+z*z), 3));
    }

    double dvydm(double x, double y, double z, double m) {
        return (-y*m*m)/(pow(sqrt(x*x+y*y+z*z), 3));;
    }

    double dvzdm(double x, double y, double z, double m) {
        return (-z*m*m)/(pow(sqrt(x*x+y*y+z*z), 3));;
    }

    void updateMatrix(Star s38, Star s55,  double m, Matrix dXdP, double G){

        Matrix dFdGm = Matrix({{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,dvxdm(s38.x, s38.y, s38.z, G)},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,dvydm(s38.x, s38.y, s38.z, G)},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,dvzdm(s38.x, s38.y, s38.z, G)},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,dvxdm(s55.x, s55.y, s55.z, G)},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,dvydm(s55.x, s55.y, s55.z, G)},
                               {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,dvzdm(s55.x, s55.y, s55.z, G)}});



//        dXdP.DebugPrint();
//        std::cout<<"\n*\n";
        Matrix dFdX = Matrix({
                {0.0f,                  0.0f,                 0.0f,  0.0f, 0.0f, 0.0f,       1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {0.0f,                  0.0f,                 0.0f,  0.0f, 0.0f, 0.0f,       0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {0.0f,                  0.0f,                 0.0f,  0.0f, 0.0f, 0.0f,       0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f},
                {0.0f,                  0.0f,                 0.0f,  0.0f, 0.0f, 0.0f,       0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f},
                {0.0f,                  0.0f,                 0.0f,  0.0f, 0.0f, 0.0f,       0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f},
                {0.0f,                  0.0f,                 0.0f,  0.0f, 0.0f, 0.0f,       0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f},
                {dvxdx(s38.x, s38.y, s38.z, m), dvxdy(s38.x, s38.y, s38.z, m), dvxdz(s38.x, s38.y, s38.z, m), dvxdx(s55.x, s55.y, s55.z, m), dvxdy(s55.x, s55.y, s55.z, m), dvxdz(s55.x, s55.y, s55.z, m), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {dvydx(s38.x, s38.y, s38.z, m), dvydy(s38.x, s38.y, s38.z, m), dvydz(s38.x, s38.y, s38.z, m), dvydx(s55.x, s55.y, s55.z, m), dvydy(s55.x, s55.y, s55.z, m), dvydz(s55.x, s55.y, s55.z, m), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {dvzdx(s38.x, s38.y, s38.z, m), dvzdy(s38.x, s38.y, s38.z, m), dvzdz(s38.x, s38.y, s38.z, m), dvzdx(s55.x, s55.y, s55.z, m), dvzdy(s55.x, s55.y, s55.z, m), dvzdz(s55.x, s55.y, s55.z, m), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {dvxdx(s38.x, s38.y, s38.z, m), dvxdy(s38.x, s38.y, s38.z, m), dvxdz(s38.x, s38.y, s38.z, m), dvxdx(s55.x, s55.y, s55.z, m), dvxdy(s55.x, s55.y, s55.z, m), dvxdz(s55.x, s55.y, s55.z, m), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {dvydx(s38.x, s38.y, s38.z, m), dvydy(s38.x, s38.y, s38.z, m), dvydz(s38.x, s38.y, s38.z, m), dvydx(s55.x, s55.y, s55.z, m), dvydy(s55.x, s55.y, s55.z, m), dvydz(s55.x, s55.y, s55.z, m), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f},
                {dvzdx(s38.x, s38.y, s38.z, m), dvzdy(s38.x, s38.y, s38.z, m), dvzdz(s38.x, s38.y, s38.z, m), dvzdx(s55.x, s55.y, s55.z, m), dvzdy(s55.x, s55.y, s55.z, m), dvzdz(s55.x, s55.y, s55.z, m), 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f}});
//        dFdX.DebugPrint();
        dXdP = dFdGm+dFdX*dXdP;
//        std::cout<<"\n=\n";


        dXdPRes = Matrix(dXdP);
//        dXdPRes.DebugPrint();
//        std::cout<<std::endl;
        //dXdPRes.DebugPrint();
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


    void save(int t){
        //dXdPRes.DebugPrint();
        allDeriv.push_back(std::pair<int, Matrix>(t, dXdPRes));
        //std::cout<<allDeriv.end()->second.GetRows()<<" "<<allDeriv.end()->second.GetCols();
        //allDeriv.end()->second.DebugPrint();
    }

    Matrix interpolate(double t){
        for(int i = 0; i < allDeriv.size(); i++){
            if(abs(allDeriv[i].first-t)<=2){
                Matrix Xt0(allDeriv[i].second);
                //Xt0.DebugPrint();
                Matrix Xt1(allDeriv[i+1].second);
                //Xt1.DebugPrint();

//                Matrix X = Xt0*(double)((double)(allDeriv[i+1].first-t)/(allDeriv[i+1].first-allDeriv[i].first))+
//                           Xt1*((double)(t-allDeriv[i].first)/(allDeriv[i+1].first-allDeriv[i].first));

                Matrix X = Xt0 + (Xt1-Xt0)*((t-allDeriv[i].first)/(allDeriv[i+1].first-allDeriv[i].first));
//                std::cout<<"X!!!(#(@)(@(!@((!(@!(@(!(@!: \n";
//                X.DebugPrint();
//                std::cout<<"\n\n\n";
                return X;
            }
        }
    }




 /*X(t) -> f(X(t))
  *
  * dxdp -> (dx/dp)' = df/dx*dxdp*/

//    double axdx(double x, double y, double z, double m){
//        return
//    }
};
#endif //FIRSTVERSION_ISOHRONDERIVATIVE_H
