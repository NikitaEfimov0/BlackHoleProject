//
// Created by Никита Ефимов on 24.05.2022.
//

#ifndef FIRSTVERSION_SOLVINGSYSTEM_H
#define FIRSTVERSION_SOLVINGSYSTEM_H
#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include "SFML/Graphics.hpp"
#include <fstream>
#include <sstream>
#include "MultMatrix/matrix.hpp"
#include "StarObject.h"
#include "Draw.h"
#include "StarStateInterpolator.h"
#include "IsohronDerivative.h"



class SolvingSystem{
    const double G = 0.01720209895;
    double mBlackHole = G*G;
    const double PI = 4*atan(1.);
    double PointOfMid = 0;
public:
    double norm(double x1, double y1, double z1, double  x2, double y2, double z2){
        return sqrt(pow((x1-x2), 2)+pow((y1-y2), 2)+pow((z1-z2), 2));
    }


    void derivative(std::vector<double>X, std::vector<double>&Xdot, IsohronDerivative *isohronDerivative, Matrix& dXdPNew, Matrix &dXdP){
        double m = mBlackHole;
        for(int i = 0; i < X.size(); i++){
            Xdot.push_back(0);
        }
        Xdot[0] = X[3];
        Xdot[1] = X[4];
        Xdot[2] = X[5];
        Xdot[3] =- X[0]*((mBlackHole)/(pow(norm(X[0], X[1], X[2], 0, 0, 0), 3)));
        Xdot[4] =- X[1]*((mBlackHole)/(pow(norm(X[0], X[1], X[2], 0, 0, 0), 3)));
        Xdot[5] =- X[2]*((mBlackHole)/(pow(norm(X[0], X[1], X[2], 0, 0, 0), 3)));
        Xdot[6] = X[9];
        Xdot[7] = X[10];
        Xdot[8] = X[11];
        Xdot[9] =- X[6]*((mBlackHole)/(pow(norm(X[6], X[7], X[8], 0, 0, 0), 3)));
        Xdot[10] =- X[7]*((mBlackHole)/(pow(norm(X[6], X[7], X[8], 0, 0, 0), 3)));
        Xdot[11] =- X[8]*((mBlackHole)/(pow(norm(X[6], X[7], X[8], 0, 0, 0), 3)));
        Xdot[12] = X[15];
        Xdot[13] = X[16];
        Xdot[14] = X[17];
        Xdot[15] =- X[12]*((mBlackHole)/(pow(norm(X[12], X[13], X[14], 0, 0, 0), 3)));
        Xdot[16] =- X[13]*((mBlackHole)/(pow(norm(X[12], X[13], X[14], 0, 0, 0), 3)));
        Xdot[17] =- X[14]*((mBlackHole)/(pow(norm(X[12], X[13], X[14], 0, 0, 0), 3)));

        isohronDerivative->updateMatrix(X[0], X[1], X[2], mBlackHole, dXdP);
        dXdPNew = Matrix(isohronDerivative->dXdPRes);
    }

    void set_tmp(std::vector<double>&tmp, std::vector<double>state, std::vector<double>k, double h, Matrix& dXdP, Matrix& kM, Matrix& tmpM)
    {
        for(int i = 0; i < state.size(); i++)
        {
            tmp[i] = (state[i]+h/2*k[i]);
        }

//    dXdP.DebugPrint();
//    std::cout<<"\n";
//    kM.DebugPrint();

        for(int i = 0; i < dXdP.GetRows(); i++){
            for(int  j = 0; j < dXdP.GetCols(); j++){
                tmpM.data[i][j] = dXdP.data[i][j]+h/2*kM.data[i][j];
            }
        }


    }


    void updateStates(std::vector<StarObject*>&stellarObjects, std::vector<double>system){
        std::vector<double>tmp;
        int iter = 0;
        for(int i = 0; i < stellarObjects.size(); i++){
            for(int j = iter; j < system.size()/stellarObjects.size()+iter; j++){
                tmp.push_back(system[j]);
            }
            iter+=6;
            stellarObjects[i]->update(tmp);
            tmp.clear();
        }
    }

    void RK4(std::vector<StarObject*>stellarObjects, IsohronDerivative *isohronDerivative, Matrix &dXdP){;
        std::vector<double>system;
        std::vector<double>mass;
        double h = 10;
        std::vector<double>tmp;
        for(int i = 0; i < stellarObjects.size(); i++){
            system.push_back(stellarObjects[i]->X());
            system.push_back(stellarObjects[i]->Y());
            system.push_back(stellarObjects[i]->Z());
            system.push_back(stellarObjects[i]->dX());
            system.push_back(stellarObjects[i]->dY());
            system.push_back(stellarObjects[i]->dZ());
            mass.push_back(stellarObjects[i]->M());
        }

        for(int i = 0; i < system.size(); i++){
            tmp.push_back(0);
        }
        std::vector<double>k1, k2, k3, k4;
        Matrix kM1 = Matrix(7, 6), kM2= Matrix(7, 6), kM3= Matrix(7, 6), kM4= Matrix(7, 6), tmpM = Matrix(7, 6);
        derivative(system, k1, isohronDerivative, kM1, dXdP);
        set_tmp(tmp, system, k1, h, dXdP, kM1, tmpM);
        derivative(tmp, k2, isohronDerivative, kM2, kM1);

        set_tmp(tmp, system, k2, h, dXdP, kM2, tmpM);
        derivative(tmp, k3, isohronDerivative, kM3, kM2);

        set_tmp(tmp, system, k3, 2*h, dXdP, kM3, tmpM);
        derivative(tmp, k4, isohronDerivative, kM4, kM3);
        //isohronDerivative->dXdPRes.DebugPrint();
        for(int i = 0; i < system.size(); i++){
            system[i]+=(h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]));
        }

        for(int i = 0; i < dXdP.GetRows(); i++){
            for(int j = 0; j < dXdP.GetCols(); j++){
                dXdP.data[i][j]+=(h/6*(kM1.data[i][j]+2*kM2.data[i][j]+2*kM3.data[i][j]+kM4.data[i][j]));
            }
        }
        isohronDerivative->dXdPRes = Matrix(dXdP);
        //isohronDerivative->dXdPRes.DebugPrint();
        updateStates(stellarObjects, system);



    }

    void projection(StarObject* object ,double OMEGA, double i ){
        Matrix R = Matrix({{object->X()},
                           {object->Y()},
                           {object->Z()}});
        Matrix V = Matrix({{object->dX()},
                           {object->dY()},
                           {object->dZ()}});
        Matrix A2 = Matrix(
                {{cos(OMEGA+PointOfMid), -cos(i)* sin(OMEGA+PointOfMid), sin(i)* sin(OMEGA+PointOfMid)},
                 {sin(OMEGA+PointOfMid), cos(i)* cos(OMEGA+PointOfMid), -sin(i)* cos(OMEGA+PointOfMid)},
                 {0, sin(i), cos(i)}}
        );

        Matrix RESULT = A2*R;
        // RESULT.DebugPrint();
        std::vector<double>s;
        s.push_back(RESULT.data[0][0]);
        s.push_back(RESULT.data[1][0]);
        s.push_back(RESULT.data[2][0]);
        RESULT = V;
        //RESULT.DebugPrint();
        s.push_back(RESULT.data[0][0]);
        s.push_back(RESULT.data[1][0]);
        s.push_back(RESULT.data[2][0]);
        object->update(s);
    }

    void initiation(std::vector<StarObject*>&s, Matrix &B, Matrix& dXdP){
        double Omega2 = (240.50*PI)/180, Omega55 = (129.9*PI)/180, Omega38 = (101.8*PI)/180;
        double i2 = (136.78*PI)/180, i38 = (166.22*PI)/180, i55 = (141.7*PI)/180;

        dXdP = Matrix({
                              {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                              {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                              {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,   0.0f},
                              {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,   0.0f},
                              {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,   0.0f},
                              {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f}});


        s.push_back(new StarObject(  120.451454,  -22.675722,       -104.524315,      -0.556251   ,      -3.6, 0.0, (14*2*pow(10, 30))));

        projection(s[s.size()-1], Omega2, i2);

        s.push_back(new StarObject(-22.146914,    207.074722,      15.702321,       -3.191719    ,   0.341359 , 0., 0));

        projection(s[s.size()-1], Omega38, i38);

        s.push_back(new StarObject(204.046539,   45.089123,      100.784233,      0.642879   ,     2.909287   ,    0.000000, 0));
        projection(s[s.size()-1], Omega55, i55);
    }


    void start(Matrix &B){
        IsohronDerivative isohronDerivative = IsohronDerivative();
        std::vector<StarObject*> system;
        Matrix dXdP = Matrix(7, 6);
        initiation(system, B, dXdP);

        StarStateInterpolator* interp = new StarStateInterpolator();
        int i = 0;
        while (i!=3000) {
            interp->addS2Data(system[0], i, 2);
            interp->addS2Data(system[1], i, 38);
            interp->addS2Data(system[2], i, 55);
            RK4(system, &isohronDerivative, dXdP);
            isohronDerivative.save(i);
            i += 2;
            if(i == 3000)
                break;
        }
        //isohronDerivative.printMatrixdXdP();
        std::cout<<"\n\n\n";
        interp->cleanLast();
    }


};


#endif //FIRSTVERSION_SOLVINGSYSTEM_H
