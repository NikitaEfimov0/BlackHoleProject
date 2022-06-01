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
#include "GaussNewton.h"
const double G = 0.01720209895;
double mBlackHole = G*G*4000000;
const double PI = 4*atan(1.);
double PointOfMid = 0;

double X, Y, Z, dX, dY, dZ;
double E(double e){
    double M = 0;
    double Z = (e*sin(M))/(sqrt(e*e+1.0-2.0*e*cos(M)));

    return M+Z-((pow(Z, 4)/6.0)*(1.0/tan(M)));
}


std::vector<double>convertToDecart(double a, double e,  double i, double w, double Ohm){
    double x = a*((cos(E(e))-e)*(cos(w)*cos(Ohm)-sin(w)*sin(Ohm)*cos(i))+((sqrt(1-pow(e, 2)))*sin(E(e)))*(-sin(w)*cos(Ohm)-cos(w)*sin(Ohm)*cos(i)));

    double y = a*((cos(E(e))-e)*(cos(w)*sin(Ohm)+sin(w)*cos(Ohm)*cos(i))+((sqrt(1-pow(e, 2)))*sin(E(e)))*(-sin(w)*sin(Ohm)+cos(w)*sin(Ohm)*cos(i)));

    double z = a*((cos(E(e))-e)*(sin(w)*sin(i))+(sqrt(1-pow(e, 2)))*sin(E(e))*(cos(w)*sin(i)));

    std::vector<double> res;
    res.push_back(x);
    res.push_back(y);
    res.push_back(z);

    return res;
}


std::vector<double> iauS2c(double theta, double phi)

{
    double cp;

    cp = cos(phi);
    return std::vector<double>({cos(theta) * cp*PI/180., sin(theta) * cp*PI/180., sin(phi)*PI/180.});

}




std::vector<double> S2pv(double theta, double phi, double r, double td, double pd, double rd)
    {
        double st, ct, sp, cp, rcp, x, y, rpd, w;


        st = sin(theta);
        ct = cos(theta);
        sp = sin(phi);
        cp = cos(phi);
        rcp = r * cp;
        x = ct;
        y = st;
        rpd = r * pd;
        w = rpd*sp - cp*rd;

        std::vector<double> res;
        res.push_back(x);
        res.push_back(y);
        res.push_back(sp);
        res.push_back(-y*td - w*ct);
        res.push_back(x*td - w*st);
        res.push_back(rpd*cp + sp*rd);
        return res;
}

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

    isohronDerivative->updateMatrix(X[0], X[1], X[2], mBlackHole, dXdP, G);
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
   // isohronDerivative->dXdPRes.DebugPrint();
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

void initiation(std::vector<StarObject*>&s, Matrix& dXdP){
    double Omega2 = (240.50*PI)/180, Omega55 = (129.9*PI)/180, Omega38 = (101.8*PI)/180;
    double i2 = (136.78*PI)/180, i38 = (166.22*PI)/180, i55 = (141.7*PI)/180;

    dXdP = Matrix({
        {1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
        {0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f,   0.0f},
        {0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f,   0.0f},
        {0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f,   0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f,   0.0f},
        {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f,   0.0f}});


    std::vector<double> s2 = convertToDecart(0.126, 0.884, 136.78, 71.36, 234.50);

    std::vector<double> s2Sph = iauS2c(0.0386, 0.0213);
    std::vector<double> s2V = iauS2c(0.0385, 0.0701);
    double r = sqrt(pow(s2Sph[0], 2)+pow(s2Sph[1], 2)+pow(s2Sph[2], 2));
    double r1 = sqrt(pow(s2V[0], 2)+pow(s2V[1], 2)+pow(s2V[2], 2));
    std::vector<double> finalS2 = S2pv(0.036, 0.0213, sqrt(pow(s2Sph[0], 2)+pow(s2Sph[1], 2)+pow(s2Sph[2], 2)), 0.0385-0.0386, 0.0701-0.0213, r1-r);

   //s.push_back(new StarObject(  120.451454,  -22.675722,       -104.524315,      -0.556251   ,      -3.6, 0.0, (14*2*pow(10, 30))));

    //s.push_back(new StarObject(finalS2[0], finalS2[1], finalS2[2], finalS2[3], finalS2[4], finalS2[5] , (14*2*pow(10, 30))));


     s.push_back(new StarObject(  s2Sph[0], s2Sph[1], s2Sph[2],      s2V[0]   ,s2V[1], s2V[2], (14*2*pow(10, 30))));

    //projection(s[s.size()-1], Omega2, i2);

    X = s[s.size()-1]->X();
    Y = s[s.size()-1]->Y();
    Z = s[s.size()-1]->Z();
    dX = s[s.size()-1]->dX();
    dY = s[s.size()-1]->dY();
    dZ = s[s.size()-1]->dZ();

    s.push_back(new StarObject(-22.146914,    207.074722,      15.702321,       -3.191719    ,   0.341359 , 0., 0));

    projection(s[s.size()-1], Omega38, i38);

    s.push_back(new StarObject(204.046539,   45.089123,      100.784233,      0.642879   ,     2.909287   ,    0.000000, 0));
    projection(s[s.size()-1], Omega55, i55);
}









int main(){
    IsohronDerivative isohronDerivative = IsohronDerivative();
    Matrix dXdP = Matrix(7, 6);
    std::vector<StarObject*> system;
    initiation(system, dXdP);
    sf::Clock loop_timer;
    Draw* draw = new Draw(system);
    StarStateInterpolator* interp = new StarStateInterpolator();
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Tides");
    sf::View view;
    view.setSize(38400-1920*12, 21600-1080*12);
    window.setView(view);
    int i = 0;

    while (/*window.isOpen()*/ i!=3000) {
        interp->addS2Data(system[0], i, 2);
        interp->addS2Data(system[1], i, 38);
        interp->addS2Data(system[2], i, 55);
//        sf::Event event;
//        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape))window.close();
//        while (window.pollEvent(event)) {
//            if (event.type == sf::Event::Closed)window.close();
//        }
//        window.clear();
//        draw->setObjects(window, system);
//        window.display();
        RK4(system, &isohronDerivative, dXdP);
        //sleep(2);
        isohronDerivative.save(i);
        i += 2;

    }
   //isohronDerivative.printMatrixdXdP();
    std::cout<<"\n\n\n";
    interp->cleanLast();
    GaussNewton gaussNewton = GaussNewton(mBlackHole/(G*G), X, Y, Z, dX, dY, dZ);
    //interp->interpolation(2004.580, 55);

   gaussNewton.findBlackHoleMass(interp, isohronDerivative);
   //gaussNewton.test();
    return 0;
}
