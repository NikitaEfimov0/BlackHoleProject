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
const double G = 0.01720209895;
const double mBlackHole = G*G*4000000;
const double PI = 4*atan(1.);
double PointOfMid = 0;


double norm(double x1, double y1, double z1, double  x2, double y2, double z2){
    return sqrt(pow((x1-x2), 2)+pow((y1-y2), 2)+pow((z1-z2), 2));
}


void derivative(std::vector<double>X, std::vector<double>&Xdot){
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
}

void set_tmp(std::vector<double>&tmp, std::vector<double>state, std::vector<double>k, double h)
{
    for(int i = 0; i < state.size(); i++)
    {
        tmp[i] = (state[i]+h/2*k[i]);
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

void RK4(std::vector<StarObject*>stellarObjects){
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
    derivative(system, k1);
    set_tmp(tmp, system, k1, h);
    derivative(tmp, k2);

    set_tmp(tmp, system, k2, h);
    derivative(tmp, k3);

    set_tmp(tmp, system, k3, 2*h);
    derivative(tmp, k4);

    for(int i = 0; i < system.size(); i++){
        system[i]+=(h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]));
    }
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

void initiation(std::vector<StarObject*>&s){
    double Omega2 = (240.50*PI)/180, Omega55 = (129.9*PI)/180, Omega38 = (101.8*PI)/180;
    double i2 = (136.78*PI)/180, i38 = (166.22*PI)/180, i55 = (141.7*PI)/180;


    s.push_back(new StarObject(  120.451454,  -22.675722,       -104.524315,      -0.556251   ,      -3.6, 0.0, (14*2*pow(10, 30))));

    projection(s[s.size()-1], Omega2, i2);

    s.push_back(new StarObject(-22.146914,    207.074722,      15.702321,       -3.191719    ,   0.341359 , 0., 0));

    projection(s[s.size()-1], Omega38, i38);

    s.push_back(new StarObject(204.046539,   45.089123,      100.784233,      0.642879   ,     2.909287   ,    0.000000, 0));
    projection(s[s.size()-1], Omega55, i55);
}









int main(){
    std::vector<StarObject*> system;
    initiation(system);
    sf::Clock loop_timer;
    Draw* draw = new Draw(system);
    StarStateInterpolator* interp = new StarStateInterpolator();
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Tides");
    sf::View view;
    view.setSize(38400-1920*12, 21600-1080*12);
    window.setView(view);
    int i = 0;

    while (window.isOpen()) {
        interp->addS2Data(system[0], i, 2);
        interp->addS2Data(system[1], i, 38);
        interp->addS2Data(system[2], i, 55);
        sf::Event event;
        if (sf::Keyboard::isKeyPressed(sf::Keyboard::Escape))window.close();
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)window.close();
        }
        window.clear();
        draw->setObjects(window, system);
        window.display();
        RK4(system);

        i += 23;
    }

    interp->cleanLast();
    interp->interpolation(50, 2);
    return 0;
}
