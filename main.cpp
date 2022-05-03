#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include "SFML/Graphics.hpp"
#include <fstream>
#include "MultMatrix/matrix.hpp"
const double G = 0.01720209895;
const double mBlackHole = G*G*4000000;
const double PI = 4*atan(1.);
double PointOfMid = 0;
class StarObject{
    double x, y, z;
    double xDot, yDot, zDot;
    double mass;
public:
    StarObject(double i, double j, double w, double iDot, double jDot, double wDot,  double m){
        x = i;
        y = j;
        z = w;
        xDot = jDot;
        yDot = iDot;
        zDot = wDot;
        mass = m;
    }
    double returnX(){return x;}
    double returnY(){return y;}
    double returnZ(){return z;}
    double returnXd(){return xDot;}
    double returnYd(){return yDot;}
    double returnZd(){return zDot;}
    double returnM(){return mass;}

    void update(std::vector<double>s){
        x = s[0];
        y = s[1];
        z = s[2];
        xDot = s[3];
        yDot = s[4];
        zDot = s[5];
    }

};

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
        system.push_back(stellarObjects[i]->returnX());
        system.push_back(stellarObjects[i]->returnY());
        system.push_back(stellarObjects[i]->returnZ());
        system.push_back(stellarObjects[i]->returnXd());
        system.push_back(stellarObjects[i]->returnYd());
        system.push_back(stellarObjects[i]->returnZd());
        mass.push_back(stellarObjects[i]->returnM());
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
//    for(int i = 0; i < system.size(); i++){
//        std::cout<<system[i]<<"\t";
//    }
//    std::cout<<"\n";


}

void projection(StarObject* object ,double OMEGA, double i ){
    Matrix R = Matrix({{object->returnX()},
                            {object->returnY()},
                            {object->returnZ()}});
    Matrix V = Matrix({{object->returnXd()},
                       {object->returnYd()},
                       {object->returnZd()}});
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
    double Omega2 = (234.50*PI)/180, Omega55 = (129.9*PI)/180, Omega38 = (101.8*PI)/180;
    double i2 = (136.78*PI)/180, i38 = (166.22*PI)/180, i55 = (141.7*PI)/180;
    double v_p=sqrt(mBlackHole*(1+0.884)/(1032.59*(1-0.884)));
    double v_p_38 = sqrt(mBlackHole*(1+0.818)/(1147.51*(1-0.818)));
    double v_p_55 = sqrt(mBlackHole*(1+0.74)/(892.32*(1-0.74)));
    double r_p=(1-0.884)*1032.59;
    double r_p_38 = (1-0.818)*1147.51;
    double r_p_55 = (1-0.74)*892.32;
//    s.push_back(new StarObject(  0, r_p, r_p*sin(PI*139.568327/180),  0, v_p, -v_p*sin(PI*139.568327/180), (14*2*pow(10, 30))));
//
//    projection(s[s.size()-1], Omega2, i2);
//
//    s.push_back(new StarObject(r_p_38*cos(PI*133/180)*cos(PI*166.22/180), r_p_38*sin(PI*133/180)*cos(PI*166.22/180), r_p_38*sin(PI*166.22/180),
//                               v_p_38*cos(PI*(133)/180)*cos(PI*166.22/180), v_p_38*sin(PI*(133)/180)*cos(PI*166.22/180), v_p_38*sin(PI*166.22/180), 0));
//
//    projection(s[s.size()-1], Omega38, i38);
//
//    s.push_back(new StarObject(r_p_55*cos(PI*(105)/180)*cos(PI*141.7/180), r_p_55*sin(PI*(105)/180)*cos(PI*141.7/180), r_p_55*sin(PI*141.7/180),
//                               v_p_55*cos(PI*(105)/180)*cos(PI*141.7/180), v_p_55*sin(PI*(105)/180)*cos(PI*141.7/180), v_p_55*sin(PI*141.7/180), 0));
//    projection(s[s.size()-1], Omega55, i55);

    s.push_back(new StarObject(  0, r_p*cos(PI*139.568327/180), r_p*sin(PI*139.568327/180),  0, v_p*cos(PI*139.568327/180), -v_p*sin(PI*139.568327/180), (14*2*pow(10, 30))));

    //projection(s[s.size()-1], Omega2, i2);

    s.push_back(new StarObject(r_p_38*cos(PI*52.96/180)*cos(PI*166.22/180), r_p_38*sin(PI*52.96/180)*cos(PI*166.22/180), r_p_38*sin(PI*166.22/180),
                               v_p_38*cos(PI*(52.96-90)/180)*cos(PI*166.22/180), v_p_38*sin(PI*(52.96-90)/180)*cos(PI*166.22/180), v_p_38*sin(PI*166.22/180), 0));

    //projection(s[s.size()-1], Omega38, i38);

    s.push_back(new StarObject(r_p_55*cos(PI*(-62.14)/180)*cos(PI*141.7/180), r_p_55*sin(PI*(-62.14)/180)*cos(PI*141.7/180), r_p_55*sin(PI*141.7/180),
                               v_p_55*cos(PI*(-62.14-90)/180)*cos(PI*141.7/180), v_p_55*sin(PI*(-62.14-90)/180)*cos(PI*141.7/180), v_p_55*sin(PI*141.7/180), 0));
    //projection(s[s.size()-1], Omega55, i55);
}

//s.push_back(new StarObject(  -793.684467559686, 278.3355419173991, 666.948420293991, 0.9699447417658503, 1.7714776457109187, 0.3239411986816966, (14*2*pow(10, 30))));



class Draw{
    std::vector<sf::CircleShape*> stObjects;
    std::vector<std::vector<sf::Vertex*>> orbits;
    int cyclicVar = 2000;
public:
    Draw(std::vector<StarObject*>obj){
        for(int i = 0; i < obj.size(); i++){;
            std::vector<sf::Vertex*> tmp;
            orbits.push_back(tmp);
        }
        for(int i = 0; i < obj.size(); i++){
            stObjects.push_back(new sf::CircleShape());
            stObjects[i]->setRadius(100);
            stObjects[i]->setOrigin(stObjects[i]->getRadius()/2, stObjects[i]->getRadius()/2);
            stObjects[i]->setPosition(obj[i]->returnX(), obj[i]->returnY());
            stObjects[i]->setFillColor(sf::Color::White);
            sf::Color orbColor(100*i, 255, 50*i);
            orbits[i].push_back(new sf::Vertex(sf::Vector2f(obj[i]->returnX(), obj[i]->returnY()), orbColor));
        }
    }

    void setObjects(sf::RenderWindow& w, std::vector<StarObject*>objects){
        if(cyclicVar==0){
            cyclicVar = 2000;
            for(int i = 0; i < orbits.size(); i++){
                orbits[i].erase(orbits[i].begin(), orbits[i].begin()+2000);
            }
        }
        else{
            cyclicVar--;
        }

        for(int i = 0; i < objects.size(); i++){
            sf::Color orbColor(100*i, 255, 50*i);
            orbits[i].push_back(new sf::Vertex(sf::Vector2f(objects[i]->returnX(), objects[i]->returnY()), orbColor));
            stObjects[i]->setPosition(objects[i]->returnX(), objects[i]->returnY());
            w.draw(*stObjects[i]);
            for(int j = 0; j < orbits[i].size(); j++){
                w.draw(orbits[i][j], 1, sf::Points);
            }
        }
    }
};


class StarStateInterpolator{
    std::ofstream toFileS2;
    std::ofstream toFileS38;
    std::ofstream toFileS55;
    std::ifstream fromFileS2;
    std::ifstream fromFileS38;
    std::ifstream fromFileS55;

    int howMuchSpaces(std::string string){
        int counter = 0;
        for(int i = 0; i < string.length(); i++){
            if(string[i] == ' ' || string[i]=='\t' || string[i] == '\n'){
                counter++;
            }
        }
        return counter;
    }

    void clean(std::ofstream& toF, std::ifstream& fromF, std::string path){
        std::string tmp;
        std::string result;
        const char del = '\0';
        while(!fromF.eof()){
            getline(fromF, tmp, del);
            if(howMuchSpaces(tmp)!=4){
                continue;
            }
            else{
                result.append(tmp);
            }
        }
        toF.close();
        toF.open(path);
        toF<<result;
    }

public:
    StarStateInterpolator(){
        toFileS2.open("../Data/S2.dat");
        toFileS38.open("../Data/S38.dat");
        toFileS55.open("../Data/S55.dat");

    }

    std::vector<double> interpolation(int t1, int t2){
        std::vector<double> coordinates;





        return coordinates;
    }


    void cleanLast(){
        fromFileS2.open("../Data/S2.dat");
        fromFileS38.open("../Data/S38.dat");
        fromFileS55.open("../Data/S55.dat");
        clean(toFileS2, fromFileS2, "../Data/S2.dat");
        clean(toFileS38, fromFileS38, "../Data/S38.dat");
        clean(toFileS55, fromFileS55, "../Data/S55.dat");

    }

    void addS2Data(StarObject* objects, int h, int n){
        switch (n) {
            case 2:
                toFileS2<<h;
                toFileS2<<"\t";
                toFileS2<<objects->returnX();
                toFileS2<<" ";
                toFileS2<<objects->returnY();
                toFileS2<<" ";
                toFileS2<<objects->returnZ();
                toFileS2<<"\n"<<'\0';
                break;
            case 38:
                toFileS38<<h;
                toFileS38<<"\t";
                toFileS38<<objects->returnX();
                toFileS38<<" ";
                toFileS38<<objects->returnY();
                toFileS38<<" ";
                toFileS38<<objects->returnZ();
                toFileS38<<"\n"<<'\0';
                break;
            case 55:
                toFileS55<<h;
                toFileS55<<"\t";
                toFileS55<<objects->returnX();
                toFileS55<<" ";
                toFileS55<<objects->returnY();
                toFileS55<<" ";
                toFileS55<<objects->returnZ();
                toFileS55<<"\n"<<'\0';
                break;
            default:
                break;

        }

    }


};



int main(){
    std::vector<StarObject*> system;
    initiation(system);
    Draw* draw = new Draw(system);
    StarStateInterpolator* interp = new StarStateInterpolator();
    sf::RenderWindow window(sf::VideoMode(1920, 1080), "Tides");
    sf::View view;
    view.setSize(38400-1920*12, 21600-1080*12);
    window.setView(view);
    int i = 0;
    //int lifeCycle = 19;
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

        i += 1;
        if (i % 500==0 && i!=0) {
            //lifeCycle--;
            interp->addS2Data(system[0], i, 2);
            interp->addS2Data(system[1], i, 38);
            interp->addS2Data(system[2], i, 55);
        }
//        if (lifeCycle == 0){
//            interp->cleanLast();
//            return 0;
//        }
    }

    interp->cleanLast();
    return 0;
}
