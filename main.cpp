#include <iostream>
#include <vector>
#include <cmath>
#include <unistd.h>
#include "SFML/Graphics.hpp"
#include <fstream>
#include <sstream>
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
    double Omega2 = (240.50*PI)/180, Omega55 = (129.9*PI)/180, Omega38 = (101.8*PI)/180;
    double i2 = (136.78*PI)/180, i38 = (166.22*PI)/180, i55 = (141.7*PI)/180;


    s.push_back(new StarObject(  120.451454,  -22.675722,       -104.524315,      -0.556251   ,      -3.6, 0.0, (14*2*pow(10, 30))));

    projection(s[s.size()-1], Omega2, i2);

    s.push_back(new StarObject(-22.146914,    207.074722,      15.702321,       -3.191719    ,   0.341359 , 0., 0));

    projection(s[s.size()-1], Omega38, i38);

    s.push_back(new StarObject(204.046539,   45.089123,      100.784233,      0.642879   ,     2.909287   ,    0.000000, 0));
    projection(s[s.size()-1], Omega55, i55);
}



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
            sf::Color orbColor(70*(i+1), 255/(i+1), 255/(i+2));
//            sf::Color orbColor(255, 255, 255);
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
            sf::Color orbColor(70*(i+1), 255/(i+1), 255/(i+2));
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
    std::vector<std::vector<double>> parseFile(int t, int n){
        std::vector<std::vector<double>>states;
        std::string tmp;
        std::string tmpNext;
        std::string delims = " ";
        std::vector<std::string>strings;
        const char del = '\0';
        double diff = 0.0;
        int b, e = 0;
        switch(n){
            case 2:
                fromFileS2.open("../Data/S2.dat");

                while(!fromFileS2.eof()){
                    getline(fromFileS2, tmp, del);
                    std::istringstream iss(tmp);
                    std::string w;
                    while (iss >> w) strings.push_back(w);

                    std::istringstream stream(strings[0]);
                    int h;
                    stream>>h;
                    diff = abs(h-t);
                    if(diff>=23){
                        std::vector<double>firstPos;
                        std::vector<double>secondPos;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            firstPos.push_back(value);
                        }
                        strings.clear();
                        getline(fromFileS2, tmpNext, del);
                        while ((b = tmpNext.find_first_not_of(delims, e)) != tmpNext.npos) {
                            e = tmpNext.find_first_of(delims, b);
                            strings.push_back(tmp.substr(b, e - b));
                            b = e;
                        }
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            secondPos.push_back(value);
                        }
                        states.push_back(firstPos);
                        states.push_back(secondPos);
                    }

                }

                break;
            case 38:
                while(!fromFileS38.eof()){
                    getline(fromFileS2, tmp, del);
                    while ((b = tmp.find_first_not_of(delims, e)) != tmp.npos) {
                        e = tmp.find_first_of(delims, b);
                        strings.push_back(tmp.substr(b, e - b));
                        b = e;
                    }

                    std::istringstream stream(strings[0]);
                    int h;
                    stream>>h;
                    diff = abs(h-t);
                    if(diff>=23){
                        std::vector<double>firstPos;
                        std::vector<double>secondPos;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            firstPos.push_back(value);
                        }
                        strings.clear();
                        getline(fromFileS2, tmpNext, del);
                        while ((b = tmpNext.find_first_not_of(delims, e)) != tmpNext.npos) {
                            e = tmpNext.find_first_of(delims, b);
                            strings.push_back(tmp.substr(b, e - b));
                            b = e;
                        }
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            secondPos.push_back(value);
                        }
                        states.push_back(firstPos);
                        states.push_back(secondPos);
                    }

                }
                break;
            case 55:
                while(!fromFileS2.eof()){
                    getline(fromFileS2, tmp, del);
                    while ((b = tmp.find_first_not_of(delims, e)) != tmp.npos) {
                        e = tmp.find_first_of(delims, b);
                        strings.push_back(tmp.substr(b, e - b));
                        b = e;
                    }

                    std::istringstream stream(strings[0]);
                    int h;
                    stream>>h;
                    diff = abs(h-t);
                    if(diff>=23){
                        std::vector<double>firstPos;
                        std::vector<double>secondPos;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            firstPos.push_back(value);
                        }
                        strings.clear();
                        getline(fromFileS2, tmpNext, del);
                        while ((b = tmpNext.find_first_not_of(delims, e)) != tmpNext.npos) {
                            e = tmpNext.find_first_of(delims, b);
                            strings.push_back(tmp.substr(b, e - b));
                            b = e;
                        }
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            secondPos.push_back(value);
                        }
                        states.push_back(firstPos);
                        states.push_back(secondPos);
                    }

                }
                break;
            default:
                break;
        }
        return states;
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
                tmp.append("\0");
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

    std::vector<double> interpolation(int t){
        std::vector<std::vector<double>> coordinatesS2 = parseFile(t, 2);
        std::vector<std::vector<double>> coordinatesS38 = parseFile(t, 38);
        std::vector<std::vector<double>> coordinatesS55 = parseFile(t, 55);
        for(int i = 0; i < coordinatesS2[0].size(); i++){
            std::cout<<coordinatesS2[0][i]<<" ";
        }





        //return coordinates;
    }


    void cleanLast(){
        fromFileS2.open("../Data/S2.dat");
        fromFileS38.open("../Data/S38.dat");
        fromFileS55.open("../Data/S55.dat");
        clean(toFileS2, fromFileS2, "../Data/S2.dat");
        clean(toFileS38, fromFileS38, "../Data/S38.dat");
        clean(toFileS55, fromFileS55, "../Data/S55.dat");
        fromFileS2.close();
        fromFileS38.close();
        fromFileS55.close();

    }


    void addS2Data(StarObject* objects, int h, int n){
            switch (n) {
                case 2:
                    toFileS2 << h;
                    toFileS2 << " ";
                    toFileS2 << objects->returnX()/8107.55245;
                    toFileS2 << " ";
                    toFileS2 << objects->returnY()/8107.55245;
                    toFileS2 << " ";
                    toFileS2 << objects->returnZ()/8107.55245;
                    toFileS2 << "\n" << '\0';
                    break;
                case 38:
                    toFileS38 << h;
                    toFileS38 << " ";
                    toFileS38 << objects->returnX()/8107.55245;
                    toFileS38 << " ";
                    toFileS38 << objects->returnY()/8107.55245;
                    toFileS38 << " ";
                    toFileS38 << objects->returnZ()/8107.55245;
                    toFileS38 << "\n" << '\0';
                    break;
                case 55:
                    toFileS55 << h;
                    toFileS55 << " ";
                    toFileS55 << objects->returnX()/8107.55245;
                    toFileS55 << " ";
                    toFileS55 << objects->returnY()/8107.55245;
                    toFileS55 << " ";
                    toFileS55 << objects->returnZ()/8107.55245;
                    toFileS55 << "\n" << '\0';
                    break;
                default:
                   return;

            }

    }


};



int main(){
    std::vector<StarObject*> system;
    initiation(system);
    float want_fps = 50;
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

//        sf::Int32 frame_duration = loop_timer.getElapsedTime().asMilliseconds();
//        sf::Int32 time_to_sleep = int(1000.f/want_fps) - frame_duration;
//        if (time_to_sleep > 0) {
//            sf::sleep(sf::milliseconds(time_to_sleep));
//        }
//        loop_timer.restart();

    }

    interp->cleanLast();
    //interp->interpolation(9);
    return 0;
}
