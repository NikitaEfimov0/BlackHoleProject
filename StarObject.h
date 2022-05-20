//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_STAROBJECT_H
#define FIRSTVERSION_STAROBJECT_H
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
    double X(){return x;}
    double Y(){return y;}
    double Z(){return z;}
    double dX(){return xDot;}
    double dY(){return yDot;}
    double dZ(){return zDot;}
    double M(){return mass;}

    void update(std::vector<double>s){
        x = s[0];
        y = s[1];
        z = s[2];
        xDot = s[3];
        yDot = s[4];
        zDot = s[5];
    }

};
#endif //FIRSTVERSION_STAROBJECT_H
