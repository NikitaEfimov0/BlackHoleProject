//
// Created by Никита Ефимов on 01.06.2022.
//

#ifndef FIRSTVERSION_TMPSTAROBJECT_H
#define FIRSTVERSION_TMPSTAROBJECT_H
struct Star{
    double x; double y; double z;
    double dx; double dy; double dz;

    Star(double i, double j, double w){
        x = i; y = j; z = w;
    }
};
#endif //FIRSTVERSION_TMPSTAROBJECT_H
