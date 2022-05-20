//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_ISOHRONDERIVATIVE_H
#define FIRSTVERSION_ISOHRONDERIVATIVE_H
#include "MultMatrix/matrix.hpp"
class IsohronDerivative{
    Matrix dXdP;
    double dFdx, dFdy, dFdz;
public:
    Matrix setIntialValues(){
        /* df = {v1, v2, v3, a1, a1, a3}
         * dx = {x, y, z, vx, vy, vz}*/
        dXdP = {{1f, 0f, 0f, 0f, 0f, 0f,   0f},
                {0f, 1f, 0f, 0f, 0f, 0f,   0f},
                {0f, 0f, 1f, 0f, 0f, 0f,   0f},
                {0f, 0f, 0f, 1f, 0f, 0f,   0f},
                {0f, 0f, 0f, 0f, 1f, 0f,   0f},
                {0f, 0f, 0f, 0f, 0f, 1f,   0f}};
        return dXdP
    }
 /*X(t) -> f(X(t))
  *
  * dxdp -> (dx/dp)' = df/dx*dxdp*/


};
#endif //FIRSTVERSION_ISOHRONDERIVATIVE_H
