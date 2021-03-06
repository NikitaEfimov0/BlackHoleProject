//
// Created by Никита Ефимов on 23.05.2022.
//

#ifndef FIRSTVERSION_GAUSSNEWTON_H
#define FIRSTVERSION_GAUSSNEWTON_H
#include <iostream>
#include <fstream>
#include "StarStateInterpolator.h"
#include "SolvingSystem.h"
#include <utility>
#include <vector>
class GaussNewton{
    std::ifstream S2Original;
    std::ifstream S38Original;
    std::ifstream S55Original;
    double BlackHoleMass;

    double prevSRA = 0;
    double prevSDEC = 0;
    std::vector<std::pair<double, std::pair<double, double>>>rAll;
    std::vector<std::pair<double, double>>varAll;
    IsohronDerivative isohronDerivative1;
    Matrix Beta = Matrix(1, 13);

    Star S38, S55;

    //StarStateInterpolator starStateInterpolator;
public:
    GaussNewton(double m, Star s38, Star s55) : S38(s38), S55(s55) {

        BlackHoleMass = m;
        Beta = Matrix({{s38.x},
                       {s38.y},
                       {s38.z},
                       {s38.dx},
                       {s38.dy},
                       {s38.dz},
                       {s55.x},
                       {s55.y},
                       {s55.z},
                       {s55.dx},
                       {s55.dy},
                       {s55.dz},
                       {BlackHoleMass}});
    }

    double getLen(double x, double y, double z){
        return sqrt(x*x+y*y+z*z);
    }

    Matrix dGdXMatr(Star s38, Star s55){
        double sign = -1;
        double dRadX, dRadY, dRadZ, dDecdX, dDecdY, dDecdZ;
        double lenS38 = getLen(s38.x, s38.y, s38.z);
        double cosF = cos(asin(s38.z / lenS38));
        dRadX = sign*s38.x*s38.z/(sqrt(1-s38.z*s38.z/(lenS38 * lenS38)) * lenS38 * lenS38 * lenS38);
        dRadY = sign*s38.x*s38.z/(sqrt(1-s38.z*s38.z/(lenS38 * lenS38)) * lenS38 * lenS38 * lenS38);
        dRadZ = (pow(s38.z, 3)-pow(s38.z, 2))/((sqrt(1-pow(s38.z, 2)/pow(lenS38, 2)) * pow(lenS38, 3)));
//
        if(s38.y / lenS38 > 0){
            sign = 1;
        }

        dDecdX = sign * (-1*s38.x*(s38.z*dRadX- s38.x * cosF / lenS38) + lenS38 * cosF) / (pow(lenS38, 2) * pow(cosF, 2) * sqrt(1 - pow(s38.x, 2) / (pow(lenS38, 2) * pow(cosF, 2))));
        dDecdY = sign * (-1*s38.x*(s38.z*dRadY- s38.y * cosF / lenS38) + lenS38 * cosF) / (pow(lenS38, 2) * pow(cosF, 2) * sqrt(1 - pow(s38.x, 2) / (pow(lenS38, 2) * pow(cosF, 2))));
        dDecdZ = sign * (-1*s38.x*(s38.z*dRadZ- s38.z * cosF / lenS38) + lenS38 * cosF) / (pow(lenS38, 2) * pow(cosF, 2) * sqrt(1 - pow(s38.x, 2) / (pow(lenS38, 2) * pow(cosF, 2))));


        double lenS55 = getLen(s55.x, s55.y, s55.z);
        double cosF55 = cos(asin(s55.z / lenS55));
       double dRadX1 = sign*s55.x*s55.z/(sqrt(1-s55.z*s55.z/(lenS55 * lenS55)) * lenS55 * lenS55 * lenS55);
        double  dRadY1 = sign*s55.x*s55.z/(sqrt(1-s55.z*s55.z/(lenS55 * lenS55)) * lenS55 * lenS55 * lenS55);
        double dRadZ1 = (pow(s55.z, 3)-pow(s55.z, 2))/((sqrt(1-pow(s55.z, 2)/pow(lenS38, 2)) * pow(lenS38, 3)));
//
        if(s55.y / lenS55 > 0){
            sign = 1;
        }

        double dDecdX1 = sign * (-1*s55.x*(s55.z*dRadX- s55.x * cosF / lenS55) + lenS55 * cosF) / (pow(lenS55, 2) * pow(cosF, 2) * sqrt(1 - pow(s55.x, 2) / (pow(lenS55, 2) * pow(cosF, 2))));
        double dDecdY1 = sign * (-1*s55.x*(s55.z*dRadY- s55.y * cosF / lenS55) + lenS55 * cosF) / (pow(lenS55, 2) * pow(cosF, 2) * sqrt(1 - pow(s55.x, 2) / (pow(lenS55, 2) * pow(cosF, 2))));
        double  dDecdZ1 = sign * (-1*s55.x*(s55.z*dRadZ- s55.z * cosF / lenS55) + lenS55 * cosF) / (pow(lenS55, 2) * pow(cosF, 2) * sqrt(1 - pow(s55.x, 2) / (pow(lenS55, 2) * pow(cosF, 2))));

        return Matrix({{dRadX, dRadY, dRadZ, dRadX1, dRadY1, dRadZ1, 0, 0, 0, 0, 0, 0},
                       {dDecdX, dDecdY, dDecdZ, dDecdX1, dDecdY1, dDecdZ1, 0, 0, 0, 0, 0, 0}});
    }


    void findBlackHoleMass(StarStateInterpolator *starStateInterpolator, IsohronDerivative isohronDerivative){
        while(true) {
            S2Original.open("../Data/S2Original.dat");
            S38Original.open("../Data/S38Original.dat");
            S55Original.open("../Data/S55Original.dat");
            isohronDerivative1 = isohronDerivative;
            double Sdec = 0, Sra = 0;
            std::string originalValuesInIMoment;
            std::vector<double> splitedValuesOfOriginals;
            std::vector<double> Gi;

            const char delim = '\n';
            double number;

            getline(S38Original, originalValuesInIMoment, delim);
            std::istringstream split(originalValuesInIMoment);

            while (split >> number)splitedValuesOfOriginals.push_back(number);

            Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 38);
            rAll.emplace_back(splitedValuesOfOriginals[0],
                    ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                            std::pair<double, double>(Gi[0], Gi[1])));

            varAll.emplace_back(pow(splitedValuesOfOriginals[3], 2), pow(splitedValuesOfOriginals[4], 2));
            while (!S38Original.eof()) {
                originalValuesInIMoment.clear();
                splitedValuesOfOriginals.clear();
                double numberInWhile;
                getline(S38Original, originalValuesInIMoment, delim);
                std::istringstream splitInWhile(originalValuesInIMoment);

                while (splitInWhile >> numberInWhile)splitedValuesOfOriginals.push_back(numberInWhile);

                Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 38);
                rAll.emplace_back(splitedValuesOfOriginals[0],
                                                                            ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                                                                               std::pair<double, double>(Gi[0],
                                                                                                         Gi[1])));
                varAll.emplace_back(pow(splitedValuesOfOriginals[3], 2), pow(splitedValuesOfOriginals[4], 2));
            }

            getline(S55Original, originalValuesInIMoment, delim);
            std::istringstream split1(originalValuesInIMoment);

            while (split1 >> number)splitedValuesOfOriginals.push_back(number);

            Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 55);
            rAll.emplace_back(splitedValuesOfOriginals[0],
                              ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                                 std::pair<double, double>(Gi[0], Gi[1])));

            varAll.emplace_back(pow(splitedValuesOfOriginals[3], 2), pow(splitedValuesOfOriginals[4], 2));
            while (!S55Original.eof()) {
                originalValuesInIMoment.clear();
                splitedValuesOfOriginals.clear();
                double numberInWhile;
                getline(S55Original, originalValuesInIMoment, delim);
                std::istringstream splitInWhile(originalValuesInIMoment);

                while (splitInWhile >> numberInWhile)splitedValuesOfOriginals.push_back(numberInWhile);

                Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 55);
                rAll.emplace_back(splitedValuesOfOriginals[0],
                                  ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                                     std::pair<double, double>(Gi[0],
                                                               Gi[1])));
                varAll.emplace_back(pow(splitedValuesOfOriginals[3], 2), pow(splitedValuesOfOriginals[4], 2));
            }
            //S2Original.close();
            for (int i = 0; i < rAll.size(); i++) {
                Sra += (pow(rAll[i].second.first, 2)) / (varAll[i].first);
                Sdec += (pow(rAll[i].second.second, 2)) / (varAll[i].second);
                //std::cout << Sra << " " << Sdec << '\n';
            }
            static int count = 0;
            if (/*(Sra-prevSRA) >= prevSRA/100*0.00001 || (Sdec-prevSDEC) >= prevSDEC/100*0.00001 */count !=5) {
               // std::cout << std::endl;
                //std::cout << Sra << " " << Sdec << '\n';
                std::cout<<Beta.data[12][0]<<"\n";
                GaussNewtonAlgorithm();
                prevSRA = Sra; prevSDEC = Sdec;
                count++;
            } else {
                std::cout << Beta.data[12][0];
                break;
            }
            S2Original.close();
            S38Original.close();
            S55Original.close();
            rAll.clear();
            varAll.clear();
        }
    }

    void GaussNewtonAlgorithm(){
        Matrix dGdX = dGdXMatr(S38, S55);
//        std::cout<<"DGDGIQNECJOOJNCQOCQ:\n\n";
//        dGdX.DebugPrint();
//        std::cout<<"\n\n\n";
        Matrix A = Matrix(13, 104);
        Matrix W = Matrix(104, 104);
        initiateW(W, varAll);
        initiateA(A, dGdX);
        //A.DebugPrint();
        Matrix At = Matrix(transpose(A));
       // At.DebugPrint();
//        std::cout<<"\n\n\n";
//        W.DebugPrint();
//        std::cout<<"\n\n\n";
        Matrix tmp = At*W;
        //tmp.DebugPrint();
        //std::cout<<"\n\n\n";

        Matrix AtWA = inversion(At*W*A);
        Matrix AWrB = AtWrBeta(A, W);
        //std::cout<<"\n\n\n";
        //AtWA.DebugPrint();
        //std::cout<<"\n\n\n";
        //AWrB.DebugPrint();
        //std::cout<<"\n\n\nmultTMP\n";
        Matrix multTMP = AtWA * AWrB;
        //multTMP.DebugPrint();
        //At.DebugPrint();
        Matrix newBeta = Matrix(Beta - (multTMP));
        //std::cout<<"\n\n";
        //newBeta.DebugPrint();
//        if(newBeta.data[6][0]<3e06 || newBeta.data[6][0]>4.8e06){
//            newBeta.data[6][0] = Beta.data[6][0];
//        }

        Beta = Matrix(newBeta);
        std::cout<<"New Beta:\n";
        Beta.DebugPrint();
//        std::cout<<"\n\n\n";
        updateAndRestart(Beta);

    }

    void updateAndRestart(Matrix &B){
        SolvingSystem solvingSystem = SolvingSystem();
//        if(B.data[6][0]<0){
//            B.data[6][0]*=-1;
//        }
//
//            B.data[0][0] *= 8107.55245;
//            B.data[1][0] *= 8107.55245;
        solvingSystem.start(B);
    }





    Matrix inversion(Matrix A)
    {
        double temp;

        Matrix E = Matrix(A.GetCols(), A.GetRows());



        for (int i = 0; i < A.GetCols(); i++)
            for (int j = 0; j < A.GetRows(); j++)
            {
                E.data[i][j] = 0.0;

                if (i == j)
                    E.data[i][j] = 1.0;
            }

        for (int k = 0; k < A.GetCols(); k++)
        {
            temp = A.data[k][k];

            for (int j = 0; j < A.GetRows(); j++)
            {
                A.data[k][j] /= temp;
                E.data[k][j] /= temp;
            }

            for (int i = k + 1; i < A.GetRows(); i++)
            {
                temp = A.data[i][k];

                for (int j = 0; j < A.GetRows(); j++)
                {
                    A.data[i][j] -= A.data[k][j] * temp;
                    E.data[i][j] -= E.data[k][j] * temp;
                }
            }
        }

        for (int k = A.GetRows() - 1; k > 0; k--)
        {
            for (int i = k - 1; i >= 0; i--)
            {
                temp = A.data[i][k];

                for (int j = 0; j < A.GetCols(); j++)
                {
                    A.data[i][j] -= A.data[k][j] * temp;
                    E.data[i][j] -= E.data[k][j] * temp;
                }
            }
        }

        for (int i = 0; i < A.GetRows(); i++)
            for (int j = 0; j < A.GetCols(); j++)
                A.data[i][j] = E.data[i][j];

        return A;
    }


    void initiateA(Matrix& A, Matrix& dGdX){
        int iter = 0;
        for(int i = 0; i < 104; i+=2){
            Matrix tmp = dGdX*isohronDerivative1.interpolate(rAll[iter].first);
            //tmp.DebugPrint();
            for(int j = 0; j < 13; j++){
                A.data[i][j] = tmp.data[0][j];
                A.data[i+1][j] = tmp.data[1][j];
            }

            iter++;
        }

    }


    Matrix transpose(Matrix m){
        Matrix result(m.GetRows(), m.GetCols());
        for (int i = 0; i < m.GetRows(); i++){
            for (int j = 0; j < m.GetCols(); j++){
                result.data[j][i] = m.data[i][j];
            }
        }

        return result;
    }

     Matrix AtWrBeta(Matrix A, Matrix W){
        Matrix rBeta = Matrix(1, rAll.size()*2);
        for(int i = 0, j = 0; i < rBeta.GetCols(); i+=2){
            rBeta.data[i][0] = rAll[j].second.first;
            rBeta.data[i+1][0] = rAll[j].second.second;
            j+=1;
        }
        return transpose(std::move(A))*W*rBeta;
    }


    void initiateW(Matrix &W, std::vector<std::pair<double, double>>vAll){
//        for(int i = 0; i < vAll.size(); i++){
//            std::cout<<"RA: "<<vAll[i].first<<" Dec: "<<vAll[i].second<<"\n";
//        }
        for(int i = 0; i < 104; i+=2){
            for(int j = 0; j < 104; j++){
                W.data[i][j] = 0.0f;
                W.data[i+1][j] = 0.0f;
            }
            W.data[i][i] = 1.0/(vAll[i/2].first);
            W.data[i+1][i+1] = 1.0/(vAll[i/2].second);
        }

    }

    std::pair<double, double> ri(std::pair<double, double>y, std::pair<double, double>g){

        double Ra = y.first-g.first;
        while ((Ra>M_PI) || (Ra<-M_PI)){
            int sign = Ra>M_PI ? -1: 1;
            Ra += sign*2*M_PI;
        }
        return std::pair<double, double>(Ra, y.second-g.second);
    }

    void test(){
        Matrix A = Matrix({{2, 0, 3, 0},{0, 1, 0, 0}, {0, 9, 1, 0}});
        Matrix B = Matrix(A);
        B = transpose(B);
        B.DebugPrint();
        std::cout<<"\n";
        Matrix C = Matrix({{2, 0, 3}, {0, 1, 0}, {0, 9, 1}});
        C = inversion(C);
        C.DebugPrint();
    }

};
#endif //FIRSTVERSION_GAUSSNEWTON_H
