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
    Matrix Beta = Matrix(1, 7);

    //StarStateInterpolator starStateInterpolator;
public:
    GaussNewton(double m, double x, double y, double z, double dx, double dy, double dz){

        BlackHoleMass = m;
        Beta = Matrix({{x},
                       {y},
                       {z},
                       {dx},
                       {dy},
                       {dz},
                       {BlackHoleMass}});
    }

    double getLen(double x, double y, double z){
        return sqrt(x*x+y*y+z*z);
    }

    Matrix dGdXMatr(double x, double y, double z){
        double sign = -1;
        double dRadX, dRadY, dRadZ, dDecdX, dDecdY, dDecdZ;
        double len = getLen(x, y, z);
        double cosF = cos(asin(z/ len));
        dRadX = sign*x*z/(sqrt(1-z*z/(len*len))*len*len*len);
        dRadY = sign*x*z/(sqrt(1-z*z/(len*len))*len*len*len);
        dRadZ = (pow(z, 3)-pow(z, 2))/((sqrt(1-pow(z, 2)/pow(len, 2))*pow(len, 3)));
//
        if(y/len>0){
            sign = 1;
        }

        dDecdX = sign*(-1*x*(z*dRadX-x*cosF/len)+len*cosF)/(pow(len, 2)*pow(cosF, 2)*sqrt(1-pow(x, 2)/(pow(len, 2)*pow(cosF, 2))));
        dDecdY = sign*(-1*x*(z*dRadY-y*cosF/len)+len*cosF)/(pow(len, 2)*pow(cosF, 2)*sqrt(1-pow(x, 2)/(pow(len, 2)*pow(cosF, 2))));
        dDecdZ = sign*(-1*x*(z*dRadZ-z*cosF/len)+len*cosF)/(pow(len, 2)*pow(cosF, 2)*sqrt(1-pow(x, 2)/(pow(len, 2)*pow(cosF, 2))));

        return Matrix({{dRadX, dRadY, dRadZ, 0, 0, 0},
                       {dDecdX, dDecdY, dDecdZ, 0, 0, 0}});
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

            getline(S2Original, originalValuesInIMoment, delim);
            std::istringstream split(originalValuesInIMoment);

            while (split >> number)splitedValuesOfOriginals.push_back(number);

            Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 2);
            rAll.emplace_back(splitedValuesOfOriginals[0],
                    ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                            std::pair<double, double>(Gi[0], Gi[1])));

            varAll.emplace_back(pow(splitedValuesOfOriginals[3], 2), pow(splitedValuesOfOriginals[4], 2));
            while (!S2Original.eof()) {
                originalValuesInIMoment.clear();
                splitedValuesOfOriginals.clear();
                double numberInWhile;
                getline(S2Original, originalValuesInIMoment, delim);
                std::istringstream splitInWhile(originalValuesInIMoment);

                while (splitInWhile >> numberInWhile)splitedValuesOfOriginals.push_back(numberInWhile);

                Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 2);
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

            if ((Sra-prevSRA) != prevSRA/100*0.1 || (Sdec-prevSDEC) != prevSDEC/100*0.1 ) {
                std::cout << std::endl;
                std::cout << Sra << " " << Sdec << '\n';
                std::cout<<Beta.data[6][0]<<"\n";
                GaussNewtonAlgorithm();
                prevSRA = Sra; prevSDEC = Sdec;
            } else {
                std::cout << Beta.data[6][0];
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
        Matrix dGdX = dGdXMatr(Beta.data[0][0], Beta.data[1][0], Beta.data[2][0]);
//        std::cout<<"DGDGIQNECJOOJNCQOCQ:\n\n";
//        dGdX.DebugPrint();
//        std::cout<<"\n\n\n";
        Matrix A = Matrix(7, 66);
        Matrix W = Matrix(66, 66);
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
        std::cout<<"\n\n";
        //newBeta.DebugPrint();
//        if(newBeta.data[6][0]<3e06 || newBeta.data[6][0]>4.8e06){
//            newBeta.data[6][0] = Beta.data[6][0];
//        }

        Beta = Matrix(newBeta);
//        std::cout<<"New Beta:\n";
//        Beta.DebugPrint();
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
        for(int i = 0; i < 66; i+=2){
            Matrix tmp = dGdX*isohronDerivative1.interpolate(rAll[iter].first);
            //tmp.DebugPrint();
            for(int j = 0; j < 7; j++){
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
        for(int i = 0; i < 66; i+=2){
            for(int j = 0; j < 66; j++){
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
