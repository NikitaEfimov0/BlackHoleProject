//
// Created by Никита Ефимов on 23.05.2022.
//

#ifndef FIRSTVERSION_GAUSSNEWTON_H
#define FIRSTVERSION_GAUSSNEWTON_H
#include <iostream>
#include <fstream>
#include "StarStateInterpolator.h"
#include "SolvingSystem.h"
#include <vector>
class GaussNewton{
    std::ifstream S2Original;
    std::ifstream S38Original;
    std::ifstream S55Original;
    double BlackHoleMass;
    std::vector<std::pair<double, std::pair<double, double>>>rAll;
    std::vector<std::pair<double, double>>varAll;
    IsohronDerivative isohronDerivative1;
    Matrix Beta = Matrix(1, 7);
    int flag = 0;
    //StarStateInterpolator starStateInterpolator;
public:
    GaussNewton(double m){


        BlackHoleMass = m;
        Beta = Matrix({{120.451454/8107.55245},
                       {-22.675722/8107.55245},
                       {-104.524315},
                       {-0.556251},
                       {-3.6},
                       {0.0},
                       {BlackHoleMass}});
    }

    std::pair<double, double> Var(double Dec, double Ra){

        double ra = Ra - (int)Ra;
        ra = 0.001 * ((int)std::trunc(ra * 1000) % 10);
        if (ra == 0){
            ra = 0.0001;
        } else {
            ra /= 2;
        }

        double dec = Dec - (int)Dec;
        dec = 0.001 * ((int)std::trunc(dec * 1000) % 10);

        if (dec == 0){
            dec = 0.0001;
        } else {
            dec /= 2;
        }

        return std::pair<double, double>(ra, dec);
    }

    bool findBlackHoleMass(StarStateInterpolator *starStateInterpolator, IsohronDerivative isohronDerivative){
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
            rAll.push_back(std::pair<double, std::pair<double, double>>(splitedValuesOfOriginals[0],
                                                                        ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                                                                           std::pair<double, double>(Gi[1], Gi[2]))));

            varAll.push_back(Var(rAll[rAll.size() - 1].second.second, rAll[rAll.size() - 1].second.first));
            while (!S2Original.eof()) {
                originalValuesInIMoment.clear();
                splitedValuesOfOriginals.clear();
                double numberInWhile;
                getline(S2Original, originalValuesInIMoment, delim);
                std::istringstream splitInWhile(originalValuesInIMoment);

                while (splitInWhile >> numberInWhile)splitedValuesOfOriginals.push_back(numberInWhile);

                Gi = starStateInterpolator->interpolation(splitedValuesOfOriginals[0], 2);
                rAll.push_back(std::pair<double, std::pair<double, double>>(splitedValuesOfOriginals[0],
                                                                            ri(std::pair<double, double>(splitedValuesOfOriginals[1],splitedValuesOfOriginals[2]),
                                                                               std::pair<double, double>(Gi[1],
                                                                                                         Gi[2]))));
                varAll.push_back(Var(rAll[rAll.size() - 1].second.second, rAll[rAll.size() - 1].second.first));
            }
            S2Original.close();
            for (int i = 0; i < rAll.size(); i++) {
                Sra += (pow(rAll[i].second.first, 2)) / (Var(rAll[i].second.second, rAll[i].second.first).first);
                Sdec += (pow(rAll[i].second.second, 2)) / (Var(rAll[i].second.second, rAll[i].second.first).second);
            }
            if (Sra != 0 || Sdec != 0) {
                std::cout << std::endl;
                std::cout << Sra << " " << Sdec << '\n';
                GaussNewtonAlgorithm();
            } else {
                std::cout << BlackHoleMass;
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
        Matrix dGdX = Matrix({{-0.0001, 0, 0, 0, 0, 0},
                              {0, -0.0001, 0, 0, 0, 0}});

        Matrix A = Matrix(7, 66);
        Matrix W = Matrix(66, 66);
        initiateW(W, varAll);
        initiateA(A, dGdX);
        Matrix At = Matrix(transpose(A));
        At.DebugPrint();
        std::cout<<"\n\n\n";
        W.DebugPrint();
        std::cout<<"\n\n\n";
        Matrix tmp = At*W;
        tmp.DebugPrint();
        std::cout<<"\n\n\n";

        Matrix AtWA = inversion(At*W*A);
        Matrix AWB = AtWrBeta(A, W);
        Matrix multTMP = AtWA*AWB;
        multTMP.DebugPrint();
        //At.DebugPrint();
        Matrix newBeta = Matrix(Beta - (AtWA*AWB));
        std::cout<<"\n\n\n";
        newBeta.DebugPrint();
        Beta = Matrix(newBeta);
        updateAndRestart(Beta);

    }

    void updateAndRestart(Matrix &B){
        SolvingSystem solvingSystem = SolvingSystem();

            B.data[0][0] *= 8107.55245;
            B.data[1][0] *= 8107.55245;
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
        for(int i = 0; i < 66; i++){
            Matrix tmp = dGdX*isohronDerivative1.interpolate(rAll[iter].first);
            for(int j = 0; j < 7; j++){
                A.data[i][j] = tmp.data[0][j];
                A.data[i+1][j] = tmp.data[1][j];
            }
            i+=2;
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
        for(int i = 0, j = 0; i < rBeta.GetCols(); i++){
            rBeta.data[i][0] = rAll[j].second.first;
            rBeta.data[i+1][0] = rAll[j].second.second;
            i+=2;
            j+=1;
        }
        return transpose(A)*W*rBeta;
    }


    void initiateW(Matrix &W, std::vector<std::pair<double, double>>vAll){
        static int count = 0;
        static int shift = 0;
        for(int i = 0; i < vAll.size(); i++){
            std::cout<<"RA: "<<vAll[i].first<<" Dec: "<<vAll[i].second<<"\n";
        }
        for(int i = 0; i < 66; i+=2){

            for(int j = 0; j < 66; j++){
                W.data[i][j] = 0.0f;
                W.data[i+1][j] = 0.0f;
            }
            W.data[i][i] = double(1.0/(vAll[i/2].first));
            W.data[i+1][i+1] = double(1.0/(vAll[i/2].second));
        }

    }

    std::pair<double, double> ri(std::pair<double, double>y, std::pair<double, double>g){
        return std::pair<double, double>(y.first-g.first, y.second-g.second);
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
