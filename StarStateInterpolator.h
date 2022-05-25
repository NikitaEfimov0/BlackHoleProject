//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_STARSTATEINTERPOLATOR_H
#define FIRSTVERSION_STARSTATEINTERPOLATOR_H
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

    std::vector<std::pair<std::vector<double>, double>> parseFile(double t, int n){
        std::vector<std::pair<std::vector<double>, double>>states;
        std::string tmp;
        std::string tmpNext;
        std::string delims = " ";
        int tn, tn1;
        std::vector<std::string>strings, strings1;
        const char del = '\n';
        double diff = 0.0;
        bool terminate = false;
        int b, e = 0;
        switch(n){
            case 2:
                fromFileS2.open("../Data/S2.dat");

                while(!fromFileS2.eof() && !terminate){
                    getline(fromFileS2, tmp, del);
                    std::istringstream iss(tmp);
                    std::string w;
                    while (iss >> w) strings.push_back(w);

                    std::istringstream stream(strings[0]);
                    int h;
                    stream>>h;
                    tn = h;
                    diff = abs(h-t);
                    if(diff<=2){
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
                        std::istringstream iss1(tmpNext);
                        std::string w1;
                        while (iss1 >> w1) strings.push_back(w1);

                        std::istringstream stream1(strings[0]);
                        int h1;
                        stream1>>h1;
                        tn1 = h1;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            secondPos.push_back(value);
                        }
                        states.push_back(std::pair(firstPos, tn));
                        states.push_back(std::pair(secondPos, tn1));
                        terminate = true;
                    }
                    else{
                        strings.clear();
                    }

                }

                fromFileS2.close();
                break;
            case 38:
                fromFileS38.open("../Data/S38.dat");

                while(!fromFileS38.eof() && !terminate){
                    getline(fromFileS38, tmp, del);
                    std::istringstream iss(tmp);
                    std::string w;
                    while (iss >> w) strings.push_back(w);

                    std::istringstream stream(strings[0]);
                    int h;
                    stream>>h;
                    tn = h;
                    diff = abs(h-t);
                    if(diff<=2){
                        std::vector<double>firstPos;
                        std::vector<double>secondPos;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            firstPos.push_back(value);
                        }
                        strings.clear();
                        getline(fromFileS38, tmpNext, del);
                        std::istringstream iss1(tmpNext);
                        std::string w1;
                        while (iss1 >> w1) strings.push_back(w1);

                        std::istringstream stream1(strings[0]);
                        int h1;
                        stream1>>h1;
                        tn1 = h1;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            secondPos.push_back(value);
                        }
                        states.push_back(std::pair(firstPos, tn));
                        states.push_back(std::pair(secondPos, tn1));
                        terminate = true;
                    }
                    else{
                        strings.clear();
                    }


                }
                fromFileS38.close();
                break;
            case 55:
                fromFileS55.open("../Data/S55.dat");

                while(!fromFileS55.eof() && !terminate){
                    getline(fromFileS55, tmp, del);
                    std::istringstream iss(tmp);
                    std::string w;
                    while (iss >> w) strings.push_back(w);

                    std::istringstream stream(strings[0]);
                    int h;
                    stream>>h;
                    tn = h;
                    diff = abs(h-t);
                    if(diff<=2){
                        std::vector<double>firstPos;
                        std::vector<double>secondPos;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            firstPos.push_back(value);
                        }
                        strings.clear();
                        getline(fromFileS55, tmpNext, del);
                        std::istringstream iss1(tmpNext);
                        std::string w1;
                        while (iss1 >> w1) strings.push_back(w1);

                        std::istringstream stream1(strings[0]);
                        int h1;
                        stream1>>h1;
                        tn1 = h1;
                        for(int i = 1; i < strings.size(); i++){
                            double value;
                            std::istringstream str(strings[i]);
                            str>>value;
                            secondPos.push_back(value);
                        }
                        states.push_back(std::pair(firstPos, tn));
                        states.push_back(std::pair(secondPos, tn1));
                        terminate = true;
                    }
                    else{
                        strings.clear();
                    }

                }
                fromFileS55.close();
                break;
            default:
                break;
        }
        return states;
    }

    void clean(std::ofstream& toF, std::ifstream& fromF, std::string path){
        std::string tmp;
        std::string result;
        const char del = '\n';
        while(!fromF.eof()){
            getline(fromF, tmp, del);
            if(howMuchSpaces(tmp)<6){
                continue;
            }
            else{
                tmp.append("\n");
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



    std::vector<double> interpolation(double t, int StarNum){

        std::vector<std::pair<std::vector<double>, double>> coordinatesS2 = parseFile(t, StarNum);
//        for(int i = 0; i < coordinatesS2.size(); i++){
//            std::cout<<"t"<<i<<" = "<<coordinatesS2[i].second<<" ";
//            for(int j = 0; j < coordinatesS2[i].first.size(); j++){
//                std::cout<<coordinatesS2[i].first[j]<<" ";
//            }
//            std::cout<<std::endl;
//        }
        Matrix Xt0({{coordinatesS2[0].first[0]}, {coordinatesS2[0].first[1]}, {coordinatesS2[0].first[2]}, {coordinatesS2[0].first[3]}, {coordinatesS2[0].first[4]}, {coordinatesS2[0].first[5]}});
        Matrix Xt1({{coordinatesS2[1].first[0]}, {coordinatesS2[1].first[1]}, {coordinatesS2[1].first[2]},{coordinatesS2[1].first[3]}, {coordinatesS2[1].first[4]}, {coordinatesS2[1].first[5]}});

        Matrix X = Xt0*(double)((double)(coordinatesS2[1].second-t)/(coordinatesS2[1].second-coordinatesS2[0].second))+
                   Xt1*((double)(t-coordinatesS2[0].second)/(coordinatesS2[1].second-coordinatesS2[0].second));
        std::vector<double> result;
        for(int i = 0; i < 6; i++){
            result.push_back(X.data[i][0]);
        }

//        std::cout<<"t = "<<t<<" ";
//        for(int i = 0; i < result.size(); i++){
//            std::cout<<result[i]<<" ";
//        }
        return result;
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
        toFileS2.close();
        toFileS38.close();
        toFileS55.close();

    }


//    double norm(double x1, double y1, double z1, double  x2, double y2, double z2){
//        return sqrt(pow((x1-x2), 2)+pow((y1-y2), 2)+pow((z1-z2), 2));
//    }

//    std::pair<double, double> toSpherical(double x, double y, double z){
//
//        double r = norm(x, y, z, 0, 0, 0);
//        double Dec = asin(z / r);
//
//        double Ra;
//
//        if (y/r > 0){
//            Ra = acos((x) / r) / cos(Dec);
//        } else {
//            Ra = 2 * M_PI - acos((x / r) / cos(Dec) );
//        }
//
//        return std::pair<double, double>(Ra, Dec);
//    }



    void addData(StarObject* objects, int h, int n){
        switch (n) {
            case 2:
                toFileS2 << h;
                toFileS2 << " ";
                toFileS2 << objects->X();
                toFileS2 << " ";
                toFileS2 << objects->Y();
                toFileS2 << " ";
                toFileS2 << objects->Z();
                toFileS2<< " ";
                toFileS2 << objects->dX();
                toFileS2<< " ";
                toFileS2 << objects->dY();
                toFileS2<< " ";
                toFileS2 << objects->dZ();
                toFileS2 << "\n";
                break;
            case 38:
                toFileS38 << h;
                toFileS38 << " ";
                toFileS38 << objects->X();
                toFileS38 << " ";
                toFileS38 << objects->Y() ;
                toFileS38 << " ";
                toFileS38 << objects->Z();
                toFileS38 << " ";
                toFileS38 << objects->dX();
                toFileS38<< " ";
                toFileS38 << objects->dY();
                toFileS38<< " ";
                toFileS38 << objects->dZ();
                toFileS38 << "\n";
                break;
            case 55:
                toFileS55 << h;
                toFileS55 << " ";
                toFileS55 << objects->X();
                toFileS55 << " ";
                toFileS55 << objects->Y();
                toFileS55 << " ";
                toFileS55 << objects->Z();
                toFileS55 << " ";
                toFileS55 << objects->dX();
                toFileS55<< " ";
                toFileS55 << objects->dY();
                toFileS55<< " ";
                toFileS55 << objects->dZ();
                toFileS55 << "\n";
                break;
            default:
                return;

        }
    }
    void addS2Data(StarObject* objects, int h, int n){
        switch (n) {
            case 2:
                toFileS2 << h;
                toFileS2 << " ";
                toFileS2 << objects->X()/8107.55245;
                toFileS2 << " ";
                toFileS2 << objects->Y() / 8107.55245;
                toFileS2 << " ";
                toFileS2 << objects->Z() / 8107.55245;
                toFileS2<< " ";
                toFileS2 << objects->dX();
                toFileS2<< " ";
                toFileS2 << objects->dY();
                toFileS2<< " ";
                toFileS2 << objects->dZ();
                toFileS2 << "\n";
                break;
            case 38:
                toFileS38 << h;
                toFileS38 << " ";
                toFileS38 << objects->X()/8107.55245;
                toFileS38 << " ";
                toFileS38 << objects->Y() / 8107.55245;
                toFileS38 << " ";
                toFileS38 << objects->Z() / 8107.55245;
                toFileS38 << " ";
                toFileS38 << objects->dX();
                toFileS38<< " ";
                toFileS38 << objects->dY();
                toFileS38<< " ";
                toFileS38 << objects->dZ();
                toFileS38 << "\n";
                break;
            case 55:
                toFileS55 << h;
                toFileS55 << " ";
                toFileS55 << objects->X()/8107.55245;
                toFileS55 << " ";
                toFileS55 << objects->Y() / 8107.55245;
                toFileS55 << " ";
                toFileS55 << objects->Z() / 8107.55245;
                toFileS55 << " ";
                toFileS38 << objects->dX();
                toFileS38<< " ";
                toFileS38 << objects->dY();
                toFileS38<< " ";
                toFileS38 << objects->dZ();
                toFileS55 << "\n";
                break;
            default:
                return;

        }



    }



};
#endif //FIRSTVERSION_STARSTATEINTERPOLATOR_H
