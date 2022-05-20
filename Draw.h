//
// Created by Никита Ефимов on 20.05.2022.
//

#ifndef FIRSTVERSION_DRAW_H
#define FIRSTVERSION_DRAW_H
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
            stObjects[i]->setPosition(obj[i]->X(), obj[i]->Y());
            stObjects[i]->setFillColor(sf::Color::White);
            sf::Color orbColor(70*(i+1), 255/(i+1), 255/(i+2));
//            sf::Color orbColor(255, 255, 255);
            orbits[i].push_back(new sf::Vertex(sf::Vector2f(obj[i]->X(), obj[i]->Y()), orbColor));
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
            orbits[i].push_back(new sf::Vertex(sf::Vector2f(objects[i]->X(), objects[i]->Y()), orbColor));
            stObjects[i]->setPosition(objects[i]->X(), objects[i]->Y());
            w.draw(*stObjects[i]);
            for(int j = 0; j < orbits[i].size(); j++){
                w.draw(orbits[i][j], 1, sf::Points);
            }
        }
    }
};
#endif //FIRSTVERSION_DRAW_H
