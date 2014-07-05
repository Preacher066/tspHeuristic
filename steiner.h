#include "graph.h"

#define APPROX 1.0

class Steiner{
public:
Steiner();
~Steiner();
void drawSteiner(FILE** svg);
bool explorer();
void driver(int root);

double T;
std::vector<std::pair<double,double> > vertices;
std::vector<bool> visited;
std::set<int>* edges;
std::vector<std::pair<int,int> > edgeList;
std::vector<bool> steiner;
std::set<std::vector<std::pair<double, double> > > tours;
std::priority_queue<std::pair<double, std::pair<int,int> > > heap;
};