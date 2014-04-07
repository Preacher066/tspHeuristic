//this is master
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "graph.h"


int SEED=10, NUM_SENSORS=50;
double MAX=1000.0;
base_generator_type generator(SEED);
boost::uniform_real<> uni_dist(0.0,MAX);
boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

FILE *svg;

extern void drawPoint(double x, double y);
extern void drawLine(double x1, double y1, double x2, double y2);
double dist(std::pair<double, double> p1, std::pair<double,double> p2);

int main() {
	std::vector< std::pair<Point,unsigned> > points;		//Sensor positions
	std::vector<std::pair<double, double> > vertices;
	//std::vector<std::pair<double, double> >& v = vertices;

	for(int i = 0; i < NUM_SENSORS; i++){
		double x = uni();
		double y = uni();
		//points.push_back( std::make_pair( Point(x,y), i ) );
		vertices.push_back(std::make_pair(x,y));
	}

	Graph g(vertices);

	printf("done man\n");

	fopen_s(&svg,"dfsTSP.svg","w");
	//g.MST();
	g.drawGraph(svg,MAX);
	fclose(svg);

	std::cin.get();
	return 0;
}