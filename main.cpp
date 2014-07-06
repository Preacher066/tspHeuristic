//this is master
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "steiner.h"


int SEED=15, NUM_SENSORS=20;
double MAX=1000.0;
base_generator_type generator(SEED);
boost::uniform_real<> uni_dist(0.0,MAX);
boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

FILE *svg;

extern void drawPoint(double x, double y);
extern void drawLine(double x1, double y1, double x2, double y2);
double dist(std::pair<double, double> p1, std::pair<double,double> p2);

int main() {
	/*fopen_s(&svg,"dfsTSP.svg","w");
	std::vector<std::pair<double, double> > vertices;
	//this is rsteiner
	for(int i = 0; i < NUM_SENSORS; i++){
	double x = uni();
	double y = uni();
	vertices.push_back(std::make_pair(x,y));
	}
	*/

	//std::vector<std::pair<double, double> > vertices;
	//FILE* vrt = fopen("Steiners\\50\\greedy0.txt","r");
	//svg = fopen("gy0.svg","w");
	//double x=0.0,y=0.0;
	//while(fscanf(vrt,"%lf %lf",&x,&y)!=EOF){
	//	std::pair<double, double> v = std::make_pair(x,y);
	//	vertices.push_back(v);
	//}
	//Graph g(vertices);
	//g.TSPCircuit(false);
	//g.drawGraph(svg, MAX);
	//fclose(svg);

	//assert(_access("STrees\\1.vrt", 06)==0 && _access("STrees\\1.edg", 06)==0);
	//assert(edg && vrt);
	//for(int i=50;i<=100;i+=10){
		//for(int j=0;j<=5;j++){
			char fname[100];
			sprintf(fname,"Steiners\\%d\\greedy%d.edg",50,2);
			Steiner s(fname);
			s.T = 1000.0;
			s.driver(0);
			char fname2[100];
			sprintf(fname2,"Steiners\\50\\mule%d.svg",50,2);
			svg = fopen(fname2,"w");
			s.drawSteiner(svg);
			fclose(svg);
		//}
	//}

	//std::cin.get();
	return 0;
}