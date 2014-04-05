//this is master
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "graph.h"

typedef boost::minstd_rand													base_generator_type;
typedef CGAL::Exact_predicates_inexact_constructions_kernel					Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>	Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>							Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>							Delaunay;
typedef Kernel::Point_2														Point;

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
		points.push_back( std::make_pair( Point(x,y), i ) );
		vertices.push_back(std::make_pair(x,y));
	}

	Graph g(vertices);

	Delaunay dt;
	dt.insert(points.begin(),points.end());

	for(Delaunay::Finite_edges_iterator it = dt.finite_edges_begin(); it != dt.finite_edges_end(); ++it){
		Delaunay::Edge e=*it;
		int i1= e.first->vertex((e.second+1)%3)->info();
		int i2= e.first->vertex((e.second+2)%3)->info();
		g.edges.push_back(std::make_pair(i1,i2));
		g.heap.push(std::make_pair(dist(vertices[i1],vertices[i2]), std::make_pair(i1,i2)));
		g.Adj[i1][i2] = true;
	}

	printf("done man\n");

	fopen_s(&svg,"dfsTSP.svg","w");
	//g.MST();
	g.drawGraph(svg,MAX);
	fclose(svg);

	std::cin.get();
	return 0;
}

double dist(std::pair<double, double> p1, std::pair<double,double> p2){
	return (sqrt(((p1.first-p2.first)*(p1.first-p2.first))+((p1.second-p2.second)*(p1.second-p2.second))));
}