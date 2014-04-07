#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include "GEOM\GeomPerfectMatching.h"
#include <vector>
#include <cstdio>
#include <queue>
#include <set>
#include <map>
#include <assert.h>

typedef boost::minstd_rand													base_generator_type;
typedef CGAL::Exact_predicates_inexact_constructions_kernel					Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>	Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>							Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>							Delaunay;
typedef Kernel::Point_2														Point;
typedef std::priority_queue<std::pair<double, std::pair<int,int> > >		EdgeHeap;

class Graph{
public:
	Graph(std::vector<std::pair<double,double> >& vertexList);
	~Graph();
	void MST();
	void dfsTSP();
	void Del();
	void PMatch();
	void drawGraph(FILE* svg, double MAX);
	int N;
	bool** Adj;
	std::vector<std::pair<double,double> > vertices;
	std::vector<std::pair<int, int> > edges;
};