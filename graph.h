#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/convex_hull_2.h> 
#include "GEOM\GeomPerfectMatching.h"
#include <vector>
#include <cstdio>
#include <queue>
#include <set>
#include <map>
#include <assert.h>
#include<iostream>
#include<fstream>

typedef boost::minstd_rand													base_generator_type;
typedef CGAL::Exact_predicates_inexact_constructions_kernel					Kernel;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned int, Kernel>	Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>							Tds;
typedef CGAL::Delaunay_triangulation_2<Kernel, Tds>							Delaunay;
typedef Kernel::Point_2														Point;
typedef std::pair<int, int>													Edge;
typedef std::priority_queue<std::pair<double, std::pair<int,int> > >		EdgeHeap;
typedef std::vector<std::pair<int, int> >									EdgeVector;
typedef std::set<std::pair<int, int> >										EdgeSet;
typedef std::vector<std::pair<double, double> >								VxVector;

class Graph{
public:
	Graph(std::vector<std::pair<double,double> >& vertexList);
	~Graph();
	void Purge();
	void MST();
	void AdjLister();
	bool isReachable(int v1, int v2);
	void dfsTSP();
	void Del();
	void PMatch();
	void OddMatch();
	void EulerCircuit(int start);
	double TSPCircuit(bool refresh);
	void drawGraph(FILE* svg, double MAX);
	void drawArrow(FILE* svg, int p1, int p2);
	int N;		//vertices.size, initialized in the constructor
	std::vector<int>* AdjList;
	std::vector<std::pair<double,double> > vertices;
	std::vector<std::pair<int, int> > edges;
	std::vector<bool> visited;
	std::vector<bool> marked;
	std::vector<int> eCircuit;
	std::vector<int> TSP;
};