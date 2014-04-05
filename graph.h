#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <vector>
#include <cstdio>
#include <queue>
#include <set>
#include <assert.h>

typedef std::priority_queue<std::pair<double, std::pair<int,int> > > EdgeHeap;

class Graph{
public:
	Graph(std::vector<std::pair<double,double> >& vertexList);
	~Graph();
	void MST();
	void dfsTSP();
	void drawGraph(FILE* svg, double MAX);
	int N;
	bool** Adj;
	std::vector<std::pair<double,double> > vertices;
	std::vector<std::pair<int, int> > edges;
	EdgeHeap heap;
};