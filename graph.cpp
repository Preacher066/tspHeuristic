#include "graph.h"

extern void drawPoint(double x, double y);
extern void drawLine(double x1, double y1, double x2, double y2);
extern FILE *svg;
extern double MAX;

int find(int* parents, int i){
if (parents[i] != i)
        parents[i] = find(parents, parents[i]);
    return parents[i];
}

void Union(int* parents, int* ranks, int x, int y){
	int xroot = find(parents, x);
    int yroot = find(parents, y);
    if (ranks[xroot] < ranks[yroot])
        parents[xroot] = yroot;
    else if (ranks[xroot] > ranks[yroot])
        parents[yroot] = xroot;
    else{
		parents[yroot] = xroot;
		ranks[xroot]++;
    }
}

double dist(std::pair<double, double> p1, std::pair<double,double> p2){
	return (sqrt(((p1.first-p2.first)*(p1.first-p2.first))+((p1.second-p2.second)*(p1.second-p2.second))));
}

Graph::Graph(std::vector<std::pair<double,double> >& vertexList){
	vertices = vertexList;
	N=vertices.size();
	Adj = new bool*[N];
	for(int i=0;i<N;i++){
		Adj[i] = new bool[N];
		for(int j=0;j<N;j++){
			Adj[i][j] = false;
		}
	}
}

Graph::~Graph(){
	edges.clear();
	vertices.clear();
	for(int i=0;i<N;i++){
		delete Adj[i];
	}
	delete Adj;
}

void Graph::Del(){
	Delaunay dt;
	std::vector< std::pair<Point,unsigned> > points;		//Sensor positions
	int i=0;
	for(std::vector<std::pair<double, double> >::iterator it = vertices.begin(); it != vertices.end(); it++){
		points.push_back( std::make_pair( Point(it->first,it->second), i++));
	}
	dt.insert(points.begin(),points.end());

	for(Delaunay::Finite_edges_iterator it = dt.finite_edges_begin(); it != dt.finite_edges_end(); ++it){
		Delaunay::Edge e=*it;
		int i1= e.first->vertex((e.second+1)%3)->info();
		int i2= e.first->vertex((e.second+2)%3)->info();
		edges.push_back(std::make_pair(i1,i2));
//		heap.push(std::make_pair(dist(vertices[i1],vertices[i2]), std::make_pair(i1,i2)));
		Adj[i1][i2] = true;
	}
}

void Graph::MST(){
	std::vector<std::pair<int, int> > temp_edges;
	EdgeHeap temp_heap;
	for(std::vector<std::pair<int, int> >::iterator it=edges.begin(); it!=edges.end(); it++){
		temp_heap.push(std::make_pair(dist(vertices[it->first],vertices[it->second]), std::make_pair(it->first,it->second)));
	}
	int *parents = new int[N];
	int *ranks = new int[N];
	for(int i=0;i<N;i++){
		parents[i] = i;
		ranks[i] = 0;
	}
	while(temp_edges.size()!=N-1){
		std::pair<int,int> e = temp_heap.top().second;
		if(find(parents,e.first) != find(parents,e.second)){
			temp_edges.push_back(e);
			Union(parents,ranks,e.first,e.second);
		}
		temp_heap.pop();
	}

	edges = temp_edges;
	//drawGraph(svg, MAX);
	return;
}

void Graph::PMatch(){
	GeomPerfectMatching gmpm(vertices.size(),2);

	//setting options
	gmpm.gpm_options.init_Delaunay =	0;
	gmpm.gpm_options.init_greedy =		0;
	gmpm.gpm_options.init_KNN =			0;
	gmpm.gpm_options.iter_max =			1;

	//sending vertices
	GeomPerfectMatching::REAL* V = new int[2*vertices.size()];
	std::map<GeomPerfectMatching::PointId,int> IndexMap;
	GeomPerfectMatching::PointId* Indices = new GeomPerfectMatching::PointId[vertices.size()];
	bool* marker = new bool[vertices.size()];
	for(unsigned int i=0;i<vertices.size();i++){
		V[2*i] = vertices[i].first;
		V[2*i+1] = vertices[i].second;
		GeomPerfectMatching::PointId d = gmpm.AddPoint(V+(2*i));
		IndexMap[d] = i;
		Indices[i] = d;
		marker[i] = false;
	}

	//sending edges
	for(unsigned int i=0;i<edges.size();i++){
		gmpm.AddInitialEdge(IndexMap[edges[i].first], IndexMap[edges[i].second]);
	}

	gmpm.SolveComplete();

	std::vector<std::pair<int, int> > temp_edges;
	for(unsigned int i=0;i<vertices.size();i++){
		if(marker[i]) continue;
		int j = IndexMap[gmpm.GetMatch(Indices[i])];
		temp_edges.push_back(std::make_pair(i,j));
		marker[j] = marker[i] = true;
	}
	edges = temp_edges;
}

void Graph::OddMatch(){
	
	//listing the odd degree vertices in oddDegree vector
	VxVector oddDegree;
	std::vector<int> matchKeeper;
	int* degParity = new int[vertices.size()];
	for(int i=0;i<vertices.size();i++){
		degParity[i] = 0;
	}
	for(EdgeVector::iterator it=edges.begin();it!=edges.end();it++){
		degParity[it->first] = 1-degParity[it->first];
		degParity[it->second] = 1-degParity[it->second];
	}
	for(int i=0;i<vertices.size();i++){
		if(degParity[i] == 1){		//odd degree vertices; ofcourse here we are assuming no duplicate edges.
			oddDegree.push_back(vertices[i]);
			matchKeeper.push_back(i);
		}
	}

	//doing the perfect matching in the odd degree vertices
	//not adding duplicate edges
	//so if a matching reuslts in two odd vertices joining which already have an edge between them, then that edge is not added
	Graph oddGraph(oddDegree);
	EdgeSet edgeSet(edges.begin(), edges.end());
	oddGraph.PMatch();
	//oddGraph.drawGraphT(svg, MAX);
	for(EdgeVector::iterator it = oddGraph.edges.begin(); it != oddGraph.edges.end(); it++){
		Edge edge = std::make_pair(matchKeeper[it->first], matchKeeper[it->second]);
		if(edgeSet.find(edge) != edgeSet.end());
		else edges.push_back(edge);
	}
}

void Graph::drawGraph(FILE* svg,double MAX){
	fprintf(svg,"<svg width=\"%.3f\" height=\"%.3f\" xmlns=\"http://www.w3.org/2000/svg\">\n", MAX, MAX+50.0);
	fprintf(svg,"<g>\n");
	fprintf(svg,"<title>Layer 1</title>\n");
	fprintf(svg,"<rect height=\"%.3f\" width=\"%.3f\" x=\"%.3f\" y=\"%.3f\" stroke-width=\"0\" stroke=\"#000000\" fill-opacity=\"0.5\" fill=\"#000000\"/>\n",MAX,MAX,0.0,0.0);	// the gray background of the whole field

	//drawing vertices
	for(std::vector<std::pair<double, double> >::iterator it = vertices.begin();it!=vertices.end();it++){
		fprintf(svg,"<circle r=\"2.4814\" cx=\"%.3f\" cy=\"%.3f\" stroke-width=\"5\" stroke=\"#000000\" fill=\"#ff0000\"/>\n",it->first,it->second);
	}
	//drawing edges
	for(std::vector<std::pair<int, int> >::iterator it = edges.begin();it!=edges.end();it++){
		fprintf(svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#ff0000\" stroke-width=\"0.5\" fill-opacity=\"1.0\" fill=\"none\"/>",vertices[it->first].first,vertices[it->first].second,vertices[it->second].first,vertices[it->second].second);
	}

	fprintf(svg,"</g>\n</svg>\n");

	return;
}

void Graph::drawGraphT(FILE* svg,double MAX){
	//drawing edges
	for(std::vector<std::pair<int, int> >::iterator it = edges.begin();it!=edges.end();it++){
		fprintf(svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#00ff00\" stroke-width=\"3\" stroke-opacity=\"0.5\" fill=\"none\"/>",vertices[it->first].first,vertices[it->first].second,vertices[it->second].first,vertices[it->second].second);
	}

	fprintf(svg,"</g>\n</svg>\n");

	return;
}