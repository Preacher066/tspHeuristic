#include "graph.h"

extern void drawPoint(double x, double y);
extern void drawLine(double x1, double y1, double x2, double y2);

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
	EdgeHeap h;
	heap=h;		//erasing heap
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
	for(std::vector<std::pair<int, int> >::iterator it; it!=edges.end(); it++){
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
	return;
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