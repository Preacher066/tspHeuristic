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

void Graph::MST(){
	std::vector<std::pair<int, int> > temp_edges;
	EdgeHeap temp_heap;

	int *parents = new int[N];
	int *ranks = new int[N];
	for(int i=0;i<N;i++){
		parents[i] = i;
		ranks[i] = 0;
	}
	while(temp_edges.size()!=N-1){
		std::pair<int,int> e = heap.top().second;
		if(find(parents,e.first) != find(parents,e.second)){
			temp_edges.push_back(e);
			temp_heap.push(heap.top());
			Union(parents,ranks,e.first,e.second);
		}
		heap.pop();
	}

	edges = temp_edges;
	heap = temp_heap;
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