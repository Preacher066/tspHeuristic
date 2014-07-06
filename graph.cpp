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
	AdjList = new std::vector<int>[N];
	for(int i=0;i<N;i++){
		marked.push_back(false);
		visited.push_back(false);
	}
}

Graph::~Graph(){
	edges.clear();
	vertices.clear();
	eCircuit.clear();
	marked.clear();
	visited.clear();
	TSP.clear();
	for(int i=0;i<N;i++){
		AdjList[i].clear();
	}
	delete [] AdjList;
}

/*
Clears the edges and other data from graph to compute TSP anew
*/
void Graph::Purge(){
	edges.clear();
	eCircuit.clear();
	TSP.clear();
	for(int i=0;i<N;i++){
		visited[i] = false;
		marked[i] = false;
	}
}

/*
	Produces the delaunay triangulation of this graph
*/
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
	}
}

/*
	Produces MST of this graph
*/
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
	return;
}

/*
	Produces Adjacency list of this graph in the form of an array of vectors
*/
void Graph::AdjLister(){
	for(EdgeVector::iterator it = edges.begin(); it!=edges.end(); it++){
		AdjList[it->first].push_back(it->second);
		AdjList[it->second].push_back(it->first);
	}
	for(int i=0; i<N; i++){
		if(AdjList[i].size()==1) marked[i]=true;
		std::sort(AdjList[i].begin(), AdjList[i].end());
	}
}

bool Graph::isReachable(int v1, int v2){			//Assumes that the visited array is clear
	visited[v1] = true;
	if(std::binary_search(AdjList[v1].begin(), AdjList[v1].end(),v2)) return true;
	for(std::vector<int>::iterator it = AdjList[v1].begin(); it!=AdjList[v1].end(); it++){
		if(!visited[*it] && isReachable(*it,v2)){
			return true;
		}
	}
	return false;
}

/*
	Produces Euler circuit of this graph
*/
void Graph::EulerCircuit(int v){
	eCircuit.push_back(v);
	if(AdjList[v].size() == 0) return;

	//finding single edge neighbors for v
	int i=0;
	while(i < AdjList[v].size()){
		int v1 = AdjList[v][i];
		if(AdjList[v1].size()==1 && marked[v1]){
			eCircuit.push_back(v1);
			eCircuit.push_back(v);
			AdjList[v1].erase(std::lower_bound(AdjList[v1].begin(),AdjList[v1].end(),v));
			AdjList[v].erase(AdjList[v].begin()+i);
		}
		else i++;
	}
	if(AdjList[v].size() == 0) return;

	//finding all non-bridges
	i=0;
	while(i<AdjList[v].size()){
		for(int j=0;j<N;j++){
			visited[j] = false;
		}
		
		int v2 = AdjList[v][i];
		//remove the edge v--i
		AdjList[v2].erase(std::lower_bound(AdjList[v2].begin(),AdjList[v2].end(),v));
		AdjList[v].erase(AdjList[v].begin()+i);

		//test whether v is still reachable from v2
		if(isReachable(v2,v)){	//if yes then EulerCircuit(v2)
			EulerCircuit(v2);
			return;
		}
		else{					//This edge is a bridge. Restore the edge and try other vertices adjacent to v.
			AdjList[v].push_back(v2);
			std::sort(AdjList[v].begin(), AdjList[v].end());
			AdjList[v2].push_back(v);
			std::sort(AdjList[v2].begin(), AdjList[v2].end());
			i++;
		}
	}

	//by this point, all of the incident edges on v are bridges
	//so pick any one and follow it

	int v3 = *(AdjList[v].begin());

	//but don't forget to erase this edge too
	AdjList[v3].erase(std::lower_bound(AdjList[v3].begin(),AdjList[v3].end(),v));
	AdjList[v].erase(AdjList[v].begin());

	EulerCircuit(v3);
	return;
}

/*
	Perfect matching algorithm from Kolmogorov's Bolssom package
	used in OddMatch function
*/
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

/*
	Used in producing Euler circuit of the graph,
	perfectly matches odd degree vertices
*/
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
	//duplicate edges are included intentionally
	Graph oddGraph(oddDegree);
	EdgeSet edgeSet(edges.begin(), edges.end());
	oddGraph.PMatch();
	//oddGraph.drawGraphT(svg, MAX);
	for(EdgeVector::iterator it = oddGraph.edges.begin(); it != oddGraph.edges.end(); it++){
		Edge e1 = std::make_pair(matchKeeper[it->first], matchKeeper[it->second]);
		/*Edge e2 = std::make_pair(e1.second, e1.first);
		if(edgeSet.find(e1) != edgeSet.end() || edgeSet.find(e2) != edgeSet.end());
		else*/ edges.push_back(e1);
	}
}

/*
	Produces TSP circuit of this graph and returns the resulting cycle's weight
*/
double Graph::TSPCircuit(bool refresh){
	if(refresh){
		delete [] AdjList;
		int N1=N;
		N=vertices.size();
		AdjList = new std::vector<int>[N];
		for(int i=0;i<N-N1;i++){
			marked.push_back(false);
			visited.push_back(false);
		}
	}
	if(N==1) return 0.0;
	if(N==2) {
		TSP.push_back(0);
		TSP.push_back(1);
		TSP.push_back(0);
		return 2*dist(vertices[0],vertices[1]);
	}
	Del();
	MST();
	OddMatch();
	AdjLister();
	EulerCircuit(0);
	double weight=0.0;
	bool* v = new bool[N];
	for(int k=0;k<N;k++){
		v[k]=false;
	}
	for(std::vector<int>::iterator it = eCircuit.begin();it!=eCircuit.end();it++){
		if(!v[*it]){TSP.push_back(*it);v[*it]=true;}
		if(TSP.size() >= 2){
			weight+=dist(vertices[*(TSP.end()-1)], vertices[*(TSP.end()-2)]);
		}
	}
	TSP.push_back(TSP[0]);
	weight+=dist(vertices[*(TSP.end()-1)], vertices[TSP[0]]);
	return weight;
}

void Graph::drawGraph(FILE* svg,double MAX){
	double offset=MAX/2.0;
	fprintf(svg,"<svg width=\"%.3f\" height=\"%.3f\" xmlns=\"http://www.w3.org/2000/svg\">\n", MAX, MAX+50.0);
	fprintf(svg,"<g>\n");
	fprintf(svg,"<title>Layer 1</title>\n");
	fprintf(svg,"<rect height=\"%.3f\" width=\"%.3f\" x=\"%.3f\" y=\"%.3f\" stroke-width=\"0\" stroke=\"#000000\" fill-opacity=\"0.5\" fill=\"#000000\"/>\n",MAX,MAX,0.0,0.0);	// the gray background of the whole field

	//drawing vertices
	int i=0;
	for(std::vector<std::pair<double, double> >::iterator it = vertices.begin();it!=vertices.end();it++,i++){
		fprintf(svg,"<circle r=\"2.5\" cx=\"%.3f\" cy=\"%.3f\" stroke-width=\"5\" stroke=\"#000000\" fill=\"#ff0000\"/>\n",offset+(it->first),offset+(it->second));
		fprintf(svg, "<text x=\"%.3f\" y=\"%.3f\" fill=\"red\"> %d </text>", offset+(it->first)+2.5,offset+(it->second)+2.5, i);
	}
	//drawing edges
	for(std::vector<std::pair<int, int> >::iterator it = edges.begin();it!=edges.end();it++){
		fprintf(svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#ff0000\" stroke-width=\"0.5\" fill-opacity=\"1.0\" fill=\"none\"/>",offset+(vertices[it->first].first),offset+(vertices[it->first].second),offset+(vertices[it->second].first),offset+(vertices[it->second].second));
	}

	for(int i = 1; i<TSP.size(); i++){
		drawArrow(svg, offset+TSP[i-1], offset+TSP[i]);
	}

	for(int j = 0; j<eCircuit.size(); j++){
		fprintf(svg, "<text x=\"%.3f\" y=\"10.0\" fill=\"red\"> _%d_ </text>", 10.0+(j*30.0), eCircuit[j]);
	}

	fprintf(svg,"</g>\n</svg>\n");

	return;
}

void Graph::drawArrow(FILE* svg, int p1, int p2){
	fprintf(svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#00ff00\" stroke-width=\"3\" stroke-opacity=\"0.5\" fill=\"none\"/>",vertices[p1].first,vertices[p1].second,vertices[p2].first,vertices[p2].second);
	return;
}