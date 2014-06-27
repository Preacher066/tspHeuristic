#include "steiner.h"

extern double MAX;
extern FILE* svg;

extern double dist(std::pair<double, double> p1, std::pair<double,double> p2);

double currT;
std::queue<int> boundary;

Steiner::Steiner(){
	std::ifstream vrt("STrees\\3.vrt");
	std::ifstream edg("STrees\\3.edg");
	std::vector<std::pair<double, double> > tvertices;		//temp vertices. After eliminating duplicates from file input, we push into vertices
	std::map<std::pair<double, double>, int> edgeMaker;
	double x=0.0,y=0.0;
	int s=0,i=0;
	while(!vrt.eof()){
		vrt >> x >> y >> s;
		std::pair<double, double> v = std::make_pair(x,y);
		tvertices.push_back(v);
		if(s!=0) steiner.push_back(true);
		else steiner.push_back(false);
		visited.push_back(false);
	}

	//eliminating duplicates (apparently there are!!)
	std::sort(tvertices.begin(),tvertices.end());
	for(int i=tvertices.size()-1;i>=1;i--){
		if(tvertices[i]!=tvertices[i-1]){
			vertices.push_back(tvertices[i]);
		}
	}
	vertices.push_back(tvertices[0]);

	for(int i=0; i<vertices.size();i++){
		edgeMaker[vertices[i]]=i;
	}

	edges = new std::set<int>[vertices.size()];
	double x1=0.0,y1=0.0;
	double x2=0.0,y2=0.0;
	while(!edg.eof()){
		edg >> x1 >> y1 >> s;
		edg >> x2 >> y2 >> s;
		std::pair<double, double> v1 = std::make_pair(x1,y1);
		std::pair<double, double> v2 = std::make_pair(x2,y2);
		int i1 = edgeMaker[v1];
		int i2 = edgeMaker[v2];
		edges[i1].insert(i2);
		edges[i2].insert(i1);
		edgeList.push_back(std::make_pair(i1,i2));
	}
}

Steiner::~Steiner(){
	for(int i=0;i<vertices.size();i++){
		edges[i].clear();
	}
	delete [] edges;
	vertices.clear();
	steiner.clear();
	edgeList.clear();
}

bool Steiner::explorer(){
	if(boundary.empty()) return false;
	int root = boundary.front();
	boundary.pop();
	int curr=0;
	std::set<int> currSet;
	double currT=0.0;
	curr=root;
	visited[root]=true;
	for (std::set<int>::iterator it = edges[curr].begin();it!=edges[curr].end();it++){
		if(!visited[*it])
			heap.push(std::make_pair(-1.0*dist(vertices[curr], vertices[*it]), std::make_pair(*it,curr)));
	}
	currSet.insert(curr);

	while (!heap.empty()){
		double nextWeight = -1.0*heap.top().first;
		currT+=nextWeight;
		if(currT>0.5*T) break;
		curr=heap.top().second.first;
		currSet.insert(curr);
		visited[curr]=true;
		heap.pop();
		for (std::set<int>::iterator it = edges[curr].begin();it!=edges[curr].end();it++){
			if(!visited[*it])
				heap.push(std::make_pair(-1.0*dist(vertices[curr], vertices[*it]), std::make_pair(*it,curr)));
		}
	}

	std::vector<std::pair<double, double> > vert,tsp;
	for(std::set<int>::iterator it = currSet.begin(); it!=currSet.end();it++){
		vert.push_back(vertices[*it]);
	}
	Graph g(vert);
	double cycleWeight=g.TSPCircuit(false);

	while (!heap.empty() && cycleWeight<T){
		curr=heap.top().second.first;
		g.vertices.push_back(vertices[curr]);
		cycleWeight=g.TSPCircuit(true);
		if(cycleWeight>T) break;
		currSet.insert(curr);
		visited[curr]=true;
		heap.pop();
		for (std::set<int>::iterator it = edges[curr].begin();it!=edges[curr].end();it++){
			if(!visited[*it])
				heap.push(std::make_pair(-1.0*dist(vertices[curr], vertices[*it]), std::make_pair(*it,curr)));
		}
	}

	for(std::vector<int>::iterator it = g.TSP.begin();it!=g.TSP.end();it++){
		tsp.push_back(g.vertices[*it]);
	}

	tours.insert(tsp);

	std::set<int> bset;
	while(!heap.empty()){
		std::pair<double, std::pair<int,int> > c = heap.top();
		bset.insert(c.second.second);
		heap.pop();
	}

	for(std::set<int>::iterator it=bset.begin();it!=bset.end();it++){
		boundary.push(*it);
	}

	if(boundary.empty()) return false;
	return true;
}

void Steiner::driver(int root){
	boundary.push(root);
	while(explorer());
	return;
}

void Steiner::drawSteiner(FILE** svg){
	*svg=fopen("new3.svg", "w");
	fprintf(*svg,"<svg width=\"%.3f\" height=\"%.3f\" xmlns=\"http://www.w3.org/2000/svg\">\n", MAX, MAX+50.0);
	fprintf(*svg,"<g>\n");
	fprintf(*svg,"<title>Layer 1</title>\n");
	fprintf(*svg,"<rect height=\"%.3f\" width=\"%.3f\" x=\"%.3f\" y=\"%.3f\" stroke-width=\"0\" stroke=\"#000000\" fill-opacity=\"0.5\" fill=\"#0000\"/>\n",MAX,MAX,0.0,0.0);	// the gray background of the whole field

	//drawing vertices
	int i=0;
	for(std::vector<std::pair<double, double> >::iterator it = vertices.begin();it!=vertices.end();it++,i++){
		if(!steiner[i]) fprintf(*svg,"<circle r=\"3\" cx=\"%.3f\" cy=\"%.3f\" fill=\"#ff0000\"/>\n",it->first,it->second);
		else fprintf(*svg,"<circle r=\"2\" cx=\"%.3f\" cy=\"%.3f\" fill=\"#000000\"/>\n",it->first,it->second);
		//fprintf(*svg, "<text x=\"%.3f\" y=\"%.3f\" style=\"fill: #ff0000; stroke: none; font-size: 11px;\"> %d </text>", it->first+2.5,it->second+2.5, i);
	}
	//drawing edges
	for(std::vector<std::pair<int, int> >::iterator it = edgeList.begin();it!=edgeList.end();it++){
		fprintf(*svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#0000ff\" stroke-width=\"0.5\" stroke-opacity=\"0.3\"/>\n",vertices[it->first].first,vertices[it->first].second,vertices[it->second].first,vertices[it->second].second);
	}
	
	int a=0;
	for(std::set<std::vector<std::pair<double, double> > >::iterator it = tours.begin();it!=tours.end();it++){
		std::vector<Point> in,out;
		std::vector<std::pair<double, double> > in1;
		in1=*it;
		
		//if(a++>7) break;
		if(in1.size()>0){
			//drawing the tour itself
			for(int i=in1.size()-1;i>=1;i--){
				//fprintf(*svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#ffff00\" stroke-width=\"1\"/>\n",in1[i].first, in1[i].second, in1[i-1].first, in1[i-1].second);
				in.push_back(Point(in1[i].first, in1[i].second));
			}
			in.push_back(Point(in1[0].first, in1[0].second));
			CGAL::convex_hull_2(in.begin(), in.end(), std::back_inserter(out));
			out.push_back(out[0]);

			//drawing the convex hull
			for(int i=out.size()-1;i>=1;i--){
				fprintf(*svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#00ff00\" stroke-width=\"2\" stroke-opacity=\"0.2\"/>\n",out[i].x(), out[i].y(), out[i-1].x(), out[i-1].y());
			}
		}
	}

	fprintf(*svg,"</g>\n</svg>\n");

	return;
}