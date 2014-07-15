#include "steiner.h"

extern double MAX;
extern FILE* svg;
extern double ASP;
double extern MSPEED;


extern double dist(std::pair<double, double> p1, std::pair<double,double> p2);

double currT;
std::queue<int> boundary;

Steiner::Steiner(char* edgName, char* dscName){
	FILE* edg = fopen(edgName,"r");
	FILE* dsc = fopen(dscName,"r");
	std::set<std::pair<double, double> > tvertices;
	std::map<std::pair<double, double>, int> SMaker;
	double x=0.0,y=0.0;
	int s1=0,s2=0,d=0;
	while(fscanf(dsc,"%lf %lf %d",&x,&y,&d)!=EOF){
		coveredCount[std::make_pair(x,y)] = d;
	}
	fclose(dsc);
	while(fscanf(edg,"%lf %lf %d",&x,&y,&s1)!=EOF){
		std::pair<double, double> v = std::make_pair(x,y);
		tvertices.insert(v);
		SMaker[v] = s1;
		visited.push_back(false);
		steiner.push_back(false);
	}

	int i=0;

	for(std::set<std::pair<double, double> >::iterator it = tvertices.begin();it!=tvertices.end();it++,i++){
		vertices.push_back(*it);
		if(SMaker[*it]==1) steiner[i]=true;
		edgeMaker[*it] = i;
	}

	fclose(edg);
	edg = fopen(edgName,"r");

	edges = new std::set<int>[vertices.size()];
	double x1=0.0,y1=0.0;
	double x2=0.0,y2=0.0;
	while(fscanf(edg,"%lf %lf %d",&x1,&y1,&s1)!=EOF){
		fscanf(edg,"%lf %lf %d",&x2,&y2,&s2);
		std::pair<double, double> v1 = std::make_pair(x1,y1);
		std::pair<double, double> v2 = std::make_pair(x2,y2);
		int i1 = edgeMaker[v1];
		int i2 = edgeMaker[v2];
		edges[i1].insert(i2);
		edges[i2].insert(i1);
		edgeList.push_back(std::make_pair(i1,i2));
	}
	fclose(edg);
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
	//std::set<int> currSet;
	std::vector<std::pair<double, double> > tvert,vert,tsp;
	double currT=0.0;
	curr=root;
	visited[root]=true;
	for (std::set<int>::iterator it = edges[curr].begin();it!=edges[curr].end();it++){
		if(!visited[*it]){
			double ed=0.0;
			if(coveredCount.count(vertices[*it]) != 0)
				ed=coveredCount[vertices[*it]];

			heap.push(std::make_pair(-1.0*(2.0*dist(vertices[curr], vertices[*it]) + ed*ASP*MSPEED),
				std::make_pair(*it,curr)));
		}
	}
	vert.push_back(vertices[curr]);

	while (!heap.empty()){
		double nextWeight = -1.0*heap.top().first;
		currT+=nextWeight;
		if(currT>=T) break;
		curr=heap.top().second.first;
		vert.push_back(vertices[curr]);
		visited[curr]=true;
		heap.pop();
		for (std::set<int>::iterator it = edges[curr].begin();it!=edges[curr].end();it++){
			if(!visited[*it]){
				double ed=0.0;
				if(coveredCount.count(vertices[*it]) != 0)
					ed=coveredCount[vertices[*it]];

				heap.push(std::make_pair(-1.0*(2.0*dist(vertices[curr], vertices[*it]) + ed*ASP*MSPEED),
					std::make_pair(*it,curr)));
			}
		}
	}

	cleanTour(vert);
	Graph g(vert);
	double cycleWeight=g.TSPCircuit(false);
	for(std::vector<std::pair<double, double> >::iterator it = vert.begin();it!=vert.end();it++){
		cycleWeight+=coveredCount[*it]*ASP*MSPEED;
	}

	while (!heap.empty() && cycleWeight<T){
		tvert=vert;
		curr=heap.top().second.first;
		tvert.push_back(vertices[curr]);
		cleanTour(tvert);
		Graph g1(tvert);
		cycleWeight=g1.TSPCircuit(false);
		for(std::vector<std::pair<double, double> >::iterator it = tvert.begin();it!=tvert.end();it++){
			cycleWeight+=coveredCount[*it]*ASP*MSPEED;
		}
		if(cycleWeight>T) break;
		vert=tvert;
		visited[curr]=true;
		heap.pop();
		for (std::set<int>::iterator it = edges[curr].begin();it!=edges[curr].end();it++){
			if(!visited[*it]){
				double ed=0.0;
				if(coveredCount.count(vertices[*it]) != 0)
					ed=coveredCount[vertices[*it]];

				heap.push(std::make_pair(-1.0*(2.0*dist(vertices[curr], vertices[*it]) + ed*ASP*MSPEED),
					std::make_pair(*it,curr)));
			}
		}
	}

	Graph g2(vert);
	g2.TSPCircuit(false);

	for(std::vector<int>::iterator it = g2.TSP.begin();it!=g2.TSP.end();it++){
		tsp.push_back(g2.vertices[*it]);
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
	double largestEdge=0.0,ce=0.0;
	for(std::vector<std::pair<int,int> >::iterator it = edgeList.begin(); it!=edgeList.end();it++){
		ce=dist(vertices[it->first],vertices[it->second]);
		if(ce > largestEdge)
			largestEdge = ce;

	}
	if(T<2*largestEdge){
		printf("\n\n**Latency not large enough**\n\n");
		return;
	}
	boundary.push(root);
	while(explorer());
return;
}

void Steiner::cleanTour(std::vector<std::pair<double,double> >& vert){
	int e=0,i=0;
	bool flag=true;
	std::vector<std::pair<double,double> > tvert;
	for(std::vector<std::pair<double,double> >::iterator it = vert.begin();it!=vert.end();it++){
		e=edgeMaker[*it];
		if(steiner[e]){
			flag=true;
			for(std::set<int>::iterator it1=edges[e].begin();it1!=edges[e].end();it1++){
				if(std::find(vert.begin(),vert.end(),vertices[*it1]) == vert.end()) {flag=false; break;}
			}
			if(!flag) tvert.push_back(*it);
		}
		else tvert.push_back(*it);
	}
	vert=tvert;
	return;
}

void Steiner::drawSteiner(FILE* svg){
	double offset=MAX/2;
	fprintf(svg,"<svg width=\"%.3f\" height=\"%.3f\" xmlns=\"http://www.w3.org/2000/svg\">\n", 2*MAX, 2*MAX+50.0);
	fprintf(svg,"<g>\n");
	fprintf(svg,"<title>Layer 1</title>\n");
	fprintf(svg,"<rect height=\"%.3f\" width=\"%.3f\" x=\"%.3f\" y=\"%.3f\" stroke-width=\"0\" stroke=\"#000000\" fill-opacity=\"0.5\" fill=\"#0000\"/>\n",2*MAX,2*MAX,0.0,0.0);	// the gray background of the whole field

	//drawing vertices
	int i=0;
	for(std::vector<std::pair<double, double> >::iterator it = vertices.begin();it!=vertices.end();it++,i++){
		if(steiner[i] == 0) fprintf(svg,"<circle r=\"3\" cx=\"%f\" cy=\"%f\" fill=\"#ff0000\"/>\n",offset+it->first,offset+it->second);
		else fprintf(svg,"<circle r=\"2\" cx=\"%f\" cy=\"%f\" fill=\"#000000\"/>\n",offset+it->first,offset+it->second);
		fprintf(svg, "<text x=\"%f\" y=\"%f\" style=\"fill: #ff0000; stroke: none; font-size: 11px;\"> %d </text>", offset+it->first+2.5,offset+it->second+2.5, i);
	}
		//drawing edges
	for(std::vector<std::pair<int, int> >::iterator it = edgeList.begin();it!=edgeList.end();it++){
		fprintf(svg,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"#0000ff\" stroke-width=\"1\" stroke-opacity=\"1\"/>\n",offset+(vertices[it->first].first),offset+(vertices[it->first].second),offset+(vertices[it->second].first),offset+(vertices[it->second].second));
	}

	for(std::set<std::vector<std::pair<double, double> > >::iterator it = tours.begin();it!=tours.end();it++){
		std::vector<Point> in,out;
		std::vector<std::pair<double, double> > in1;
		in1=*it;

		if(in1.size()>0){
			//drawing the tour itself
			for(int i=in1.size()-1;i>=1;i--){
				fprintf(svg,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"#ffff00\" stroke-width=\"1\"/>\n",offset+(in1[i].first), offset+(in1[i].second), offset+(in1[i-1].first), offset+(in1[i-1].second));
				in.push_back(Point(in1[i].first, in1[i].second));
			}
			in.push_back(Point(in1[0].first, in1[0].second));
			CGAL::convex_hull_2(in.begin(), in.end(), std::back_inserter(out));
			out.push_back(out[0]);

			
		}
	}
	fprintf(svg, "<text x=\"500.0\" y=\"500.0\" style=\"fill: #ff0000; stroke: none; font-size: 40px;\"> %d </text>",tours.size());
	fprintf(svg,"</g>\n</svg>\n");

	return;
}