#include <vector>
#include <cstdio>

extern FILE *svg;

void drawPoint(double x, double y){
	fprintf(svg,"<circle r=\"2.4814\" cx=\"%.3f\" cy=\"%.3f\" stroke-width=\"5\" stroke=\"#000000\" fill=\"#ff0000\"/>\n",x,y);
	return;
}

void drawLine(double x1, double y1, double x2, double y2){
	fprintf(svg,"<line x1=\"%.3f\" y1=\"%.3f\" x2=\"%.3f\" y2=\"%.3f\" stroke=\"#ff0000\" stroke-width=\"0.5\" fill-opacity=\"1.0\" fill=\"none\"/>",x1,y1,x2,y2);
	return;
}

void drawPointsPair(std::pair<std::pair<double,double>, std::pair<double,double> > C){
	drawLine(C.first.first, C.first.second, C.second.first, C.second.second);
	fprintf(svg,"<circle r=\"0.3\" cx=\"%.3f\" cy=\"%.3f\" stroke-width=\"1\" stroke=\"#ff0000\" fill=\"#ff0000\"/>\n",C.first.first,C.first.second);
	fprintf(svg,"<circle r=\"0.3\" cx=\"%.3f\" cy=\"%.3f\" stroke-width=\"1\" stroke=\"#ff0000\" fill=\"#ff0000\"/>\n",C.second.first,C.second.second);
	return;
}
