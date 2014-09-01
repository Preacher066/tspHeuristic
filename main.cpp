#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include "steiner.h"


int SEED=15, NUM_SENSORS=20;
double MAX=1000.0;
double ASP=5.0;
double MSPEED=10.0;
base_generator_type generator(SEED);
boost::uniform_real<> uni_dist(0.0,MAX);
boost::variate_generator<base_generator_type&, boost::uniform_real<> > uni(generator, uni_dist);

FILE *svg;

extern void drawPoint(double x, double y);
extern void drawLine(double x1, double y1, double x2, double y2);
double dist(std::pair<double, double> p1, std::pair<double,double> p2);

int main() {

	/**************************************for generating res.txt file in every folder******************************************************/
	//for(int i=50;i<=100;i+=10){
	//	char fname2[100];
	//	sprintf(fname2,"Steiners\\%d\\PResults\\res.txt",i);
	//	FILE* res = fopen(fname2,"w");
	//	for(double T=1000.0;T<=18000.0;T+=100.0){
	//		std::vector<int> ms;
	//		ASP = ((T/MSPEED)*4*8)/(13.4*1024);
	//		for(int j=0;j<=5;j++){
	//			char fname[50];
	//			sprintf(fname,"Steiners\\%d\\Steiner trees\\greedy%d.edg",i,j);
	//			char fname3[50];
	//			sprintf(fname3,"Steiners\\%d\\Disc placement\\greedy%d.txt",i,j);
	//			Steiner s(fname,fname3);
	//			s.T = T;
	//			s.driver(0);
	//			ms.push_back(s.tours.size());
	//		
	//			/*char fname2[100];
	//			sprintf(fname2,"Steiners\\%d\\mule%d.svg",i,j);
	//			svg = fopen(fname2,"w");
	//			s.drawSteiner(svg);
	//			fclose(svg);*/
	//		}
	//		printf("%d: %f\n",i,T);
	//		fprintf(res,"%.3f: ",T);
	//		int total=0;
	//		for(std::vector<int>::iterator it1=ms.begin();it1!=ms.end();it1++){
	//			fprintf(res, "%d ",*it1);
	//			total+=*it1;
	//		}
	//		fprintf(res,"\n");
	//		if(total==6) break;
	//	}
	//	fclose(res);
	//}

	/********************************************************************************************************************/

	/***for generating res_avg.txt and res_lat.txt files in every folder using res.txt files in them (generated using above code)********/

	FILE* res;
	FILE* res_avg;
	FILE* res_lat;
	for(int i=50;i<=100;i+=10){
		double record[30];
		int count[30];
		for(int l=0;l<30;l++) {count[l]=0;record[l]=0.0;}
		char fname[100];
		sprintf(fname,"Steiners\\%d\\PResults\\res.txt",i);
		char fname2[100];
		sprintf(fname2,"Steiners\\%d\\PResults\\res_avg.txt",i);
		char fname3[100];
		sprintf(fname3,"Steiners\\%d\\PResults\\res_lat.txt",i);
		res = fopen(fname,"r");
		res_avg = fopen(fname2,"w");
		res_lat = fopen(fname3,"w");
		double lat=0.0,avg=0.0;
		int m[6];
		while(fscanf(res,"%lf: %d %d %d %d %d %d", &lat, m, m+1, m+2, m+3, m+4, m+5)!=EOF){
			avg=m[0]+m[1]+m[2]+m[3]+m[4]+m[5];
			avg=avg/6.0;
			fprintf(res_avg,"%.3f: %.3f\n",lat,avg);
			record[m[0]]+=lat;
			count[m[0]]++;
			record[m[1]]+=lat;
			count[m[1]]++;
			record[m[2]]+=lat;
			count[m[2]]++;
			record[m[3]]+=lat;
			count[m[3]]++;
			record[m[4]]+=lat;
			count[m[4]]++;
			record[m[5]]+=lat;
			count[m[5]]++;
		}
		for(int k=25;k>=1;k--){
			if(count[k])
				fprintf(res_lat,"%d: %.3f\n",k, (double)(record[k]/count[k]));
		}
		fclose(res);
		fclose(res_avg);
		fclose(res_lat);
	}

	/********************************************************************************************************************/

	//std::cin.get();
	return 0;
}