#include <iomanip>
#include <boost/unordered_map.hpp>
//#include <boost/foreach.hpp>
#include <utility>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <map>
#include <stdlib.h>
#include <cmath>
#include <ctime>
#include <sys/time.h>                // for gettimeofday()
//#include "cpptimer.h"

#define infidouble std::numeric_limits<double>::infinity()

typedef std::vector<double> Point;
typedef std::vector<unsigned int> IDs;

unsigned int MAX_C_ITERATION = 5;


struct cluster {
	cluster(): count(0) { }
	std::vector<double> center;
	unsigned int count;
};

struct vertex{
	std::vector<double> point;
	unsigned int best_cluster;
	double best_distance;
};

typedef std::vector<vertex> Vertexes;
typedef std::vector<cluster> Clusters;

// helper function to compute distance between points
double sqr_distance(const std::vector<double>& a,
		const std::vector<double>& b) {
	//ASSERT_EQ(a.size(), b.size());
	double total = 0;
	for (unsigned int i = 0;i < a.size(); ++i) {
		double d = a[i] - b[i];
		total += d * d;
	}
	return sqrt(total);
}


// helper function to add two vectors
std::vector<double>& plus_equal_vector(std::vector<double>& a,
		const std::vector<double>& b) {
	//ASSERT_EQ(a.size(), b.size());
	for (unsigned int i = 0;i < a.size(); ++i) {
		a[i] += b[i];
	}
	return a;
}

// helper function to subs two vectors
std::vector<double>& subs_equal_vector(std::vector<double>& a,
		const std::vector<double>& b) {
	//ASSERT_EQ(a.size(), b.size());
	for (unsigned int i = 0;i < a.size(); ++i) {
		a[i] -= b[i];
	}
	return a;
}

// helper function to scale a vector vectors
std::vector<double>& scale_vector(std::vector<double>& a, double d) {
	for (unsigned int i = 0;i < a.size(); ++i) {
		a[i] *= d;
	}
	return a;
}

/*
// Dump a point
std::ostream& operator << (std::ostream& os, Point& p){

BOOST_FOREACH(Point::value_type d, p){ os << d << " "; }    
return os;
}
 */

/*
// Dump a Set of Query IDs
std::ostream& operator << (std::ostream& os, IDs & sp){
BOOST_FOREACH(IDs::value_type pid, sp){ os << pid << " " ;}
return os;
};
 */

bool  getassignment(std::vector<vertex> &points, std::vector<cluster> &clusters)
{
	bool somechanged = false;
	unsigned int distancecount = 0;
	double di, dbest;
	unsigned int cbest;
	unsigned int num_points = points.size();
	unsigned int num_clusters = clusters.size();
	for (unsigned int i = 0; i < num_points;  ++i) {
		dbest = infidouble ;
		cbest = int(-1);
		for (unsigned int j = 0; j < num_clusters; ++j) {
			di = sqr_distance(points[i].point, clusters[j].center);
			distancecount++;

			if(dbest > di){
				dbest = di;
				cbest = j;
			}
		}
		points[i].best_distance = dbest;
		if(cbest != points[i].best_cluster){
			points[i].best_cluster = cbest;
			somechanged = true;
		}
	}
	//std::cout << "distancecount " << distancecount << std::endl;
	return somechanged;
}

void updatecenter(std::vector<vertex> &points, std::vector<cluster> &clusters)
{
	unsigned int dimensions = points[0].point.size(); 
	unsigned int num_points = points.size();
	unsigned int num_clusters = clusters.size();
	std::vector<double> emptycenter;
	emptycenter.resize(dimensions);
	for(unsigned int i = 0; i < dimensions; i++){
		emptycenter[i] = 0;
	}

	for (unsigned int i = 0; i < num_clusters; ++i) {
		clusters[i].center = emptycenter;
		clusters[i].count = 0;
	}

	for (unsigned int i = 0; i < num_points;  ++i) {
		plus_equal_vector(clusters[points[i].best_cluster].center,points[i].point);
		clusters[points[i].best_cluster].count  +=1; 
	}

	for (unsigned int i = 0; i < num_clusters; ++i) {
		// std::cout << "clusters["<<i<<"].count " <<clusters[i].count <<std::endl;
		double d = clusters[i].count;
		scale_vector(clusters[i].center, 1.0 /d);
	}
}


class VertexesSpace{
	unsigned int num_points;
	unsigned int num_dimensions;
	Vertexes points;
	int option;
	public:
	VertexesSpace(unsigned int num_points, unsigned int num_dimensions, unsigned int option, std::string filename)  : num_points(num_points), num_dimensions(num_dimensions), option(option)
	{ 
		if(option == 0){//randomly generated inputs
			for (unsigned int i=0; i < num_points; i++){
				vertex p;
				for (unsigned int d=0 ; d < num_dimensions; d++)
				{ p.point.push_back( rand() % 10 ); }
				points.push_back(p);
				//std::cout << "pid[" << i << "]= ("<< p.point << ")" <<std::endl;;
			}
		}
		//option1: readfromfile each line represents a vector
		else if(option == 1){
			std::ifstream ifs(filename.c_str());
			std::string s;
			for (unsigned int i=0; i < num_points; i++){
				vertex p;
				std::getline( ifs, s );
				std::istringstream iss(s);
				copy( std::istream_iterator<double>( iss ), std::istream_iterator<double>(),std::back_inserter(p.point));
				points.push_back(p);
				//std::cout << "pid[" << i << "]= ("<< p.point << ")" <<std::endl;
			}
		}
	}
	inline const unsigned int getNumVertexes() const {return num_points;}
	inline const unsigned int getNumDimensions() const {return num_dimensions;}
	inline vertex& getPointatindex(unsigned int pid) { return points[pid];}
	inline Vertexes& getPoints() { return points;}
};


class RepsSpace{
	VertexesSpace&  ps;
	Clusters  ss;
	unsigned int num_clusters;

	public:
	RepsSpace(VertexesSpace & ps, unsigned int num_clusters)
		: ps(ps), num_clusters(num_clusters)
	{
		//step 1: initialize Rep points to the first k vertexes.
		cluster p;
		for (unsigned int i = 0; i < num_clusters; i++){
			p.center = ps.getPointatindex(i).point;
			//std::cout << "center[" << i << "]= ("<< p.center << ")" <<std::endl;
			ss.push_back(p);
		}
	}
	inline const unsigned int getNumClusters() const {return num_clusters;}
	inline cluster& getClusteratindex(unsigned int pid) { return ss[pid];}
	inline Clusters& getClusters() { return ss;}
	inline vertex& getPointatindex(unsigned int pid) { return ps.getPointatindex(pid);}
	inline Vertexes& getPoints() { return ps.getPoints();}
};

void kmeans_standard(RepsSpace rs){
	bool somechanged = true;
	unsigned int iteration_count = 0;
	Vertexes points = rs.getPoints();
	Clusters clusters = rs.getClusters();



	while(somechanged && iteration_count <= 100) {
		//std::cout << "Enter while loop" <<std::endl;
		somechanged = getassignment(points, clusters);
		updatecenter(points, clusters);
		iteration_count++;
	}
	/*
	   std::cout << "itertions " << iteration_count << std::endl;
	   for(unsigned int i = 0; i < points.size(); i++){
	   std::cout << "point[" << i << "]= ("<< points[i].point << ").bestcenter "<< points[i].best_cluster << " Distance: " << points[i].best_distance << std::endl;
	   }
	 * */

	//
	double dist = 0;
	double tmp;
	for (unsigned int i = 0; i < points.size(); i++) {
		tmp = sqr_distance(points[i].point, clusters[points[i].best_cluster].center);
		dist += tmp;
	}
    std::cout << "K, error, iter num\n";
    std::cout << rs.getNumClusters() << ", " << dist << ", " << iteration_count << std::endl; 
}


int main(int argc, char** argv) {
	if(argc == 1)
	{
		std::cout<<"Usage: [argument]"<<std::endl;
		std::cout<<"Number of points "<<std::endl;
		std::cout<<"Dimenstions of the point "<<std::endl;
		std::cout<<"Randomly generately the graph or read from file? 0 or 1"<<std::endl;
		std::cout<<"filename"<<std::endl;
		std::cout<<"Number of clustes: k "<<std::endl;
		exit(1);
	}

	unsigned int num_points =  std::atoi(argv[1]);
	unsigned int num_dimensions =  std::atoi(argv[2]);
	unsigned int option =  std::atoi(argv[3]);
	std::string filename = argv[4];
	unsigned int num_clusters = std::atoi(argv[5]);

	VertexesSpace ps(num_points, num_dimensions, option, filename);
	RepsSpace rs(ps, num_clusters);
	kmeans_standard(rs);
}


