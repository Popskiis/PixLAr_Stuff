/////////////////////////////////////////////////////////////////////
//////////////////CLUSTERING ALGORITHM BEGINS////////////////////////
/////////////////////////////////////////////////////////////////////

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "TFile.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TGraph.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TRandom.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TH3F.h"
#include "TProfile2D.h"
#include "TChain.h"
#include "TStyle.h"
#include "TString.h"
#include "TVector3.h"
#include <list>
#include <fstream>
#include <sstream>
#include <string>

using namespace std;

#define PI 3.14159
#define euler 2.71828
#define cal_factor 0.0000385

// make a delta ray class and store all of the points here. Also make a function which will find the closest point in the muon track and 
//		add it to the vector???

struct Point {
	float x;
	float y;
	float z;
	float q;
};

struct PCAResults {
	TVector3 centroid;
	pair<TVector3,TVector3> endPoints;
	float length;
	TVector3 eVals;
	vector<TVector3> eVecs;
};

struct TrkPoint{
	int c;
	int cid;
	int pass;
	double x;
	double y;
	double z;
	double q;
};


struct Cluster_pt{
	double rn;
	double ev;
	double cl;
	double x;
	double y;
	double z;
	double q;
	double px;
	double py;
	double pz;
	double ord;
	double closeness;
	double vertex;
	double sp_vertex;
	double true_vertex;
	double is_ord;
};


typedef vector<TrkPoint> track_def;
typedef vector<track_def> superTrack;
typedef vector<Point> PointCloud;
typedef vector<Cluster_pt> Cluster;
void LoadPointCloud(PointCloud &points, const track_def &ord_trk);
PCAResults DoPCA(const PointCloud &points);
double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2);
vector<double> Unit_Vec(double x1,double y1,double z1);
void sameCluster( TrkPoint &p1, TrkPoint &p2 ); 
vector<double> Unit_Vec_NO(double x1,double x2,double y1,double y2,double z1,double z2);
double dotProdFunc(double x1,double x2,double y1,double y2,double z1,double z2);
vector<double> CrossProd(double x1,double x2,double y1,double y2,double z1,double z2);
void setCID( track_def &eliminate, int newCID );


struct by_y { 
	bool operator()(TrkPoint const &a, TrkPoint const &b) { 
		if(a.y == b.y) return a.x > b.x;
		else return a.y > b.y;
	}
};

class newDelta {
	public: 
		superTrack finalDeltas;
		int numPoints = 0;
		
		void storeDelta( track_def delta ) {
			track_def tempStorage;
			for (int j = 0; j < delta.size(); j++)
			{
				//~ cout << "CID: " << delta[j].cid << endl;
				if ( delta[j].cid == 5 )
				{
					++numPoints;
					tempStorage.push_back( delta[j] );
				}
			}
			
			if ( tempStorage.size() != 0)
			{
				finalDeltas.push_back( tempStorage );
			}
			
		};
		
		
		void LOG() {
			cout << "\n==========\nBEGIN LOG\n==========\n\n"
			<< "Final number of passing delta rays: " << finalDeltas.size() << endl
			<< "Total number of points: " << numPoints << endl << endl;
			
				
			for (int iRay = 0; iRay < finalDeltas.size(); iRay++)
			{
				cout << "Ray " << iRay+1 << ": \n"
				<< "  " << finalDeltas[iRay].size() << " points "<< endl;
			}
			
			cout << "\n==========\nEND LOG\n==========\n\n";

		};
		
};



int main()
{
//Parameters
//

double alpha = 0.3;
double ext_alpha = 2.;
double max_phi = acos(0.97);
double back_ang = -0.4; // Angle of back cone 113 deg -> 1.98 Rad
double eps_ang = 0.314159; // epsilon_phi paramter
double eps_dist= 1.5; // epsilon_r parameter
double weight;
unsigned ang_points = 8; // Points in window usied for ordering algorithm

TFile* f = TFile::Open("test.root");	
	if (f == 0) {
		std::cout << "Could not open the requested file.\n";
		return 0;
	}


TFile* g = TFile::Open("deltaRays.root", "RECREATE");
	if(g == 0) {
		std::cout << "Could not open the requested file.\n";
		return 0;
	}

	double tempX[100000], tempY[100000], tempZ[100000], tempID[100000], tempQ[100000], isDelta[100000];		int tempEntry, nHits;

	TTree newTree("newTree", "newTree");

	newTree.Branch("nHits", &nHits, "nHits/I");
	newTree.Branch("sp_x", tempX, "sp_x[nHits]/D");
	newTree.Branch("sp_y", tempY, "sp_y[nHits]/D");
	newTree.Branch("sp_z", tempZ, "sp_z[nHits]/D");
	newTree.Branch("cid", tempID, "cid[nHits]/D");
	//newTree.Branch("sp_charge", tempQ, "sp_charge[nHits]/D");
	newTree.Branch("entry", &tempEntry, "entry/I");
	newTree.Branch("isDelta", isDelta, "isDelta[nHits]/D");

	TTreeReader myReader("treeTest", f);
	TTreeReaderValue<int> getHits(myReader, "nHits");
	TTreeReaderArray<double> getX(myReader, "sp_x");
	TTreeReaderArray<double> getY(myReader, "sp_y");
	TTreeReaderArray<double> getZ(myReader, "sp_z");
// 	TTreeReaderArray<double> getQ(myReader, "sp_charge");
	TTreeReaderValue<int> grabEntry(myReader, "entry");

// myReader.SetEntry( 8 );
while( myReader.Next() ) {
	if( *grabEntry != 315 ) continue;

	cout << "Entry num: " << *grabEntry << endl;
	cout << "Total points: " << getX.GetSize() << endl;
track_def trk;


// store all values from our entry into object trk
for (int iPoint = 0; iPoint < getX.GetSize(); iPoint++) {
	TrkPoint tempPoint;
	tempPoint.c = *grabEntry;
	tempPoint.x = getX[iPoint];
	tempPoint.y = getY[iPoint];
	tempPoint.q = 0;
	tempPoint.z = getZ[iPoint];
	trk.push_back(tempPoint);
}


//Sort track in descending y value
std::sort(trk.begin(), trk.end(), by_y());
int trk_size = trk.size();
track_def points_left; // points yet to be clustered
track_def points_gd; // points left out of cluster
track_def ord_trk; // points clustered
ord_trk.push_back(trk[0]); //Highest y value point is the first point in the oredered track

//Store points being tested by clustering algorithm
for (int i = 1; i < trk_size; ++i){
	TrkPoint tempPoint;
	tempPoint.c = trk[i].c;
	tempPoint.x = trk[i].x;
	tempPoint.y = trk[i].y;
	tempPoint.z = trk[i].z;
	tempPoint.q = trk[i].q;
	points_left.push_back(tempPoint);
}
int pl_size = points_left.size();
std::sort(points_left.begin(), points_left.end(), by_y());


//Clustering algorithm parameters
double old_dist = 10000000.;
int hi_weight_at = -1;
double dist;
double low_ord_y = 10000000.;
double old_weight  = -1.;
double vertex_res,sp_vertex_res;
double num_close_unclustered = 0.;
double bottom_dist;
double shortest_dist = 1000000;


//Variable needed by algorithm
track_def ang_chunk;
std::vector<double> vec_pca_ang;
std::vector<double> vec_ob_temp;
double pca_dotP;
std::vector<double> pca_vec;
std::vector<double> cone_vec;
double cone_dotP, dotP_low_dist;
double flip;
TrkPoint VertexPoint;
TrkPoint sp_VertexPoint;
double far_ucp; 
double far_ucp_ratio; 


//Boolean parameters
bool cone_test_fail = true;
bool ongoing_cone_test = true;
bool closest_point_found = true;
bool closest_point_found_cone = false;
bool closest_point_found_volume = false;
bool closest_point_start_found = false;
bool print_vals = false;
bool closest_point_start_backcone = false;
bool close_ucp_found = false;
bool closest_point_found_short = false;
bool true_vertex_found = false;


//Start of clustering
while(pl_size != 0){
	ang_chunk.clear();
	vec_pca_ang.clear();
	hi_weight_at = -1;
	old_weight = -1.;


	if(ord_trk.size() > ang_points){// PCA will be used when the first ang_points are clustered
		for(unsigned p = ord_trk.size() - ang_points; p < ord_trk.size(); ++p){ // load points for PCA
			TrkPoint ang_tempPoint;
			ang_tempPoint.c = ord_trk[p].c;
			ang_tempPoint.x = ord_trk[p].x;
			ang_tempPoint.y = ord_trk[p].y;
			ang_tempPoint.z = ord_trk[p].z;
			ang_tempPoint.q = ord_trk[p].q;
			ang_chunk.push_back(ang_tempPoint);	
		}


		PointCloud ang_pointcloud;
		PCAResults ang_results;
		LoadPointCloud(ang_pointcloud, ang_chunk); 
		ang_results = DoPCA(ang_pointcloud); // Do PCA
	
		// Load PC eigenvector (unit vector)
		vec_pca_ang.push_back(ang_results.eVecs[0](0));
		vec_pca_ang.push_back(ang_results.eVecs[0](1));
		vec_pca_ang.push_back(ang_results.eVecs[0](2));


		//Make unit vector from first and last point used in PCA calculations := pca_vec
		pca_vec = Unit_Vec_NO(ang_chunk.at(0).x,ang_chunk.back().x,ang_chunk.at(0).y,ang_chunk.back().y,ang_chunk.at(0).z,ang_chunk.back().z);
		
		//Dot Product of PC eigenvector and pca_vec.
		pca_dotP = dotProdFunc(pca_vec[0],vec_pca_ang[0],pca_vec[1],vec_pca_ang[1],pca_vec[2],vec_pca_ang[2]);


		//Decide the orientation of eigenvector to be in direction of pca_vec
		if (pca_dotP > -pca_dotP){
			flip = 1.;
		}else{
			flip = -1.;
		}
		//flip -> flip = -1. ; no flip -> flip = -1.
		//Flip eigenvector
		vec_pca_ang[0] = flip * vec_pca_ang[0];
		vec_pca_ang[1] = flip * vec_pca_ang[1];
		vec_pca_ang[2] = flip * vec_pca_ang[2];
	}else{
		pca_vec.push_back(0);
		pca_vec.push_back(0);
		pca_vec.push_back(0);
		vec_pca_ang.push_back(0);
		vec_pca_ang.push_back(0);
		vec_pca_ang.push_back(0);
	}


	//Boolean reset
	cone_test_fail = true;
	ongoing_cone_test = true;
	closest_point_found_cone = false;
	closest_point_found_volume = false;
	closest_point_start_found = false;
	closest_point_start_backcone = false;
	closest_point_found = false;
	closest_point_found_short = false;


	dotP_low_dist = 0;


	//Cluster first ang_points based on distance alone
	if (ord_trk.size() <= ang_points){
		for (int j = 0; j < pl_size; ++j){
			cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
			cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
			cone_dotP = acos(cone_dotP);
			dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
			weight = pow(euler, -dist/eps_dist) * pow(euler, -cone_dotP/eps_ang);
			if(weight > old_weight){
				old_weight = weight;
				hi_weight_at = j;
				closest_point_start_found = true;
				closest_point_found = true;
			}
		}
	}
	//Priority #1: find point in cone of radius
	if (ord_trk.size() > ang_points){	
		for (int j = 0; j < pl_size; ++j){
			cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
			cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
			cone_dotP = acos(cone_dotP);
			dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
			weight = pow(euler, -dist/eps_dist) * pow(euler, -cone_dotP/eps_ang);
			if(cone_dotP < max_phi){
				if(dist < ext_alpha){
					if(weight > old_weight){
						old_weight = weight;
						dotP_low_dist = cone_dotP;
						hi_weight_at = j;
						closest_point_found_cone = true;
						closest_point_found = true;
					}
				}
			}
		}
	}
	//Priority #2: find point in sphere of radius alpha, excluding back cone
	if (closest_point_found == false && ord_trk.size() > ang_points){
		for (int j = 0; j < pl_size; ++j){
			cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
			cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
			cone_dotP = acos(cone_dotP);
			dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
			weight = pow(euler, -dist/eps_dist) * pow(euler, -cone_dotP/eps_ang);
			if(dist < alpha){
				if(cone_dotP < back_ang){	
					if(weight > old_weight){
						old_weight = weight;
						hi_weight_at = j;
						closest_point_found_volume = true;
						closest_point_found = true;
					}
				}
			}
		}
	}
	//find point with highest weight. DOES NOT CLUSTER THIS POINT
	if(closest_point_found == false || ord_trk.size() > ang_points){
		for (int j = 0; j < pl_size; ++j){
			cone_vec = Unit_Vec_NO(ord_trk.back().x,points_left[j].x,ord_trk.back().y,points_left[j].y,ord_trk.back().z,points_left[j].z);
			cone_dotP = dotProdFunc(cone_vec[0],vec_pca_ang[0],cone_vec[1],vec_pca_ang[1],cone_vec[2],vec_pca_ang[2]);
			cone_dotP = acos(cone_dotP);
			dist = Pythagoras(points_left[j].x,ord_trk.back().x,points_left[j].y,ord_trk.back().y,points_left[j].z,ord_trk.back().z);
			weight = pow(euler, -dist/eps_dist) * pow(euler, -cone_dotP/eps_ang);
			if(isnan(weight)) weight = 0.;	
			if (weight > old_weight){
				old_weight = weight;
				old_dist = dist;
				hi_weight_at = j;
			}
		}
	}
	// Store point with highest weight
	TrkPoint tempPoint;
	tempPoint.c = points_left[hi_weight_at].c;
	tempPoint.x = points_left[hi_weight_at].x;
	tempPoint.y = points_left[hi_weight_at].y;
	tempPoint.z = points_left[hi_weight_at].z;
	tempPoint.q = points_left[hi_weight_at].q;
	if (closest_point_found){// if closest point found, store in ord_trk object
		ord_trk.push_back(tempPoint);
		if (tempPoint.y < low_ord_y){
			low_ord_y = tempPoint.y;
		}
		old_dist = 10000000;
		old_weight = -1.;
		points_left.erase(points_left.begin() + hi_weight_at);
		pl_size = points_left.size();
	}else{ // else, store for future reclustering
		points_gd.push_back(tempPoint);
		old_dist = 10000000;
		old_weight = -1.;
		points_left.erase (points_left.begin() + hi_weight_at);
		pl_size = points_left.size();
	}
	if (pl_size == 0) break; // end loop if there are no more points to cluster
}


track_def final_muon = ord_trk;
track_def delta_points;
//////////////////////////////
// RECLUSTERING ALGORITHM
//////////////////////////////

double threshold_for_cylinder = 0.4; //cm

TrkPoint new_trk_point;
TrkPoint unc_trk_point;
double dist_cluster, low_dist_cluster;
double x1_diff, y1_diff, z1_diff;
double x2_diff, y2_diff, z2_diff;
double Rx,Ry,Rz,Rq;
double Ax,Ay,Az,Bx,By,Bz,Px,Py,Pz,C,Vx,Vy,Vz;
double normAB, normAV, normVP;
double dotProd_AVAB;
double Pq;
int closest_line_seg; 
bool closest_pt_to_line_found;
vector<double> CrossP;
double closeness;
int i_low_dist;
// cout << "\n SIZE:  " << points_gd.size() << endl;
for(int i = 0; i < points_gd.size(); ++i){
	// cout << "\n INCREMENT:  " << i << endl;
	Cluster_pt temp_pt;
	low_dist_cluster = 100000.;
	closest_pt_to_line_found = false;
	Px = points_gd.at(i).x;
	Py = points_gd.at(i).y;
	Pz = points_gd.at(i).z;
	// Pq = points_gd.at(i).q;
// cout << "\n SIZE:  " << ord_trk.size() << endl;
	for(int j = 0; j < ord_trk.size() - 1; ++j){
		Ax = ord_trk.at(j).x;
		Bx = ord_trk.at(j+1).x;
		Ay = ord_trk.at(j).y;
		By = ord_trk.at(j+1).y;
		Az = ord_trk.at(j).z;
		Bz = ord_trk.at(j+1).z;
							
		C = (Px-Ax)*(Bx-Ax) + (Py-Ay)*(By-Ay) + (Pz-Az)*(Bz-Az);
		C = C/(pow((Bx-Ax),2.0) + pow((By-Ay),2.0) + pow((Bz-Az),2.0));
		
		Vx = Ax + C*(Bx-Ax);
		Vy = Ay + C*(By-Ay);
		Vz = Az + C*(Bz-Az);

		normAB = Pythagoras((Bx-Ax),0.0,(By-Ay),0.0,(Bz-Az),0.0);
		normAV = Pythagoras((Vx-Ax),0.0,(Vy-Ay),0.0,(Vz-Az),0.0);
		dotProd_AVAB = (Vx-Ax)*(Bx-Ax) + (Vy-Ay)*(By-Ay) + (Vz-Az)*(Bz-Az);
		dotProd_AVAB = dotProd_AVAB/(normAB*normAV);	
		normVP = Pythagoras(Px,Vx,Py,Vy,Pz,Vz);
		
		if(dotProd_AVAB == -1) continue;
		if(normAV > normAB) continue;
		if(normVP < low_dist_cluster){
			low_dist_cluster = normVP;
			i_low_dist = i;
		}	
	}

	if(low_dist_cluster < threshold_for_cylinder){
		final_muon.push_back(points_gd.at(i_low_dist));
	}else{
		delta_points.push_back(points_gd.at(i_low_dist));
	}					
// printf("\n HERE \n");
}



// FINISHED CLUSTERING POINTS
/////////////////////////////////////////////////////////


printf( " Muon First point: ( %f, %f, %f ) 		Muon Last Point: ( %f, %f, %f )\n", ord_trk[0].x, ord_trk[0].y, ord_trk[0].z,  ord_trk[ord_trk.size() - 1].x,  ord_trk[ord_trk.size() - 1].y,  ord_trk[ord_trk.size() - 1].z );

//////////////////////////////////////////////////////////////////////////////////////////////////
// RECLUSTERING DELTA RAY CANDIDATES & DOING CUTS ///////////////////////////////////////////////
// //////////////////////////////////////////////////////////////////////////////////////////////


// walk through all points in the delta candidates and initialize cid
	for (int kk = 0; kk < delta_points.size(); kk++) {
		delta_points[kk].cid = kk;	
	}

// now we can cluster points close to eachother
	double max_distance = 1.0;
	for (int clusterPoint = 0; clusterPoint < delta_points.size(); clusterPoint++) {
		for (int otherPoint = 0; otherPoint < delta_points.size(); otherPoint++) {
			
			// initial clustering thresh was 0.8 initially
			if ( Pythagoras( delta_points[clusterPoint].x,  delta_points[otherPoint].x,   delta_points[clusterPoint].y,   delta_points[otherPoint].y,   delta_points[clusterPoint].z,   delta_points[otherPoint].z ) < max_distance ){
				sameCluster( delta_points[clusterPoint], delta_points[otherPoint] );			
			}
		}
	}
	


			superTrack deltaClusters;

//			loop through all possible cluster IDs			
			for (int cc = 0; cc < delta_points.size(); ++cc)
			{
//				make a tempVec for every cluster ID
				track_def myTempVec;
//				loop through all points in delta_points vector
				for (int iPoint = 0; iPoint < delta_points.size(); ++iPoint)
				{
//					push back all points into my temporary vector
					if ( delta_points[iPoint].cid == cc )
					{
						myTempVec.push_back(delta_points[iPoint]);	
					}
				}
//				if TempVec is not filled with any points, don't store it in deltaClusters
				if( myTempVec.size() != 0 ) {
					deltaClusters.push_back(myTempVec);
				}
			}

	std::cout << "\n\n~~> Number of initial clusters: " << deltaClusters.size() << "\n\n";


superTrack tempCompare = deltaClusters;
int passingDeltaRays = 0;
// do cuts for everything now in the delta ray vector
	for (int i = 0; i < deltaClusters.size(); i++) {
		bool done = false;
		// if the cluster is too small we omit as candidate
		if ( deltaClusters[i].size() <= 3 ) {
			setCID( deltaClusters[i], 10);
		}else{
		// loop all other clusters
		for (int j = 0; j < deltaClusters.size(); j++) {
			if(done) break;
			// if not same cluster
			if ( deltaClusters[i][0].cid != tempCompare[j][0].cid ) {
				for (int thisPoint = 0; thisPoint < deltaClusters[i].size(); thisPoint++) {
					if(done) break;
					// if within 5cm of either end of the track, cut it
					if ( deltaClusters[i][thisPoint].y > (ord_trk[0].y - 5) || deltaClusters[i][thisPoint].y < (ord_trk[ ord_trk.size() - 1 ].y + 5) )  {
						//~ deltaClusters[i][thisPoint].cid = 10;
						//~ continue;
						setCID( deltaClusters[i], 10 );
						done = true;
						break;
					}
					
					for (int otherPoint = 0; otherPoint < deltaClusters[j].size(); otherPoint++) {
						// if any points are within 5cm of another cluster then we cut it
						if (   Pythagoras( deltaClusters[i][thisPoint].x,  deltaClusters[j][otherPoint].x,   deltaClusters[i][thisPoint].y,   deltaClusters[j][otherPoint].y,   deltaClusters[i][thisPoint].z,   deltaClusters[j][otherPoint].z ) <= 5.0 ) {
							//~ deltaClusters[i][thisPoint].cid = 10;
							//~ continue;
							setCID( deltaClusters[i], 10 );
							done = true;
							break;
						}
						// if we pass the last two if statements, then we have a good candidate
						//~ deltaClusters[i][thisPoint].cid = 5;
						cout << "Number of passing delta rays: " << ++passingDeltaRays << "\n";
						setCID( deltaClusters[i], 5 );
						done = true;
						break;
					}
				}
			}				//somehow make all points cid = 10 and leave when I have a cut
		}					
	}
}

cout << "delta points: " << delta_points.size() << endl;
cout << "muon: " << final_muon.size() << endl;


newDelta DRays;
for (int i = 0; i < deltaClusters.size(); i++)
{
	DRays.storeDelta( deltaClusters[i] );
}
DRays.LOG();

//////////////////////////////////////////////////////////////////////////////////////////////////////
// WRITING TTREE AND WHATNOT ///////////////////////////////////////////////////////////////////////
// /////////////////////////////////////////////////////////////////////////////////////////////////
nHits = *getHits;
tempEntry = *grabEntry;

int pointshift = 0;
for (int i = 0; i < deltaClusters.size(); i++) {
	for (int iPoint = 0; iPoint < deltaClusters[i].size(); ++iPoint)
	{
		tempX[pointshift] = deltaClusters[i][iPoint].x;
		tempY[pointshift] = deltaClusters[i][iPoint].y;		
		tempZ[pointshift] = deltaClusters[i][iPoint].z;
		tempID[pointshift] = deltaClusters[i][iPoint].cid;
		isDelta[pointshift] = 1;
		++pointshift;
	}
}


	for (int iPoint = 0; iPoint < final_muon.size(); ++iPoint)
	{
		tempX[iPoint + pointshift] = final_muon[iPoint].x;
		tempY[iPoint + pointshift] = final_muon[iPoint].y;		
		tempZ[iPoint + pointshift] = final_muon[iPoint].z;
		final_muon[iPoint].cid = 0;
		tempID[iPoint + pointshift] = final_muon[iPoint].cid;
		isDelta[iPoint + pointshift] = 0;
	}

}
newTree.Fill();
g->Write();
g->Close();
f->Close();
}



////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////START OF FUNCTIONS//////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////


void sameCluster( TrkPoint &p1, TrkPoint &p2 ) {

	if (p1.cid < p2.cid)
	{
		p2.cid = p1.cid;
	}

	else p1.cid = p2.cid;
}


double Pythagoras(double x1,double x2,double y1,double y2,double z1,double z2){
	double dist;
	dist = sqrt(pow(x2-x1,2.) + pow(y2-y1,2.) + pow(z2-z1,2.));
	return dist;
}


vector<double> Unit_Vec(double x1,double y1,double z1){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,0.0,y1,0.0,z1,0.0);
	v.push_back(x1/norm);
	v.push_back(y1/norm);
	v.push_back(z1/norm);
	return v;
}


vector<double> Unit_Vec_NO(double x1,double x2,double y1,double y2,double z1,double z2){
	std::vector<double> v;
	double norm;
	norm = Pythagoras(x1,x2,y1,y2,z1,z2);
	v.push_back((x2-x1)/norm);
	v.push_back((y2-y1)/norm);
	v.push_back((z2-z1)/norm);
	return v;
}


double dotProdFunc(double x1,double x2,double y1,double y2,double z1,double z2){
	double dotP;
	dotP = x1*x2 + y1*y2 + z1*z2;
	return dotP;
}


vector<double> CrossProd(double x1,double x2,double y1,double y2,double z1,double z2){
	vector<double> v;
	double x,y,z;
	x = y1*z2 - z1*y2;
	y = x1*z2 - x2*z1;
	y = -y;
	z = x1*y2 - x2*y1;	
	v.push_back(x);
	v.push_back(y);
	v.push_back(z);	
	return v;
}


void LoadPointCloud(PointCloud &points, const track_def &ord_trk) {
	for (int i = 0; i < ord_trk.size(); ++i){
		Point tempPoint;
		tempPoint.x = ord_trk.at(i).x;
		tempPoint.y = ord_trk.at(i).y;
		tempPoint.z = ord_trk.at(i).z;
		tempPoint.q = ord_trk.at(i).q;
		points.push_back(tempPoint);


	}
	return;
}
	
PCAResults DoPCA(const PointCloud &points) {
	TVector3 outputCentroid;
	pair<TVector3,TVector3> outputEndPoints;
	float outputLength;
	TVector3 outputEigenValues;
	vector<TVector3> outputEigenVecs;
	float meanPosition[3] = {0., 0., 0.};
	unsigned int nThreeDHits = 0;
	for (unsigned int i = 0; i < points.size(); i++) {
		meanPosition[0] += points[i].x;
		meanPosition[1] += points[i].y;
		meanPosition[2] += points[i].z;
		++nThreeDHits;
	}
	if (nThreeDHits == 0) {
		PCAResults results;
		return results; 
	}
	const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
	meanPosition[0] /= nThreeDHitsAsFloat;
	meanPosition[1] /= nThreeDHitsAsFloat;
	meanPosition[2] /= nThreeDHitsAsFloat;
	outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);
	float xi2 = 0.0;
	float xiyi = 0.0;
	float xizi = 0.0;
	float yi2 = 0.0;
	float yizi = 0.0;
	float zi2 = 0.0;
	float weightSum = 0.0;
	for (unsigned int i = 0; i < points.size(); i++) {
		const float weight(1.);
		const float x((points[i].x - meanPosition[0]) * weight);
		const float y((points[i].y - meanPosition[1]) * weight);
		const float z((points[i].z - meanPosition[2]) * weight);
		xi2	+= x * x;
		xiyi += x * y;
		xizi += x * z;
		yi2	+= y * y;
		yizi += y * z;
		zi2	+= z * z;
		weightSum += weight * weight;
	}


	Eigen::Matrix3f sig;


	sig <<  xi2, xiyi, xizi,
			xiyi, yi2, yizi,
			xizi, yizi, zi2;


	sig *= 1.0 / weightSum;


	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);


	typedef std::pair<float,size_t> EigenValColPair;
	typedef std::vector<EigenValColPair> EigenValColVector;


	EigenValColVector eigenValColVector;
	const auto &resultEigenMat(eigenMat.eigenvalues());
	eigenValColVector.emplace_back(resultEigenMat(0), 0);
	eigenValColVector.emplace_back(resultEigenMat(1), 1);
	eigenValColVector.emplace_back(resultEigenMat(2), 2);


	std::sort(eigenValColVector.begin(), eigenValColVector.end(), [](const EigenValColPair &left, const EigenValColPair &right){return left.first > right.first;} );


	outputEigenValues = TVector3(eigenValColVector.at(0).first, eigenValColVector.at(1).first, eigenValColVector.at(2).first);


	const Eigen::Matrix3f &eigenVecs(eigenMat.eigenvectors());


	for (const EigenValColPair &pair : eigenValColVector) {
		outputEigenVecs.emplace_back(eigenVecs(0, pair.second), eigenVecs(1, pair.second), eigenVecs(2, pair.second));
	}


	PCAResults results;


	Eigen::ParametrizedLine<float,3> priAxis(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)),Eigen::Vector3f(outputEigenVecs[0](0),outputEigenVecs[0](1),outputEigenVecs[0](2)));


	Eigen::Vector3f endPoint1(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));
	Eigen::Vector3f endPoint2(Eigen::Vector3f(outputCentroid(0),outputCentroid(1),outputCentroid(2)));


	Eigen::Vector3f testPoint;
	Eigen::Vector3f projTestPoint;
	float maxDist1 = -1.0;
	float maxDist2 = -1.0;
	float dist;
	float dotP;
	for (unsigned int i = 0; i < points.size(); i++) {
		testPoint = Eigen::Vector3f(points[i].x,points[i].y,points[i].z);
		projTestPoint = priAxis.projection(testPoint);
		dist = sqrt(pow(projTestPoint(0)-outputCentroid(0),2.0)+pow(projTestPoint(1)-outputCentroid(1),2.0)+pow(projTestPoint(2)-outputCentroid(2),2.0));
		dotP = (projTestPoint(0)-outputCentroid(0))*outputEigenVecs[0](0) + (projTestPoint(1)-outputCentroid(1))*outputEigenVecs[0](1) + (projTestPoint(2)-outputCentroid(2))*outputEigenVecs[0](2);




		if ((dotP < 0.0) && (dist > maxDist1)) {
			endPoint1 = projTestPoint;
			maxDist1 = dist;
		}
		else if ((dotP > 0.0) && (dist > maxDist2)) {
			endPoint2 = projTestPoint;
			maxDist2 = dist;
		}
	}
	outputEndPoints.first = TVector3(endPoint1(0),endPoint1(1),endPoint1(2));
	outputEndPoints.second = TVector3(endPoint2(0),endPoint2(1),endPoint2(2));
	outputLength = sqrt(pow(endPoint2(0)-endPoint1(0),2.0)+pow(endPoint2(1)-endPoint1(1),2.0)+pow(endPoint2(2)-endPoint1(2),2.0));
	results.centroid = outputCentroid;
	results.endPoints = outputEndPoints;
	results.length = outputLength;
	results.eVals = outputEigenValues;
	results.eVecs = outputEigenVecs;
	return results;
}

void setCID( track_def &eliminate, int newCID )
{
	for (int i = 0; i < eliminate.size(); i++)
	{
		eliminate[i].cid = newCID;
	}
	
}
