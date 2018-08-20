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

struct Point
{
	double x;
	double y;
	double z;
	double q;
	double cluster_id;
	double entry;
	double event;
};


struct PCAResults {
  TVector3 centroid;
  pair<TVector3,TVector3> endPoints;
  float length;
  TVector3 eVals;
  vector<TVector3> eVecs;
};



typedef std::vector<Point> pointCluster;
typedef std::vector<pointCluster> familyClusters;
void findBestClusterMoveAllPoints( const Point &moveCent, const familyClusters &otherClusters, const pointCluster &myCluster, double &newNdz, double &newMdy, int &bestCID, int &tbestCID, int &eliminate, int &connector, double &distanceCheck );
void moveCluster( const pointCluster &moveMe, const Point &theCentroids, double MdY, double NdZ);
void mergeFriends( pointCluster &mergeMe, pointCluster &bestFriend, familyClusters &myClusterFamily );
void mergeClusters( familyClusters &myClusterFamily, pointCluster &mergeMe, double threshold );
void mergeToClosest( familyClusters &myClusterFamily, pointCluster &mergeMe );
void pushNewCentroid( pointCluster &newMergerCl );
double changed_sq_eucl_dist( Point p1, Point p2, double NdZ, double MdY );
PCAResults DoPCA(const pointCluster &points);
pointCluster tempMerge( pointCluster moveMe, familyClusters wholeTrack, double NdZ, double MdY );
// pointCluster tempMerge1( pointCluster moveMe, familyClusters wholeTrack, double NdZ, double MdY );
double movedVal( double oldVal, double shifterVal );
double getAvgX( pointCluster cl );
double getAvgY( pointCluster cl );
double getAvgZ( pointCluster cl );
double square( double x );
double getMinSquare( Point moveCent, Point closestCent, double minDistance );
void sameCluster( Point &p1, Point &p2 );
double eucl_dist( Point p1, Point p2 );
void snapshot( familyClusters myClusterFamily, int movenum );
void adjustTrack( familyClusters &myClusterFamily, const pointCluster oldPoints );
familyClusters doReclustering( pointCluster track, bool &allFinished );


struct compare
{
	bool operator()(const pointCluster &first, const pointCluster &second) {
		return first.size() < second.size();
	}	
};


// ===================================================================================================================================
// =========================================== MAIN FUNCTION =========================================================================
// ===================================================================================================================================


int main()
{
	TStopwatch myWatch;
	myWatch.Start();
	
	// change this to choose whether u want pictures of each cluster move or not
	bool wantPictures = false;

//	Open and read TTree from anaTree file
	// original file is tracks.root
	TFile* f = TFile::Open("tracks.root");
		if(f == 0) {
			std::cout << "Could not open the requested file.\n";
			return -1;
		}

//	double values used in for loop later to write TTree
	double tempX[100000], tempY[100000], tempZ[100000], tempQ[100000], tempID[100000];		int tempEvent, tempEntry, nHits;
	double centX, centY, centZ;		int centID;



	TFile* g = TFile::Open("test.root", "RECREATE");
		if(g == 0) {
			std::cout << "Could not open the requested file.\n";
			return -1;
		}

		TTree treeTest("treeTest", "test tree");

		treeTest.Branch("nHits", &nHits, "nHits/I");
		treeTest.Branch("sp_x", tempX, "sp_x[nHits]/D");
		treeTest.Branch("sp_y", tempY, "sp_y[nHits]/D");
		treeTest.Branch("sp_z", tempZ, "sp_z[nHits]/D");
		treeTest.Branch("sp_charge", tempQ, "sp_charge[nHits]/D");
		treeTest.Branch("cluster_id", tempID, "cluster_id[nHits]/D");
		treeTest.Branch("entry", &tempEntry, "entry/I");
		treeTest.Branch("event", &tempEvent, "event/I");



	//	Gain access to data containers with TTreeReader and friends
		
		TTreeReader myReader("tree", f);
		TTreeReaderValue<int> getHits(myReader, "nHits");
		TTreeReaderArray<double> grabYPos(myReader, "sp_y");
		TTreeReaderArray<double> grabXPos(myReader, "sp_x");
		TTreeReaderArray<double> grabZPos(myReader, "sp_z");
		TTreeReaderArray<int> getQ(myReader, "sp_charge");
		TTreeReaderValue<int> grabEntry(myReader, "entry");
		TTreeReaderValue<int> grabEvent(myReader, "event");



while( myReader.Next() ) {
	if(*grabEntry != 315) continue;

		pointCluster allPoints, oldPoints;


//			for all points in the event, store as an element in the allPoints vector. Give each point a unique
//				cluster ID to be adjusted later.
			for(int iPoint = 0; iPoint < grabYPos.GetSize(); ++iPoint) {
				Point tempPoint;


				tempPoint.x = grabXPos[iPoint];
				tempPoint.y = grabYPos[iPoint];		
				tempPoint.z = grabZPos[iPoint];		
				tempPoint.cluster_id = iPoint;
				tempPoint.entry = *grabEntry;
				tempPoint.event = *grabEvent;
				allPoints.push_back(tempPoint);
				oldPoints.push_back(tempPoint);
			}


//================================================ Setting clusterIDs here ===========================================================


//			threshold for how close a point must be to get clustered
// 			originally thresh was 1.0 when testing on 4430016 / 690
// 			then, we had it a 0.8 and it worked better

//			iterate through all points on track
			for (int iPoint = 0; iPoint < grabXPos.GetSize(); ++iPoint) {
//				then for each point compare euclidian distance to all other points
				double threshold = 0.6;
				for (int closestPoint = 0; closestPoint < grabXPos.GetSize(); ++closestPoint)
				{
//					find first point within threshold and group them together in the same cluster				
//						& keep going until we have the closest point
					if ( eucl_dist( allPoints[iPoint], allPoints[closestPoint] ) < threshold )
					{
						sameCluster( allPoints[iPoint], allPoints[closestPoint] );
						threshold = eucl_dist( allPoints[iPoint], allPoints[closestPoint] ); 
					}
				}
			}




//========================================= Making a vector with cluster of points as elements ================================================


			familyClusters myClusterFamily;

//			loop through all possible cluster IDs			
			for (int cid = 0; cid < grabXPos.GetSize(); ++cid)
			{
//				make a tempVec for every cluster ID
				pointCluster myTempVec;
//				loop through all points in allPoints vector
				for (int iPoint = 0; iPoint < grabXPos.GetSize(); ++iPoint)
				{
//					push back all points into my temporary vector
					if ( allPoints[iPoint].cluster_id == cid )
					{
						myTempVec.push_back(allPoints[iPoint]);	
					}
				}
//				if TempVec is not filled with any points, don't store it in myClusterFamily
				if( myTempVec.size() != 0 ) {
					myClusterFamily.push_back(myTempVec);
				}
			}


//======================================= Finding centroids and sorting clusters ================================================


	//	push back centroids into each cluster -> we know the last element will always be the centroid
	//	-> later we remove them by using pop_back to delete last element
	for (int iCluster = 0; iCluster < myClusterFamily.size(); ++iCluster)
	{
		Point centroidTemp;
		centroidTemp.x = getAvgX( myClusterFamily[iCluster] );
		centroidTemp.y = getAvgY( myClusterFamily[iCluster] );
		centroidTemp.z = getAvgZ( myClusterFamily[iCluster] );
		centroidTemp.cluster_id = myClusterFamily[iCluster][0].cluster_id;
		myClusterFamily[iCluster].push_back(centroidTemp);
	}


// 	sort clusters by smallest to largest
	compare c;
	std::sort( myClusterFamily.begin(), myClusterFamily.end(), c );
//	now we can feed the list into the next part -> finding the best way to move each cluster



//======================================= Assign randomID number to clusterID ====================================================


//		make mersenne twister engine and distribution selection	
		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_int_distribution<> dist(0,1000);


		int randID;		// our random number for ID change!
		// loop through all clusters in family
		for (int iCluster = 0; iCluster < myClusterFamily.size(); ++iCluster)
		{
			randID = dist(gen);	// call the mersenne twister engine inside of our distribution and set to randID
		//	loop through all clusterIDs and change them to random
			for (int iPoint = 0; iPoint < myClusterFamily[iCluster].size(); ++iPoint)
			{
				myClusterFamily[iCluster][iPoint].cluster_id = randID;
			}
		}




//======================================== Finding best move for each cluster =====================================================


//	vector of points pointCluster created here to store all of our centroids
	pointCluster theCentroids;
	familyClusters storageVector;


				std::cout << "\n==========================================================================================" <<
				"=============================\n 			STARTING ENTRY: " << *grabEntry << "		\n======================================================"
				<< "===============================================================================\n\n";



	int movenum = 0;
	int moveThis = 0;
	bool done = false;
	bool firstTime = true;
	bool allFinished = false;
	
	
//	we just keep doing the code here until we are done!
	while( !done )
	{
		++movenum;
		std::cout << "\nCurrent number of clusters: " << myClusterFamily.size() << "\n";

		// make a snapshot
		if (wantPictures) {
			snapshot( myClusterFamily, movenum );
		}

//			use these for shifting each cluster associated with moveCent by bestNdZ and bestNdY
		double bestNdz, bestMdy;
		int bestCID; 		int eliminate = 0;
		int tbestCID = -1;
		int connector = 0;
		double distanceCheck = 0;
//			make temp point moveCent = last element of smallest cluster = centroid for that cluster
		Point moveCent = myClusterFamily[moveThis][myClusterFamily[moveThis].size() - 1]; 

//			if none of them want to move -> then we are done
		if ( myClusterFamily.size() == 1 )
		{
			if (firstTime) {
				allFinished = true;
				firstTime = false;
			}
			familyClusters newClusterFamily;
			//recluster again to see if we have more than one cluster
			// newClusterFamily = doReclustering( myClusterFamily[0], allFinished );

			if (allFinished) {
				std::cout << "\n========================================================================================================" <<
				"=============================\n 						CLUSTERING DONE 		\n======================================================"
				<< "===============================================================================\n\n";
				done = true; // we are now done!
				break;
			}else myClusterFamily = newClusterFamily;
		}

//			find the best Mdy and Ndz values for moving any one of the clusters around to minimize distances
		findBestClusterMoveAllPoints( moveCent, myClusterFamily, myClusterFamily[moveThis], bestNdz, bestMdy, bestCID, tbestCID, eliminate, connector, distanceCheck );
		printf( " \nBest CID: %d \n", bestCID );
//				add Mdy and Ndz to all points within new smallest cluster
		for (int iPoint = 0; iPoint < myClusterFamily[0].size(); ++iPoint)
		{
			myClusterFamily[moveThis][iPoint].y += bestMdy;
			myClusterFamily[moveThis][iPoint].z += bestNdz;
		}	
		mergeFriends( myClusterFamily[moveThis], myClusterFamily[bestCID], myClusterFamily );
		std::sort( myClusterFamily.begin(), myClusterFamily.end(), c ); 
	}
	++movenum;

	// remove all centroids now
	for (int killCentroid = 0; killCentroid < myClusterFamily.size(); ++killCentroid)
		{
			myClusterFamily[killCentroid].pop_back();
		}	

	for (int i = 0; i < storageVector.size(); i++) {
		myClusterFamily.push_back(storageVector[i]);	
	}

	//make a snapshot
	if (wantPictures) {
		snapshot( myClusterFamily, movenum );
	}

std::cout << "\n ADJUSTING TRACK...\n";
adjustTrack( myClusterFamily, oldPoints );
std::cout << "\n==================================== DONE ====================================\n\n";



//======================================= Fill our TTree with our new Data! ========================================================
			nHits = *getHits;
			tempEntry = *grabEntry;
			tempEvent = *grabEvent;
			int pointShift = 0;
//			filling TTree branches here using for loop and temp xyz values.	
			for (int iCluster = 0; iCluster < myClusterFamily.size(); ++iCluster)
			{
				for (int iPoint = 0; iPoint < myClusterFamily[iCluster].size(); ++iPoint)
				{
					tempX[pointShift] = myClusterFamily[iCluster][iPoint].x;
					tempY[pointShift] = myClusterFamily[iCluster][iPoint].y;		
					tempZ[pointShift] = myClusterFamily[iCluster][iPoint].z;
					tempQ[pointShift] = myClusterFamily[iCluster][iPoint].q;
					tempID[pointShift] = myClusterFamily[iCluster][iPoint].cluster_id;	
					++pointShift;
				}
			}

			treeTest.Fill();
}

//======================================== output stuff ====================================================================



g->Write();
myWatch.Stop();
myWatch.Print();
g->Close();		
f->Close();

return 0;
}



// ====================================== Functions list ====================================================



familyClusters doReclustering( pointCluster track, bool &allFinished )
{
	bool moreThanPoint = false;
	familyClusters newClusterFamily;
	track.pop_back();

	for (int points = 0; points < track.size(); points++) {
		track[points].cluster_id = points;
	}

//	iterate through all points on track
	for (int iPoint = 0; iPoint < track.size(); ++iPoint) {
//		then for each point compare euclidian distance to all other points
		double threshold = 0.8;
		for (int closestPoint = 0; closestPoint < track.size(); ++closestPoint)
		{
//			find first point within threshold and group them together in the same cluster				
//				& keep going until we have the closest point
			if ( eucl_dist( track[iPoint], track[closestPoint] ) < threshold )
			{
				sameCluster( track[iPoint], track[closestPoint] );
				threshold = eucl_dist( track[iPoint], track[closestPoint] ); 
			}
		}
	}

//	loop through all possible cluster IDs			
	for (int cid = 0; cid < track.size(); ++cid)
	{
//		make a tempVec for every cluster ID
		pointCluster myTempVec;
//		loop through all points in track vector
		for (int iPoint = 0; iPoint < track.size(); ++iPoint)
		{
//			push back all points into my temporary vector
			if ( track[iPoint].cluster_id == cid )
			{
				myTempVec.push_back(track[iPoint]);	
			}
		}
//		if TempVec is not filled with any points, don't store it in myClusterFamily
		if( myTempVec.size() != 0 ) {
			newClusterFamily.push_back(myTempVec);
		}
	}
std::cout << "\n  size:  " << newClusterFamily.size() << std::endl;

	//	push back centroids into each cluster -> we know the last element will always be the centroid
	//	-> later we remove them by using pop_back to delete last element
	for (int iCluster = 0; iCluster < newClusterFamily.size(); ++iCluster)
	{
		if (newClusterFamily[iCluster].size() > 2) moreThanPoint = true;

		Point centroidTemp;
		centroidTemp.x = getAvgX( newClusterFamily[iCluster] );
		centroidTemp.y = getAvgY( newClusterFamily[iCluster] );
		centroidTemp.z = getAvgZ( newClusterFamily[iCluster] );
		centroidTemp.cluster_id = newClusterFamily[iCluster][0].cluster_id;
		newClusterFamily[iCluster].push_back(centroidTemp);
	}


	// 	sort clusters by smallest to largest
	compare comp;
	std::sort( newClusterFamily.begin(), newClusterFamily.end(), comp );
	if (newClusterFamily.size() > 1 && moreThanPoint) allFinished = false;
	return newClusterFamily;
}


void adjustTrack( familyClusters &myClusterFamily, const pointCluster oldPoints )
{
	double bestLinearity = 0;
	double newLinearity;
	PCAResults newResults;
	double dZ = 4.5;
	double dY = 2.4;
	int bestN = 0;
	int bestM = 0;

	for (int N = -2; N < 2; N++) {
		for (int M = -2; M < 2; M++) {
			pointCluster amIGood = oldPoints; 
			for (int pushTrack = 0; pushTrack < myClusterFamily[0].size(); pushTrack++) {
				Point tempPoint;
				tempPoint.x = myClusterFamily[0][pushTrack].x;
				tempPoint.y = movedVal( myClusterFamily[0][pushTrack].y, M*dY );
				tempPoint.z = movedVal( myClusterFamily[0][pushTrack].z, N*dZ );
				amIGood.push_back( tempPoint );
			} 
			newResults = DoPCA( amIGood );
			newLinearity = ( newResults.eVals(0) ) / ( newResults.eVals(0) + newResults.eVals(1) + newResults.eVals(2) );
			if ( newLinearity > bestLinearity ) {
				std::cout << "\n linearity of new track:  " << newLinearity << std::endl; 
				bestLinearity = newLinearity;
				bestN = N;
				bestM = M;
			}
		}
	}
	for (int pnt = 0; pnt < myClusterFamily[0].size(); pnt++) {
		myClusterFamily[0][pnt].y += bestM*dY;
		myClusterFamily[0][pnt].z += bestN*dZ;
	}
}




void snapshot( familyClusters myClusterFamily, int movenum )
{

	// zx image
	TH2D* h2 = new TH2D("h2","h2", 2000, 0, 54, 2000, 6, 90);
	TCanvas *c1 = new TCanvas("c1","c1");
	c1->cd();
		for (int i = 0; i < myClusterFamily.size(); i++) {
			for (int j = 0; j < ( myClusterFamily[i].size() - 1 ); j++) {
				h2->Fill(myClusterFamily[i][j].x, myClusterFamily[i][j].z );	
			}
		}
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	h2->Draw();
	c1->SaveAs(Form("myPics/move%d_zx.png", movenum));
	delete h2;
	delete c1;

	// yx image
	TH2D* h3 = new TH2D("h3","h3", 2000, 0, 54, 2000, -20, -20);
	TCanvas *c2 = new TCanvas("c2","c2");
	c2->cd();
		for (int i = 0; i < myClusterFamily.size(); i++) {
			for (int j = 0; j < ( myClusterFamily[i].size() - 1 ); j++) {
				h3->Fill(myClusterFamily[i][j].x, myClusterFamily[i][j].y );	
			}
		}
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	h3->Draw();
	c2->SaveAs(Form("myPics/move%d_yx.png", movenum));
	delete h3;
	delete c2;        

}


void mergeFriends( pointCluster &mergeMe, pointCluster &bestFriend, familyClusters &myClusterFamily )
{
	
//				first we remove the centroids of both clusters...
				mergeMe.pop_back();	bestFriend.pop_back();
						
//				then we give all points in mergeMe the same cluster ID as the first point in bestfriend 
//					and we move all points in mergeMe into the same cluster vector as best friend
				for (int mergePoint = 0; mergePoint < mergeMe.size(); ++mergePoint)
				{
//					give all points the same CID first	
					mergeMe[mergePoint].cluster_id = bestFriend[0].cluster_id;	
//					printf("\nMergeMe CID: %f 	best Friend CID: %f\n", mergeMe[mergePoint].cluster_id,  bestFriend[0].cluster_id);	
//					now we put that point into place at the end of the cluster which was closest					
					bestFriend.push_back( mergeMe[mergePoint] );
				}

//				now we can find a new centroid and push it back into the new cluster
					pushNewCentroid( bestFriend );
//				then we can remove mergeMe entirely from myClusterFamily because it is unused now
				myClusterFamily.erase( myClusterFamily.begin() );
				return;
}



void mergeClusters( familyClusters &myClusterFamily, pointCluster &mergeMe, double threshold )
{

//	first, we need to get access to all clusters in myClusterFamily
	for (int iCluster = 0; iCluster < myClusterFamily.size(); ++iCluster)
	{
//		now we need access to all of the points within the other clusters
		for (int iPoint = 0; iPoint < myClusterFamily[iCluster].size(); ++iPoint)
		{
//			then, find the first point that is close enough to any point inside of cluster mergeMe and make them the same cluster
//			so long as it is a point not in the same cluster as mergeMe!			
			for (int kk = 0; kk < mergeMe.size(); ++kk)
			{
				if ( eucl_dist( myClusterFamily[iCluster][iPoint], mergeMe[kk] ) < threshold ) 
				{

					if ( mergeMe[0].cluster_id != myClusterFamily[iCluster][iPoint].cluster_id )
					{

		//				first we remove the centroids of both clusters...
						mergeMe.pop_back();		myClusterFamily[iCluster].pop_back();
						
		//				then we give all points in mergeMe the same cluster ID as the first point within threshold	
		//					and we move all points in mergeMe into the same cluster vector as the first point within threshold
						for (int mergePoint = 0; mergePoint < mergeMe.size(); ++mergePoint)
						{
		//					give all points the same CID first	
							mergeMe[mergePoint].cluster_id = myClusterFamily[iCluster][iPoint].cluster_id;	
		//					now we put that point into place at the end of the cluster which was closest					
							myClusterFamily[iCluster].push_back( mergeMe[mergePoint] );
						}

		//				now we can find a new centroid and push it back into the new cluster
							pushNewCentroid( myClusterFamily[iCluster] );
		//				then we can remove mergeMe entirely from myClusterFamily because it is unused now
						myClusterFamily.erase( myClusterFamily.begin() );
						return;
					}		
				}
			}	
		}
	}
} 



void mergeToClosest( familyClusters &myClusterFamily, pointCluster &mergeMe )
{
	double tempDist = 0;
	double smallestDistance = 100000000; // temporary value for initialization
	int bestIndex;

//	first, we need to get access to all clusters in myClusterFamily
	for (int iCluster = 0; iCluster < myClusterFamily.size(); ++iCluster)
	{
//		make sure cluster ID is not the same
		if( myClusterFamily[iCluster][0].cluster_id != mergeMe[0].cluster_id )
		{

//		now we need access to all of the points within the other clusters
			for (int iPoint = 0; iPoint < myClusterFamily[iCluster].size(); ++iPoint)
			{

//				then, we find the CLOSEST POINT WE CAN to any point inside of cluster mergeMe and make them the same cluster
				for (int kk = 0; kk < mergeMe.size(); ++kk)
				{
					tempDist = eucl_dist( myClusterFamily[iCluster][iPoint], mergeMe[kk] );

					if ( tempDist < smallestDistance )
					{
						smallestDistance = tempDist;
						bestIndex = iCluster;
					}
				}	
			}
		}
	}



//				first we remove the centroids of both clusters...
				mergeMe.pop_back();		myClusterFamily[bestIndex].pop_back();
					
//				then we give all points in mergeMe the same cluster ID as the first point within threshold	
//					and we move all points in mergeMe into the same cluster vector as the first point within threshold
				for (int mergePoint = 0; mergePoint < mergeMe.size(); ++mergePoint)
				{
//					give all points the same CID first	
					mergeMe[mergePoint].cluster_id = myClusterFamily[bestIndex][0].cluster_id;	
//					now we put that point into place at the end of the cluster which was closest					
					myClusterFamily[bestIndex].push_back( mergeMe[mergePoint] );
				}

//				now we can find a new centroid and push it back into the new cluster
				pushNewCentroid( myClusterFamily[bestIndex] );
//				then we can remove mergeMe entirely from myClusterFamily because it is unused now
				myClusterFamily.erase( myClusterFamily.begin() );
				return;		
}



void pushNewCentroid( pointCluster &newMergerCl )
{
		Point centroidTemp;
		centroidTemp.x = getAvgX( newMergerCl );
		centroidTemp.y = getAvgY( newMergerCl );
		centroidTemp.z = getAvgZ( newMergerCl );
		centroidTemp.cluster_id = newMergerCl[0].cluster_id;
		newMergerCl.push_back(centroidTemp);
}



double changed_sq_eucl_dist( Point p1, Point p2, double NdZ, double MdY )
{
	double distMagnitude = 0;
	double dx, dy, dz;

	dx = p2.x - p1.x;
	dy = p2.y - p1.y - MdY;
	dz = p2.z - p1.z - NdZ;

	distMagnitude +=  sqrt( square( dx ) + square( dy ) + square( dz ) );

	return distMagnitude;
}



double eucl_dist(Point p1, Point p2)
{
	double distMagnitude = 0;
	distMagnitude += sqrt( square( p2.x - p1.x ) + square( p2.y - p1.y ) + square( p2.z - p1.z ) );

	return distMagnitude;
}


double getMinSquare( Point moveCent, Point closestCent, double minDistance ) {
//	ignore the same centroid that we are on because the distance between it and itself will always be zero
	if (closestCent.cluster_id != moveCent.cluster_id)
	{
		if ( square( eucl_dist( moveCent, closestCent ) ) < minDistance )
		{
			minDistance = square( eucl_dist( moveCent, closestCent ) );
		}
	}
	return minDistance;
}


double square(double x) {

	return x*x;
}


double getAvgX( pointCluster cl ) {
	double sum = 0;		int numPoints = 0;
	for (int iPoint = 0; iPoint < cl.size(); ++iPoint)
	{
		sum += cl[iPoint].x;
		++numPoints;
	}
	return sum/numPoints;
}


double getAvgY( pointCluster cl ) {
	double sum = 0;		int numPoints = 0;
	for (int iPoint = 0; iPoint < cl.size(); ++iPoint)
	{
		sum += cl[iPoint].y;
		++numPoints;
	}
	return sum/numPoints;
}


double getAvgZ( pointCluster cl ) {
	double sum = 0;		int numPoints = 0;
	for (int iPoint = 0; iPoint < cl.size(); ++iPoint)
	{
		sum += cl[iPoint].z;
		++numPoints;
	}
	return sum/numPoints;
}


void sameCluster( Point &p1, Point &p2 ) {

	if (p1.cluster_id < p2.cluster_id)
	{
		p2.cluster_id = p1.cluster_id;
	}

	else p1.cluster_id = p2.cluster_id;
}



pointCluster tempMerge( pointCluster moveMe, pointCluster bestCluster, double NdZ, double MdY )
{
	pointCluster tempCluster;
// 	first we need to remove all centroids
	moveMe.pop_back(); 	bestCluster.pop_back();

// 	 	loop through all points in bestCluster and add them to tempCluster
		for (int i = 0; i < bestCluster.size(); i++) {
			tempCluster.push_back( bestCluster[i] );
		}
		
// 		loop through all points in moveMe and move them then add them to tempCluster
		for (int i = 0; i < moveMe.size(); i++) {
			moveMe[i].y += MdY; 	moveMe[i].z += NdZ;
			tempCluster.push_back( moveMe[i] );
		}

		return tempCluster;
}


double movedVal( double oldVal, double shifterVal )
{
	return oldVal += shifterVal;	
}


void findBestClusterMoveAllPoints( const Point &moveCent, const familyClusters &otherClusters, const pointCluster &myCluster, double &newNdz, double &newMdy, int &bestCID, int &tbestCID, int &eliminate, int &connector, double &distanceCheck)
{
// ======= Adjustable parameters =======================
	double threshold = 0.6;
	// worked well on 1cm & also was at 1.2
	double holdIt = 1.0;
	int bigClusters = 8;
// =====================================================

	double tempDist = 0.0;
	double anotherTempDist;
	double smallestDistance = 1000000000.0;  // temporary value for initialization later
	double tsmallestDistance = 1000000000.0;  // another temporary value for initialization later
	int bestN = 99; 		int bestM = 99;
	int tempCID = -1;
	PCAResults results;
//	N is for z 	&  M is for y
	double dZ = 4.5;
	double dY = 2.4;
	double linearity;
	double bestLinearity = 0;

// 	the edges of the board in Y are always the same for both boards
	double BoardminY = -18; 	double BoardmaxY = 18;
	double BoardminZ = 12; 		double BoardmaxZ = 86;

// 7 z 14 y
PCAResults endPoints = DoPCA( myCluster );
// printf( "End Points:  ( %f, %f, %f )    &&    ( %f, %f, %f ) \n", endPoints.endPoints.first(0), endPoints.endPoints.first(1), endPoints.endPoints.first(2), endPoints.endPoints.second(0), endPoints.endPoints.second(1), endPoints.endPoints.second(2) );


//	loop through each combination of moving 0, +/- 1, and +/- 2 ROI in y and z direction for pair of centroids passed in function
	for (int N = -16; N <= 16; ++N)
	{
		// was -15 to 15
		for (int M = -30; M <= 30; ++M)
		{
// 			if the centroid of the cluster we move is outside of the detector -> then we continue
			if ( movedVal( moveCent.z, N*dZ ) < 6 || movedVal( moveCent.z, N*dZ ) > 90 || movedVal( moveCent.y, M*dY ) < -24 || movedVal( moveCent.y, M*dY ) > 24  ||( movedVal(moveCent.z, N*dZ) < 54.5 && movedVal(moveCent.z, N*dZ) > 50 && movedVal(moveCent.y, M*dY) < -13.2 ) ) {
				continue;
			}


//			loop through all available clusters
			for (int iCent = 0; iCent < otherClusters.size(); ++iCent)
			{
//				as long as CID is not the same as what I want to move then we continue
				if ( otherClusters[iCent][0].cluster_id != moveCent.cluster_id )
				{
					
// 					make a simpler name for the centroid of the other cluster
					// Point OtherCentroid = otherClusters[iCent][ otherClusters[iCent].size() - 1 ];

// 					store second best cluster choice
					/* if( bestCID != tempCID [> && linearity > 0.90 <]) {
					 *         tbestCID = tempCID;
					 *         tempCID = bestCID;
					 *         tsmallestDistance = anotherTempDist;
					 *         anotherTempDist = smallestDistance;
					 * } */

						results = DoPCA( tempMerge( myCluster, otherClusters[iCent], N*dZ, M*dY ) );
						linearity = ( results.eVals(0) ) / ( results.eVals(0) + results.eVals(1) + results.eVals(2) ); 


/*                                         find where the centroid is located -> less than or greater than 50cm in z
 *                                         if ( moveCent.z > 50 ) {
 *                                                 set the min and maximum values for the edges of the detector if greater than 50cm in z
 *                                                 BoardminZ = 50; 	BoardmaxZ = 86;
 *                                                         if we are in a dead ROI -> continue to next cluster
 *                                                         if ( moveCent.z < 54.5 && moveCent.y < -13.2 && moveCent.y > -18 ){
 *                                                                 continue;					
 *                                                         }
 *                                         }
 * 
 *                                         if the centroid of the closest cluster is not within the same board as the cluster we want to move, then we move to the next centroid
 *                                         if ( (otherClusters.size() > 4) && (OtherCentroid.z < BoardminZ || OtherCentroid.z > BoardmaxZ || OtherCentroid.y < BoardminY || OtherCentroid.y > BoardmaxY) ) {
 *                                                 continue;
 *                                         } */

//					loop all other points in otherclusters
					for (int iPoint = 0; iPoint < (otherClusters[iCent].size() - 1); ++iPoint)
					{
//						now we loop through all other points in the same cluster and moveCent and find the closest point to the closest point						
						for( int centPoints = 0; centPoints < (myCluster.size() - 1); ++centPoints )
						{
							tempDist = changed_sq_eucl_dist( myCluster[centPoints], otherClusters[ iCent ][ iPoint ], N*dZ, M*dY );		
	//						then, if our temp dist is smaller than our smallest of the small distances, AND the principle eigen value is bigger  
	//							store it and our N and M values to loop again!
							/* if(movenum == 103){
							 *               printf("Nval: %d  Mval:  %d   Coords: (%f, %f, %f)    CID: %f    distance: %f\n", N,M,otherClusters[iCent][iPoint].x, otherClusters[iCent][iPoint].y, otherClusters[iCent][iPoint].z, otherClusters[iCent][iPoint].cluster_id, tempDist);
							 * } */
							if ( (tempDist < smallestDistance) )
							{
								if ( tempDist < threshold && linearity > 0.90 ) {
									if ( linearity > bestLinearity ) {
										smallestDistance = tempDist;
										bestN = N;
										bestM = M;
										bestCID = iCent;
										bestLinearity = linearity;
									//	printf(" Coords: (%f, %f, %f) 		CID: %f   2CID: %f 		best distance: %f\n", otherClusters[iCent][iPoint].x, otherClusters[iCent][iPoint].y, otherClusters[iCent][iPoint].z, otherClusters[iCent][iPoint].cluster_id, 2bestCID, smallestDistance);
										// printf( " \nBest CID: %d, 	Runner up: %d   &&    DIST: %f \n", bestCID, tbestCID, tsmallestDistance);
										std::cout << "linearity:  " << (results.eVals(0)) / ( results.eVals(0) + results.eVals(1) + results.eVals(2) ) << "\n";
									}
								}
					
								else {
									smallestDistance = tempDist;
									bestCID = iCent;
									bestN = N;
									bestM = M;
								//	printf(" Coords: (%f, %f, %f) 		CID: %f   2CID: %f 		best distance: %f\n", otherClusters[iCent][iPoint].x, otherClusters[iCent][iPoint].y, otherClusters[iCent][iPoint].z, otherClusters[iCent][iPoint].cluster_id, 2bestCID, smallestDistance);
									// printf( " Best CID: %d, 	Runner up: %d   &&    DIST: %f \n", bestCID, tbestCID, tsmallestDistance);
								}
							}	
						}
					}		
				}
			}
		}
	}

	
	printf( "\n ~~> Dist: %f \n", smallestDistance );

	/* if ( smallestDistance > holdIt [> || (bestN && bestM) == 99) <] && otherClusters[0].size() < bigClusters) {
	 *         eliminate = 1;
	 * }     */

	
	/* if ( tsmallestDistance < 0.5 && tbestCID != bestCID && tbestCID != -1) {
	 *         printf( " \nBest CID: %d, 	Runner up: %d   &&    DIST: %f \n", bestCID, tbestCID, tsmallestDistance);
	 *         // put a flag here that we can use back in the main loop for knowing if we want to merge a second cluster
	 *         connector = 1;
	 * }else{ printf("We couldnt find shit yo!\n"); } */


//	store our new best values for everything here to be used back in the for loop in the main function
distanceCheck = smallestDistance;
newNdz = bestN*dZ;
newMdy = bestM*dY;
}




PCAResults DoPCA(const pointCluster &points) {

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
    return results; // FAIL FROM NO INPUT POINTS
  }

  const float nThreeDHitsAsFloat(static_cast<float>(nThreeDHits));
  meanPosition[0] /= nThreeDHitsAsFloat;
  meanPosition[1] /= nThreeDHitsAsFloat;
  meanPosition[2] /= nThreeDHitsAsFloat;
  outputCentroid = TVector3(meanPosition[0], meanPosition[1], meanPosition[2]);

  // Define elements of our covariance matrix
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

      xi2  += x * x;
      xiyi += x * y;
      xizi += x * z;
      yi2  += y * y;
      yizi += y * z;
      zi2  += z * z;
      weightSum += weight * weight;
  }

  // Using Eigen package
  Eigen::Matrix3f sig;

  sig << xi2, xiyi, xizi,
         xiyi, yi2, yizi,
         xizi, yizi, zi2;

  //if (std::fabs(weightSum) < std::numeric_limits<float>::epsilon())
  //{
  //    std::cout << "PCAShowerParticleBuildingAlgorithm::RunPCA - The total weight of three dimensional hits is 0!" << std::endl;
  //    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
  //}

  sig *= 1.0 / weightSum;

  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigenMat(sig);

  //if (eigenMat.info() != Eigen::ComputationInfo::Success)
  //{
  //    std::cout << "PCAShowerParticleBuildingAlgorithm::RunPCA - PCA decompose failure, number of three D hits = " << nThreeDHits << std::endl;
  //    throw StatusCodeException(STATUS_CODE_INVALID_PARAMETER);
  //}

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

//  could replace these projTestPoint with just TestPoint    
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
