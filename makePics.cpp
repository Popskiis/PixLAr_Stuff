#include <iostream>
#include "TNtuple.h"
#include "TTree.h"
#include "TFile.h"

void makePics()
{
	TFile* f = TFile::Open("test.root");	
		if (f == 0) {
			std::cout << "Could not open the requested file.\n";
			return;
		}

TTreeReader myReader("treeTest", f);
TTree* drawTree = myReader.GetTree();
TTreeReaderArray<double> getX(myReader, "sp_x");
TTreeReaderArray<double> getY(myReader, "sp_y");
TTreeReaderArray<double> getZ(myReader, "sp_z");
//TTreeReaderArray<int> getQ(myReader, "sp_charge");
TTreeReaderValue<int> getEntry(myReader, "entry");


	while(myReader.Next())
	{
		int tempEntry = *getEntry;
		TCanvas *c1 = new TCanvas("c1","c1");
		c1->cd();
		drawTree->Draw("sp_z:sp_x", Form("entry == %d", tempEntry));
		c1->SaveAs(Form("myPics/track_zx_%d.png", tempEntry));

		TCanvas *c2 = new TCanvas("c2","c2");
		c2->cd();
		drawTree->Draw("sp_y:sp_x", Form("entry == %d", tempEntry));
		c2->SaveAs(Form("myPics/track_yx_%d.png", tempEntry));
	}

f->Close();
}

