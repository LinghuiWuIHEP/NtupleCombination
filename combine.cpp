#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TKey.h"
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph2D.h"
#include "TPolyLine3D.h"
#include "TDirectory.h"
#include "TDirectoryFile.h"

#include "RecEvent.h"
#include "ImpactPoints.h"
#include "XCET.h"

using namespace std;

int main(int argc, char* argv[]){
     int runNo = 0;
     if(argc > 1){
	  sscanf(argv[1], "%d", &runNo);
     } else{
	  cout << "too few arguments" << endl;
	  return -1;
     }

     char pathHG[200];
     sprintf(pathHG, "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/ntuples/v13");
     // sprintf(pathHG, "/afs/cern.ch/work/l/lwu/public/ntuple");
     char fnameHG[200];
     sprintf(fnameHG, "%s/ntuple_%d.root", pathHG, runNo);

     char pathAH[200];
     sprintf(pathAH, "/eos/cms/store/group/dpg_hgcal/tb_hgcal/2018/cern_h2_october/offline_analysis/AHCAL_ntuples/v3");
     // sprintf(pathAH, "/afs/cern.ch/work/l/lwu/public/ntuple");
     char fnameAH[200];
     sprintf(fnameAH, "%s/Run_%d.root", pathAH, runNo);

     // open HGCAL and AHCAL root files and get the trees
     TFile fHG(fnameHG);
     TDirectory* dirRechit = (TDirectory*)fHG.Get("rechitntupler");
     TTree* trRechit = (TTree*)dirRechit->FindObjectAny("hits");
     int nHgEvents = trRechit->GetEntries();
     cout << "Read file: " << fnameHG << endl;
     cout << "rechit entries: " << nHgEvents << endl;

     TDirectory* dirTrk = (TDirectory*)fHG.Get("trackimpactntupler");
     TTree* trTrk = (TTree*)dirTrk->FindObjectAny("impactPoints");
     cout << "track entries: " << trTrk->GetEntries() << endl;

     TDirectory* dirXcet = (TDirectory*)fHG.Get("XCETntupler");
     TTree* trXcet = (TTree*)dirXcet->FindObjectAny("XCET");
     cout << "Xcet entries: " << trXcet->GetEntries() << endl;

     TFile fAH(fnameAH);
     TTree* trAh = (TTree*)fAH.Get("bigtree");
     int nAhEvents = trAh->GetEntries();
     cout << "Read file: " << fnameAH << endl;
     cout << "Ah entries: " << nAhEvents << endl;

     int nEvents = (nHgEvents>nAhEvents)? nAhEvents : nHgEvents;

     // set output file
     char fnameOut[200];
     sprintf(fnameOut, "fullNtuple_%d.root", runNo);
     TFile outFile(fnameOut, "recreate");
     outFile.cd();

     // initialize classes to read the trees and set new trees which will be output
     RecEvent* rec = new RecEvent(trRechit, trAh);
     TTree* newTrHits = new TTree("hits", "hits");
     rec->SetNewTree(newTrHits);

     ImpactPoints* trk = new ImpactPoints(trTrk);
     TTree* newTrTrk = new TTree("impactPoints", "impactPoints");
     trk->SetNewTree(newTrTrk);

     XCET* pXcet = new XCET(trXcet);
     TTree* newTrXcet = new TTree("XCET", "XCET");
     pXcet->SetNewTree(newTrXcet);

     for(int i=0; i<nEvents; i++){
	  rec->GetEntry(i);
	  trk->GetEntry(i);
	  pXcet->GetEntry(i);
	  UInt_t evtRec = rec->event;

	  if(0==(evtRec%1000)) cout << "event " << evtRec << endl;

	  // removing the last extra events in HGCAL
	  if(i>=nAhEvents) break;

	  rec->FillNewTree();
	  trk->FillNewTree();
	  pXcet->FillNewTree();
     }

     // set directory in the output root file
     outFile.cd();
     TDirectoryFile* newDirRechit = new TDirectoryFile("rechitntupler", "rechitntupler");
     newDirRechit->Add(newTrHits);

     TDirectoryFile* newDirTrk = new TDirectoryFile("trackimpactntupler", "trackimpactntupler");
     newDirTrk->Add(newTrTrk);

     TDirectoryFile* newDirXcet = new TDirectoryFile("XCETntupler", "XCETntupler");
     newDirXcet->Add(newTrXcet);

     // output
     outFile.cd();
     newDirRechit->Write();
     newDirTrk->Write();
     newDirXcet->Write();
     outFile.Close();
     cout << "full ntuple was written" << endl;

     return 0;
}

