//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed May 15 11:08:25 2019 by ROOT version 6.10/05
// from TTree bigtree/bigtree
// found on file: Run_700.root
//////////////////////////////////////////////////////////

#ifndef AhTree_h
#define AhTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class AhTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *newTr;   //!pointer to the analyzed TTree or TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           eventNumber;
   Long64_t        eventTime;
   Int_t           ahc_iEvt;
   Int_t           ahc_nHits;
   Int_t           ahc_nLayers;
   Float_t         ahc_energySum;
   Float_t         ahc_energyDensity;
   Float_t         ahc_radius;
   Float_t         ahc_radiusEw;
   Float_t         ahc_cogX;
   Float_t         ahc_cogY;
   Float_t         ahc_cogZ;
   Float_t         ahc_frac25;
   Float_t         ahc_energySum5Layer;
   Int_t           ahc_nHits5Layer;
   Float_t         ahc_cogX5Layer;
   Float_t         ahc_cogY5Layer;
   Float_t         ahc_cogZ5Layer;
   vector<float>   *ahc_energyPerLayer;
   vector<int>     *ahc_nHitsPerLayer;
   vector<float>   *ahc_energyPerLayer_err;
   vector<int>     *ahc_cellSize;
   vector<float>   *ahc_cogXPerLayer;
   vector<float>   *ahc_cogYPerLayer;
   vector<float>   *ahc_radiusPerLayer;
   vector<float>   *ahc_radiusEwPerLayer;
   vector<int>     *ahc_hitCellID;
   vector<int>     *ahc_hitI;
   vector<int>     *ahc_hitJ;
   vector<int>     *ahc_hitK;
   vector<float>   *ahc_hitEnergy;
   vector<float>   *ahc_hitTime;
   vector<int>     *ahc_hitType;
   vector<float>   *ahc_hitRadius;
   vector<float>   *ahc_hitEnergyDensity;
   vector<float>   *ahc_hitX;
   vector<float>   *ahc_hitY;
   vector<float>   *ahc_hitZ;
   ULong64_t       lda_trigTime;
   Int_t           event_BXID;

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_eventTime;   //!
   TBranch        *b_ahc_iEvt;   //!
   TBranch        *b_ahc_nHits;   //!
   TBranch        *b_ahc_nLayers;   //!
   TBranch        *b_ahc_energySum;   //!
   TBranch        *b_ahc_energyDensity;   //!
   TBranch        *b_ahc_radius;   //!
   TBranch        *b_ahc_radiusEw;   //!
   TBranch        *b_ahc_cogX;   //!
   TBranch        *b_ahc_cogY;   //!
   TBranch        *b_ahc_cogZ;   //!
   TBranch        *b_ahc_frac25;   //!
   TBranch        *b_ahc_energySum5Layer;   //!
   TBranch        *b_ahc_nHits5Layer;   //!
   TBranch        *b_ahc_cogX5Layer;   //!
   TBranch        *b_ahc_cogY5Layer;   //!
   TBranch        *b_ahc_cogZ5Layer;   //!
   TBranch        *b_ahc_energyPerLayer;   //!
   TBranch        *b_ahc_nHitsPerLayer;   //!
   TBranch        *b_ahc_energyPerLayer_err;   //!
   TBranch        *b_ahc_cellSize;   //!
   TBranch        *b_ahc_cogXPerLayer;   //!
   TBranch        *b_ahc_cogYPerLayer;   //!
   TBranch        *b_ahc_radiusPerLayer;   //!
   TBranch        *b_ahc_radiusEwPerLayer;   //!
   TBranch        *b_ahc_hitCellID;   //!
   TBranch        *b_ahc_hitI;   //!
   TBranch        *b_ahc_hitJ;   //!
   TBranch        *b_ahc_hitK;   //!
   TBranch        *b_ahc_hitEnergy;   //!
   TBranch        *b_ahc_hitTime;   //!
   TBranch        *b_ahc_hitType;   //!
   TBranch        *b_ahc_hitRadius;   //!
   TBranch        *b_ahc_hitEnergyDensity;   //!
   TBranch        *b_ahc_hitX;   //!
   TBranch        *b_ahc_hitY;   //!
   TBranch        *b_ahc_hitZ;   //!
   TBranch        *b_lda_trigTime;   //!
   TBranch        *b_event_BXID;   //!

   AhTree(TTree *tree=0);
   virtual ~AhTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     SetNewTree(TTree* tr=0);
   virtual void     FillNewTree();
};

#endif

#ifdef AhTree_cxx
AhTree::AhTree(TTree *tree) : fChain(0), newTr(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   /* if (tree == 0) { */
   /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Run_700.root"); */
   /*    if (!f || !f->IsOpen()) { */
   /*       f = new TFile("Run_700.root"); */
   /*    } */
   /*    f->GetObject("bigtree",tree); */

   /* } */
   Init(tree);
}

AhTree::~AhTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AhTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AhTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AhTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   ahc_energyPerLayer = 0;
   ahc_nHitsPerLayer = 0;
   ahc_energyPerLayer_err = 0;
   ahc_cellSize = 0;
   ahc_cogXPerLayer = 0;
   ahc_cogYPerLayer = 0;
   ahc_radiusPerLayer = 0;
   ahc_radiusEwPerLayer = 0;
   ahc_hitCellID = 0;
   ahc_hitI = 0;
   ahc_hitJ = 0;
   ahc_hitK = 0;
   ahc_hitEnergy = 0;
   ahc_hitTime = 0;
   ahc_hitType = 0;
   ahc_hitRadius = 0;
   ahc_hitEnergyDensity = 0;
   ahc_hitX = 0;
   ahc_hitY = 0;
   ahc_hitZ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("eventTime", &eventTime, &b_eventTime);
   fChain->SetBranchAddress("ahc_iEvt", &ahc_iEvt, &b_ahc_iEvt);
   fChain->SetBranchAddress("ahc_nHits", &ahc_nHits, &b_ahc_nHits);
   fChain->SetBranchAddress("ahc_nLayers", &ahc_nLayers, &b_ahc_nLayers);
   fChain->SetBranchAddress("ahc_energySum", &ahc_energySum, &b_ahc_energySum);
   fChain->SetBranchAddress("ahc_energyDensity", &ahc_energyDensity, &b_ahc_energyDensity);
   fChain->SetBranchAddress("ahc_radius", &ahc_radius, &b_ahc_radius);
   fChain->SetBranchAddress("ahc_radiusEw", &ahc_radiusEw, &b_ahc_radiusEw);
   fChain->SetBranchAddress("ahc_cogX", &ahc_cogX, &b_ahc_cogX);
   fChain->SetBranchAddress("ahc_cogY", &ahc_cogY, &b_ahc_cogY);
   fChain->SetBranchAddress("ahc_cogZ", &ahc_cogZ, &b_ahc_cogZ);
   fChain->SetBranchAddress("ahc_frac25", &ahc_frac25, &b_ahc_frac25);
   fChain->SetBranchAddress("ahc_energySum5Layer", &ahc_energySum5Layer, &b_ahc_energySum5Layer);
   fChain->SetBranchAddress("ahc_nHits5Layer", &ahc_nHits5Layer, &b_ahc_nHits5Layer);
   fChain->SetBranchAddress("ahc_cogX5Layer", &ahc_cogX5Layer, &b_ahc_cogX5Layer);
   fChain->SetBranchAddress("ahc_cogY5Layer", &ahc_cogY5Layer, &b_ahc_cogY5Layer);
   fChain->SetBranchAddress("ahc_cogZ5Layer", &ahc_cogZ5Layer, &b_ahc_cogZ5Layer);
   fChain->SetBranchAddress("ahc_energyPerLayer", &ahc_energyPerLayer, &b_ahc_energyPerLayer);
   fChain->SetBranchAddress("ahc_nHitsPerLayer", &ahc_nHitsPerLayer, &b_ahc_nHitsPerLayer);
   fChain->SetBranchAddress("ahc_energyPerLayer_err", &ahc_energyPerLayer_err, &b_ahc_energyPerLayer_err);
   fChain->SetBranchAddress("ahc_cellSize", &ahc_cellSize, &b_ahc_cellSize);
   fChain->SetBranchAddress("ahc_cogXPerLayer", &ahc_cogXPerLayer, &b_ahc_cogXPerLayer);
   fChain->SetBranchAddress("ahc_cogYPerLayer", &ahc_cogYPerLayer, &b_ahc_cogYPerLayer);
   fChain->SetBranchAddress("ahc_radiusPerLayer", &ahc_radiusPerLayer, &b_ahc_radiusPerLayer);
   fChain->SetBranchAddress("ahc_radiusEwPerLayer", &ahc_radiusEwPerLayer, &b_ahc_radiusEwPerLayer);
   fChain->SetBranchAddress("ahc_hitCellID", &ahc_hitCellID, &b_ahc_hitCellID);
   fChain->SetBranchAddress("ahc_hitI", &ahc_hitI, &b_ahc_hitI);
   fChain->SetBranchAddress("ahc_hitJ", &ahc_hitJ, &b_ahc_hitJ);
   fChain->SetBranchAddress("ahc_hitK", &ahc_hitK, &b_ahc_hitK);
   fChain->SetBranchAddress("ahc_hitEnergy", &ahc_hitEnergy, &b_ahc_hitEnergy);
   fChain->SetBranchAddress("ahc_hitTime", &ahc_hitTime, &b_ahc_hitTime);
   fChain->SetBranchAddress("ahc_hitType", &ahc_hitType, &b_ahc_hitType);
   fChain->SetBranchAddress("ahc_hitRadius", &ahc_hitRadius, &b_ahc_hitRadius);
   fChain->SetBranchAddress("ahc_hitEnergyDensity", &ahc_hitEnergyDensity, &b_ahc_hitEnergyDensity);
   fChain->SetBranchAddress("ahc_hitX", &ahc_hitX, &b_ahc_hitX);
   fChain->SetBranchAddress("ahc_hitY", &ahc_hitY, &b_ahc_hitY);
   fChain->SetBranchAddress("ahc_hitZ", &ahc_hitZ, &b_ahc_hitZ);
   fChain->SetBranchAddress("lda_trigTime", &lda_trigTime, &b_lda_trigTime);
   fChain->SetBranchAddress("event_BXID", &event_BXID, &b_event_BXID);
   Notify();
}

Bool_t AhTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AhTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AhTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void AhTree::SetNewTree(TTree* tr){
     newTr = tr;
     newTr->Branch("runNumber", &runNumber, "runNumber/I");
     newTr->Branch("eventNumber", &eventNumber, "eventNumber/I");
     newTr->Branch("eventTime", &eventTime, "eventTime/L");
     newTr->Branch("ahc_iEvt", &ahc_iEvt, "ahc_iEvt/I");
     newTr->Branch("ahc_nHits", &ahc_nHits, "ahc_nHits/I");
     newTr->Branch("ahc_nLayers", &ahc_nLayers, "ahc_nLayers/I");
     newTr->Branch("ahc_energySum", &ahc_energySum, "ahc_energySum/F");
     newTr->Branch("ahc_energyDensity", &ahc_energyDensity, "ahc_energyDensity/F");
     newTr->Branch("ahc_radius", &ahc_radius, "ahc_radius/F");
     newTr->Branch("ahc_radiusEw", &ahc_radiusEw, "ahc_radiusEw/F");
     newTr->Branch("ahc_cogX", &ahc_cogX, "ahc_cogX/F");
     newTr->Branch("ahc_cogY", &ahc_cogY, "ahc_cogY/F");
     newTr->Branch("ahc_cogZ", &ahc_cogZ, "ahc_cogZ/F");
     newTr->Branch("ahc_frac25", &ahc_frac25, "ahc_frac25/F");
     newTr->Branch("ahc_energySum5Layer", &ahc_energySum5Layer, "ahc_energySum5Layer/F");
     newTr->Branch("ahc_nHits5Layer", &ahc_nHits5Layer, "ahc_nHits5Layer/I");
     newTr->Branch("ahc_cogX5Layer", &ahc_cogX5Layer, "ahc_cogX5Layer/F");
     newTr->Branch("ahc_cogY5Layer", &ahc_cogY5Layer, "ahc_cogY5Layer/F");
     newTr->Branch("ahc_cogZ5Layer", &ahc_cogZ5Layer, "ahc_cogZ5Layer/F");
     newTr->Branch("ahc_energyPerLayer", ahc_energyPerLayer);
     newTr->Branch("ahc_nHitsPerLayer", ahc_nHitsPerLayer);
     newTr->Branch("ahc_energyPerLayer_err", ahc_energyPerLayer_err);
     newTr->Branch("ahc_cellSize", ahc_cellSize);
     newTr->Branch("ahc_cogXPerLayer", ahc_cogXPerLayer);
     newTr->Branch("ahc_cogYPerLayer", ahc_cogYPerLayer);
     newTr->Branch("ahc_radiusPerLayer", ahc_radiusPerLayer);
     newTr->Branch("ahc_radiusEwPerLayer", ahc_radiusEwPerLayer);
     newTr->Branch("ahc_hitCellID", ahc_hitCellID);
     newTr->Branch("ahc_hitI", ahc_hitI);
     newTr->Branch("ahc_hitJ", ahc_hitJ);
     newTr->Branch("ahc_hitK", ahc_hitK);
     newTr->Branch("ahc_hitEnergy", ahc_hitEnergy);
     newTr->Branch("ahc_hitTime", ahc_hitTime);
     newTr->Branch("ahc_hitType", ahc_hitType);
     newTr->Branch("ahc_hitRadius", ahc_hitRadius);
     newTr->Branch("ahc_hitEnergyDensity", ahc_hitEnergyDensity);
     newTr->Branch("ahc_hitX", ahc_hitX);
     newTr->Branch("ahc_hitY", ahc_hitY);
     newTr->Branch("ahc_hitZ", ahc_hitZ);
     newTr->Branch("lda_trigTime", lda_trigTime);
     newTr->Branch("event_BXID", event_BXID);
     newTr->Branch("lda_trigTime", &lda_trigTime, "lda_trigTime/l");
     newTr->Branch("event_BXID", &event_BXID, "event_BXID/I");
}
void AhTree::FillNewTree(){
     newTr->Fill();
}
#endif // #ifdef AhTree_cxx
