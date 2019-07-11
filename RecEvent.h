//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Apr 29 15:22:05 2019 by ROOT version 6.10/05
// from TTree hits/HGC rechits
// found on file: ntuple_700.root
//////////////////////////////////////////////////////////

#ifndef RecEvent_h
#define RecEvent_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>

using namespace std;

class RecEvent {
public :
   TTree          *fChain1;   //!pointer to the analyzed TTree or TChain
   TTree          *fChain2;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent1; //!current Tree number in a TChain
   Int_t           fCurrent2; //!current Tree number in a TChain
   TTree          *newTr;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types of tree1
   UInt_t          event;
   ULong64_t       trigger_timestamp;
   UInt_t          run;
   Int_t           pdgID;
   Float_t         beamEnergy;
   Float_t         trueBeamEnergy;
   Float_t         energyLostEE;
   Float_t         energyLostFH;
   Float_t         energyLostBH;
   Float_t         energyLostBeam;
   Float_t         energyLostOutside;
   UInt_t          NRechits;
   vector<unsigned int> *rechit_detid;
   vector<unsigned int> *rechit_module;
   vector<unsigned int> *rechit_layer;
   vector<unsigned int> *rechit_chip;
   vector<unsigned int> *rechit_channel;
   vector<unsigned int> *rechit_type;
   vector<float>   *rechit_x;
   vector<float>   *rechit_y;
   vector<float>   *rechit_z;
   vector<short>   *rechit_iu;
   vector<short>   *rechit_iv;
   vector<short>   *rechit_iU;
   vector<short>   *rechit_iV;
   vector<float>   *rechit_energy;
   vector<float>   *rechit_energy_noHG;
   vector<float>   *rechit_amplitudeHigh;
   vector<float>   *rechit_amplitudeLow;
   vector<bool>    *rechit_hg_goodFit;
   vector<bool>    *rechit_lg_goodFit;
   vector<bool>    *rechit_hg_saturated;
   vector<bool>    *rechit_lg_saturated;
   vector<bool>    *rechit_fully_calibrated;
   vector<bool>    *rechit_noise_flag;
   vector<float>   *rechit_TS2High;
   vector<float>   *rechit_TS2Low;
   vector<float>   *rechit_TS3High;
   vector<float>   *rechit_TS3Low;
   vector<unsigned short> *rechit_Tot;
   vector<short>   *rechit_toa_calib_flag;
   vector<short>   *rechit_toaFall_flag;
   vector<short>   *rechit_toaRise_flag;
   vector<float>   *rechit_toaFall_norm;
   vector<float>   *rechit_toaRise_norm;
   vector<float>   *rechit_toaFall_time;
   vector<float>   *rechit_toaRise_time;
   vector<float>   *rechit_toaFall_corr_time;
   vector<float>   *rechit_toaRise_corr_time;
   vector<float>   *rechit_calib_time_toaFall;
   vector<float>   *rechit_calib_time_toaRise;
   vector<unsigned short> *rechit_toaRise;
   vector<unsigned short> *rechit_toaFall;
   vector<float>   *rechit_timeMaxHG;
   vector<float>   *rechit_timeMaxLG;

   // Declaration of leaf types of tree2
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


   // List of branches of tree1
   TBranch        *b_event;   //!
   TBranch        *b_trigger_timestamp;   //!
   TBranch        *b_run;   //!
   TBranch        *b_pdgID;   //!
   TBranch        *b_beamEnergy;   //!
   TBranch        *b_trueBeamEnergy;   //!
   TBranch        *b_energyLostEE;   //!
   TBranch        *b_energyLostFH;   //!
   TBranch        *b_energyLostBH;   //!
   TBranch        *b_energyLostBeam;   //!
   TBranch        *b_energyLostOutside;   //!
   TBranch        *b_NRechits;   //!
   TBranch        *b_rechit_detid;   //!
   TBranch        *b_rechit_module;   //!
   TBranch        *b_rechit_layer;   //!
   TBranch        *b_rechit_chip;   //!
   TBranch        *b_rechit_channel;   //!
   TBranch        *b_rechit_type;   //!
   TBranch        *b_rechit_x;   //!
   TBranch        *b_rechit_y;   //!
   TBranch        *b_rechit_z;   //!
   TBranch        *b_rechit_iu;   //!
   TBranch        *b_rechit_iv;   //!
   TBranch        *b_rechit_iU;   //!
   TBranch        *b_rechit_iV;   //!
   TBranch        *b_rechit_energy;   //!
   TBranch        *b_rechit_energy_noHG;   //!
   TBranch        *b_rechit_amplitudeHigh;   //!
   TBranch        *b_rechit_amplitudeLow;   //!
   TBranch        *b_rechit_hg_goodFit;   //!
   TBranch        *b_rechit_lg_goodFit;   //!
   TBranch        *b_rechit_hg_saturated;   //!
   TBranch        *b_rechit_lg_saturated;   //!
   TBranch        *b_rechit_fully_calibrated;   //!
   TBranch        *b_rechit_noise_flag;   //!
   TBranch        *b_rechit_TS2High;   //!
   TBranch        *b_rechit_TS2Low;   //!
   TBranch        *b_rechit_TS3High;   //!
   TBranch        *b_rechit_TS3Low;   //!
   TBranch        *b_rechit_Tot;   //!
   TBranch        *b_rechit_toa_calib_flag;   //!
   TBranch        *b_rechit_toaFall_flag;   //!
   TBranch        *b_rechit_toaRise_flag;   //!
   TBranch        *b_rechit_toaFall_norm;   //!
   TBranch        *b_rechit_toaRise_norm;   //!
   TBranch        *b_rechit_toaFall_time;   //!
   TBranch        *b_rechit_toaRise_time;   //!
   TBranch        *b_rechit_toaFall_corr_time;   //!
   TBranch        *b_rechit_toaRise_corr_time;   //!
   TBranch        *b_rechit_calib_time_toaFall;   //!
   TBranch        *b_rechit_calib_time_toaRise;   //!
   TBranch        *b_rechit_toaRise;   //!
   TBranch        *b_rechit_toaFall;   //!
   TBranch        *b_rechit_timeMaxHG;   //!
   TBranch        *b_rechit_timeMaxLG;   //!

   // List of branches of tree2
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


   /* RecEvent(char* fname, TChain *tree=0); */
   RecEvent(TTree *tree1=0, TTree *tree2=0);
   virtual ~RecEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetTree1Entries();
   virtual Int_t    GetTree2Entries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree, TTree *tree2);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     SetNewTree(TTree* tr=0);
   virtual void     FillNewTree();
};

#endif

#ifdef RecEvent_cxx
RecEvent::RecEvent(TTree *tree1, TTree* tree2) : fChain1(0), fChain2(0), newTr(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   /* if (tree == 0) { */
   /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(fname); */
   /*    if (!f || !f->IsOpen()) { */
   /*       f = new TFile(fname); */
   /*    } */
   /*    char dirName[300]; */
   /*    sprintf(dirName, "%s:/rechitntupler", fname); */
   /*    TDirectory * dir = (TDirectory*)f->Get(dirName); */
   /*    dir->GetObject("hits",tree); */

   /* } */
     Init(tree1, tree2);
}

RecEvent::~RecEvent()
{
   if (!fChain1) return;
   delete fChain1->GetCurrentFile();
   if (!fChain2) return;
   delete fChain2->GetCurrentFile();
}

Int_t RecEvent::GetTree1Entries(){
     return fChain1->GetEntries();
}

Int_t RecEvent::GetTree2Entries(){
     return fChain2->GetEntries();
}

Int_t RecEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain1) return 0;
   if (!fChain2) return 0;
   fChain1->GetEntry(entry);
   fChain2->GetEntry(entry);
   return 1;
}
Long64_t RecEvent::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain1) return -5;
   Long64_t centry1 = fChain1->LoadTree(entry);
   Long64_t centry2 = fChain2->LoadTree(entry);
   if (centry1 < 0) return centry1;
   if (centry2 < 0) return centry2;
   if (fChain1->GetTreeNumber() != fCurrent1) {
      fCurrent1 = fChain1->GetTreeNumber();
      Notify();
   }
   return centry1;
}

void RecEvent::Init(TTree *tree1, TTree *tree2)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer of tree1
   rechit_detid = 0;
   rechit_module = 0;
   rechit_layer = 0;
   rechit_chip = 0;
   rechit_channel = 0;
   rechit_type = 0;
   rechit_x = 0;
   rechit_y = 0;
   rechit_z = 0;
   rechit_iu = 0;
   rechit_iv = 0;
   rechit_iU = 0;
   rechit_iV = 0;
   rechit_energy = 0;
   rechit_energy_noHG = 0;
   rechit_amplitudeHigh = 0;
   rechit_amplitudeLow = 0;
   rechit_hg_goodFit = 0;
   rechit_lg_goodFit = 0;
   rechit_hg_saturated = 0;
   rechit_lg_saturated = 0;
   rechit_fully_calibrated = 0;
   rechit_noise_flag = 0;
   rechit_TS2High = 0;
   rechit_TS2Low = 0;
   rechit_TS3High = 0;
   rechit_TS3Low = 0;
   rechit_Tot = 0;
   rechit_toa_calib_flag = 0;
   rechit_toaFall_flag = 0;
   rechit_toaRise_flag = 0;
   rechit_toaFall_norm = 0;
   rechit_toaRise_norm = 0;
   rechit_toaFall_time = 0;
   rechit_toaRise_time = 0;
   rechit_toaFall_corr_time = 0;
   rechit_toaRise_corr_time = 0;
   rechit_calib_time_toaFall = 0;
   rechit_calib_time_toaRise = 0;
   rechit_toaRise = 0;
   rechit_toaFall = 0;
   rechit_timeMaxHG = 0;
   rechit_timeMaxLG = 0;
   // Set branch addresses and branch pointers
   if (!tree1) return;
   fChain1 = tree1;
   fCurrent1 = -1;
   fChain1->SetMakeClass(1);

   fChain1->SetBranchAddress("event", &event, &b_event);
   fChain1->SetBranchAddress("trigger_timestamp", &trigger_timestamp, &b_trigger_timestamp);
   fChain1->SetBranchAddress("run", &run, &b_run);
   fChain1->SetBranchAddress("pdgID", &pdgID, &b_pdgID);
   fChain1->SetBranchAddress("beamEnergy", &beamEnergy, &b_beamEnergy);
   fChain1->SetBranchAddress("trueBeamEnergy", &trueBeamEnergy, &b_trueBeamEnergy);
   fChain1->SetBranchAddress("energyLostEE", &energyLostEE, &b_energyLostEE);
   fChain1->SetBranchAddress("energyLostFH", &energyLostFH, &b_energyLostFH);
   fChain1->SetBranchAddress("energyLostBH", &energyLostBH, &b_energyLostBH);
   fChain1->SetBranchAddress("energyLostBeam", &energyLostBeam, &b_energyLostBeam);
   fChain1->SetBranchAddress("energyLostOutside", &energyLostOutside, &b_energyLostOutside);
   fChain1->SetBranchAddress("NRechits", &NRechits, &b_NRechits);
   fChain1->SetBranchAddress("rechit_detid", &rechit_detid, &b_rechit_detid);
   fChain1->SetBranchAddress("rechit_module", &rechit_module, &b_rechit_module);
   fChain1->SetBranchAddress("rechit_layer", &rechit_layer, &b_rechit_layer);
   fChain1->SetBranchAddress("rechit_chip", &rechit_chip, &b_rechit_chip);
   fChain1->SetBranchAddress("rechit_channel", &rechit_channel, &b_rechit_channel);
   fChain1->SetBranchAddress("rechit_type", &rechit_type, &b_rechit_type);
   fChain1->SetBranchAddress("rechit_x", &rechit_x, &b_rechit_x);
   fChain1->SetBranchAddress("rechit_y", &rechit_y, &b_rechit_y);
   fChain1->SetBranchAddress("rechit_z", &rechit_z, &b_rechit_z);
   fChain1->SetBranchAddress("rechit_iu", &rechit_iu, &b_rechit_iu);
   fChain1->SetBranchAddress("rechit_iv", &rechit_iv, &b_rechit_iv);
   fChain1->SetBranchAddress("rechit_iU", &rechit_iU, &b_rechit_iU);
   fChain1->SetBranchAddress("rechit_iV", &rechit_iV, &b_rechit_iV);
   fChain1->SetBranchAddress("rechit_energy", &rechit_energy, &b_rechit_energy);
   fChain1->SetBranchAddress("rechit_energy_noHG", &rechit_energy_noHG, &b_rechit_energy_noHG);
   fChain1->SetBranchAddress("rechit_amplitudeHigh", &rechit_amplitudeHigh, &b_rechit_amplitudeHigh);
   fChain1->SetBranchAddress("rechit_amplitudeLow", &rechit_amplitudeLow, &b_rechit_amplitudeLow);
   fChain1->SetBranchAddress("rechit_hg_goodFit", &rechit_hg_goodFit, &b_rechit_hg_goodFit);
   fChain1->SetBranchAddress("rechit_lg_goodFit", &rechit_lg_goodFit, &b_rechit_lg_goodFit);
   fChain1->SetBranchAddress("rechit_hg_saturated", &rechit_hg_saturated, &b_rechit_hg_saturated);
   fChain1->SetBranchAddress("rechit_lg_saturated", &rechit_lg_saturated, &b_rechit_lg_saturated);
   fChain1->SetBranchAddress("rechit_fully_calibrated", &rechit_fully_calibrated, &b_rechit_fully_calibrated);
   fChain1->SetBranchAddress("rechit_noise_flag", &rechit_noise_flag, &b_rechit_noise_flag);
   fChain1->SetBranchAddress("rechit_TS2High", &rechit_TS2High, &b_rechit_TS2High);
   fChain1->SetBranchAddress("rechit_TS2Low", &rechit_TS2Low, &b_rechit_TS2Low);
   fChain1->SetBranchAddress("rechit_TS3High", &rechit_TS3High, &b_rechit_TS3High);
   fChain1->SetBranchAddress("rechit_TS3Low", &rechit_TS3Low, &b_rechit_TS3Low);
   fChain1->SetBranchAddress("rechit_Tot", &rechit_Tot, &b_rechit_Tot);
   fChain1->SetBranchAddress("rechit_toa_calib_flag", &rechit_toa_calib_flag, &b_rechit_toa_calib_flag);
   fChain1->SetBranchAddress("rechit_toaFall_flag", &rechit_toaFall_flag, &b_rechit_toaFall_flag);
   fChain1->SetBranchAddress("rechit_toaRise_flag", &rechit_toaRise_flag, &b_rechit_toaRise_flag);
   fChain1->SetBranchAddress("rechit_toaFall_norm", &rechit_toaFall_norm, &b_rechit_toaFall_norm);
   fChain1->SetBranchAddress("rechit_toaRise_norm", &rechit_toaRise_norm, &b_rechit_toaRise_norm);
   fChain1->SetBranchAddress("rechit_toaFall_time", &rechit_toaFall_time, &b_rechit_toaFall_time);
   fChain1->SetBranchAddress("rechit_toaRise_time", &rechit_toaRise_time, &b_rechit_toaRise_time);
   fChain1->SetBranchAddress("rechit_toaFall_corr_time", &rechit_toaFall_corr_time, &b_rechit_toaFall_corr_time);
   fChain1->SetBranchAddress("rechit_toaRise_corr_time", &rechit_toaRise_corr_time, &b_rechit_toaRise_corr_time);
   fChain1->SetBranchAddress("rechit_calib_time_toaFall", &rechit_calib_time_toaFall, &b_rechit_calib_time_toaFall);
   fChain1->SetBranchAddress("rechit_calib_time_toaRise", &rechit_calib_time_toaRise, &b_rechit_calib_time_toaRise);
   fChain1->SetBranchAddress("rechit_toaRise", &rechit_toaRise, &b_rechit_toaRise);
   fChain1->SetBranchAddress("rechit_toaFall", &rechit_toaFall, &b_rechit_toaFall);
   fChain1->SetBranchAddress("rechit_timeMaxHG", &rechit_timeMaxHG, &b_rechit_timeMaxHG);
   fChain1->SetBranchAddress("rechit_timeMaxLG", &rechit_timeMaxLG, &b_rechit_timeMaxLG);

   // Set object pointer of tree2
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
   // Set branch addresses and branch pointers of tree2
   if (!tree2) return;
   fChain2 = tree2;
   fCurrent2 = -1;
   fChain2->SetMakeClass(1);

   fChain2->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain2->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain2->SetBranchAddress("eventTime", &eventTime, &b_eventTime);
   fChain2->SetBranchAddress("ahc_iEvt", &ahc_iEvt, &b_ahc_iEvt);
   fChain2->SetBranchAddress("ahc_nHits", &ahc_nHits, &b_ahc_nHits);
   fChain2->SetBranchAddress("ahc_nLayers", &ahc_nLayers, &b_ahc_nLayers);
   fChain2->SetBranchAddress("ahc_energySum", &ahc_energySum, &b_ahc_energySum);
   fChain2->SetBranchAddress("ahc_energyDensity", &ahc_energyDensity, &b_ahc_energyDensity);
   fChain2->SetBranchAddress("ahc_radius", &ahc_radius, &b_ahc_radius);
   fChain2->SetBranchAddress("ahc_radiusEw", &ahc_radiusEw, &b_ahc_radiusEw);
   fChain2->SetBranchAddress("ahc_cogX", &ahc_cogX, &b_ahc_cogX);
   fChain2->SetBranchAddress("ahc_cogY", &ahc_cogY, &b_ahc_cogY);
   fChain2->SetBranchAddress("ahc_cogZ", &ahc_cogZ, &b_ahc_cogZ);
   fChain2->SetBranchAddress("ahc_frac25", &ahc_frac25, &b_ahc_frac25);
   fChain2->SetBranchAddress("ahc_energySum5Layer", &ahc_energySum5Layer, &b_ahc_energySum5Layer);
   fChain2->SetBranchAddress("ahc_nHits5Layer", &ahc_nHits5Layer, &b_ahc_nHits5Layer);
   fChain2->SetBranchAddress("ahc_cogX5Layer", &ahc_cogX5Layer, &b_ahc_cogX5Layer);
   fChain2->SetBranchAddress("ahc_cogY5Layer", &ahc_cogY5Layer, &b_ahc_cogY5Layer);
   fChain2->SetBranchAddress("ahc_cogZ5Layer", &ahc_cogZ5Layer, &b_ahc_cogZ5Layer);
   fChain2->SetBranchAddress("ahc_energyPerLayer", &ahc_energyPerLayer, &b_ahc_energyPerLayer);
   fChain2->SetBranchAddress("ahc_nHitsPerLayer", &ahc_nHitsPerLayer, &b_ahc_nHitsPerLayer);
   fChain2->SetBranchAddress("ahc_energyPerLayer_err", &ahc_energyPerLayer_err, &b_ahc_energyPerLayer_err);
   fChain2->SetBranchAddress("ahc_cellSize", &ahc_cellSize, &b_ahc_cellSize);
   fChain2->SetBranchAddress("ahc_cogXPerLayer", &ahc_cogXPerLayer, &b_ahc_cogXPerLayer);
   fChain2->SetBranchAddress("ahc_cogYPerLayer", &ahc_cogYPerLayer, &b_ahc_cogYPerLayer);
   fChain2->SetBranchAddress("ahc_radiusPerLayer", &ahc_radiusPerLayer, &b_ahc_radiusPerLayer);
   fChain2->SetBranchAddress("ahc_radiusEwPerLayer", &ahc_radiusEwPerLayer, &b_ahc_radiusEwPerLayer);
   fChain2->SetBranchAddress("ahc_hitCellID", &ahc_hitCellID, &b_ahc_hitCellID);
   fChain2->SetBranchAddress("ahc_hitI", &ahc_hitI, &b_ahc_hitI);
   fChain2->SetBranchAddress("ahc_hitJ", &ahc_hitJ, &b_ahc_hitJ);
   fChain2->SetBranchAddress("ahc_hitK", &ahc_hitK, &b_ahc_hitK);
   fChain2->SetBranchAddress("ahc_hitEnergy", &ahc_hitEnergy, &b_ahc_hitEnergy);
   fChain2->SetBranchAddress("ahc_hitTime", &ahc_hitTime, &b_ahc_hitTime);
   fChain2->SetBranchAddress("ahc_hitType", &ahc_hitType, &b_ahc_hitType);
   fChain2->SetBranchAddress("ahc_hitRadius", &ahc_hitRadius, &b_ahc_hitRadius);
   fChain2->SetBranchAddress("ahc_hitEnergyDensity", &ahc_hitEnergyDensity, &b_ahc_hitEnergyDensity);
   fChain2->SetBranchAddress("ahc_hitX", &ahc_hitX, &b_ahc_hitX);
   fChain2->SetBranchAddress("ahc_hitY", &ahc_hitY, &b_ahc_hitY);
   fChain2->SetBranchAddress("ahc_hitZ", &ahc_hitZ, &b_ahc_hitZ);
   fChain2->SetBranchAddress("lda_trigTime", &lda_trigTime, &b_lda_trigTime);
   fChain2->SetBranchAddress("event_BXID", &event_BXID, &b_event_BXID);


   Notify();
}

Bool_t RecEvent::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void RecEvent::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain1) return;
   fChain1->Show(entry);
   if (!fChain2) return;
   fChain2->Show(entry);

}
Int_t RecEvent::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void RecEvent::SetNewTree(TTree* tr){
     newTr = tr;
     newTr->Branch("event", &event, "event/i");
     newTr->Branch("trigger_timestamp", &trigger_timestamp, "trigger_timestamp/l");
     newTr->Branch("run", &run, "run/i");
     newTr->Branch("pdgID", &pdgID, "pdgID/I");
     newTr->Branch("beamEnergy", &beamEnergy, "beamEnergy/F");
     newTr->Branch("trueBeamEnergy", &trueBeamEnergy, "trueBeamEnergy/F");
     newTr->Branch("energyLostEE", &energyLostEE, "energyLostEE/F");
     newTr->Branch("energyLostFH", &energyLostFH, "energyLostFH/F");
     newTr->Branch("energyLostBH", &energyLostBH, "energyLostBH/F");
     newTr->Branch("energyLostBeam", &energyLostBeam, "energyLostBeam/F");
     newTr->Branch("energyLostOutside", &energyLostOutside, "energyLostOutside/F");
     newTr->Branch("NRechits", &NRechits, "NRechits/i");
     newTr->Branch("rechit_detid", rechit_detid);
     newTr->Branch("rechit_module", rechit_module);
     newTr->Branch("rechit_layer", rechit_layer);
     newTr->Branch("rechit_chip", rechit_chip);
     newTr->Branch("rechit_channel", rechit_channel);
     newTr->Branch("rechit_type", rechit_type);
     newTr->Branch("rechit_x", rechit_x);
     newTr->Branch("rechit_y", rechit_y);
     newTr->Branch("rechit_z", rechit_z);
     newTr->Branch("rechit_iu", rechit_iu);
     newTr->Branch("rechit_iv", rechit_iv);
     newTr->Branch("rechit_iU", rechit_iU);
     newTr->Branch("rechit_iV", rechit_iV);
     newTr->Branch("rechit_energy", rechit_energy);
     newTr->Branch("rechit_energy_noHG", rechit_energy_noHG);
     newTr->Branch("rechit_amplitudeHigh", rechit_amplitudeHigh);
     newTr->Branch("rechit_amplitudeLow", rechit_amplitudeLow);
     newTr->Branch("rechit_hg_goodFit", rechit_hg_goodFit);
     newTr->Branch("rechit_lg_goodFit", rechit_lg_goodFit);
     newTr->Branch("rechit_hg_saturated", rechit_hg_saturated);
     newTr->Branch("rechit_lg_saturated", rechit_lg_saturated);
     newTr->Branch("rechit_fully_calibrated", rechit_fully_calibrated);
     newTr->Branch("rechit_noise_flag", rechit_noise_flag);
     newTr->Branch("rechit_TS2High", rechit_TS2High);
     newTr->Branch("rechit_TS2Low", rechit_TS2Low);
     newTr->Branch("rechit_TS3High", rechit_TS3High);
     newTr->Branch("rechit_TS3Low", rechit_TS3Low);
     newTr->Branch("rechit_Tot", rechit_Tot);
     newTr->Branch("rechit_toa_calib_flag", rechit_toa_calib_flag);
     newTr->Branch("rechit_toaFall_flag", rechit_toaFall_flag);
     newTr->Branch("rechit_toaRise_flag", rechit_toaRise_flag);
     newTr->Branch("rechit_toaFall_norm", rechit_toaFall_norm);
     newTr->Branch("rechit_toaRise_norm", rechit_toaRise_norm);
     newTr->Branch("rechit_toaFall_time", rechit_toaFall_time);
     newTr->Branch("rechit_toaRise_time", rechit_toaRise_time);
     newTr->Branch("rechit_toaFall_corr_time", rechit_toaFall_corr_time);
     newTr->Branch("rechit_toaRise_corr_time", rechit_toaRise_corr_time);
     newTr->Branch("rechit_calib_time_toaFall", rechit_calib_time_toaFall);
     newTr->Branch("rechit_calib_time_toaRise", rechit_calib_time_toaRise);
     newTr->Branch("rechit_toaRise", rechit_toaRise);
     newTr->Branch("rechit_toaFall", rechit_toaFall);
     newTr->Branch("rechit_timeMaxHG", rechit_timeMaxHG);
     newTr->Branch("rechit_timeMaxLG", rechit_timeMaxLG);

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
     newTr->Branch("lda_trigTime", &lda_trigTime, "lda_trigTime/l");
     newTr->Branch("event_BXID", &event_BXID, "event_BXID/I");
}
void RecEvent::FillNewTree(){
     newTr->Fill();
}

#endif // #ifdef RecEvent_cxx
