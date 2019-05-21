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
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *newTr;
// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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

   // List of branches
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

   /* RecEvent(char* fname, TChain *tree=0); */
   RecEvent(TTree *tree=0);
   virtual ~RecEvent();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Int_t    GetEntries();
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual void     SetNewTree(TTree* tr=0);
   virtual void     FillNewTree();
};

#endif

#ifdef RecEvent_cxx
/* RecEvent::RecEvent(char* fname, TChain *tree) : fChain(0)  */
RecEvent::RecEvent(TTree *tree) : fChain(0), newTr(0)
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
   Init(tree);
}

RecEvent::~RecEvent()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t RecEvent::GetEntries(){
     return fChain->GetEntries();
}

Int_t RecEvent::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t RecEvent::LoadTree(Long64_t entry)
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

void RecEvent::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
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
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("trigger_timestamp", &trigger_timestamp, &b_trigger_timestamp);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("pdgID", &pdgID, &b_pdgID);
   fChain->SetBranchAddress("beamEnergy", &beamEnergy, &b_beamEnergy);
   fChain->SetBranchAddress("trueBeamEnergy", &trueBeamEnergy, &b_trueBeamEnergy);
   fChain->SetBranchAddress("energyLostEE", &energyLostEE, &b_energyLostEE);
   fChain->SetBranchAddress("energyLostFH", &energyLostFH, &b_energyLostFH);
   fChain->SetBranchAddress("energyLostBH", &energyLostBH, &b_energyLostBH);
   fChain->SetBranchAddress("energyLostBeam", &energyLostBeam, &b_energyLostBeam);
   fChain->SetBranchAddress("energyLostOutside", &energyLostOutside, &b_energyLostOutside);
   fChain->SetBranchAddress("NRechits", &NRechits, &b_NRechits);
   fChain->SetBranchAddress("rechit_detid", &rechit_detid, &b_rechit_detid);
   fChain->SetBranchAddress("rechit_module", &rechit_module, &b_rechit_module);
   fChain->SetBranchAddress("rechit_layer", &rechit_layer, &b_rechit_layer);
   fChain->SetBranchAddress("rechit_chip", &rechit_chip, &b_rechit_chip);
   fChain->SetBranchAddress("rechit_channel", &rechit_channel, &b_rechit_channel);
   fChain->SetBranchAddress("rechit_type", &rechit_type, &b_rechit_type);
   fChain->SetBranchAddress("rechit_x", &rechit_x, &b_rechit_x);
   fChain->SetBranchAddress("rechit_y", &rechit_y, &b_rechit_y);
   fChain->SetBranchAddress("rechit_z", &rechit_z, &b_rechit_z);
   fChain->SetBranchAddress("rechit_iu", &rechit_iu, &b_rechit_iu);
   fChain->SetBranchAddress("rechit_iv", &rechit_iv, &b_rechit_iv);
   fChain->SetBranchAddress("rechit_iU", &rechit_iU, &b_rechit_iU);
   fChain->SetBranchAddress("rechit_iV", &rechit_iV, &b_rechit_iV);
   fChain->SetBranchAddress("rechit_energy", &rechit_energy, &b_rechit_energy);
   fChain->SetBranchAddress("rechit_energy_noHG", &rechit_energy_noHG, &b_rechit_energy_noHG);
   fChain->SetBranchAddress("rechit_amplitudeHigh", &rechit_amplitudeHigh, &b_rechit_amplitudeHigh);
   fChain->SetBranchAddress("rechit_amplitudeLow", &rechit_amplitudeLow, &b_rechit_amplitudeLow);
   fChain->SetBranchAddress("rechit_hg_goodFit", &rechit_hg_goodFit, &b_rechit_hg_goodFit);
   fChain->SetBranchAddress("rechit_lg_goodFit", &rechit_lg_goodFit, &b_rechit_lg_goodFit);
   fChain->SetBranchAddress("rechit_hg_saturated", &rechit_hg_saturated, &b_rechit_hg_saturated);
   fChain->SetBranchAddress("rechit_lg_saturated", &rechit_lg_saturated, &b_rechit_lg_saturated);
   fChain->SetBranchAddress("rechit_fully_calibrated", &rechit_fully_calibrated, &b_rechit_fully_calibrated);
   fChain->SetBranchAddress("rechit_noise_flag", &rechit_noise_flag, &b_rechit_noise_flag);
   fChain->SetBranchAddress("rechit_TS2High", &rechit_TS2High, &b_rechit_TS2High);
   fChain->SetBranchAddress("rechit_TS2Low", &rechit_TS2Low, &b_rechit_TS2Low);
   fChain->SetBranchAddress("rechit_TS3High", &rechit_TS3High, &b_rechit_TS3High);
   fChain->SetBranchAddress("rechit_TS3Low", &rechit_TS3Low, &b_rechit_TS3Low);
   fChain->SetBranchAddress("rechit_Tot", &rechit_Tot, &b_rechit_Tot);
   fChain->SetBranchAddress("rechit_toa_calib_flag", &rechit_toa_calib_flag, &b_rechit_toa_calib_flag);
   fChain->SetBranchAddress("rechit_toaFall_flag", &rechit_toaFall_flag, &b_rechit_toaFall_flag);
   fChain->SetBranchAddress("rechit_toaRise_flag", &rechit_toaRise_flag, &b_rechit_toaRise_flag);
   fChain->SetBranchAddress("rechit_toaFall_norm", &rechit_toaFall_norm, &b_rechit_toaFall_norm);
   fChain->SetBranchAddress("rechit_toaRise_norm", &rechit_toaRise_norm, &b_rechit_toaRise_norm);
   fChain->SetBranchAddress("rechit_toaFall_time", &rechit_toaFall_time, &b_rechit_toaFall_time);
   fChain->SetBranchAddress("rechit_toaRise_time", &rechit_toaRise_time, &b_rechit_toaRise_time);
   fChain->SetBranchAddress("rechit_toaFall_corr_time", &rechit_toaFall_corr_time, &b_rechit_toaFall_corr_time);
   fChain->SetBranchAddress("rechit_toaRise_corr_time", &rechit_toaRise_corr_time, &b_rechit_toaRise_corr_time);
   fChain->SetBranchAddress("rechit_calib_time_toaFall", &rechit_calib_time_toaFall, &b_rechit_calib_time_toaFall);
   fChain->SetBranchAddress("rechit_calib_time_toaRise", &rechit_calib_time_toaRise, &b_rechit_calib_time_toaRise);
   fChain->SetBranchAddress("rechit_toaRise", &rechit_toaRise, &b_rechit_toaRise);
   fChain->SetBranchAddress("rechit_toaFall", &rechit_toaFall, &b_rechit_toaFall);
   fChain->SetBranchAddress("rechit_timeMaxHG", &rechit_timeMaxHG, &b_rechit_timeMaxHG);
   fChain->SetBranchAddress("rechit_timeMaxLG", &rechit_timeMaxLG, &b_rechit_timeMaxLG);
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
   if (!fChain) return;
   fChain->Show(entry);
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
}
void RecEvent::FillNewTree(){
     newTr->Fill();
}

#endif // #ifdef RecEvent_cxx
