//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 13 17:33:43 2019 by ROOT version 6.10/05
// from TTree XCET/XCET
// found on file: ntuple_700.root
//////////////////////////////////////////////////////////

#ifndef XCET_h
#define XCET_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class XCET {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree          *newTr;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          event;
   UInt_t          run;
   Short_t         XCET_021507_signal;
   Short_t         XCET_021523_signal;
   Short_t         scintillator_coincidences;
   Short_t         scintillator_vetos;

   // List of branches
   TBranch        *b_event;   //!
   TBranch        *b_run;   //!
   TBranch        *b_XCET_021507_signal;   //!
   TBranch        *b_XCET_021523_signal;   //!
   TBranch        *b_scintillator_coincidences;   //!
   TBranch        *b_scintillator_vetos;   //!

   XCET(TTree *tree=0);
   virtual ~XCET();
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

#ifdef XCET_cxx
XCET::XCET(TTree *tree) : fChain(0), newTr(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   /* if (tree == 0) { */
   /*    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ntuple_700.root"); */
   /*    if (!f || !f->IsOpen()) { */
   /*       f = new TFile("ntuple_700.root"); */
   /*    } */
   /*    TDirectory * dir = (TDirectory*)f->Get("ntuple_700.root:/XCETntupler"); */
   /*    dir->GetObject("XCET",tree); */

   /* } */
   Init(tree);
}

XCET::~XCET()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t XCET::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t XCET::LoadTree(Long64_t entry)
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

void XCET::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("XCET_021507_signal", &XCET_021507_signal, &b_XCET_021507_signal);
   fChain->SetBranchAddress("XCET_021523_signal", &XCET_021523_signal, &b_XCET_021523_signal);
   fChain->SetBranchAddress("scintillator_coincidences", &scintillator_coincidences, &b_scintillator_coincidences);
   fChain->SetBranchAddress("scintillator_vetos", &scintillator_vetos, &b_scintillator_vetos);
   Notify();
}

Bool_t XCET::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void XCET::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t XCET::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
void XCET::SetNewTree(TTree* tr){
     newTr = tr;
     newTr->Branch("event", &event, "event/i");
     newTr->Branch("run", &run, "run/i");
     newTr->Branch("XCET_021507_signal", &XCET_021507_signal, "XCET_021507_signal/S");
     newTr->Branch("XCET_021523_signal", &XCET_021523_signal, "XCET_021523_signal/S");
     newTr->Branch("scintillator_coincidences", &scintillator_coincidences, "scintillator_coincidences/S");
     newTr->Branch("scintillator_vetos", &scintillator_vetos, "scintillator_vetos/S");
}
void XCET::FillNewTree(){
     newTr->Fill();
}
#endif // #ifdef XCET_cxx
