#define RecEvent_cxx
#include "RecEvent.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void RecEvent::Loop()
{
//   In a ROOT session, you can do:
//      root> .L RecEvent.C
//      root> RecEvent t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain1 == 0) return;
   if (fChain2 == 0) return;

   Long64_t nentries1 = fChain1->GetEntriesFast();
   Long64_t nentries2 = fChain2->GetEntriesFast();

   Long64_t nentries = (nentries1>nentries2)? nentries2:nentries1;

   Long64_t nbytes1 = 0, nb1 = 0;
   Long64_t nbytes2 = 0, nb2 = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb1 = fChain1->GetEntry(jentry);   nbytes1 += nb1;
      nb2 = fChain2->GetEntry(jentry);   nbytes2 += nb2;
      // if (Cut(ientry) < 0) continue;
   }
}
