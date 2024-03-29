//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 14 20:50:38 2014 by ROOT version 5.34/18
// from TTree rawData/Raw Data Tree
// found on file: carbon_triumf_09-13_t.root
//////////////////////////////////////////////////////////

#ifndef PosCal_h
#define PosCal_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class PosCal {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UChar_t         si_mul;
   UChar_t         si_det[12];   //[si_mul]
   UChar_t         si_quad[12];   //[si_mul]
   Int_t           si_ch_e[12];   //[si_mul]
   Float_t         si_cal_e[12];   //[si_mul]
   Int_t           si_ch_t[12];   //[si_mul]
   UChar_t         pc_mul;
   UChar_t         pc_wire[8];   //[pc_mul]
   Int_t           pc_ch_right_e[8];   //[pc_mul]
   Int_t           pc_ch_left_e[8];   //[pc_mul]
   Int_t           pc_ch_right_t[8];   //[pc_mul]
   Int_t           pc_ch_left_t[8];   //[pc_mul]
   Int_t           ic_ch;

   // List of branches
   TBranch        *b_si_mul;   //!
   TBranch        *b_si_det;   //!
   TBranch        *b_si_quad;   //!
   TBranch        *b_si_ch_e;   //!
   TBranch        *b_si_cal_e;   //!
   TBranch        *b_si_ch_t;   //!
   TBranch        *b_pc_mul;   //!
   TBranch        *b_pc_wire;   //!
   TBranch        *b_pc_ch_right_e;   //!
   TBranch        *b_pc_ch_left_e;   //!
   TBranch        *b_pc_ch_right_t;   //!
   TBranch        *b_pc_ch_left_t;   //!
   TBranch        *b_ic_ch;   //!

   PosCal(TTree *tree=0);
   virtual ~PosCal();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

};

#endif

#ifdef PosCal_cxx
PosCal::PosCal(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("carbon_triumf_09-13_t.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("carbon_triumf_09-13_t.root");
      }
      f->GetObject("rawData",tree);

   }
   Init(tree);
}

PosCal::~PosCal()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t PosCal::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t PosCal::LoadTree(Long64_t entry)
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

void PosCal::Init(TTree *tree)
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

   fChain->SetBranchAddress("si_mul", &si_mul, &b_si_mul);
   fChain->SetBranchAddress("si_det", si_det, &b_si_det);
   fChain->SetBranchAddress("si_quad", si_quad, &b_si_quad);
   fChain->SetBranchAddress("si_ch_e", si_ch_e, &b_si_ch_e);
   fChain->SetBranchAddress("si_cal_e", si_cal_e, &b_si_cal_e);
   fChain->SetBranchAddress("si_ch_t", si_ch_t, &b_si_ch_t);
   fChain->SetBranchAddress("pc_mul", &pc_mul, &b_pc_mul);
   fChain->SetBranchAddress("pc_wire", pc_wire, &b_pc_wire);
   fChain->SetBranchAddress("pc_ch_right_e", pc_ch_right_e, &b_pc_ch_right_e);
   fChain->SetBranchAddress("pc_ch_left_e", pc_ch_left_e, &b_pc_ch_left_e);
   fChain->SetBranchAddress("pc_ch_right_t", pc_ch_right_t, &b_pc_ch_right_t);
   fChain->SetBranchAddress("pc_ch_left_t", pc_ch_left_t, &b_pc_ch_left_t);
   fChain->SetBranchAddress("ic_ch", &ic_ch, &b_ic_ch);
   Notify();
}

Bool_t PosCal::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void PosCal::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t PosCal::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef PosCal_cxx
