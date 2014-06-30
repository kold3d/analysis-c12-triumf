//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun 27 09:51:28 2014 by ROOT version 5.34/18
// from TTree energyAngle/Events tagged with reconstructed energy and angle.
// found on file: energy_angle.root
//////////////////////////////////////////////////////////

#ifndef Spectra_h
#define Spectra_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
class EnergyLoss;

typedef struct PCBoundEntry {
  Float_t cmEnergy;
  Float_t LeftPCBound;
  Float_t RightPCBound;
} PCBoundEntry;
typedef std::map<int,std::vector<PCBoundEntry> > PCBoundTable;

typedef struct Vect3d {
  Float_t x;
  Float_t y;
  Float_t z;
} Point3d;

class Vec3d {
public:
 Vec3d(Point3d start, Point3d stop) {
    origin.x = start.x;
    origin.y = start.y;
    origin.z = start.z;
    mag = sqrt((start.x-stop.x)*(start.x-stop.x)+
	       (start.y-stop.y)*(start.y-stop.y)+
	       (start.z-stop.z)*(start.z-stop.z));
    x=(start.x-stop.x)/mag;
    y=(start.y-stop.y)/mag;
    z=(start.z-stop.z)/mag;
  };
  Point3d PointIntersectXPlane(Float_t xp) {
   Point3d point = {xp+origin.x,xp/x*y+origin.y,xp/x*z+origin.z};
   return point;
  };
  Point3d PointIntersectYPlane(Float_t yp) {
    Point3d point = {yp/y*x+origin.x,yp+origin.y,yp/y*z+origin.z};
    return point;
  };
  Point3d PointIntersectZPlane(Float_t zp) {
    Point3d point = {zp/z*x+origin.x,zp/z*y+origin.y,zp+origin.z};
    return point;
  };
  Point3d origin;
  Float_t mag;
  Float_t x;
  Float_t y;
  Float_t z;
};

typedef std::pair<Float_t,Float_t> Boundary;

class Spectra {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         measured_energy;
   Int_t           quadrant;
   Int_t           detector;
   Float_t         cm_energy[2];
   Float_t         lab_angle[2];
   Float_t         position[2];
   Float_t         sum_dE[2];
   Int_t           wire[2];

   // List of branches
   TBranch        *b_measured_energy;   //!
   TBranch        *b_quadrant;   //!
   TBranch        *b_detector;   //!
   TBranch        *b_cm_energy;   //!
   TBranch        *b_lab_angle;   //!
   TBranch        *b_position;   //!
   TBranch        *b_sum_dE;   //!
   TBranch        *b_wire;   //!

   Spectra(TTree *tree=0);
   virtual ~Spectra();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   static void      CalcPCBoundTable();
   static void      ReadPCBoundTable();
   static std::pair<Float_t,Float_t> LookupPCBound(Int_t,Float_t);

 private:
   static Float_t pressure;
   static Float_t temperature;
   static EnergyLoss* proton;
   static EnergyLoss* carbon;
   static Float_t density;
   static Float_t m1;
   static Float_t m2;
   static Float_t beam_energy;

   static void InitParameters();
   void DivideTargetThickness(TH1F*);
   void EstimateSolidAngleNorm(TH1F*,Int_t);
   void CalcSolidAngleNorm(TH1F*,Int_t);
   std::pair<Int_t,Float_t> CalcPCCell(Float_t,Float_t,Float_t,Float_t);
   std::pair<Float_t,Float_t> CalcPCBoundary(Int_t,Float_t);

   static PCBoundTable pctable;
};

#endif

#ifdef Spectra_cxx
Spectra::Spectra(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("energy_angle.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("energy_angle.root");
      }
      f->GetObject("energyAngle",tree);

   }
   Init(tree);
   InitParameters();
}

Spectra::~Spectra()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Spectra::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Spectra::LoadTree(Long64_t entry)
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

void Spectra::Init(TTree *tree)
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

   fChain->SetBranchAddress("measured_energy", &measured_energy, &b_measured_energy);
   fChain->SetBranchAddress("quadrant", &quadrant, &b_quadrant);
   fChain->SetBranchAddress("detector", &detector, &b_detector);
   fChain->SetBranchAddress("cm_energy", cm_energy, &b_cm_energy);
   fChain->SetBranchAddress("lab_angle", lab_angle, &b_lab_angle);
   fChain->SetBranchAddress("position", position, &b_position);
   fChain->SetBranchAddress("sum_dE", sum_dE, &b_sum_dE);
   fChain->SetBranchAddress("wire", wire, &b_wire);
   Notify();
}

Bool_t Spectra::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Spectra::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Spectra::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Spectra_cxx
