#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>

#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TChain.h"
#include "TMath.h"
#include "TString.h"

#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/samples.h"

const int nTotBins = 59;
class BaseHistogram
    {
      
    public:
      void BookHistgram(const char *, const int&, const char *);
      TFile *oFile;

      TH1D *h_cutFlow_aux,*h_cutFlow, *h_cutFlow_misc;    
      TH1D *hNtops;
      TH1D *hNbjets;
      TH1D *hNjets;
      TH1D *hMT2;
      TH1D *hMET;
      TH1D *hHT;
      TH1D *hMHT;
      TH2D *h2_nB_nT,*h2_nB_nT_baseline, *h2_nB_nT_loose, *h2_MT2_vs_met,*h2_MT2_vs_met_loose, *h2_MT2_vs_met_baseline;
      TH2D *h2_nB_nJ,*h2_nB_nJ_loose,*h2_nB_nJ_baseline;
      TH2D *h2_nT_nJ,*h2_nT_nJ_loose,*h2_nT_nJ_baseline;
      TH1D *h_mtw_loose, h_mtw_baseline;
      TH1D *h_searchBinYields;
    };

void BaseHistogram::BookHistgram(const char *outFileName, const int& filerun, const char *spec)
    {
      TString filename(outFileName);
      TString index(std::to_string(filerun));
      TString specT(spec);
      filename += "_"+specT+"_"+index+".root";
      oFile = new TFile(filename, "recreate");

      TString sampleNameT(outFileName);
      h_cutFlow = new TH1D(sampleNameT+"_h1_cutFlow", sampleNameT+": cut flow table", 20, 0, 20);
      h_cutFlow_aux = new TH1D(sampleNameT+"_h1_cutFlow_aux", sampleNameT+": more cut flow table", 20, 0, 20);
      h_cutFlow_misc = new TH1D(sampleNameT+"_h1_cutFlow_misc", sampleNameT+": more cut flow table", 20, 0, 20);
      h_cutFlow->SetCanExtend(TH1::kAllAxes); h_cutFlow->Sumw2();
      h_cutFlow_aux->SetCanExtend(TH1::kAllAxes); h_cutFlow_aux->Sumw2();
      h_cutFlow_misc->SetCanExtend(TH1::kAllAxes); h_cutFlow_misc->Sumw2();

      hNtops = new TH1D("hNtops", "No. of top;N_{t};Events", 5, 0, 5);
      hNbjets = new TH1D("hNbjets", "No. of bjets;N_{b};Events", 5, 0, 5);
      hNjets = new TH1D("hNjets", "No. of jets;N_{jets};Events", 14, 0, 14);

      hMT2 = new TH1D("hMT2", "MT2; M_{T2} GeV; Events" , 100, 0., 1000.);
      hMET = new TH1D("hMET", "Missing E_{T}; E_{T}^{miss} GeV; Events" , 100, 0., 1000.);
      hHT = new TH1D("hHT", "HT; HT GeV; Events" , 100, 0., 1000.);
      hMHT = new TH1D("hMHT", "MHT; H_{T}^{miss} GeV; Events" , 100, 0., 1000.);

      h2_nB_nT = new TH2D("h2_nB_nT", "bjets_vs_nTops;N_{b};N_{t}", 5, 0, 5, 5, 0, 5);
      h2_nB_nT->Sumw2();
      hNtops->Sumw2();
      hNbjets->Sumw2();
      hNjets->Sumw2();
      hMT2->Sumw2();
      hMET->Sumw2();
      hHT->Sumw2();
      hMHT->Sumw2();

      h_searchBinYields = new TH1D("h_searchBinYields", "search bin yields", nTotBins, 0, nTotBins);
      h_searchBinYields->Sumw2();

    }


void FillDouble(TH1D *hist, const double &a, const double &w)
    {
      int nbin = hist->GetNbinsX();
      double low = hist->GetBinLowEdge(nbin);
      double high = hist->GetBinLowEdge(nbin + 1);
      double copy = a;
      if(copy >= high) copy = low;
      hist->Fill(copy, w);
    }

void FillInt(TH1D *hist, const int &a, const double &w)
    {
      int nbin = hist->GetNbinsX();
      int low = (int)hist->GetBinLowEdge(nbin);
      int high = (int)hist->GetBinLowEdge(nbin + 1);
      int copy = a;
      if(copy >= high) copy = low;
      hist->Fill(copy, w);
    }

template <typename T>
void Fill2D(TH2 *hist, const T& a, const T& b, const double &w){
  int nbinx = hist->GetNbinsX();
  int nbiny = hist->GetNbinsY();
  T lowx = hist->GetXaxis()->GetBinLowEdge(nbinx);
  T highx = hist->GetXaxis()->GetBinLowEdge(nbinx + 1);
  T lowy = hist->GetYaxis()->GetBinLowEdge(nbiny);
  T highy = hist->GetYaxis()->GetBinLowEdge(nbiny + 1);
  T copyx = a;
  if(copyx >= highx) copyx = lowx;
  T copyy = b;
  if(copyy >= highy) copyy = lowy;
  hist->Fill(copyx, copyy, w);
}
