#include "DataMC.h"

#include <iostream>
#include <algorithm>
#include <cstring>
#include <string>
#include <map>
#include <cmath>
#include <set>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>

#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/searchBins.h"
#include "SusyAnaTools/Tools/customize.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/baselineDef.h"

#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TString.h"
#include "TVector2.h"
#include "TMath.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TChain.h"



const bool isblind = false;
const bool doSingleMuonCS = false;
const bool doInvDphi = false;
const bool usegenmet = false;


const char  *spec1 = "DataMC";
const std::string spec = "MVA";
using namespace std;
// === Main Function ===================================================
int main(int argc, char* argv[]) {
  if (argc < 5)
    {
      std::cerr <<"Please give 5 arguments "<<"SubsampleName"<< " MaxEvent"<<" Startfile"<<" No. of Files to run"<<std::endl;
      std::cerr <<" Valid configurations are " << std::endl;
      std::cerr <<" ./Closure TTbarInc 1000 0 1" << std::endl;
      return -1;
    }

  // Manage arguments
  const char *subSampleName = argv[1];
  const char *Maxevent = argv[2];
  const  char *Stratfile = argv[3];
  const  char *Filerun = argv[4];
  

  const  int startfile = std::atoi(Stratfile);
  const int filerun = std::atoi(Filerun);
  int maxEvent = std::atoi(Maxevent);


  // Prepare file list and finalize it
  TChain *fChain = 0;
  BaseHistogram bh;
  int stfile = int(startfile/filerun);
  bh.BookHistgram(subSampleName, stfile, spec1);
  
  const string condor =  (argc == 6) ? argv[5]: "";
  
  AnaSamples::SampleSet ss = condor.empty()? AnaSamples::SampleSet():AnaSamples::SampleSet(argv[5]);
  AnaSamples::SampleCollection sc(ss);
                                   
  double scaleMC = 1.;                                                                              
  if(ss[subSampleName] != ss.null())                                                                             
    {                                                                                                               
      fChain = new TChain(ss[subSampleName].treePath.c_str());                                                           
      ss[subSampleName].addFilesToChain(fChain, startfile, filerun);

      scaleMC = ss[subSampleName].getWeight();
    } 


  bool isData = false;
  TString sampleString(subSampleName);
  if(sampleString.Contains("Data")){scaleMC = 1.0; isData = true;}

  NTupleReader *tr = 0;
  if( isData )  tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames_DataOnly);
  else {tr = new NTupleReader(fChain, AnaConsts::activatedBranchNames);}
  



  //BaselineVessel blvT3(tr, "Type3");
  //BaselineVessel blvICHEP(tr, "ICHEP");
  BaselineVessel blvMVA(*tr, "MVA");

  //blvT3.SetupTopTagger(false);
  //blvICHEP.SetupTopTagger(true, "Example_Legacy_TopTagger.cfg");
  blvMVA.SetupTopTagger(true,"Example_TopTagger.cfg" );

  //tr.registerFunction(blvT3);
  //tr.registerFunction(blvICHEP);
  tr->registerFunction(blvMVA);

  //Searchbin                                                                                                                                                                      
  SearchBins SB("SB_59_2016");

  // --- Analyse events --------------------------------------------
  std::cout<<"First loop begin: "<<std::endl;
  int entries = tr->getNEntries();
  int entryToProcess = -1;

  std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<scaleMC<<std::endl; 
  cout<<"maxevent: "<<maxEvent<<endl;
  while(tr->getNextEvent())
    {

      if(maxEvent>=0 && tr->getEvtNum() > maxEvent ) break;
      if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (entries>=10 && (tr->getEvtNum()-1)%(entries/10) == 0) ) std::cout<<
							     "\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
      
      // Internal evtWeight in the sample: default is 1.0 execept for MC samples with intrinsic weight, e.g., QCD flat sample.
      double iniWeight = tr->getVar<double>("evtWeight");
      double puWeight = 1.0; // currently set to be 1.0
      
      double stored_weight = sampleString.Contains("Data") ? 1 : tr->getVar<double>("stored_weight");
      int sign_of_stored_weight = (stored_weight > 0) ? 1 : ((stored_weight < 0) ? -1 : 0);
      
      if( sign_of_stored_weight == 1 ) bh.h_cutFlow_aux->Fill("posSign", 1);
      if( sign_of_stored_weight == -1 ) bh.h_cutFlow_aux->Fill("negSign", 1);
      if( sign_of_stored_weight == 0 ) bh.h_cutFlow_aux->Fill("zero", 1);
      
      double evtWeight = iniWeight >=0 ? iniWeight * puWeight * sign_of_stored_weight : iniWeight * puWeight;
      scaleMC = evtWeight*scaleMC;

      // Get branches out directly from what stored in the tree
      const unsigned int & run = tr->getVar<unsigned int>("run"); 
      const unsigned int & lumi = tr->getVar<unsigned int>("lumi"); 
      const unsigned int & event = tr->getVar<unsigned int>("event");
      
      if( isData && isblind && (!doSingleMuonCS || !doInvDphi) ){
	if( run > 274240 ) continue;
      }
      
      const double & genmet_tmp = tr->getVar<double>("genmet");
      const double & genmetphi_tmp = tr->getVar<double>("genmetphi");
      const double genmet = (&genmet_tmp) != nullptr ? tr->getVar<double>("genmet") : 0;
      const double genmetphi = (&genmetphi_tmp) != nullptr ? tr->getVar<double>("genmetphi") : -999;
      TLorentzVector genmetLVec;
      if( (&genmet_tmp) != nullptr ){
	genmetLVec.SetPtEtaPhiM(tr->getVar<double>("genmet"), 0, tr->getVar<double>("genmetphi"), 0);
      }
      
      const double met = ( usegenmet && (&genmet_tmp) != nullptr )? genmet : tr->getVar<double>("met");
      const double metphi = ( usegenmet && (&genmet_tmp) != nullptr )? genmetphi : tr->getVar<double>("metphi");
      TLorentzVector metLVec; metLVec.SetPtEtaPhiM(met, 0, metphi, 0);
      
      const int nbJets = tr->getVar<int>("cntCSVS"+spec), nTops = tr->getVar<int>("nTopCandSortedCnt"+spec );
      
      const double MT2 = tr->getVar<double>("best_had_brJet_MT2"+spec );
      const double HT = tr->getVar<double>("HT"+spec );
      const int nJets = tr->getVar<int>("cntNJetsPt30Eta24"+spec );
      const TLorentzVector mhtLVec = AnaFunctions::calcMHT(tr->getVec<TLorentzVector>("jetsLVec"), AnaConsts::pt30Arr);
      const double MHT  = mhtLVec.Pt();
      
      const std::vector<double> & dPhiVec = tr->getVec<double>("dPhiVec" + spec);
      
      const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
      const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger");
     
      if(isData)
	{
	  bool foundTrigger = false;
	  for(unsigned it=0; it<TriggerNames.size(); it++)
	    {
	      // CS trigger
	      if( sampleString.Contains("SingleMuon") )
		{
		  if( TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_v") != string::npos
		      || TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT400_v")!= string::npos
		      || TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT600_v")!= string::npos 
		      || TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_v")!= string::npos 
		      ||TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT400_PFMET50_v")!= string::npos)
		    {
		      //if( TriggerNames[it].find("HLT_Mu15_IsoVVVL_PFHT350_v") != string::npos){
		      if( PassTrigger[it] ) foundTrigger = true;
		    }
		}

	      // Search Trigger
	      if( sampleString.Contains("HTMHT") )
		{
		  /*
		    if(    TriggerNames[it].find("HLT_PFMET170_NoiseCleaned_v") != std::string::npos
		    || TriggerNames[it].find("HLT_PFMET170_JetIdCleaned_v") != std::string::npos 
		    || TriggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v") != std::string::npos
		    || TriggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v") != std::string::npos
		    || TriggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v") != std::string::npos
		    || TriggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v") != std::string::npos
		    || TriggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v") != std::string::npos
		    || TriggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v") != std::string::npos
		    ){
		  */
		  /* if(    TriggerNames[it].find("HLT_PFHT350_PFMET100_JetIdCleaned_v") != std::string::npos
                     || TriggerNames[it].find("HLT_PFHT350_PFMET100_NoiseCleaned_v") != std::string::npos
                     || TriggerNames[it].find("HLT_PFHT350_PFMET100_v") != std::string::npos ){
	      */
		  if( TriggerNames[it].find("HLT_PFHT300_PFMET100_v") != std::string::npos )
		    {

		      if( PassTrigger[it] ) foundTrigger = true;
		    }
		}
	      
	      
	      
	    }
	  if( !foundTrigger ) continue;
	}

      
      const std::vector<TLorentzVector> & muonsLVec = tr->getVec<TLorentzVector>("muonsLVec");
      const std::vector<double> & muonsRelIso = tr->getVec<double>("muonsRelIso");
      const std::vector<double> & muonsMiniIso = tr->getVec<double>("muonsMiniIso");
      const std::vector<double> & muonsMtw = tr->getVec<double>("muonsMtw");
      const std::vector<int> & muonsFlagMedium = tr->getVec<int>("muonsFlagMedium");
      
      const std::vector<TLorentzVector> & elesLVec = tr->getVec<TLorentzVector>("elesLVec");
      const std::vector<double> & elesRelIso = tr->getVec<double>("elesRelIso");
      const std::vector<double> & elesMiniIso = tr->getVec<double>("elesMiniIso");
      const std::vector<double> & elesMtw = tr->getVec<double>("elesMtw");
      
      const std::vector<int> W_tau_prongsVec = sampleString.Contains("Data")? std::vector<int>() : tr->getVec<int>("W_tau_prongsVec");
      const std::vector<int> W_tau_emuVec = sampleString.Contains("Data")? std::vector<int>() : tr->getVec<int>("W_tau_emuVec");
      const std::vector<int> W_emuVec = sampleString.Contains("Data")? std::vector<int>() : tr->getVec<int>("W_emuVec");
      
      std::vector<int> emuVec_merge;
      emuVec_merge.reserve( W_emuVec.size() + W_tau_emuVec.size() ); 
      emuVec_merge.insert( emuVec_merge.end(), W_emuVec.begin(), W_emuVec.end() );
      emuVec_merge.insert( emuVec_merge.end(), W_tau_emuVec.begin(), W_tau_emuVec.end() );
      
      // Get branches out from further computation in the baselineDef.h -> note that the spec string needs to be added
      const std::vector<TLorentzVector> & jetsLVec_forTagger = tr->getVec<TLorentzVector>("jetsLVec_forTagger" + spec);
      const std::vector<double> & recoJetsBtag_forTagger = tr->getVec<double>("recoJetsBtag_forTagger" + spec);
      
      
      
      const bool passBaselineNoLepVeto = tr->getVar<bool>("passBaselineNoLepVeto" +spec);
      const bool passBaseline = tr->getVar<bool>("passBaseline" + spec);
      
      const int searchBinIdx = SB.find_Binning_Index(nbJets, nTops, MT2, met);
      
      int cnt_eleTop =0, cnt_muTop =0, cnt_taumuTop =0, cnt_taueleTop =0, cnt_tauhadTop =0, cnt_allhadTop =0;
      
      if (!passBaseline) continue;
      
      FillDouble(bh.hMT2, MT2,scaleMC);
      FillDouble(bh.hMET, met, scaleMC);
      FillDouble(bh.hHT, HT, scaleMC);
      FillDouble(bh.hMHT, MHT, scaleMC);
      
      
      FillInt(bh.hNjets,nJets, scaleMC);
      FillInt(bh.hNtops,nTops, scaleMC);
      FillInt(bh.hNbjets,nbJets, scaleMC);

      if(searchBinIdx != -1) bh.h_searchBinYields->Fill(searchBinIdx, scaleMC);
      
      Fill2D(bh.h2_nB_nT, nbJets, nTops, scaleMC);
      //cout << "nbJets  nTops   MT2   HT  :: " << nbJets <<  "\t" <<  nTops << "\t" << MT2 << "\t" << HT << endl;
    }
  
  bh.hMT2->Write();
  bh.hMET->Write();
  bh.hHT->Write();
  bh.hMHT->Write();
  bh.hNbjets->Write();
  bh.hNjets->Write();
  bh.hNtops->Write();
  bh.h2_nB_nT->Write();
  bh.h_searchBinYields->Write();
  bh.oFile->Close();
  
  fChain->Reset();
  return 0;
}
