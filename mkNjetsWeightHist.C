#include <iostream>
#include <string>
#include <cstring>
#include <sstream>

#include "SusyAnaTools/Tools/samples.h"
#include "SusyAnaTools/Tools/NTupleReader.h"
#include "SusyAnaTools/Tools/baselineDef.h"
#include "SusyAnaTools/Tools/customize.h"

#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TChain.h"

void FillInt(TH1D *hist, const int &a, const double &w);
const std::string spec = "njets_reweigh";
const int njetBins = 7;
const double  jetBins[] = {4, 5, 6, 7, 8, 9, 10, 20};
// === Main Function ===================================================
int main(int argc, char* argv[]) 
   {
     if (argc < 5)
       {
	 std::cerr <<"Please give 5 arguments "<<"SubsampleName"<< " MaxEvent"<<" No. of Files to run"<<" StartFile"<<std::endl;
	 std::cerr <<" Valid configurations are " << std::endl;
	 std::cerr <<" ./Closure TTbarInc 1000 0 1" << std::endl;
	 return -1;
       }


     // Manage arguments
     const char *subSampleName = argv[1];
     const char *Maxevent = argv[2];
     const  char *Stratfile = argv[4];
     const  char *Filerun = argv[3];
     
     
     const  int startfile = std::atoi(Stratfile);
     const int filerun = std::atoi(Filerun);
     int maxEvent = std::atoi(Maxevent);
     
     
     // Prepare file list and finalize it
     TChain *fChain = 0;
     int stFile = int(startfile/filerun);
     std::stringstream ssss;
     ssss << stFile;
     std::string str = ssss.str();

     const std::string condor =  (argc == 6) ? argv[5]: "";
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

     BaselineVessel blvMVA(*tr, spec);
     blvMVA.SetupTopTagger(true,"Example_TopTagger.cfg");
     tr->registerFunction(blvMVA);  

   
     // --- Analyse events --------------------------------------------
     std::cout<<"First loop begin: "<<std::endl;
     int entries = tr->getNEntries();
     std::cout<<"\nentries : "<<entries<<"\t MC Scale: "<<scaleMC<<std::endl; 
     std::cout<<"maxevent: "<<maxEvent<<std::endl;


     //Create Histograms and files.
     TH1D* h1_njets = NULL;
     TFile* infile = new TFile(sampleString+ "_nJets_Reweigh_"+str+".root", "RECREATE");
     if(isData)
       {
	 h1_njets = new TH1D("h1_nJets_Data", "Number of Jets Data;N_{jets};Events",njetBins, jetBins );
       }
     else
       {
	 h1_njets = new TH1D("h1_nJets_MC", "Number of Jets Data;N_{jets};Events",njetBins, jetBins );
       }


     while(tr->getNextEvent())
       {
	 
	 if(maxEvent>=0 && tr->getEvtNum() > maxEvent ) break;
	 if( tr->getEvtNum()-1 == 0 || tr->getEvtNum() == entries || (entries>=10 && (tr->getEvtNum()-1)%(entries/10) == 0) )
	   {
	     std::cout<<"\n   Processing the "<<tr->getEvtNum()-1<<"th event ..."<<std::endl;
	   }

	 const bool passBaseline = tr->getVar<bool>("passBaseline" + spec);
	 if(!passBaseline) continue;

	 // Internal evtWeight in the sample: default is 1.0 execept for MC samples with intrinsic weight, e.g., QCD flat sample.
	 double iniWeight = tr->getVar<double>("evtWeight");
	 double puWeight = 1.0; // currently set to be 1.0
	 double stored_weight = sampleString.Contains("Data") ? 1 : tr->getVar<double>("stored_weight");
	 int sign_of_stored_weight = (stored_weight > 0) ? 1 : ((stored_weight < 0) ? -1 : 0);
	 double evtWeight = iniWeight >=0 ? iniWeight * puWeight * sign_of_stored_weight : iniWeight * puWeight;
	 scaleMC = evtWeight*scaleMC;
	 const unsigned int & run = tr->getVar<unsigned int>("run");
	 if(isData)
	   {
	     if ( run <= 278820 || run >= 279931) continue;
	   } 

	 const std::vector<std::string> & TriggerNames = tr->getVec<std::string>("TriggerNames");
	 const std::vector<int> & PassTrigger = tr->getVec<int>("PassTrigger"); 
	 if(isData)
	   {
	     bool foundTrigger = false;
	     for(unsigned it=0; it<TriggerNames.size(); it++)
	       {
		 // Search Trigger
		 if( sampleString.Contains("MET") )
		   {
		     if(TriggerNames[it].find("HLT_PFMET170_NoiseCleaned_v") != std::string::npos
			|| TriggerNames[it].find("HLT_PFMET170_JetIdCleaned_v") != std::string::npos 
			|| TriggerNames[it].find("HLT_PFMET100_PFMHT100_IDTight_v") != std::string::npos
			|| TriggerNames[it].find("HLT_PFMET110_PFMHT110_IDTight_v") != std::string::npos
			|| TriggerNames[it].find("HLT_PFMET120_PFMHT120_IDTight_v") != std::string::npos
			|| TriggerNames[it].find("HLT_PFMET130_PFMHT130_IDTight_v") != std::string::npos
			|| TriggerNames[it].find("HLT_PFMET140_PFMHT140_IDTight_v") != std::string::npos
			|| TriggerNames[it].find("HLT_PFMET150_PFMHT150_IDTight_v") != std::string::npos
			)
		       {
			 if( PassTrigger[it] ) foundTrigger = true;                                                                           
		       }                                           
		   }
	       }
	     if( !foundTrigger ) continue;
	   }



	 const int nJets = tr->getVar<int>("cntNJetsPt30Eta24"+spec );	 

	 FillInt(h1_njets, nJets, scaleMC);

       }
     h1_njets->Write();
     infile->Close();
     fChain->Reset();
     return 0;
			  

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
