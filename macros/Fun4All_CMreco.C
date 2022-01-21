#include <GlobalVariables.C>

#include <TPCGemGainCalb.h>

#include <G4Setup_sPHENIX.C>
#include <G4_Bbc.C>
#include <G4_CaloTrigger.C>
#include <G4_DSTReader.C>
#include <G4_Global.C>
#include <G4_HIJetReco.C>
#include <G4_Input.C>
#include <G4_Jets.C>
#include <G4_ParticleFlow.C>
#include <G4_Production.C>
#include <G4_TopoClusterReco.C>
#include <G4_User.C>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <G4_Tracking.C>


#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)

int Fun4All_CMreco(
		   	  const int process = 1,
			  const int nEvents  = 1,
			  const string &embed_input_str1 = "DST_TRACKS_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000001-0",
			  const string &embed_input_str2 = "DST_TRKR_CLUSTER_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000001-0",
			  const string &embed_input_str3 = "DST_TRUTH_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000001-0",
			  const string &embed_input_str4 = "DST_TRKR_G4HIT_sHijing_0_20fm_50kHz_bkg_0_20fm-0000000001-0",
			  const string &outputFile_str = "CM_HitsDistTree", //changed
			  const int skip = 0,
			  const string &outdir = ".")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(1);

  // construct the input and output file names
  char num_field[500];
  sprintf(num_field,"%04d.root", process);
  string numin = num_field;
  string embed_infile1 = embed_input_str1+numin;
  string embed_infile2 = embed_input_str2+numin;
  string embed_infile3 = embed_input_str3+numin;
  string embed_infile4 = embed_input_str4+numin;
  string outputFile = outputFile_str+numin;

  cout << "input file1: " << embed_infile1 << endl;
  cout << "outputFile: " << outputFile << endl;

 // just if we set some flags somewhere in this macro
  recoConsts *rc = recoConsts::instance();


  Enable::TRACKING_EVAL = true;



  
  Input::READHITS = true;
  Input::VERBOSITY = 10;
  INPUTREADHITS::filename[0] = embed_infile1;
  INPUTREADHITS::filename[1] = embed_infile2;
  INPUTREADHITS::filename[2] = embed_infile3;
  INPUTREADHITS::filename[3] = embed_infile4;

  InputInit();

  // Reads in the file
  InputManagers();

  TPCGemGainCalb *track_readback = new TPCGemGainCalb();
  track_readback->Verbosity(0);
  track_readback->set_process(process);
  se->registerSubsystem(track_readback);


  string outfile = "mDistHitsTree_Event";
  outfile += std::to_string(process);
  outfile += ".root";

  if(Enable::TRACKING_EVAL)
    Tracking_Eval(outfile);

  //se->skip(skip);
  se->run(nEvents);
  se->PrintTimer();
  //-----
  // Exit
  //-----

  se->End();
  std::cout << "All done" << std::endl;
  delete se;

return 0;
}
