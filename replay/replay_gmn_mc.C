#include <iostream>
#include <unordered_map> 

#include "TSystem.h"
#include "THaGlobals.h"
#include "TString.h"
#include "TFile.h"
#include "TList.h"
#include "TObject.h"

#include "THaEvData.h"
#include "THaRun.h"
#include "THaAnalyzer.h"
#include "THaVarList.h"
#include "THaInterface.h"

#include "SBSBigBite.h"
#include "SBSBBShower.h"
#include "SBSBBTotalShower.h"
#include "SBSGRINCH.h"
#include "SBSEArm.h"
#include "SBSHCal.h"
#include "SBSGEMSpectrometerTracker.h"
#include "SBSTimingHodoscope.h"

#include "THaGoldenTrack.h"
#include "THaPrimaryKine.h"
#include "THaRunParameters.h"

#include "SBSSimDecoder.h"

TDatime get_datime(uint sbsconfig)
/* Returns TDatime for a given SBS configuration */
{
  std::unordered_map<uint,TDatime> m = {{4,  "2021-10-21 07:09:00"},
					{7,  "2021-11-13 22:06:00"},
					{11, "2021-11-25 12:17:00"},
					{14, "2022-01-12 11:55:00"},
					{8,  "2022-01-22 10:16:00"},
					{9,  "2022-02-02 01:55:00"}};
  if (m.find(sbsconfig)==m.end()) 
    throw std::invalid_argument("Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  return m[sbsconfig];
}

void replay_gmn_mc(const char* filebase, uint sbsconfig, uint nev = -1, TString experiment="gmn")
{
  SBSBigBite* bigbite = new SBSBigBite("bb", "BigBite spectrometer" );
  //bigbite->AddDetector( new SBSBBShower("ps", "BigBite preshower") );
  //bigbite->AddDetector( new SBSBBShower("sh", "BigBite shower") );
  bigbite->AddDetector( new SBSBBTotalShower("ts", "sh", "ps", "BigBite shower") );
  bigbite->AddDetector( new SBSGRINCH("grinch", "GRINCH PID") );
  //bigbite->AddDetector( new SBSGenericDetector("grinch", "GRINCH PID") );
  bigbite->AddDetector( new SBSTimingHodoscope("hodo", "timing hodo") );
  //bigbite->AddDetector( new SBSGEMSpectrometerTracker("gem", "GEM tracker") );
  bigbite->AddDetector( new SBSGEMSpectrometerTracker("gem", "GEM tracker") );
  gHaApps->Add(bigbite);
  
  SBSEArm *harm = new SBSEArm("sbs","Hadron Arm with HCal");
  harm->AddDetector( new SBSHCal("hcal","HCAL") );
  gHaApps->Add(harm);

  // bigbite->SetDebug(5);
  //harm->SetDebug(2);

  THaAnalyzer* analyzer = new THaAnalyzer;
  
  THaInterface::SetDecoder( SBSSimDecoder::Class() );
  
  TString run_file = Form("%s.root", filebase);
  if(std::getenv("DATA_DIR")){
    run_file = Form("%s/%s.root", std::string(std::getenv("DATA_DIR")).c_str(), filebase);
  }
  
  if( gSystem->AccessPathName(run_file) ) {
    Error("replay.C", "Input file does not exist: %s", run_file.Data() );
    exit(-1);
  }
  
  THaRunBase *run = new SBSSimFile(run_file.Data(), "gmn", "");
  run->SetFirstEvent(0);

  //cout << "Number of events to replay (-1=all)? ";
  //if( nev > 0 )
  //run->SetFirstEvent(110);
  run->SetLastEvent(nev);
  
  run->SetDataRequired(0);

  run->SetDate(get_datime(sbsconfig));
  //run->SetDate(TDatime());
 

  gHaPhysics->Add( new THaGoldenTrack( "BB.gold", "BigBite golden track", "bb" ));
  gHaPhysics->Add( new THaPrimaryKine( "e.kine", "electron kinematics", "bb", 0.0, 0.938272 ));
  
  TString out_dir = gSystem->Getenv("OUT_DIR");
  if( out_dir.IsNull() )
    out_dir = ".";
  TString out_file = out_dir + "/" + Form("replayed_%s.root", filebase);
  
  analyzer->SetOutFile( out_file.Data() );
  cout << "output file " << out_file.Data() << " set up " << endl; 
  // File to record cuts accounting information
  analyzer->SetSummaryFile(out_file.ReplaceAll(".root",".log")); // optional

  // Change the cratemap to point to the sim one
  analyzer->SetCrateMapFileName("sbssim_cratemap");

  cout << "sim crate map setup " << endl;

  TString prefix = gSystem->Getenv("SBS_REPLAY");
  prefix += "/replay/";
  
  TString odef_filename = "replay_gmn_mc.odef";
  odef_filename.Prepend( prefix );
  analyzer->SetOdefFile( odef_filename );
  
  //added cut list in order to have 
  TString cdef_filename = "replay_gmn_mc.cdef";
  cdef_filename.Prepend( prefix );
  analyzer->SetCutFile( cdef_filename );

  //analyzer->SetCutFile( "replay_gmn.cdef" );
  //analyzer->SetOdefFile( "replay_gmn_mc.odef" );
  
  cout << "cut file and out file processed " << endl;
  
  analyzer->SetVerbosity(2);  // write cut summary to stdout
  analyzer->EnableBenchmarks();

  run->Init();
  run->Print();

  run->GetParameters()->Print();
  
  cout << "about to process " << endl;
  analyzer->Process(run);

  run->GetParameters()->Print();
  
  // Clean up

  analyzer->Close();
  delete analyzer;
  //gHaCuts->Clear();
  gHaVars->Clear();
  gHaPhysics->Delete();
  gHaApps->Delete();
 
  
}

int main(int argc, char *argv[])
{
  new THaInterface( "The Hall A analyzer", &argc, argv, 0, 0, 1 );
  //const 
  string filebase; 
  uint nev = -1;
  if(argc<2 || argc>3){
    cout << "Usage: replay_gmn filebase(char*) nev(uint)" << endl;
    return -1;
  }
  filebase = argv[1];
  if(argc==3) nev = atoi(argv[2]);

  replay_gmn_mc(filebase.c_str(), nev);
  return 0;
}
