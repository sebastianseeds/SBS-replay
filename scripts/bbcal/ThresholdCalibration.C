#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
//#include "gmn_tree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TEventList.h"
#include <iostream>
#include "TF1.h"

//Simple macro to loop on gmn ROOT trees and determine a rough constant of proportionality between trigger sum amplitudes and calibrated preshower and shower energies:
void ThresholdCalibration( const char *rootfilename, const char *outputfilename="BBthresholdcalibtemp.root", double Atrigmin=0.2, double t0trig=176.5,double t0PS=-3.5, double t0SH=-1.0, double tcut=10.0){
  
  TChain *C = new TChain("T");

  C->Add(rootfilename);

  TFile *fout = new TFile(outputfilename,"RECREATE");

  //approximately 5 MeV/mV --> 5 GeV/V
  TH1D *hAtime_trig = new TH1D("hAtime", "BBCAL; trigger ADC time", 500, 0, 500 );
  TH1D *hAtime_PS = new TH1D("hAtime_PS", "BBCAL; Preshower ADC time", 250, -125, 125 );
  TH1D *hAtime_SH = new TH1D("hAtime_SH", "BBCAL; Shower ADC time", 250, -125, 125 );
  TH2D *hAtime_trig_vs_PS = new TH2D("hAtime_trig_vs_PS", "BBCAL; PS ADC time (ns); trigger ADC time (ns)", 250, -125, 125, 500,0,500);
  TH2D *hAtime_trig_vs_SH = new TH2D("hAtime_trig_vs_SH", "BBCAL; SH ADC time (ns); trigger ADC time (ns)", 250, -125, 125, 500,0,500);
  TH2D *hAtime_PS_vs_SH = new TH2D("hAtime_PS_vs_SH", "BBCAL; SH ADC time (ns); PS ADC time (ns)", 250, -125, 125, 250, -125, 125 );
  TH1D *htdiff_trig_PS = new TH1D("htdiff_trig_PS", "BBCAL; #Deltat(trig - PS) (ns);",250,-50,50);
  TH1D *htdiff_trig_SH = new TH1D("htdiff_trig_SH", "BBCAL; #Deltat(trig - SH) (ns);",250,-50,50);
  TH1D *htdiff_PS_SH = new TH1D("htdiff_PS_SH", "BBCAL; #Deltat(PS-SH) (ns);",250,-50,50);
  
  TH1D *hratio = new TH1D("hratio","BBCAL; E_{PS+SH}/Max trigger amplitude (GeV/Volt);", 250, 0.0, 12.5 );

  TH1D *hAtrigmax = new TH1D("hAtrigmax","BBCAL; Max trigger amplitude (V);", 125, Atrigmin, Atrigmin+0.75);
  TH1D *hEtot = new TH1D("hEtot", "BBCAL; E_{PS+SH} (GeV);", 250, 0.0, 5.0 );
  TH2D *hEtot_vs_Atrig = new TH2D("hEtot_vs_Atrig", "BBCAL; Max trigger amplitude (V); E_{PS+SH} (GeV)", 250, Atrigmin, Atrigmin+0.75, 250, 0.0, 5.0 );
  TH2D *hratio_vs_Atrig = new TH2D("hratio_vs_Atrig", "BBCAL; Max trigger amplitude (V); E_{PS+SH}/Max trigger amplitude (GeV/Volt) ", 250, Atrigmin, Atrigmin+0.75, 250, 0.0,12.5 );
  
  long nevent=0;

  Double_t Eps;
  Double_t Esh;
  Double_t Atime_ps;
  Double_t Atime_sh;
  
  Int_t Nhits_bbtrig;
  Double_t Atime_bbtrig[25];
  Double_t Aamp_bbtrig[25];

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.ps.atimeblk",1);
  C->SetBranchStatus("bb.sh.atimeblk",1);
  C->SetBranchStatus("Ndata.bb.bbtrig.a_amp_p",1);
  C->SetBranchStatus("bb.bbtrig.a_amp_p",1);
  C->SetBranchStatus("bb.bbtrig.a_time",1);

  C->SetBranchAddress("bb.ps.e",&Eps);
  C->SetBranchAddress("bb.sh.e",&Esh);
  C->SetBranchAddress("bb.sh.atimeblk",&Atime_sh);
  C->SetBranchAddress("bb.ps.atimeblk",&Atime_ps);
  C->SetBranchAddress("Ndata.bb.bbtrig.a_amp_p",&Nhits_bbtrig);
  C->SetBranchAddress("bb.bbtrig.a_amp_p",Aamp_bbtrig);
  C->SetBranchAddress("bb.bbtrig.a_time",Atime_bbtrig);
  
  while( C->GetEntry( nevent++ ) ){
    if( nevent%1000 == 0 ) cout << nevent << endl;
    //choose trigger sum with maximum amplitude:

    int itrigmax=-1;
    double Atrigmax = 0.0;
    
    for( int itrig=0; itrig<Nhits_bbtrig; itrig++ ){
      if( itrig == 0 || Aamp_bbtrig[itrig] > Atrigmax ){
	itrigmax = itrig;
	Atrigmax = Aamp_bbtrig[itrig];
      }
    }

    if( itrigmax >= 0 && Atrigmax/1000. > Atrigmin ){
      hAtime_trig->Fill( Atime_bbtrig[itrigmax] );
      hAtime_trig_vs_PS->Fill( Atime_ps, Atime_bbtrig[itrigmax] );
      hAtime_trig_vs_SH->Fill( Atime_sh, Atime_bbtrig[itrigmax] );
      hAtime_PS_vs_SH->Fill( Atime_sh, Atime_ps );
      hAtime_PS->Fill( Atime_ps );
      hAtime_SH->Fill( Atime_sh );

      double Etot = Eps + Esh;

      //for amplitude ratios, cut on timing:
      double ttrigtemp = Atime_bbtrig[itrigmax] - t0trig;
      double tPStemp = Atime_ps - t0PS;
      double tSHtemp = Atime_sh - t0SH;

      htdiff_trig_PS->Fill( ttrigtemp-tPStemp );
      htdiff_trig_SH->Fill( ttrigtemp-tSHtemp );
      htdiff_PS_SH->Fill( tPStemp-tSHtemp );

      if( fabs( ttrigtemp - tPStemp ) <= tcut && fabs( ttrigtemp - tSHtemp ) <= tcut &&
	  fabs( tPStemp - tSHtemp ) <= tcut ){
	hEtot_vs_Atrig->Fill( Atrigmax/1000., Etot );
	hEtot->Fill( Etot );
	hratio->Fill( Etot/(Atrigmax/1000.) );

	hratio_vs_Atrig->Fill( Atrigmax/1000.0, Etot/(Atrigmax/1000.) );
	
	hAtrigmax->Fill( Atrigmax/1000.0 );
      }
    }
    
    
  }

  //A model for the trigger amplitude distribution is a rising exponential "leading edge" function convoluted with a falling exponential:

  // N(x) = N0 * (1 - exp(-(x-thr)/width)) * exp(-slope*x) 
  
  TF1 *fitfunc_Atrigmax = new TF1("fitfunc_Atrigmax", "(1.-exp(([1]-x)/[2]))*exp([0]-[3]*x)",Atrigmin, Atrigmin+0.75);

  double threshold_guess = hAtrigmax->GetBinCenter(hAtrigmax->GetMaximumBin());

  int binlow = hAtrigmax->GetMaximumBin();
  double ymax = hAtrigmax->GetBinContent(binlow);
  
  while( hAtrigmax->GetBinContent(binlow--) >= 0.7 * ymax && binlow >= 1 ){};

  double xlow = hAtrigmax->GetBinCenter(binlow);

  cout << "xlow = " << xlow << endl;
  
  double threshold_width = 0.01;
  double eslope = 23.4;
  double N0 = 17.0;

  hAtrigmax->Fit("expo","S","",threshold_guess + 5.*threshold_width, threshold_guess+5.*threshold_width + 0.5);

  double Const = ( (TF1*) (hAtrigmax->GetListOfFunctions()->FindObject("expo") ) )->GetParameter("Constant");
  double Slope = ( (TF1*) (hAtrigmax->GetListOfFunctions()->FindObject("expo") ) )->GetParameter("Slope");
  
  
  double startpar[] = {Const, threshold_guess, threshold_width, -Slope};
  
  fitfunc_Atrigmax->SetParameters( startpar );

  hAtrigmax->Draw();
  fitfunc_Atrigmax->SetNpx(1000);
  fitfunc_Atrigmax->Draw("CSAME");

  
  
  hAtrigmax->Fit(fitfunc_Atrigmax,"S","",xlow,threshold_guess+0.6);

  double effthreshold = fitfunc_Atrigmax->GetParameter(1);
  double dthr = fitfunc_Atrigmax->GetParError(1);

  double width = fitfunc_Atrigmax->GetParameter(2);
  
  //Plot total energy over Atrig ratio for events with trigger amplitude well above calibrated threshold:

  int binx = hratio_vs_Atrig->GetXaxis()->FindBin( effthreshold + 4.5*width );

  //double xlow = hratio_vs_Atrig->GetXaxis()->GetBinCenter( binx );

  int nbinsx = hratio_vs_Atrig->GetNbinsX();
  
  TH1D *hratio_projy = hratio_vs_Atrig->ProjectionY("hratio_Atrig_py",binx,nbinsx);

  double mean = hratio_projy->GetMean();
  double rms = hratio_projy->GetRMS();
  
  hratio_projy->Fit("gaus","S","",mean-2.5*rms, mean+1.5*rms);

  double Ethreshold = ( (TF1*) hratio_projy->GetListOfFunctions()->FindObject("gaus") )->GetParameter("Mean") * effthreshold;

  cout << "Effective energy threshold = " << Ethreshold << " GeV" << endl;
  
  
  fout->Write();
  
  

}
