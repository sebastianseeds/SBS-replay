#include "TH1D.h"
#include "TF1.h"
#include "TChain.h"
#include "TTree.h"
#include "TCut.h"
#include "TCanvas.h"

TH1D *hW2elastic;

double fitfunc( double *x, double *par ){
  double W2 = x[0];

  double Norm = par[0];
  // Let's use (up to) 6th order polynomial for the background:
  double bg = 0.0;
  for( int i = 0; i<7; i++ ){
    bg += par[i+1]*pow(W2,i);
  }

  double signal = Norm * hW2elastic->Interpolate(W2);

  return signal + bg; 
  
}

void Fit_W2_inclusive( TH1D *hW2all, TH1D *hW2elastic_temp, int order_bg=4, double W2min=-1.0, double W2max=2.0 ){
  TF1 *FitFunc = new TF1( "FitFunc", fitfunc, -1, 3, 8 );

  FitFunc->SetNpx(1000);
  
  if( order_bg > 6 ) order_bg = 6;

  hW2elastic = new TH1D( *hW2elastic_temp );
  
  if( hW2all != nullptr && hW2elastic != nullptr ){
    double startpar[] = { 1.0, 0, 0, 0, 0, 0, 0, 0 };

    FitFunc->SetParameters(startpar);

   
    
    int ipar=0;
    while( ipar <= order_bg ){
      ipar++;
    }

    while( ipar <= 6 ){
      FitFunc->FixParameter( ipar+1, 0.0 );
      ipar++;
    }

    hW2all->Fit( FitFunc, "", "", W2min, W2max );

    hW2all->SetLineColor(1);
    hW2all->SetLineWidth(2);
    
    hW2all->DrawCopy();
    //hW2elastic->Scale( FitFunc->GetParameter(0) );

    TH1D *hW2elastic_copy = new TH1D(*hW2elastic); 
    //    hW2elastic->Draw("SAME");

    hW2elastic_copy->SetName("hW2elastic_copy");

    hW2elastic_copy->Scale( FitFunc->GetParameter(0) );
    
    TH1D *hW2_bg = new TH1D( *hW2all );

    hW2_bg->SetName("hW2_bg");

    hW2_bg->Add( hW2all, hW2elastic_copy, 1., -1. );

    hW2_bg->SetMarkerStyle(21);
    hW2_bg->SetMarkerColor(4);
    hW2_bg->SetLineColor(4);
    
    hW2_bg->DrawCopy("SAME");

    hW2elastic_copy->SetMarkerStyle(20);
    hW2elastic_copy->SetMarkerColor(2);
    hW2elastic_copy->SetLineColor(2);
    
    hW2elastic_copy->DrawCopy("same");

    *hW2elastic_temp = *hW2elastic_copy;
    
  }


}

void GetElasticCounts( const char *rootfilename, double W2min=0.4, double W2max = 1.4, TCut cut="", TCut fiducialcut="", double W2min_fit=0, double W2max_fit=1.4 ){
  TChain *C = new TChain("Tout");

  C->Add(rootfilename);

  TH1D *hW2_nocut = new TH1D("hW2_nocut", ";W^{2} (GeV^{2});", 100, -1.0,3.0 );
  TH1D *hW2_HCALcut = new TH1D("hW2_HCALcut", ";W^{2} (GeV^{2});", 100, -1.0, 3.0 );
  TH1D *hW2_anticut = new TH1D("hW2_anticut", ";W^{2} (GeV^{2});", 100, -1.0, 3.0 );

  C->Project("hW2_nocut", "W2", fiducialcut);
  C->Project("hW2_HCALcut", "W2", cut&&fiducialcut);
  C->Project("hW2_anticut", "W2", !cut&&fiducialcut);

  TCanvas *c1 = new TCanvas("c1","c1",1200,750);
  c1->Divide(2,1,.001,.001);

  c1->cd(1);
  
  Fit_W2_inclusive( hW2_nocut, hW2_HCALcut, 6, W2min_fit, W2max_fit );

  double Counts,dCounts;

  Counts = hW2_HCALcut->IntegralAndError( hW2_HCALcut->FindBin(W2min), hW2_HCALcut->FindBin(W2max), dCounts );

  cout << "Elastic yield (total with fiducial cut) = " << Counts << " +/- " << dCounts << endl;

  double CountsAcut,dCountsAcut;

  c1->cd(2);
  
  Fit_W2_inclusive( hW2_anticut, hW2_HCALcut, 6, W2min_fit, W2max_fit );

  CountsAcut = hW2_HCALcut->IntegralAndError( hW2_HCALcut->FindBin(W2min), hW2_HCALcut->FindBin(W2max), dCountsAcut );

  cout << "Elastic yield (passed fiducial cut, failed HCAL cut) = " << CountsAcut << " +/- " << dCountsAcut << endl;

  double efficiency = 1.-CountsAcut/Counts;
  double defficiency = sqrt(efficiency*(1.-efficiency)/Counts);
  
  cout << "HCAL effective proton efficiency = " << efficiency << " +/- " << defficiency << endl;
}
