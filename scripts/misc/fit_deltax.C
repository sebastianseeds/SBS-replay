#include "TH1F.h"
#include "TF1.h"
#include "TAxis.h"
#include "TStyle.h"

void fit_deltax(TH1F *hdx, double Np, double Nn, double mup=-0.8, double mun=0.0, double sigmap=0.1, double sigman=0.1){

  gStyle->SetOptFit();
  
  double xmin = hdx->GetXaxis()->GetXmin();
  double xmax = hdx->GetXaxis()->GetXmin();

  TF1 *fitfunc = new TF1("fitfunc", "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-0.5*pow((x-[4])/[5],2))+[6]+[7]*x+[8]*x*x+[9]*x*x*x+[10]*pow(x,4)",-3,2);

  fitfunc->SetNpx(1000);
  
  double startpar[11] = { Np, mup, sigmap, Nn, mun, sigman, 0, 0, 0, 0, 0 };

  fitfunc->SetParameters(startpar);

  if( Nn <= 0.1 ){
    fitfunc->FixParameter(3,0.0);
    fitfunc->FixParameter(4,0.0);
    fitfunc->FixParameter(5,1.0);
  }
  
  hdx->Fit(fitfunc,"S");

  TF1 *ngaus = new TF1( "ngaus", "[0]*exp(-0.5*pow((x-[1])/[2],2))",-3,2);
  TF1 *pgaus = new TF1( "pgaus", "[0]*exp(-0.5*pow((x-[1])/[2],2))",-3,2);

  for( int i=0; i<=2; i++ ){
    ngaus->SetParameter( i, fitfunc->GetParameter(i+3) );
    pgaus->SetParameter( i, fitfunc->GetParameter(i) );
  }

  double Nproton, Nneutron;
  
  double xminp = mup - 10.0*fitfunc->GetParameter( 2 );
  double xmaxp = mup + 10.0*fitfunc->GetParameter( 2 );
  
  //double Nproton = pgaus->Integral(xminp,xmaxp) / hdx->GetBinWidth(10);

  double xminn = mun - 10.0*fitfunc->GetParameter( 5 );
  double xmaxn = mun + 10.0*fitfunc->GetParameter( 5 );
  
  Nproton = pgaus->Integral(xminp,xmaxp) / hdx->GetBinWidth(10);
  Nneutron = ngaus->Integral(xminn, xmaxn) / hdx->GetBinWidth(10);
  
  cout << "Number of neutrons = " << Nneutron << " +/- " << sqrt(Nneutron) << endl;
  cout << "Number of protons = " << Nproton << " +/- " << sqrt(Nproton) << endl;

  double ratio = Nneutron/Nproton;
  double dratio = fabs(ratio) * sqrt( 1.0/Nneutron + 1.0/Nproton );

  cout << "n/p ratio = " << ratio << " +/- " << dratio << endl;

  double FFratio = sqrt(ratio) * 2.793/1.913;
  double dFFratio = dratio / (2.0*FFratio) * 2.793/1.913;

  cout << "GMn/GMp (neglect GEp, GEn) = " << FFratio << " +/- " << dFFratio << endl;
  
  //y = x^(1/2), dy = 1/(2y) dx

}
