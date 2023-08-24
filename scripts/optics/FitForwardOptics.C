#include "TChain.h"
#include "TTree.h"
#include "TEventList.h"
#include "TCut.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TCut.h"
#include "TEventList.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDecompSVD.h"

void FitForwardOptics( const char *configfilename, const char *outfilename="foptics_bb_temp.root" ){

  TChain *C = new TChain("T");

  ifstream infile(configfilename);
 
  TString currentline;

  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline);
    }
  }

  TCut cut = "";
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      cut += currentline;
    }
  }

  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",cut);

  cout << "Number of events passing global cut = " << elist->GetN() << endl;

  int order = 4;
  infile >> order;

  cout << "order = " << order << endl;
  
  vector<int> xtar_expon;
  //vector<int> xptar_expon;
  //vector<int> yptar_expon
  
  int nparams = 0;
  //The purpose of this macro is only to fit forward optics matrix elements:
  for( int i=0; i<=order; i++){
    for( int j=0; j<=order-i; j++ ){
      for( int k=0; k<=order-i-j; k++ ){
	for( int l=0; l<=order-i-j-k; l++ ){
	  for( int m=0; m<=order-i-j-k-l; m++ ){
	    xtar_expon.push_back( i );
	    nparams++;
	  }
	}
      }
    }
  }

  cout << "xtar_expon.size() = " << xtar_expon.size();
  cout << "n params = " << nparams << endl;
  
  TMatrixD Mforward( nparams, nparams );
  TVectorD b_xfp( nparams );
  TVectorD b_yfp( nparams );
  TVectorD b_thfp( nparams );
  TVectorD b_phfp( nparams );
  
  for( int i=0; i<nparams; i++ ){
    b_xfp(i) = 0.0;
    b_yfp(i) = 0.0;
    b_thfp(i) = 0.0;
    b_phfp(i) = 0.0;
    for( int j=0; j<nparams; j++ ){
      Mforward( i, j ) = 0.0;
    }
  }
  
  // double fpmin[4], fpmax[4];    //ranges for plotting of fp variables: order is: x, y, xp, yp
  // double tgtmin[5], tgtmax[5];  //ranges for plotting of tgt variables: order is: xp, yp, y, p, x

  // // for(int i=0; i<4; i++){
  // //   infile >> fpmin[i] >> fpmax[i];
  // // }
  // // for(int i=0; i<5; i++){
  // //   infile >> tgtmin[i] >> tgtmax[i];
  // // }

  int xtar_flag=1; //include and fit xtar-dependent terms in the expansion or not? (also interpret flag as maximum order for xtar-dependent terms)
  infile >> xtar_flag;

  int nbins = 250;

  TFile *fout = new TFile(outfilename,"RECREATE");
  
   //Histograms of "recon" versions of fp track parameters:
  TH1D *hxfpdiff = new TH1D("hxfpdiff", "", nbins, -0.005, 0.005);
  TH1D *hyfpdiff = new TH1D("hyfpdiff", "", nbins, -0.005, 0.005);
  TH1D *hthfpdiff = new TH1D("hthfpdiff", "", nbins, -0.001, 0.001 );
  TH1D *hphfpdiff = new TH1D("hphfpdiff", "", nbins, -0.001, 0.001 );

  TH2D *hdxfp_xfp = new TH2D("hdxfp_xfp", ";x_{fp};x_{fp}(fit-true)", nbins, -0.75, 0.75, 250, -0.005, 0.005 );
  TH2D *hdxfp_yfp = new TH2D("hdxfp_yfp", ";y_{fp};x_{fp}(fit-true)", nbins, -0.2, 0.2, 250, -0.005, 0.005 );
  TH2D *hdxfp_thfp = new TH2D("hdxfp_thfp", ";#theta_{fp};x_{fp}(fit-true)", nbins, -0.5, 0.5, 250, -0.005, 0.005 );
  TH2D *hdxfp_phfp = new TH2D("hdxfp_phfp", ";#phi_{fp};x_{fp}(fit-true)", nbins, -0.15, 0.15, 250, -0.005, 0.005 );

  TH2D *hdyfp_xfp = new TH2D("hdyfp_xfp", ";x_{fp};y_{fp}(fit-true)", nbins, -0.75, 0.75, 250, -0.005, 0.005 );
  TH2D *hdyfp_yfp = new TH2D("hdyfp_yfp", ";y_{fp};y_{fp}(fit-true)", nbins, -0.2, 0.2, 250, -0.005, 0.005 );
  TH2D *hdyfp_thfp = new TH2D("hdyfp_thfp", ";#theta_{fp};y_{fp}(fit-true)", nbins, -0.5, 0.5, 250, -0.005, 0.005 );
  TH2D *hdyfp_phfp = new TH2D("hdyfp_phfp", ";#phi_{fp};y_{fp}(fit-true)", nbins, -0.15, 0.15, 250, -0.005, 0.005 );


  TH2D *hdthfp_xfp = new TH2D("hdthfp_xfp", ";x_{fp};#theta_{fp}(fit-true)", nbins, -0.75, 0.75, 250, -0.002, 0.002);
  TH2D *hdthfp_yfp = new TH2D("hdthfp_yfp", ";y_{fp};#theta_{fp}(fit-true)", nbins, -0.2, 0.2, 250, -0.002, 0.002 );
  TH2D *hdthfp_thfp = new TH2D("hdthfp_thfp", ";#theta_{fp};#theta_{fp}(fit-true)", nbins, -0.5, 0.5, 250, -0.002, 0.002 );
  TH2D *hdthfp_phfp = new TH2D("hdthfp_phfp", ";#phi_{fp};#theta_{fp}(fit-true)", nbins, -0.15, 0.15, 250, -0.002, 0.002 );

  TH2D *hdphfp_xfp = new TH2D("hdphfp_xfp", ";x_{fp};#phi_{fp}(fit-true)", nbins, -0.75, 0.75, 250, -0.002, 0.002);
  TH2D *hdphfp_yfp = new TH2D("hdphfp_yfp", ";y_{fp};#phi_{fp}(fit-true)", nbins, -0.2, 0.2, 250, -0.002, 0.002 );
  TH2D *hdphfp_thfp = new TH2D("hdphfp_thfp", ";#theta_{fp};#phi_{fp}(fit-true)", nbins, -0.5, 0.5, 250, -0.002, 0.002 );
  TH2D *hdphfp_phfp = new TH2D("hdphfp_phfp", ";#phi_{fp};#phi_{fp}(fit-true)", nbins, -0.15, 0.15, 250, -0.002, 0.002 );
  
  //C->Add(rootfilename);

 
  
  C->SetBranchStatus("*",0);
  
  //BigBite track variables:
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vx",1);
  C->SetBranchStatus("bb.tr.vy",1);
  C->SetBranchStatus("bb.tr.vz",1);
  
  C->SetBranchStatus("bb.tr.r_x",1);
  C->SetBranchStatus("bb.tr.r_y",1);
  C->SetBranchStatus("bb.tr.r_th",1);
  C->SetBranchStatus("bb.tr.r_ph",1);
  
  C->SetBranchStatus("bb.tr.tg_x",1);
  C->SetBranchStatus("bb.tr.tg_y",1);
  C->SetBranchStatus("bb.tr.tg_th",1);
  C->SetBranchStatus("bb.tr.tg_ph",1);

  C->SetBranchStatus("bb.gem.track.ntrack",1);
  C->SetBranchStatus("bb.gem.track.nhits",1);
  C->SetBranchStatus("bb.gem.track.chi2ndf",1);
  
  //Shower and preshower variables:
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);

  //C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("bb.x_bcp",1);
  C->SetBranchStatus("bb.y_bcp",1);
  C->SetBranchStatus("bb.z_bcp",1);

  C->SetBranchStatus("bb.x_fcp",1);
  C->SetBranchStatus("bb.y_fcp",1);
  C->SetBranchStatus("bb.z_fcp",1);
  
  int MAXNTRACKS=100;
  
  double ntracks;
  double p[MAXNTRACKS],px[MAXNTRACKS],py[MAXNTRACKS],pz[MAXNTRACKS];
  double vx[MAXNTRACKS],vy[MAXNTRACKS],vz[MAXNTRACKS];

  double xfp[MAXNTRACKS],yfp[MAXNTRACKS],thfp[MAXNTRACKS],phfp[MAXNTRACKS];
  double xtgt[MAXNTRACKS],ytgt[MAXNTRACKS],thtgt[MAXNTRACKS],phtgt[MAXNTRACKS];

  double rxfp[MAXNTRACKS],ryfp[MAXNTRACKS],rthfp[MAXNTRACKS],rphfp[MAXNTRACKS];

  double EPS,ESH,xSH,ySH,xPS;
  double x_bcp,y_bcp,z_bcp,x_fcp,y_fcp,z_fcp;

  C->SetBranchAddress("bb.gem.track.ntrack",&ntracks);
  C->SetBranchAddress("bb.tr.p",p);
  C->SetBranchAddress("bb.tr.px",px);
  C->SetBranchAddress("bb.tr.py",py);
  C->SetBranchAddress("bb.tr.pz",pz);
  C->SetBranchAddress("bb.tr.vx",vx);
  C->SetBranchAddress("bb.tr.vy",vy);
  C->SetBranchAddress("bb.tr.vz",vz);

  //Focal-plane track variables (use "rotated" versions):
  C->SetBranchAddress("bb.tr.r_x",rxfp);
  C->SetBranchAddress("bb.tr.r_y",ryfp);
  C->SetBranchAddress("bb.tr.r_th",rthfp);
  C->SetBranchAddress("bb.tr.r_ph",rphfp);

  //Focal-plane track variables (use "unrotated" versions):
  C->SetBranchAddress("bb.tr.x",xfp);
  C->SetBranchAddress("bb.tr.y",yfp);
  C->SetBranchAddress("bb.tr.th",thfp);
  C->SetBranchAddress("bb.tr.ph",phfp);

  //Target track variables (other than momentum):
  C->SetBranchAddress("bb.tr.tg_x",xtgt);
  C->SetBranchAddress("bb.tr.tg_y",ytgt);
  C->SetBranchAddress("bb.tr.tg_th",thtgt);
  C->SetBranchAddress("bb.tr.tg_ph",phtgt);

  C->SetBranchAddress("bb.ps.e",&EPS);
  C->SetBranchAddress("bb.sh.e",&ESH);
  C->SetBranchAddress("bb.ps.x",&xPS);
  C->SetBranchAddress("bb.sh.x",&xSH);
  C->SetBranchAddress("bb.sh.y",&ySH);
  
  long nevent=0;
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;

    if( int(ntracks) == 1 ){
      
      
      vector<double> term(nparams);
      int ipar=0;
      for(int i=0; i<=order; i++){
	for(int j=0; j<=order-i; j++){
	  for(int k=0; k<=order-i-j; k++){
	    for(int l=0; l<=order-i-j-k; l++){
	      for(int m=0; m<=order-i-j-k-l; m++){
		term[ipar] = pow( thtgt[0], m ) * pow( phtgt[0], l ) * pow( ytgt[0], k ) * pow( 1.0/p[0], j ) * pow( xtgt[0], i );

		b_xfp[ipar] += term[ipar] * rxfp[0];
		b_yfp[ipar] += term[ipar] * ryfp[0];
		b_thfp[ipar] += term[ipar] * rthfp[0];
		b_phfp[ipar] += term[ipar] * rphfp[0];

		ipar++;
	      }
	    }
	  }
	}
      }

      for( ipar=0; ipar<nparams; ipar++ ){
	for( int jpar=0; jpar<nparams; jpar++ ){
	  Mforward(ipar,jpar) += term[ipar]*term[jpar];
	}
      }
    }
  }

  cout << "xtar_flag, order = " << xtar_flag << ", " << order << endl;
  
  if( xtar_flag != order ){ //then zero-out all xtar-dependent terms
    for( int ipar=0; ipar<nparams; ipar++ ){
      cout << "ipar, xtar_expon = " << ipar << ", " << xtar_expon[ipar] << endl;
      if( xtar_expon[ipar] > 0 ){
	Mforward(ipar,ipar) = 1.0;
	b_xfp(ipar) = 0.0;
	b_yfp(ipar) = 0.0;
	b_thfp(ipar) = 0.0;
	b_phfp(ipar) = 0.0;
      
	for( int jpar=0; jpar<nparams; jpar++ ){
	  if( jpar != ipar ){
	    Mforward(ipar,jpar) = 0.0;
	  }
	}
      }
    }
  }

  // Solving for forward matrix elements
  TDecompSVD A_xfp(Mforward);
  TDecompSVD A_yfp(Mforward);
  TDecompSVD A_thfp(Mforward);
  TDecompSVD A_phfp(Mforward);

  A_xfp.Solve(b_xfp);
  A_yfp.Solve(b_yfp);
  A_thfp.Solve(b_thfp);
  A_phfp.Solve(b_phfp);

  TString fopticsfilename = outfilename;
  fopticsfilename.ReplaceAll(".root",".dat");
  
  ofstream fopticsfile( fopticsfilename.Data() );
  fopticsfile << nparams << endl;

  int ipar=0;
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    TString ccoeff;
	    ccoeff.Form( " %15.8g ", b_xfp(ipar) );
	    fopticsfile << ccoeff;
	    ccoeff.Form( " %15.8g ", b_yfp(ipar) );
	    fopticsfile << ccoeff;
	    ccoeff.Form( " %15.8g ", b_thfp(ipar) );
	    fopticsfile << ccoeff;
	    ccoeff.Form( " %15.8g ", b_phfp(ipar) );
	    fopticsfile << ccoeff;

	    TString cexpon;
	    cexpon.Form( " %d %d %d %d %d", m, l, k, j, i );
	    fopticsfile << cexpon << endl;

	    ipar++;
	  }
	}
      }
    }
  }

  nevent=0;
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;

    if( int(ntracks) == 1 ){
      double xfp_recon = 0.0, yfp_recon = 0.0, thfp_recon = 0.0, phfp_recon = 0.0;

      int ipar=0;
      for(int i=0; i<=order; i++){
	for(int j=0; j<=order-i; j++){
	  for(int k=0; k<=order-i-j; k++){
	    for(int l=0; l<=order-i-j-k; l++){
	      for(int m=0; m<=order-i-j-k-l; m++){
		double term = pow( thtgt[0], m ) * pow( phtgt[0], l ) * pow( ytgt[0], k ) * pow( 1.0/p[0], j ) * pow( xtgt[0], i );
		
		xfp_recon += term * b_xfp(ipar);
		yfp_recon += term * b_yfp(ipar);
		thfp_recon += term * b_thfp(ipar);
		phfp_recon += term * b_phfp(ipar);
		
		ipar++;
	      }
	    }
	  }
	}
      }

      hxfpdiff->Fill( rxfp[0]-xfp_recon );
      hyfpdiff->Fill( ryfp[0]-yfp_recon );
      hthfpdiff->Fill( rthfp[0]-thfp_recon );
      hphfpdiff->Fill( rphfp[0]-phfp_recon );

      hdxfp_xfp->Fill( rxfp[0], rxfp[0]-xfp_recon );
      hdyfp_xfp->Fill( rxfp[0], ryfp[0]-yfp_recon );
      hdthfp_xfp->Fill( rxfp[0], rthfp[0]-thfp_recon );
      hdphfp_xfp->Fill( rxfp[0], rphfp[0]-phfp_recon );

      hdxfp_yfp->Fill( ryfp[0], rxfp[0]-xfp_recon );
      hdyfp_yfp->Fill( ryfp[0], ryfp[0]-yfp_recon );
      hdthfp_yfp->Fill( ryfp[0], rthfp[0]-thfp_recon );
      hdphfp_yfp->Fill( ryfp[0], rphfp[0]-phfp_recon );

      hdxfp_thfp->Fill( rthfp[0], rxfp[0]-xfp_recon );
      hdyfp_thfp->Fill( rthfp[0], ryfp[0]-yfp_recon );
      hdthfp_thfp->Fill( rthfp[0], rthfp[0]-thfp_recon );
      hdphfp_thfp->Fill( rthfp[0], rphfp[0]-phfp_recon );
      
      hdxfp_phfp->Fill( rphfp[0], rxfp[0]-xfp_recon );
      hdyfp_phfp->Fill( rphfp[0], ryfp[0]-yfp_recon );
      hdthfp_phfp->Fill( rphfp[0], rthfp[0]-thfp_recon );
      hdphfp_phfp->Fill( rphfp[0], rphfp[0]-phfp_recon );
    }
  }

  fout->Write();
  
}
