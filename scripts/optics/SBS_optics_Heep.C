#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TTreeFormula.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDecompSVD.h"
#include <iostream>
#include <fstream>
#include "TMatrixD.h"
#include "TVectorD.h"
#include <vector>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TFile.h"

const double PI = TMath::Pi();
const double Mp = 0.938272;

void SBS_optics_Heep( const char *configfilename, const char *outfilename="optics_SBS_Heep_temp.root" ){
  ifstream configfile(configfilename);

  if( configfile ){
    TString currentline;

    TChain *C = new TChain("T");
    
    while ( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C->Add(currentline.Data());
      }
    }

    TCut globalcut = "";
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut += currentline.Data();
      }
    }

    TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);

    double ebeam=4.291, bbtheta=29.5, sbstheta=34.7, bbdist=1.63, sbsdist=2.8, hcaldist=17.0;
    double dx0=-1.0, dy0=0.0, dxsigma=0.2, dysigma=0.1;
    int usedxdycut = 0;
    double dt0=0.0, dtsigma=2.0;
    int usedtcut = 0;

    int optics_order = 2; 

    double W2min = 0.5, W2max = 1.3;
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig") ){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");

	if( tokens->GetEntries() >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	  if( skey == "W2min" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    W2min = sval.Atof();
	  }

	  if( skey == "W2max" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    W2max = sval.Atof();
	  }
	  
	  if( skey == "ebeam" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    ebeam = sval.Atof();
	  }

	  if( skey == "bbtheta" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    bbtheta = sval.Atof() * PI/180.0;
	  }

	  if( skey == "sbstheta" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    sbstheta = sval.Atof() * PI/180.0;
	  }

	  if( skey == "bbdist" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    bbdist = sval.Atof();
	  }

	  if( skey == "sbsdist" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    sbsdist = sval.Atof();
	  }

	  if( skey == "hcaldist" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    hcaldist = sval.Atof();
	  }

	  if( skey == "dx0" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    dx0 = sval.Atof();
	  }

	  if( skey == "dy0" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    dy0 = sval.Atof();
	  }

	  if( skey == "dxsigma" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    dxsigma = sval.Atof();
	  }

	  if( skey == "dysigma" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    dysigma = sval.Atof();
	  }

	  if( skey == "usedxdycut" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    usedxdycut = sval.Atoi();
	  }

	  if( skey == "dt0" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    dt0 = sval.Atof();
	  }

	  if( skey == "dtsigma" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    dtsigma = sval.Atof();
	  }

	  if( skey == "usedtcut" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    usedtcut = sval.Atoi();
	  }

	  if( skey == "optics_order" ){
	    TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	    optics_order = sval.Atoi();
	  }
	  
	  
	}
	
	
      }
    }

    int nterms = 0;

    vector<double> xtar_expon;
    
    for( int i=0; i<=optics_order; i++ ){
      for( int j=0; j<=optics_order-i; j++ ){
	for( int k=0; k<=optics_order-i-j; k++ ){
	  for( int l=0; l<=optics_order-i-j-k; l++ ){
	    for( int m=0; m<=optics_order-i-j-k-l; m++ ){
	      nterms++;
	      xtar_expon.push_back( i );
	    }
	  }
	}
      }
    }

    TMatrixD Moptics( nterms, nterms );
    TVectorD b_ytar( nterms );
    TVectorD b_thtar( nterms );
    TVectorD b_phtar( nterms );
    TVectorD b_pth( nterms );

    for( int i=0; i<nterms; i++ ){
      for( int j=0; j<nterms; j++ ){
	Moptics( i, j ) = 0.0;
      }
      b_ytar( i ) = 0.0;
      b_thtar( i ) = 0.0;
      b_phtar( i ) = 0.0;
      b_pth( i ) = 0.0;
    }
  
    
    TFile *fout = new TFile(outfilename,"RECREATE");

    C->SetBranchStatus("*",0);

    //BigBite Tracking branches:
    C->SetBranchStatus("bb.tr.n",1);
    C->SetBranchStatus("bb.tr.vz",1);
    C->SetBranchStatus("bb.tr.px",1);
    C->SetBranchStatus("bb.tr.py",1);
    C->SetBranchStatus("bb.tr.pz",1);
    C->SetBranchStatus("bb.tr.p",1);
    C->SetBranchStatus("bb.tr.tg_th",1);
    C->SetBranchStatus("bb.tr.tg_ph",1);
    C->SetBranchStatus("bb.tr.tg_y",1);
    C->SetBranchStatus("bb.tr.tg_x",1);
    C->SetBranchStatus("bb.gem.track.nhits",1);
    C->SetBranchStatus("bb.gem.track.chi2ndf",1);
    
    const int MAXNTRACKS = 100;
    
    double BB_ntrack;
    double BB_chi2ndf[MAXNTRACKS];
    double BB_vz[MAXNTRACKS];
    double BB_px[MAXNTRACKS], BB_py[MAXNTRACKS], BB_pz[MAXNTRACKS];
    double BB_p[MAXNTRACKS];
    double BB_thtar[MAXNTRACKS], BB_phtar[MAXNTRACKS], BB_xtar[MAXNTRACKS], BB_ytar[MAXNTRACKS];
    double BB_tracknhits[MAXNTRACKS];
    C->SetBranchAddress("bb.tr.n",&BB_ntrack);
    C->SetBranchAddress("bb.tr.vz", BB_vz);
    C->SetBranchAddress("bb.tr.px", BB_px);
    C->SetBranchAddress("bb.tr.py", BB_py);
    C->SetBranchAddress("bb.tr.pz", BB_pz);
    C->SetBranchAddress("bb.tr.p",BB_p);
    C->SetBranchAddress("bb.tr.tg_th",BB_thtar);
    C->SetBranchAddress("bb.tr.tg_ph",BB_phtar);
    C->SetBranchAddress("bb.tr.tg_x",BB_xtar);
    C->SetBranchAddress("bb.tr.tg_y",BB_ytar);
    C->SetBranchAddress("bb.gem.track.nhits",BB_tracknhits);
    C->SetBranchAddress("bb.gem.track.chi2ndf", BB_chi2ndf);
    
    C->SetBranchStatus("e.kine.W2",1);

    //don't need to set the branch address for W2 (we aren't going to use it directly)

    //SBS Tracking branches:
    C->SetBranchStatus("sbs.tr.n",1);
    C->SetBranchStatus("sbs.tr.vz",1);
    C->SetBranchStatus("sbs.tr.px",1);
    C->SetBranchStatus("sbs.tr.py",1);
    C->SetBranchStatus("sbs.tr.pz",1);
    C->SetBranchStatus("sbs.tr.p",1);
    C->SetBranchStatus("sbs.tr.r_x",1);
    C->SetBranchStatus("sbs.tr.r_y",1);
    C->SetBranchStatus("sbs.tr.r_th",1);
    C->SetBranchStatus("sbs.tr.r_ph",1);
    C->SetBranchStatus("sbs.tr.x",1);
    C->SetBranchStatus("sbs.tr.y",1);
    C->SetBranchStatus("sbs.tr.th",1);
    C->SetBranchStatus("sbs.tr.ph",1);
    C->SetBranchStatus("sbs.tr.tg_x",1);
    C->SetBranchStatus("sbs.tr.tg_y",1);
    C->SetBranchStatus("sbs.tr.tg_th",1);
    C->SetBranchStatus("sbs.tr.tg_ph",1);
    C->SetBranchStatus("sbs.gem.track.nhits",1);
    C->SetBranchStatus("sbs.gem.track.chi2ndf",1);

    double SBS_ntrack;
    double SBS_vz[MAXNTRACKS];
    double SBS_px[MAXNTRACKS], SBS_py[MAXNTRACKS], SBS_pz[MAXNTRACKS];
    double SBS_p[MAXNTRACKS];
    double SBS_rx[MAXNTRACKS], SBS_ry[MAXNTRACKS], SBS_rth[MAXNTRACKS], SBS_rph[MAXNTRACKS];
    double SBS_x[MAXNTRACKS], SBS_y[MAXNTRACKS], SBS_th[MAXNTRACKS], SBS_ph[MAXNTRACKS];
    double SBS_xtar[MAXNTRACKS], SBS_ytar[MAXNTRACKS], SBS_thtar[MAXNTRACKS], SBS_phtar[MAXNTRACKS];
    double SBS_tracknhits[MAXNTRACKS], SBS_chi2ndf[MAXNTRACKS];

    C->SetBranchAddress("sbs.tr.n",&SBS_ntrack);
    C->SetBranchAddress("sbs.tr.vz",SBS_vz);
    C->SetBranchAddress("sbs.tr.px",SBS_px);
    C->SetBranchAddress("sbs.tr.py",SBS_py);
    C->SetBranchAddress("sbs.tr.pz",SBS_pz);
    C->SetBranchAddress("sbs.tr.p",SBS_p);
    C->SetBranchAddress("sbs.tr.r_x", SBS_rx);
    C->SetBranchAddress("sbs.tr.r_y", SBS_ry);
    C->SetBranchAddress("sbs.tr.r_th", SBS_rth);
    C->SetBranchAddress("sbs.tr.r_ph", SBS_rph);
    C->SetBranchAddress("sbs.tr.x", SBS_x);
    C->SetBranchAddress("sbs.tr.y", SBS_y);
    C->SetBranchAddress("sbs.tr.th", SBS_th);
    C->SetBranchAddress("sbs.tr.ph", SBS_ph);
    C->SetBranchAddress("sbs.tr.tg_x", SBS_xtar);
    C->SetBranchAddress("sbs.tr.tg_y", SBS_ytar);
    C->SetBranchAddress("sbs.tr.tg_th", SBS_thtar);
    C->SetBranchAddress("sbs.tr.tg_ph", SBS_phtar);
    C->SetBranchAddress("sbs.gem.track.nhits", SBS_tracknhits);
    C->SetBranchAddress("sbs.gem.track.chi2ndf", SBS_chi2ndf);
    
    
    //BigBite calorimeter branches:
    C->SetBranchStatus("bb.ps.e",1);
    C->SetBranchStatus("bb.sh.e",1);
    C->SetBranchStatus("bb.ps.atimeblk",1);
    C->SetBranchStatus("bb.sh.atimeblk",1);

    double EPS, ESH, TPS, TSH;

    C->SetBranchAddress("bb.ps.e",&EPS);
    C->SetBranchAddress("bb.sh.e",&ESH);
    C->SetBranchAddress("bb.ps.atimeblk",&TPS);
    C->SetBranchAddress("bb.sh.atimeblk",&TSH);
    
    //SBS HCAL branches:
    C->SetBranchStatus("sbs.hcal.nclus",1);
    C->SetBranchStatus("sbs.hcal.x",1);
    C->SetBranchStatus("sbs.hcal.y",1);
    C->SetBranchStatus("sbs.hcal.e",1);
    C->SetBranchStatus("sbs.hcal.atimeblk",1);

    double xHCAL, yHCAL, EHCAL, THCAL;
    C->SetBranchAddress("sbs.hcal.x",&xHCAL);
    C->SetBranchAddress("sbs.hcal.y",&yHCAL);
    C->SetBranchAddress("sbs.hcal.e",&EHCAL);
    C->SetBranchAddress("sbs.hcal.atimeblk",&THCAL);

    TVector3 BB_zaxis( sin(bbtheta), 0, cos(bbtheta) );
    TVector3 BB_xaxis( 0, -1, 0 );
    TVector3 BB_yaxis = BB_zaxis.Cross( BB_xaxis ).Unit();

    TVector3 SBS_zaxis( -sin(sbstheta), 0, cos(sbstheta) );
    TVector3 SBS_xaxis( 0, -1, 0 );
    TVector3 SBS_yaxis = SBS_zaxis.Cross( SBS_xaxis ).Unit();

    TVector3 HCAL_center = hcaldist * SBS_zaxis;

    TH1D *hdeltat = new TH1D("hdeltat",";t_{HCAL}-t_{SH} (ns);", 600,-150,150);
    TH2D *hdxdy = new TH2D("hdxdy",";y_{HCAL}-y_{BB} (m); x_{HCAL}-x_{BB} (m)", 250,-1.25,1.25, 500,-5,0);
    TH1D *hdthtar = new TH1D("hdthtar", ";#theta_{tgt} (SBS) - #theta_{tgt} (BB) (rad);", 150, -0.1,0.1);
    TH1D *hdphtar = new TH1D("hdphtar", ";#phi_{tgt} (SBS) - #phi_{tgt} (BB) (rad);", 150, -0.1, 0.1 );
    TH1D *hdytar = new TH1D("hdytar", ";y_{tgt} (SBS) - y_{tgt} (BB) (m);", 150, -0.1, 0.1 );

    TH1D *hdpBBall = new TH1D("hdpBBall", ";p_{BB}/p_{el}(#theta_{BB})-1;", 150,-0.15,0.15);
    TH1D *hdpBBcut = new TH1D("hdpBBcut", ";p_{BB}/p_{el}(#theta_{BB}-1;", 150,-0.15,0.15);
    
    TH1D *hW2all = new TH1D("hW2all", "All events;W^{2} (GeV^{2});", 150, -1, 3 );
    TH1D *hW2cut = new TH1D("hW2cut", "Elastic;W^{2} (GeV^{2});", 150, -1, 3 );

    TTree *Tout = new TTree("Tout","SBS H(e,e'p) optics tree");

    //SBS quantities as reconstructed:
    //Focal-plane quantities for SBS. These are the important ones for optics:
    double T_xSBS, T_ySBS, T_thSBS, T_phSBS;
    //Target quantities for SBS:
    double T_xtarSBS, T_ytarSBS, T_thtarSBS, T_phtarSBS;
    double T_pSBS;
    double T_vzSBS;
    double T_pthetaSBS, T_pphiSBS;
    //SBS target quantities as predicted by BB:
    double T_eytarSBS, T_ethtarSBS, T_ephtarSBS;
    //physics quantities as reconstructed by BB:
    double T_etheta, T_ephi, T_ep, T_ptheta_expect, T_pphi_expect, T_ppexpect;
    double T_Q2, T_W2;
    double T_vzBB;
    double T_dx, T_dy, T_dt;
    double T_thetabendSBS;
    double T_ethetabendSBS;

    Tout->Branch( "xSBS", &T_xSBS );
    Tout->Branch( "ySBS", &T_ySBS );
    Tout->Branch( "thSBS", &T_thSBS );
    Tout->Branch( "phSBS", &T_phSBS );
    Tout->Branch( "xtarSBS", &T_xtarSBS );
    Tout->Branch( "ytarSBS", &T_ytarSBS );
    Tout->Branch( "thtarSBS", &T_thtarSBS );
    Tout->Branch( "phtarSBS", &T_phtarSBS );
    Tout->Branch( "pSBS", &T_pSBS );
    Tout->Branch( "vzSBS", &T_vzSBS );
    Tout->Branch( "pthetaSBS", &T_pthetaSBS );
    Tout->Branch( "pphiSBS", &T_pphiSBS );
    Tout->Branch( "thetabendSBS", &T_thetabendSBS );
    Tout->Branch( "ytarexpectSBS", &T_eytarSBS );
    Tout->Branch( "thtarexpectSBS", &T_ethtarSBS );
    Tout->Branch( "phtarexpectSBS", &T_ephtarSBS );
    Tout->Branch( "thetabendSBS_expect", &T_ethetabendSBS );
    Tout->Branch( "pexpectSBS", &T_ppexpect );
    Tout->Branch( "epBB", &T_ep );
    Tout->Branch( "ethetaBB", &T_etheta );
    Tout->Branch( "ephiBB", &T_ephi );
    Tout->Branch( "pthetaBB", &T_ptheta_expect );
    Tout->Branch( "pphiBB", &T_pphi_expect );
    Tout->Branch( "Q2", &T_Q2 );
    Tout->Branch( "vzBB", &T_vzBB );
    Tout->Branch( "W2", &T_W2 );
    Tout->Branch( "deltax", &T_dx );
    Tout->Branch( "deltay", &T_dy );
    Tout->Branch( "deltat", &T_dt );
    int T_HCALcut;
    Tout->Branch( "HCALcut", &T_HCALcut );
    
    
    
    int treenum=0, currenttreenum=0;
    
    long nevent=0;
    while( C->GetEntry( nevent++ ) ){
      if( nevent % 1000 == 0 ) cout << nevent << endl;

      currenttreenum = C->GetTreeNumber();

      if( nevent == 1 || currenttreenum != treenum ){
	GlobalCut->UpdateFormulaLeaves();
	treenum = currenttreenum;
      }

      bool passedglobal = GlobalCut->EvalInstance(0) != 0;

      if( passedglobal && int(BB_ntrack) >= 1 && int(SBS_ntrack)>=1 ){
	TVector3 ep( BB_px[0], BB_py[0], BB_pz[0] );
	TVector3 ephat = ep.Unit();

	TVector3 vertexBB( 0, 0, BB_vz[0] );

	TLorentzVector ep4vect( ep, BB_p[0] );
	TLorentzVector Ptarg(0,0,0,Mp);
	TLorentzVector Pbeam(0,0,ebeam,ebeam);
	TLorentzVector q4vect = Pbeam - ep4vect;

	double W2recon = (Ptarg + q4vect).M2();
	double Q2recon = -q4vect.M2();
	
	double etheta = ephat.Theta();
	double ephi = ephat.Phi();

	double eprime_etheta = ebeam/(1.+ebeam/Mp*(1.-cos(etheta)));
	double nu_etheta = ebeam - eprime_etheta;
	double pp_etheta = sqrt(pow(nu_etheta,2)+2.0*Mp*nu_etheta);
	double ptheta_etheta = acos( (ebeam-eprime_etheta*cos(etheta))/pp_etheta );
	double pphi_etheta = ephi + PI;

	TVector3 pphat_expect_global( sin(ptheta_etheta)*cos(pphi_etheta),
				      sin(ptheta_etheta)*sin(pphi_etheta),
				      cos(ptheta_etheta) );

	TVector3 pphat_expect_SBS( pphat_expect_global.Dot( SBS_xaxis ),
				   pphat_expect_global.Dot( SBS_yaxis ),
				   pphat_expect_global.Dot( SBS_zaxis ) );

	TVector3 pphat_SBS( SBS_thtar[0], SBS_phtar[0], 1.0 );
	pphat_SBS = pphat_SBS.Unit();

	TVector3 pphat_fp_SBS( SBS_rth[0], SBS_rph[0], 1.0 );
	pphat_fp_SBS = pphat_fp_SBS.Unit();
	
	TVector3 pphat_global = pphat_SBS.X() * SBS_xaxis + pphat_SBS.Y() * SBS_yaxis +
	  pphat_SBS.Z() * SBS_zaxis;

	double ptheta_SBS = pphat_global.Theta();
	double pphi_SBS = pphat_global.Phi();

	if( pphi_SBS < 0.0 ) pphi_SBS += 2.0*PI;
	
	double sbsthtar_expect = pphat_expect_SBS.X()/pphat_expect_SBS.Z();
	double sbsphtar_expect = pphat_expect_SBS.Y()/pphat_expect_SBS.Z();

	double sbsytar_expect = vertexBB.Z() * (sin(sbstheta) - sbsphtar_expect*cos(sbstheta));

	double sintersect_hcal = ( HCAL_center - vertexBB ).Dot( SBS_zaxis ) / (pphat_expect_global.Dot( SBS_zaxis ) );

	TVector3 HCAL_intersect = vertexBB + pphat_expect_global * sintersect_hcal;

	double xHCAL_expect = (HCAL_intersect - HCAL_center).Dot( SBS_xaxis );
	double yHCAL_expect = (HCAL_intersect - HCAL_center).Dot( SBS_yaxis );

	bool passedxy = sqrt(pow( (xHCAL-xHCAL_expect-dx0)/dxsigma, 2 ) + pow( (yHCAL-yHCAL_expect-dy0)/dysigma, 2 ) ) <= 3.5;
	bool passedtime = abs( THCAL-TSH-dt0 ) <= 3.5*dtsigma;

	bool passedHCAL = (passedxy || usedxdycut == 0 ) && (passedtime || usedtcut == 0);

	if( passedHCAL ){
	  //hW2cut->Fill( W2recon );
	  //hdpBBcut->Fill( BB_p[0]/eprime_etheta-1. );
	}

	T_thetabendSBS = acos( pphat_SBS.Dot( pphat_fp_SBS ) );
	T_ethetabendSBS = acos( pphat_expect_SBS.Dot( pphat_fp_SBS ) );
	
	
	if( W2recon >= W2min && W2recon <= W2max ){
	  //hdeltat->Fill( THCAL - TSH );
	  //hdxdy->Fill( yHCAL - yHCAL_expect, xHCAL - xHCAL_expect );

	  if( passedHCAL ){
	    //hdthtar->Fill( SBS_thtar[0]-sbsthtar_expect );
	    //hdphtar->Fill( SBS_phtar[0]-sbsphtar_expect );
	    //hdytar->Fill( SBS_ytar[0]-sbsytar_expect );

	    vector<double> term( nterms, 0.0 );

	    int ipar=0;

	    double xtar = 0.0;
	    //Increment optics matrix elements
	    for( int i=0; i<=optics_order; i++ ){
	      for( int j=0; j<=optics_order-i; j++ ){
		for( int k=0; k<=optics_order-i-j; k++ ){
		  for( int l=0; l<=optics_order-i-j-k; l++ ){
		    for( int m=0; m<=optics_order-i-j-k-l; m++ ){
		      term[ipar] = pow( SBS_rx[0], m ) * pow( SBS_ry[0], l ) * pow( SBS_rth[0], k ) * pow( SBS_rph[0], j ) * pow( xtar, i );

		      b_ytar( ipar ) += sbsytar_expect * term[ipar];
		      b_thtar( ipar ) += sbsthtar_expect * term[ipar];
		      b_phtar( ipar ) += sbsphtar_expect * term[ipar];
		      b_pth( ipar ) += pp_etheta * T_ethetabendSBS * term[ipar];

		      ipar++;
		    }
		  }
		}
	      }
	    }

	    for( ipar = 0; ipar<nterms; ipar++ ){
	      for( int jpar=0; jpar<nterms; jpar++ ){
		Moptics( ipar, jpar ) += term[ipar]*term[jpar];
	      }
	    }

		      
	  }
	}
	//	hdpBBall->Fill( BB_p[0]/eprime_etheta-1. );
	//hW2all->Fill( W2recon );
      

	T_xSBS = SBS_rx[0];
	T_ySBS = SBS_ry[0];
	T_thSBS = SBS_rth[0];
	T_phSBS = SBS_rph[0];

	T_xtarSBS = 0.0;
	T_ytarSBS = SBS_ytar[0];
	T_thtarSBS = SBS_thtar[0];
	T_phtarSBS = SBS_phtar[0];
	T_pSBS = SBS_p[0];
	//T_vzSBS = //SBS_vz[0];
	T_vzSBS = SBS_ytar[0]/(sin(sbstheta)-cos(sbstheta)*SBS_phtar[0]);
	T_pthetaSBS = ptheta_SBS;
	T_pphiSBS = pphi_SBS;
	T_eytarSBS = sbsytar_expect;
	T_ethtarSBS = sbsthtar_expect;
	T_ephtarSBS = sbsphtar_expect;
	T_ppexpect = pp_etheta;
	T_etheta = etheta;
	T_ephi = ephi;
	T_ptheta_expect = ptheta_etheta;
	T_pphi_expect = pphi_etheta;
	T_Q2 = Q2recon;
	T_vzBB = BB_vz[0];
	T_W2 = W2recon;
	T_dx = xHCAL - xHCAL_expect;
	T_dy = yHCAL - yHCAL_expect;
	T_dt = THCAL-TSH;
	T_HCALcut = passedHCAL;
	//T_thetabendSBS = SBS_thtar[0] - SBS_rth[0];
	//T_ethetabendSBS = sbsthtar_expect - SBS_rth[0];

	
	
	T_ep = BB_p[0];
	
	//Tout->Fill();
      }
    }

    for( int ipar=0; ipar<nterms; ipar++ ){
      if( xtar_expon[ipar] > 0 ){
	Moptics(ipar,ipar)=1.0;
	b_ytar(ipar) = 0.0;
	b_thtar(ipar) = 0.0;
	b_phtar(ipar) = 0.0;
	b_pth(ipar) = 0.0;
      }
      for( int jpar=0; jpar<nterms; jpar++ ){
	if( jpar != ipar && (xtar_expon[ipar]>0 || xtar_expon[jpar]>0) ){
	  Moptics(ipar,jpar) = 0.0;
	  Moptics(jpar,ipar) = 0.0;
	}
      }

    }
    
    TDecompSVD A_ytar( Moptics );
    TDecompSVD A_thtar( Moptics );
    TDecompSVD A_phtar( Moptics );
    TDecompSVD A_pth( Moptics );

    A_ytar.Solve( b_ytar );
    A_thtar.Solve( b_thtar );
    A_phtar.Solve( b_phtar );
    A_pth.Solve( b_pth );

    nevent = 0;

    treenum = 0;
    currenttreenum = 0;

    double T_ytarSBS_fit, T_thtarSBS_fit, T_phtarSBS_fit, T_pSBS_fit, T_vzSBS_fit;
    double T_thetabendSBS_fit;

    Tout->Branch( "ytarSBSfit", &T_ytarSBS_fit );
    Tout->Branch( "thtarSBSfit", &T_thtarSBS_fit );
    Tout->Branch( "phtarSBSfit", &T_phtarSBS_fit );
    Tout->Branch( "pSBSfit", &T_pSBS_fit );
    Tout->Branch( "vzSBSfit", &T_vzSBS_fit );
    Tout->Branch( "thetabendSBSfit", &T_thetabendSBS_fit );
    
    //loop over the tree again; this time fill the output tree including reconstructed quantities AFTER the fit:
    
    while( C->GetEntry( nevent++ ) ){
      if( nevent % 1000 == 0 ) cout << nevent << endl;

      currenttreenum = C->GetTreeNumber();

      if( nevent == 1 || currenttreenum != treenum ){
	GlobalCut->UpdateFormulaLeaves();
	treenum = currenttreenum;
      }

      bool passedglobal = GlobalCut->EvalInstance(0) != 0;

      if( passedglobal && int(BB_ntrack) >= 1 && int(SBS_ntrack)>=1 ){
	TVector3 ep( BB_px[0], BB_py[0], BB_pz[0] );
	TVector3 ephat = ep.Unit();

	TVector3 vertexBB( 0, 0, BB_vz[0] );

	TLorentzVector ep4vect( ep, BB_p[0] );
	TLorentzVector Ptarg(0,0,0,Mp);
	TLorentzVector Pbeam(0,0,ebeam,ebeam);
	TLorentzVector q4vect = Pbeam - ep4vect;

	double W2recon = (Ptarg + q4vect).M2();
	double Q2recon = -q4vect.M2();
	
	double etheta = ephat.Theta();
	double ephi = ephat.Phi();

	double eprime_etheta = ebeam/(1.+ebeam/Mp*(1.-cos(etheta)));
	double nu_etheta = ebeam - eprime_etheta;
	double pp_etheta = sqrt(pow(nu_etheta,2)+2.0*Mp*nu_etheta);
	double ptheta_etheta = acos( (ebeam-eprime_etheta*cos(etheta))/pp_etheta );
	double pphi_etheta = ephi + PI;

	TVector3 pphat_expect_global( sin(ptheta_etheta)*cos(pphi_etheta),
				      sin(ptheta_etheta)*sin(pphi_etheta),
				      cos(ptheta_etheta) );

	TVector3 pphat_expect_SBS( pphat_expect_global.Dot( SBS_xaxis ),
				   pphat_expect_global.Dot( SBS_yaxis ),
				   pphat_expect_global.Dot( SBS_zaxis ) );

	TVector3 pphat_SBS( SBS_thtar[0], SBS_phtar[0], 1.0 );
	pphat_SBS = pphat_SBS.Unit();

	TVector3 pphat_fp_SBS( SBS_rth[0], SBS_rph[0], 1.0 );
	pphat_fp_SBS = pphat_fp_SBS.Unit();
	
	TVector3 pphat_global = pphat_SBS.X() * SBS_xaxis + pphat_SBS.Y() * SBS_yaxis +
	  pphat_SBS.Z() * SBS_zaxis;

	double ptheta_SBS = pphat_global.Theta();
	double pphi_SBS = pphat_global.Phi();

	if( pphi_SBS < 0.0 ) pphi_SBS += 2.0*PI;
	
	double sbsthtar_expect = pphat_expect_SBS.X()/pphat_expect_SBS.Z();
	double sbsphtar_expect = pphat_expect_SBS.Y()/pphat_expect_SBS.Z();

	double sbsytar_expect = vertexBB.Z() * (sin(sbstheta) - sbsphtar_expect*cos(sbstheta));

	double sintersect_hcal = ( HCAL_center - vertexBB ).Dot( SBS_zaxis ) / (pphat_expect_global.Dot( SBS_zaxis ) );

	TVector3 HCAL_intersect = vertexBB + pphat_expect_global * sintersect_hcal;

	double xHCAL_expect = (HCAL_intersect - HCAL_center).Dot( SBS_xaxis );
	double yHCAL_expect = (HCAL_intersect - HCAL_center).Dot( SBS_yaxis );

	bool passedxy = sqrt(pow( (xHCAL-xHCAL_expect-dx0)/dxsigma, 2 ) + pow( (yHCAL-yHCAL_expect-dy0)/dysigma, 2 ) ) <= 3.5;
	bool passedtime = abs( THCAL-TSH-dt0 ) <= 3.5*dtsigma;

	bool passedHCAL = (passedxy || usedxdycut == 0 ) && (passedtime || usedtcut == 0);

	if( passedHCAL ){
	  hW2cut->Fill( W2recon );
	  hdpBBcut->Fill( BB_p[0]/eprime_etheta-1. );
	}

	T_thetabendSBS = acos( pphat_SBS.Dot( pphat_fp_SBS ) );
	T_ethetabendSBS = acos( pphat_expect_SBS.Dot( pphat_fp_SBS ) );

	//Reconstruct optics based on fit results:
	T_ytarSBS_fit = 0.0;
	T_thtarSBS_fit = 0.0;
	T_phtarSBS_fit = 0.0;
	double pthSBS_fit = 0.0;

	int ipar = 0;

	vector<double> term(nterms,0.0);
	
	double xtar = 0.0;
	for( int i=0; i<=optics_order; i++ ){
	  for( int j=0; j<=optics_order-i; j++ ){
	    for( int k=0; k<=optics_order-i-j; k++ ){
	      for( int l=0; l<=optics_order-i-j-k; l++ ){
		for( int m=0; m<=optics_order-i-j-k-l; m++ ){
		  term[ipar] = pow( SBS_rx[0], m ) * pow( SBS_ry[0], l ) * pow( SBS_rth[0], k ) * pow( SBS_rph[0], j ) * pow( xtar, i );
		  T_ytarSBS_fit += term[ipar]*b_ytar(ipar);
		  T_thtarSBS_fit += term[ipar]*b_thtar(ipar);
		  T_phtarSBS_fit += term[ipar]*b_phtar(ipar);
		  pthSBS_fit += term[ipar]*b_pth(ipar);
		  
		  ipar++;
		}
	      }
	    }
	  }
	}

	TVector3 pphat_SBS_fit( T_thtarSBS_fit, T_phtarSBS_fit, 1.0 );
	pphat_SBS_fit = pphat_SBS_fit.Unit();

	T_thetabendSBS_fit = acos( pphat_SBS_fit.Dot( pphat_fp_SBS ) );

	T_pSBS_fit = pthSBS_fit / T_thetabendSBS_fit;
	T_vzSBS_fit = T_ytarSBS_fit / ( sin(sbstheta) - cos(sbstheta)*T_phtarSBS_fit );
	
	if( W2recon >= W2min && W2recon <= W2max ){
	  hdeltat->Fill( THCAL - TSH );
	  hdxdy->Fill( yHCAL - yHCAL_expect, xHCAL - xHCAL_expect );

	  if( passedHCAL ){
	    hdthtar->Fill( SBS_thtar[0]-sbsthtar_expect );
	    hdphtar->Fill( SBS_phtar[0]-sbsphtar_expect );
	    hdytar->Fill( SBS_ytar[0]-sbsytar_expect );
	      
	  }
	}
	hdpBBall->Fill( BB_p[0]/eprime_etheta-1. );
	hW2all->Fill( W2recon );
      

	T_xSBS = SBS_rx[0];
	T_ySBS = SBS_ry[0];
	T_thSBS = SBS_rth[0];
	T_phSBS = SBS_rph[0];

	T_xtarSBS = 0.0;
	T_ytarSBS = SBS_ytar[0];
	T_thtarSBS = SBS_thtar[0];
	T_phtarSBS = SBS_phtar[0];
	T_pSBS = SBS_p[0];
	//T_vzSBS = //SBS_vz[0];
	T_vzSBS = SBS_ytar[0]/(sin(sbstheta)-cos(sbstheta)*SBS_phtar[0]);
	T_pthetaSBS = ptheta_SBS;
	T_pphiSBS = pphi_SBS;
	T_eytarSBS = sbsytar_expect;
	T_ethtarSBS = sbsthtar_expect;
	T_ephtarSBS = sbsphtar_expect;
	T_ppexpect = pp_etheta;
	T_etheta = etheta;
	T_ephi = ephi;
	T_ptheta_expect = ptheta_etheta;
	T_pphi_expect = pphi_etheta;
	T_Q2 = Q2recon;
	T_vzBB = BB_vz[0];
	T_W2 = W2recon;
	T_dx = xHCAL - xHCAL_expect;
	T_dy = yHCAL - yHCAL_expect;
	T_dt = THCAL-TSH;
	T_HCALcut = passedHCAL;
	//T_thetabendSBS = SBS_thtar[0] - SBS_rth[0];
	//T_ethetabendSBS = sbsthtar_expect - SBS_rth[0];

	
	
	T_ep = BB_p[0];
	
	Tout->Fill();
      }
    }
    
    fout->Write();
  }
}
