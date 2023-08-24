#include "TVector3.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TRandom3.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TCutG.h"
#include "TEventList.h"
#include "TMinuit.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLine.h"
#include "TText.h"
#include "TSystem.h"
#include "TF1.h"
#include "TEllipse.h"
#include <iostream>
#include <fstream>

// const int nxsieve = 13;
// const int nysieve = 7;

const double PI = TMath::Pi();

//Example code to align the GEM stack relative to the spectrometer coordinate system and determine distance from target to sieve plate using
//known dimensions and hole spacings of sieve and known dimensions and internal alignments of GEMs:
//We are going to be lazy and use MINUIT even though this problem could easily be linearized:

const double sigma_pos = 0.001; //1 mm w/multiple scattering (wild guess)
const double sigma_slope = sigma_pos/1.6; //1 mm/1.6 m ~= 6e-4 gives rough order of magnitude for slope uncertainty, of course multiple scattering might invalidate that.

double GEMX0, GEMY0, GEMZ0;
double GEMtheta, GEMphi;
double sigmaX0=sigma_pos, sigmaY0=sigma_pos, sigmaZ0=sigma_pos; //widths for "penalty terms"
double sigmatheta=0.0001, sigmaphi=0.001;

//these are user-configurable:
double sigma_xclose = 0.007;
double sigma_yclose = 0.007;
double sigma_zclose = 0.001;

int NTRACKS;

//These arrays will be filled by Holly's code to select tracks going through the individual sieve holes:
vector<double> XTRACK,YTRACK,XPTRACK,YPTRACK; //Measured track parameters in GEM internal coordinate system
//vector<double> XSIEVE,YSIEVE; //known x and y center positions of the holes the tracks went through.
//We assume the single foil position defines the global origin and acts as a point source of straight-line rays coming through the sieve holes into the GEMS

// This version doesn't rely on the sieve slit, it relies on minimizing the sum of squared distances of closest approach between the tracks and the origin:
void CHI2_FCN( int &npar, double *gin, double &f, double *par, int flag ){
  double chi2 = 0.0;

  //Parameters are:
  // (x0,y0,z0): coordinates of the GEM origin relative to the target center in spectrometer coordinates
  // (theta,phi): polar and azimuthal angles of the Z axis of the GEM stack relative to the spectrometer central ray:
  TVector3 GEMPOS(par[0], par[1], par[2] );
  
  //double ZSIEVE = par[3];
  double thetaGEM = par[3];
  double phiGEM = par[4]; 

  //Define unit vectors along the GEM internal axes with respect to the global system:
  TVector3 GEM_zaxis( sin(thetaGEM)*cos(phiGEM), sin(thetaGEM)*sin(phiGEM), cos(thetaGEM) );

  TVector3 Global_yaxis(0,1,0);
  //Not clear if we will need these:
  TVector3 Global_xaxis(1,0,0);
  TVector3 Global_zaxis(0,0,1); 
  
  TVector3 GEM_xaxis = (Global_yaxis.Cross( GEM_zaxis )).Unit();
  TVector3 GEM_yaxis = (GEM_zaxis.Cross( GEM_xaxis )).Unit();
  
  chi2 = 0.0;

  //The distance of closest approach between the track and the origin is found from:
  // sclose^2 = (rGEM + s * nhatGEM )^2
  // sclose^2 = rGEM^2 + 2 s * nhatGEM dot rGEM + s^2
  // dsclose^2/ds = 2 nhatGEM dot rGEM + 2s = 0
  // smin = -rGEM dot nhatGEM 
  
  for( int i=0; i<NTRACKS; i++ ){
    //Track position at sieve plate: XSIEVE and YSIEVE have to be determined by Holly's code:
    //TVector3 TrackSievePos( XSIEVE[i], YSIEVE[i], ZSIEVE );

    //Track direction 
    //TVector3 TrackDirGlobal = TrackSievePos.Unit();

    TVector3 TrackPosLocal( XTRACK[i], YTRACK[i], 0.0 );
    TVector3 TrackDirLocal( XPTRACK[i], YPTRACK[i], 1.0 );

    TrackDirLocal = TrackDirLocal.Unit();
    
    TVector3 TrackPosGlobal = GEMPOS + TrackPosLocal.X() * GEM_xaxis + TrackPosLocal.Y() * GEM_yaxis + TrackPosLocal.Z() * GEM_zaxis;
    TVector3 TrackDirGlobal = TrackDirLocal.X() * GEM_xaxis + TrackDirLocal.Y() * GEM_yaxis + TrackDirLocal.Z() * GEM_zaxis; 
    
    //Now compute track intersection with plane of first GEM:
    //Since all tracks start at origin, 
    //Equation is: ( s * TrackDir - GEMPOS ) dot GEM_zaxis = 0
    // OR: s = (GEMPOS dot GEM_zaxis)/(TrackDir dot GEM_zaxis)
    
    //  double sintersect = GEMPOS.Dot( GEM_zaxis )/(TrackDirGlobal.Dot( GEM_zaxis ) );
    double sintersect = -TrackPosGlobal.Dot( TrackDirGlobal ); 
    
    //TVector3 TrackIntersect_FirstGEMplane = sintersect * TrackDirGlobal;

    TVector3 TrackOriginProjection = TrackPosGlobal + sintersect * TrackDirGlobal;

    double sclose2 = TrackOriginProjection.Mag2();
    
    // double xfp_expect = (TrackIntersect_FirstGEMplane - GEMPOS).Dot( GEM_xaxis );
    // double yfp_expect = (TrackIntersect_FirstGEMplane - GEMPOS).Dot( GEM_yaxis );
    // double xpfp_expect = TrackDirGlobal.Dot( GEM_xaxis )/TrackDirGlobal.Dot( GEM_zaxis );
    // double ypfp_expect = TrackDirGlobal.Dot( GEM_yaxis )/TrackDirGlobal.Dot( GEM_zaxis );

    //chi2 += ( pow( xfp_expect - XTRACK[i], 2 ) + pow( yfp_expect - YTRACK[i], 2 ) ) / (sigma_pos*sigma_pos) +
    //( pow( xpfp_expect - XPTRACK[i], 2 ) + pow( ypfp_expect - YPTRACK[i], 2 ) ) / (sigma_slope*sigma_slope);
    
    //chi2 += sclose2/pow(0.002,2);
    double xclose = TrackOriginProjection.X();
    double yclose = TrackOriginProjection.Y();
    double zclose = TrackOriginProjection.Z();
    
    //Now we should also project tracks to the first GEM plane ASSUMING they start at the origin, and compare the predicted track parameters to the measured ones:
    TVector3 Origin(0,0,0);

    TVector3 TrackDirOriginFixed = (TrackPosGlobal - Origin).Unit();

    //Calculate intersection with first GEM plane under fixed origin assumption:
    // GEMzaxis dot ( Origin + Direction * s - GEMPOS ) = 0
    // --> s = (GEMPOS - Origin) dot GEMzaxi / Direction dot GEMzaxis
    sintersect = (GEMPOS - Origin).Dot( GEM_zaxis ) / (TrackDirOriginFixed.Dot( GEM_zaxis ) );

    TVector3 TrackIntersect_FirstGEMplane = Origin + sintersect * TrackDirOriginFixed;

    double xfp_expect = (TrackIntersect_FirstGEMplane - GEMPOS).Dot( GEM_xaxis );
    double yfp_expect = (TrackIntersect_FirstGEMplane - GEMPOS).Dot( GEM_yaxis );
    double xpfp_expect = TrackDirOriginFixed.Dot( GEM_xaxis )/TrackDirOriginFixed.Dot( GEM_zaxis );
    double ypfp_expect = TrackDirOriginFixed.Dot( GEM_yaxis )/TrackDirOriginFixed.Dot( GEM_zaxis );
    

    //chi2 += pow( xclose/sigma_xclose,2) + pow( yclose/sigma_yclose,2) + pow( zclose/sigma_zclose,2);
    //chi2 += sclose2 / pow( 0.007,2) ;

    chi2 += pow( (xfp_expect - XTRACK[i])/sigma_pos, 2 ) +
      pow( (yfp_expect - YTRACK[i])/sigma_pos, 2 ) +
      pow( (xpfp_expect - XPTRACK[i])/sigma_slope, 2 ) +
      pow( (ypfp_expect - YPTRACK[i])/sigma_slope, 2 );
  }

  //double penaltyterm = pow( (par[0]-GEMX0)/sigmaX0, 2 ) + pow( (par[1]-GEMY0)/sigmaY0, 2 ) + pow( (par[2]-GEMZ0)/sigmaZ0,2) + pow( (par[3]-GEMtheta)/sigmatheta,2) + pow( (par[4]-GEMphi)/sigmaphi,2);

  //chi2 += penaltyterm;
  
  f = chi2;
}

void AlignZeroField( const char *configfilename ){

  TChain *C = new TChain("T");
  
  ifstream configfile(configfilename);

  TString prefix = "bb.gem";

  //These are in spectrometer transport coordinates for an assumed magnet distance of 1.85 m
  double GEMX0=-0.1637, GEMY0=0.0, GEMZ0=2.958;
  double ZSIEVE=1.85-13.37*2.54/100.0; // = 1.51 m
  double XOFFSIEVE = 0.0;
  double YOFFSIEVE = 0.0; 
  double GEMtheta = 10.0*PI/180.0; //this will become the pitch angle
  double GEMphi = 180.0*PI/180.0; 

  TCut globalcut = "";

  double sclosemax = 0.1; 
  
  if( configfile ){
    TString currentline;
    
    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endlist") ){
      if( !currentline.BeginsWith("#") ){
	C->Add(currentline.Data());
      }
    }

    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endconfig")){
      if( !currentline.BeginsWith("#") ){
	TObjArray *tokens = currentline.Tokenize(" ");
	
	int ntokens = tokens->GetEntries();

	if( ntokens >= 2 ){
	  TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

	  if( skey == "prefix" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    prefix = stemp;
	  }

	  if( skey == "GEMX0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMX0 = stemp.Atof();
	  }
	  if( skey == "GEMY0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMY0 = stemp.Atof();
	  }

	  if( skey == "GEMZ0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMZ0 = stemp.Atof();
	  }

	  if( skey == "sigmaX0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sigmaX0 = stemp.Atof();
	  }

	  if( skey == "sigmaY0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sigmaY0 = stemp.Atof();
	  }

	  if( skey == "sigmaZ0" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sigmaZ0 = stemp.Atof();
	  }

	  if( skey == "sigmatheta" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sigmatheta = stemp.Atof() * PI/180.0;
	  }

	  if( skey == "sigmaphi" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sigmaphi = stemp.Atof() * PI/180.0;
	  }
	  // if( skey == "X0SIEVE" ){
	  //   TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   XOFFSIEVE = stemp.Atof();
	  // }

	  // if( skey == "Y0SIEVE" ){
	  //   TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   YOFFSIEVE = stemp.Atof();
	  // }

	  // if( skey == "ZSIEVE" ){
	  //   TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	  //   ZSIEVE = stemp.Atof();
	  // }

	  if( skey == "GEMtheta" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMtheta = stemp.Atof()*PI/180.0;
	  }

	  if( skey == "GEMphi" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    GEMphi = stemp.Atof()*PI/180.0;
	  }

	  if( skey == "sclosemax" ){
	    TString stemp = ( (TObjString*) (*tokens)[1] )->GetString();
	    sclosemax = stemp.Atof();
	  }
	  
	}
      }
    }

    while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
      if( !currentline.BeginsWith("#") ){
	globalcut += currentline;
      }
    }
  } else {
    return;
  }

  TEventList *elist = new TEventList("elist");
  
  C->Draw(">>elist",globalcut);
  
  cout << "Number of events passing global cut = " << elist->GetN() << endl;

  //TODO: need to declare variables and set tree branch addresses

  UInt_t MAXNTRACKS = 50000;
  
  double ntracks, besttrack;

  //We'll use the GEM track variables directly rather than the BigBite ones to avoid confusion initially.
  //The GEM track variables are guaranteed to be given in the internal GEM coordinates:
  vector<double> tracknhits(MAXNTRACKS);
  vector<double> trackChi2NDF(MAXNTRACKS);
  vector<double> trackX(MAXNTRACKS);
  vector<double> trackY(MAXNTRACKS);
  vector<double> trackXp(MAXNTRACKS);
  vector<double> trackYp(MAXNTRACKS); 

  C->SetBranchStatus("*",0);

  TString branchname;
  
  C->SetBranchStatus( branchname.Format("%s.track.ntrack",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.chi2ndf",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.x",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.y",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.xp",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.yp",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.nhits",prefix.Data() ), 1 );
  C->SetBranchStatus( branchname.Format("%s.track.besttrack",prefix.Data() ), 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  //main event loop here:
  C->SetBranchAddress( branchname.Format("%s.track.ntrack",prefix.Data() ), &ntracks );
  C->SetBranchAddress( branchname.Format("%s.track.chi2ndf",prefix.Data() ), &(trackChi2NDF[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.x",prefix.Data() ), &(trackX[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.y",prefix.Data() ), &(trackY[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.xp",prefix.Data() ), &(trackXp[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.yp",prefix.Data() ), &(trackYp[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.nhits",prefix.Data() ), &(tracknhits[0]) );
  C->SetBranchAddress( branchname.Format("%s.track.besttrack",prefix.Data() ), &besttrack );
  
  //TObjArray HList(0);
  TFile *fout = new TFile("ZeroFieldAlign_results_new.root","RECREATE");

  TH2D *hxyfp = new TH2D("hxyfp"," ; y_{fp} (m) ; x_{fp} (m)",100,-0.3,0.3,250,-1.1,1.1);
  // TH2D *hxysieve = new TH2D("hxysieve", " ; y_{sieve} (m) ; x_{sieve} (m)", 200,-0.2,0.2,200,-0.4,0.4); HList.Add(hxysieve);
  // TH2F *hYpFpYFp = new TH2F("hYpFpYFp"," ; Ypfp ; Yfp",100,-.3,.3,100,-0.3,0.3); HList.Add(hYpFpYFp);
  // TH2F *hXpFpXFp = new TH2F("hXpFpXFp"," ; Xpfp ; Xfp",100,-.7,.7,100,-0.7,0.7); HList.Add(hXpFpXFp);
  // TH2D *hxyfp_sel = new TH2D("hxyfp_sel","with sieve cuts; y_{fp} (m) ; x_{fp} (m)",100,-0.3,0.3,250,-1.1,1.1); HList.Add(hxyfp);
  // TH2D *hxysieve_sel = new TH2D("hxysieve_sel", "with sieve cuts; y_{sieve} (m) ; x_{sieve} (m)", 200,-0.4,0.4,200,-0.4,0.4); HList.Add(hxysieve);
  // TH2F *hYpFpYFp_sel = new TH2F("hYpFpYFp_sel","with sieve cuts; Ypfp ; Yfp",100,-.3,.3,100,-0.3,0.3); HList.Add(hYpFpYFp);
  // TH2F *hXpFpXFp_sel = new TH2F("hXpFpXFp_sel","with sieve cuts; Xpfp ; Xfp",100,-.7,.7,100,-0.7,0.7); HList.Add(hXpFpXFp);

  TH1D *hytar_old = new TH1D("hytar_old",";y_{tar} (m);", 250,-0.25,0.25);
  TH1D *hxtar_old = new TH1D("hxtar_old",";x_{tar} (m);", 250,-0.25,0.25);

  TH1D *hsclose_old = new TH1D("hsclose_old",";s_{close} (m);", 250, 0.0, 0.1 );
  TH1D *hsclose_new = new TH1D("hsclose_new",";s_{close} (m);", 250, 0.0, 0.1 );

  TH1D *hxclose_old = new TH1D("hxclose_old",";x_{close} (m);", 250, -0.1, 0.1 );
  TH1D *hyclose_old = new TH1D("hyclose_old", ";y_{close} (m);", 250, -0.1, 0.1 );
  TH1D *hzclose_old = new TH1D("hzclose_old", ";z_{close} (m);", 250, -0.02, 0.02 );

  TH1D *hxclose_new = new TH1D("hxclose_new",";x_{close} (m);", 250, -0.1, 0.1 );
  TH1D *hyclose_new = new TH1D("hyclose_new", ";y_{close} (m);", 250, -0.1, 0.1 );
  TH1D *hzclose_new = new TH1D("hzclose_new", ";z_{close} (m);", 250, -0.02, 0.02 );
  
  long nevent=0;

  TVector3 Global_yaxis(0,1,0);
  //Not clear if we will need these:
  TVector3 Global_xaxis(1,0,0);
  TVector3 Global_zaxis(0,0,1);

  TVector3 GEMzaxis( sin(GEMtheta)*cos(GEMphi), sin(GEMtheta)*sin(GEMphi), cos(GEMtheta) );
  TVector3 GEMxaxis = (Global_yaxis.Cross(GEMzaxis)).Unit();
  TVector3 GEMyaxis = (GEMzaxis.Cross(GEMxaxis)).Unit();
  TVector3 GEMorigin( GEMX0, GEMY0, GEMZ0 );
  //TVector3 SieveOrigin(0,0,ZSIEVE);

  cout << endl << "GEM z axis in global coordinates:" << endl;
  GEMzaxis.Print();

  cout << endl << "GEM x axis in global coordinates:" << endl;
  GEMxaxis.Print();

  cout << endl << "GEM y axis in global coordinates:" << endl;
  GEMyaxis.Print();

  cout << endl << "GEM origin in global coordinates:" << endl;
  GEMorigin.Print();

  // cout << endl << "Sieve slit origin in global coordinates: " << endl;
  // SieveOrigin.Print();

  // TString outCutFile;
  // outCutFile=Form("sieve_cut.root");
  // cout << "outCutFile =  " << outCutFile << endl;
  // TFile fcut(outCutFile,"RECREATE");
  // TString hname_cut;

  //Here are the desired sieve positions:
  // vector <Double_t> xs_cent{-(0.3+0.0492)+0.0493/cos(18.*PI/180.),
  //     -(0.3+0.0492)+(0.0493+0.0492)/cos(18.*PI/180.),
  //     -(0.3+0.0492)+0.1493/cos(9.*PI/180.),
  //     -(0.3+0.0492)+(0.1493+0.0492)/cos(9.*PI/180.),
  //     -(0.3+0.0492)+(0.1493+0.0492*2.)/cos(9.*PI/180.),
  //     -0.0492,
  //     0.0,
  //     0.0492,
  //     0.3+0.0492-(0.1493+0.0492*2.)/cos(9.*PI/180.),
  //     0.3+0.0492-(0.1493+0.0492)/cos(9.*PI/180.),
  //     0.3+0.0492-0.1493/cos(9.*PI/180.),
  //     0.3+0.0492-(0.0493+0.0492)/cos(18.*PI/180),
  //     0.3+0.0492-0.0493/cos(18.*PI/180.)};
   
  //Add in any XSIEVE offset defined by the user:
  // for( int ix=0; ix<nxsieve; ix++ ){
  //   xs_cent[ix] += XOFFSIEVE;
  // }

  // vector <Double_t> ys_cent;
  // for (Int_t nys=0;nys<nysieve;nys++) {
  //   Double_t pos=nys*0.0381-0.0381*3;//old sieve
  //   ys_cent.push_back(pos + YOFFSIEVE);
  // }
  
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){

    if( nevent % 1000 == 0 ) cout << "nevent = " << nevent << endl;
    
    //do stuff: grab track info from tree, apply cuts, fill the track arrays defined near the top of this macro (see XTRACK, YTRACK, XPTRACK, YPTRACK, XSIEVE, YSIEVE above)
    int itr = int(besttrack);
    if( besttrack >= 0 && ntracks <= MAXNTRACKS ){
      hxyfp->Fill( trackY[itr], trackX[itr] );
      // hYpFpYFp->Fill(trackYp[itr],trackY[itr]);
      // hXpFpXFp->Fill(trackXp[itr],trackX[itr]);
  
      TVector3 TrackPos_local( trackX[itr], trackY[itr], 0.0 );
      TVector3 TrackDir_local( trackXp[itr], trackYp[itr], 1.0 );
      TrackDir_local = TrackDir_local.Unit();
      
      TVector3 TrackDir_global = TrackDir_local.X() * GEMxaxis +
	TrackDir_local.Y() * GEMyaxis +
	TrackDir_local.Z() * GEMzaxis;

      TVector3 TrackPos_global = GEMorigin + TrackPos_local.X() * GEMxaxis + TrackPos_local.Y() * GEMyaxis + TrackPos_local.Z() * GEMzaxis; 
      
      //Now compute intersection of track ray in global coordinates with sieve slit:
      // (trackpos + s * trackdir - sieveorigin) dot globalzaxis = 0
      // --> s * trackdir dot globalzaxis = (sieveorigin - trackpos) dot globalzaxis
      // double sintersect = (SieveOrigin - TrackPos_global).Dot( Global_zaxis ) /
      // 	( TrackDir_global.Dot( Global_zaxis ) );
      // TVector3 TrackSievePos = TrackPos_global + sintersect * TrackDir_global;

      // hxysieve->Fill( TrackSievePos.Y(), TrackSievePos.X() );
      //cout<<"ysieve: "<<TrackSievePos.Y()<<" xsieve: "<<TrackSievePos.X()<<endl;
      double sintersect_target = -TrackPos_global.Dot( Global_zaxis ) / (TrackDir_global.Dot( Global_zaxis ) );

      TVector3 TrackTargPos = TrackPos_global + sintersect_target * TrackDir_global;

      //Compute closest approach to origin (apply a loose cut): 
      TVector3 TrackClosestApproach = TrackPos_global - TrackPos_global.Dot( TrackDir_global ) * TrackDir_global;
      
      hytar_old->Fill( TrackTargPos.Y() );
      hxtar_old->Fill( TrackTargPos.X() );

      double sclose = TrackClosestApproach.Mag();

      cout << "sclose = " << sclose << endl;
      
      hsclose_old->Fill( sclose );

      hxclose_old->Fill( TrackClosestApproach.X() );
      hyclose_old->Fill( TrackClosestApproach.Y() );
      hzclose_old->Fill( TrackClosestApproach.Z() );

      if( sclose <= sclosemax ){
      
	XTRACK.push_back(trackX[itr]);
	YTRACK.push_back(trackY[itr]);
	XPTRACK.push_back(trackXp[itr]);
	YPTRACK.push_back(trackYp[itr]);
	//XSIEVE.push_back(xs_cent[ix_found]);
	//YSIEVE.push_back(ys_cent[iy_found]);
	NTRACKS++;
      }
      
    }
  }//end while event

  //do the thing
  TMinuit *FitZeroField = new TMinuit( 5 );
  
  FitZeroField->SetFCN( CHI2_FCN );
  
  int ierflg = 0;
  
  FitZeroField->mnparm( 0, "GEMX0", GEMX0, 0.01,0,0,ierflg ); //guesstimate 1 cm as initial position accuracy
  FitZeroField->mnparm( 1, "GEMY0", GEMY0, 0.03,0,0,ierflg ); //start with 3-cm uncertainty for y position
  FitZeroField->mnparm( 2, "GEMZ0", GEMZ0, 0.01,0,0,ierflg );
  FitZeroField->mnparm( 3, "GEMtheta", GEMtheta, 0.3*PI/180.0, 0, 0, ierflg ); //guesstimate 0.3 degrees as initial angular accuracy
  FitZeroField->mnparm( 4, "GEMphi", GEMphi, 0.3*PI/180.0, 0, 0, ierflg );
    
  double arglist[10];
  arglist[0]=1;
  FitZeroField->mnexcm("SET ERR",arglist,1,ierflg);
  
  arglist[0] = 5000;
  arglist[1] = 1.;

  //Test whether function changes parameter values:
  cout << "BEFORE fit: " << endl;
  cout << "GEM (X0, Y0, Z0) = (" << GEMX0 << ", " << GEMY0 << ", " << GEMZ0 << ")" << endl;
  cout << "GEM (theta, phi) = (" << GEMtheta * 180.0/PI << ", " << GEMphi * 180.0/PI << ")" << endl;
  
  //Uncomment this line when we are actually set up to do a fit:   
  FitZeroField->mnexcm("MIGRAD",arglist,2,ierflg);

  cout << "AFTER fit: " << endl;
  cout << "GEM (X0, Y0, Z0) = (" << GEMX0 << ", " << GEMY0 << ", " << GEMZ0 << ")" << endl;
  cout << "GEM (theta, phi) = (" << GEMtheta * 180.0/PI << ", " << GEMphi * 180.0/PI << ")" << endl;
  
  //AFTER the fit, 
  
  //TODO: grab parameters, write them out to file. Profit.
  double fitpar[6];
  double fitparerr[6];

  ofstream outfile("ZeroFieldAlign_results_new.txt");

  TString outline;
  
  for (int ii=0;ii<6;ii++){
  //fitpar[ii]=FitZeroField->GetParameter(ii);fitparerr[ii]=FitZeroField->GetParameter(ii);?????
    FitZeroField->GetParameter(ii,fitpar[ii],fitparerr[ii]);
  }
  outline.Form("GEMX0 = %18.9g +/- %18.9g",fitpar[0], fitparerr[0]);
  outfile << outline << endl;
  cout << outline << endl;

 outline.Form("GEMY0 = %18.9g +/- %18.9g",fitpar[1], fitparerr[1]);
  outfile << outline << endl;
  cout << outline << endl;

  outline.Form("GEMZ0 = %18.9g +/- %18.9g",fitpar[2], fitparerr[2]);
  outfile << outline << endl;
  cout << outline << endl;

  outline.Form("GEMtheta = %18.9g +/- %18.9g",fitpar[3]*TMath::RadToDeg(), fitparerr[3]*TMath::RadToDeg());
  outfile << outline << endl;
  cout << outline << endl;

  outline.Form("GEMphi = %18.9g +/- %18.9g",fitpar[4]*TMath::RadToDeg(), fitparerr[4]*TMath::RadToDeg());
  outfile << outline << endl;
  cout << outline << endl;

  GEMtheta = fitpar[3];
  GEMphi = fitpar[4];
  GEMX0 = fitpar[0];
  GEMY0 = fitpar[1];
  GEMZ0 = fitpar[2];

  GEMzaxis.SetXYZ( sin(GEMtheta)*cos(GEMphi), sin(GEMtheta)*sin(GEMphi), cos(GEMtheta) );
  GEMxaxis = (Global_yaxis.Cross(GEMzaxis)).Unit();
  GEMyaxis = (GEMzaxis.Cross(GEMxaxis)).Unit();
  GEMorigin.SetXYZ( GEMX0, GEMY0, GEMZ0 );
  
  for( int i=0; i<NTRACKS; i++ ){
    //Track position at sieve plate: XSIEVE and YSIEVE have to be determined by Holly's code:
    //TVector3 TrackSievePos( XSIEVE[i], YSIEVE[i], ZSIEVE );
    
    //Track direction 
    //TVector3 TrackDirGlobal = TrackSievePos.Unit();
    
    TVector3 TrackPosLocal( XTRACK[i], YTRACK[i], 0.0 );
    TVector3 TrackDirLocal( XPTRACK[i], YPTRACK[i], 1.0 );
    
    TrackDirLocal = TrackDirLocal.Unit();
    
    TVector3 TrackPosGlobal = GEMorigin + TrackPosLocal.X() * GEMxaxis + TrackPosLocal.Y() * GEMyaxis + TrackPosLocal.Z() * GEMzaxis;
    TVector3 TrackDirGlobal = TrackDirLocal.X() * GEMxaxis + TrackDirLocal.Y() * GEMyaxis + TrackDirLocal.Z() * GEMzaxis; 
    
    //Now compute track intersection with plane of first GEM:
    //Since all tracks start at origin, 
    //Equation is: ( s * TrackDir - GEMPOS ) dot GEM_zaxis = 0
    // OR: s = (GEMPOS dot GEM_zaxis)/(TrackDir dot GEM_zaxis)
    
    //  double sintersect = GEMPOS.Dot( GEM_zaxis )/(TrackDirGlobal.Dot( GEM_zaxis ) );
    double sintersect = -TrackPosGlobal.Dot( TrackDirGlobal ); 
    
    //TVector3 TrackIntersect_FirstGEMplane = sintersect * TrackDirGlobal;

    TVector3 TrackOriginProjection = TrackPosGlobal + sintersect * TrackDirGlobal;

    double sclose2 = TrackOriginProjection.Mag2();

    hsclose_new->Fill( TrackOriginProjection.Mag() );

    hxclose_new->Fill( TrackOriginProjection.X() );
    hyclose_new->Fill( TrackOriginProjection.Y() );
    hzclose_new->Fill( TrackOriginProjection.Z() );
  }
  
  fout->Write();
 
}//end program
