#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
//#include "TEventList.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include <iostream>
#include "TTreeFormula.h"

void ElasticQuickAndDirty( const char *rootfilename, double Ebeam=4.291, double bbtheta=29.5, double sbstheta=34.7, double sbsdist=2.8, double hcaldist=17.0, double sbsfieldscale=1.0){

  sbstheta *= TMath::Pi()/180.0;
  bbtheta *= TMath::Pi()/180.0;

  TChain *C = new TChain("T");

  C->Add(rootfilename); 

  TCut globalcut = "bb.ps.e>0.2&&sbs.hcal.nclus>0&&bb.tr.n>=1&&abs(bb.etot_over_p-1.)<0.25";
  
  //TEventList *elist = new TEventList("elist","");
  
  //  C->Draw(">>elist",globalcut);

  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",globalcut,C);

  int MAXNTRACKS=10;
  
  //variables we need are BigBite track px,py,pz and sbs hcal x, y, e

  double ntrack;

  double epx[MAXNTRACKS];
  double epy[MAXNTRACKS];
  double epz[MAXNTRACKS];
  double ep[MAXNTRACKS];

  double vx[MAXNTRACKS];
  double vy[MAXNTRACKS];
  double vz[MAXNTRACKS];

  int MAXTRIGTDCHITS = 10;
  
  double bb_tdctrig_tdc[MAXTRIGTDCHITS];
  double bb_tdctrig_elemID[MAXTRIGTDCHITS];
  
  double xhcal,yhcal,ehcal, EPS, ESH, atimehcal, atimebbsh, atimebbps;

  C->SetBranchStatus("*",0);

  C->SetBranchStatus("bb.tdctrig.tdc",1);
  C->SetBranchStatus("bb.tdctrig.tdcelemID",1);
  
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);

  C->SetBranchStatus("bb.sh.atimeblk",1);
  C->SetBranchStatus("bb.ps.atimeblk",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.etot_over_p",1);

  C->SetBranchStatus("sbs.hcal.nclus",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);

  C->SetBranchAddress("bb.tdctrig.tdc",bb_tdctrig_tdc);
  C->SetBranchAddress("bb.tdctrig.tdcelemID",bb_tdctrig_elemID);
  
  C->SetBranchAddress("bb.tr.n",&ntrack);
  C->SetBranchAddress("bb.tr.vz",vz);
  C->SetBranchAddress("bb.tr.px",epx);
  C->SetBranchAddress("bb.tr.py",epy);
  C->SetBranchAddress("bb.tr.pz",epz);
  C->SetBranchAddress("bb.tr.p",ep);
  C->SetBranchAddress("sbs.hcal.x",&xhcal);
  C->SetBranchAddress("sbs.hcal.y",&yhcal);
  C->SetBranchAddress("sbs.hcal.e",&ehcal);

  TLorentzVector Pbeam(0,0,Ebeam,Ebeam);
  TLorentzVector Ptarg(0,0,0,0.5*(0.938272+0.939565));

  double Mp = 0.938272;

  long nevent=0;

  double W2min = 0.88-0.4;
  double W2max = 0.88+0.4;

  TFile *fout = new TFile("elastic_temp.root","RECREATE");

  TH2D *hdxdy_all = new TH2D("hdxdy_all","All events;#Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wcut = new TH2D("hdxdy_Wcut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wanticut = new TH2D("hdxdy_Wanticut","|W^{2}-0.88|<0.4;Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);

  TH1D *hW2_all = new TH1D("hW2_all","All events; W^{2} (GeV^{2});", 250, -1, 4 );

  TH1D *hEHCAL_all = new TH1D("hEHCAL_all","All events; HCAL energy sum (GeV);",250,0,0.5);
  TH1D *hEHCAL_Wcut = new TH1D("hEHCAL_Wcut","|W^{2}-0.88|<0.4;HCAL energy sum (GeV);",250,0,0.5);

  TH2D *hdxdy_all_anglesonly = new TH2D("hdxdy_all_anglesonly","All events; #Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  TH2D *hdxdy_Wcut_anglesonly = new TH2D("hdxdy_Wcut_anglesonly","|W^{2}-0.88|<0.4; #Deltay (m);#Deltax (m)",125,-2,2,125,-4,6);
  
  TVector3 hcal_origin( -hcaldist*sin(sbstheta), 0, hcaldist*cos(sbstheta) );

  TVector3 hcal_zaxis = hcal_origin.Unit();
  TVector3 hcal_xaxis(0,-1,0);
  TVector3 hcal_yaxis = hcal_zaxis.Cross( hcal_xaxis ).Unit();

  TTree *Tout = new TTree("Tout", "BB-HCAL elastic selection");

  double dx_4vect, dy_4vect, dx_angles, dy_angles;
  double W2_out;
  double Q2_out;
  double xHCAL_out, yHCAL_out;
  double ep_out, etheta_out, ephi_out;
  double EHCAL_out;
  double vz_out;
  
  Tout->Branch( "W2", &W2_out, "W2/D");
  Tout->Branch( "Q2", &Q2_out, "Q2/D");
  Tout->Branch( "dx_4vect", &dx_4vect, "dx_4vect/D");
  Tout->Branch( "dy_4vect", &dy_4vect, "dy_4vect/D");
  Tout->Branch( "dx", &dx_angles, "dx/D" );
  Tout->Branch( "dy", &dy_angles, "dy/D" );
  Tout->Branch( "xHCAL", &xHCAL_out, "xHCAL/D");
  Tout->Branch( "yHCAL", &yHCAL_out, "yHCAL/D");
  Tout->Branch( "ep", &ep_out, "ep/D");
  Tout->Branch( "EHCAL", &EHCAL_out, "EHCAL/D");
  Tout->Branch( "etheta", &etheta_out, "etheta/D");
  Tout->Branch( "ephi", &ephi_out, "ephi/D");
  Tout->Branch( "vz", &vz_out, "vz/D");

  int treenum=0, currenttreenum=0;
  
  while( C->GetEntry(nevent++) ){
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    }

    if( nevent % 10000 == 0 ) cout << nevent << endl;

    bool goodevent = GlobalCut->EvalInstance(0) != 0;

    //cout << "goodevent = " << goodevent << endl;

    if( ntrack >= 1.0 && goodevent ){
      TLorentzVector kprime( epx[0], epy[0], epz[0], ep[0] );
      TLorentzVector q = Pbeam - kprime;

      TVector3 qdir = q.Vect().Unit();

      TVector3 vertex(0,0,vz[0]);

      double sintersect = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir.Dot( hcal_zaxis );

      TVector3 hcal_intersect = vertex + sintersect * qdir; 

      double xhcal_expect = hcal_intersect.Dot( hcal_xaxis );
      double yhcal_expect = hcal_intersect.Dot( hcal_yaxis );

      double etheta = acos( kprime.Pz()/kprime.P() );
      double ephi = atan2( kprime.Py(), kprime.Px() );

      double Eprime_eth = Ebeam/ (1.0 + Ebeam/Mp * (1.0 - cos(etheta) ) );
      double nu = Ebeam - Eprime_eth;
      double pp_eth = sqrt(pow(nu,2)+2.*Mp*nu);

      double pth_eth = acos( (Ebeam-Eprime_eth*cos(etheta))/pp_eth );
      double pphi_eph = ephi + TMath::Pi();
      
      TVector3 qdir_angles( sin(pth_eth)*cos(pphi_eph), sin(pth_eth)*sin(pphi_eph), cos(pth_eth) );

      double sint_angle = (hcal_origin-vertex).Dot( hcal_zaxis )/qdir_angles.Dot( hcal_zaxis );

      TVector3 hcal_intersect_angles = vertex + sintersect * qdir_angles;

      double xhcal_expect_angles = hcal_intersect_angles.Dot( hcal_xaxis );
      double yhcal_expect_angles = hcal_intersect_angles.Dot( hcal_yaxis );
      
      hdxdy_all->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
      double W2recon = (Ptarg + q).M2();

      hW2_all->Fill( W2recon );

      hEHCAL_all->Fill( ehcal );

      W2_out = W2recon;
      Q2_out = -q.M2();

      dx_4vect = xhcal - xhcal_expect;
      dy_4vect = yhcal - yhcal_expect;
      dx_angles = xhcal - xhcal_expect_angles;
      dy_angles = yhcal - yhcal_expect_angles;
      
      if( W2recon > W2min && W2recon < W2max ){
	hdxdy_Wcut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
	hEHCAL_Wcut->Fill( ehcal );
      } else {
	hdxdy_Wanticut->Fill( yhcal - yhcal_expect, xhcal - xhcal_expect );
      }

      xHCAL_out = xhcal;
      yHCAL_out = yhcal;

      ep_out = ep[0];
      etheta_out = etheta;
      ephi_out = ephi;

      EHCAL_out = ehcal;
      vz_out = vz[0];
      
      Tout->Fill();
    }
  }

  cout << "Event loop finished" << endl;
  
  //  elist->Delete(); 
  fout->Write();

  cout << "Output file closed" << endl;

  GlobalCut->Delete();
} 
