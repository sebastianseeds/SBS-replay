/*
  Script to plot various correlation plots involving BBCAL.
  ------
  P. Datta <pdbforce@jlab.org> Created 10/25/22 
*/
#include "overlay.C"

void bbcal_corr_plots() {

  TPad *padptr = (TPad*) gPad;
  padptr->Divide(3,2);

  // plot E/p
  padptr->cd(1);
  overlay("hetot_over_p", "hEoverP_electron", "E/p (all tracks)", "E/p (electrons)");

  // plot total cluster energy
  padptr->cd(2);
  TH1F* hHEclus_cor_sh_ps; hHEclus_cor_sh_ps = (TH1F*)gDirectory->Get("hHEclus_cor_sh_ps");
  hHEclus_cor_sh_ps->Draw("colz");

  // plot sh counts
  padptr->cd(3);
  TH1F* hps_shPlusps_e; hps_shPlusps_e = (TH1F*)gDirectory->Get("hps_shPlusps_e");
  hps_shPlusps_e->Draw("colz");

  // plot coincidence time (HCAL - BBCAL)
  TH1F *h_coin = new TH1F("h_bbh_trigT_diff","HCal-BBCal Trigger Time Diff. (ns)",200,320,720);
  // TH2F *h2_bbcal_hcal_corr = new TH2F("h2_bbh_corr","BBCal-HCal Cluster Correlation; BB Shower Rows; HCal Rows",27,1,28,24,1,25);
  TH2F *h2_trigsum_shrow = new TH2F("h2_trigsum_shrow","BBCal Trigger Sums vs Shower rows; Shower rows; Trigger Sums",27,0,27,25,0,25);
  TH2F *h2_trigsum_shid = new TH2F("h2_trigsum_shid","BBCal Trigger Sums ADC vs Shower ID; Shower ID; Trigger Sums (ADC)",189,0,189,25,0,25);
  TH2F *h2_trigsumTDC_shid = new TH2F("h2_trigsum_shid","BBCal Trigger Sums TDC vs Shower ID; Shower ID; Trigger Sums (TDC)",189,0,189,25,0,25);

  // Declare trees
  TTree *T = (TTree*) gDirectory->Get("T");
  T->SetBranchStatus("*",0);
  Double_t nclusSH;      T->SetBranchStatus("bb.sh.nclus",1);     T->SetBranchAddress("bb.sh.nclus", &nclusSH);
  Double_t rblkSH;       T->SetBranchStatus("bb.sh.rowblk",1);    T->SetBranchAddress("bb.sh.rowblk", &rblkSH);
  Double_t idblkSH;      T->SetBranchStatus("bb.sh.idblk",1);     T->SetBranchAddress("bb.sh.idblk", &idblkSH);
  Double_t nclusHCAL;    T->SetBranchStatus("sbs.hcal.nclus",1);  T->SetBranchAddress("sbs.hcal.nclus", &nclusHCAL);
  Double_t rblkHCAL;     T->SetBranchStatus("sbs.hcal.rowblk",1); T->SetBranchAddress("sbs.hcal.rowblk", &rblkHCAL);
  Double_t bbtrigelID[25]; T->SetBranchStatus("bb.bbtrig.adcelemID",1);  T->SetBranchAddress("bb.bbtrig.adcelemID", &bbtrigelID);
  Int_t NbbtrigelID;     T->SetBranchStatus("Ndata.bb.bbtrig.adcelemID",1);  T->SetBranchAddress("Ndata.bb.bbtrig.adcelemID", &NbbtrigelID);
  Double_t BBtrigsumTDC[25]; T->SetBranchStatus("bb.bbtrig.tdc",1);  T->SetBranchAddress("bb.bbtrig.tdc", &BBtrigsumTDC);
  Int_t NBBtrigsumTDCelemID;  T->SetBranchStatus("Ndata.bb.bbtrig.tdcelemID",1); T->SetBranchAddress("Ndata.bb.bbtrig.tdcelemID", &NBBtrigsumTDCelemID);
  Double_t BBtrigsumTDCelemID[25]; T->SetBranchStatus("bb.bbtrig.tdcelemID",1);   T->SetBranchAddress("bb.bbtrig.tdcelemID", &BBtrigsumTDCelemID);
  Double_t tdcBBtrig[6]; T->SetBranchStatus("bb.tdctrig.tdc",1);  T->SetBranchAddress("bb.tdctrig.tdc", &tdcBBtrig);
  Int_t NtdcelemBBtrig;  T->SetBranchStatus("Ndata.bb.tdctrig.tdcelemID",1); T->SetBranchAddress("Ndata.bb.tdctrig.tdcelemID", &NtdcelemBBtrig);
  Double_t tdcelemBBtrig[6]; T->SetBranchStatus("bb.tdctrig.tdcelemID",1);   T->SetBranchAddress("bb.tdctrig.tdcelemID", &tdcelemBBtrig);

  // Loop through events
  Long64_t nevents = T->GetEntries();
  for (Long64_t nevent=0; nevent<nevents; nevent++){
    T->GetEntry(nevent);

    if(nclusSH==0 || nclusHCAL==0) continue;

    Double_t bbcal_time=0., hcal_time=0.;
    for (Int_t ihit=0; ihit<NtdcelemBBtrig; ihit++) {
      if(tdcelemBBtrig[ihit]==5) bbcal_time=tdcBBtrig[ihit];
      if(tdcelemBBtrig[ihit]==0) hcal_time=tdcBBtrig[ihit];
    }

    Double_t diff = hcal_time - bbcal_time; 
    h_coin->Fill(diff);

    // trigger sum vs sh row
    if(fabs(diff-510.)<20.){
      //h2_bbcal_hcal_corr->Fill(rblkSH+1, rblkHCAL+1);`
      for (Int_t ihit=0; ihit<NbbtrigelID; ihit++) {
	h2_trigsum_shrow->Fill(rblkSH, bbtrigelID[ihit]);
	h2_trigsum_shid->Fill(idblkSH, bbtrigelID[ihit]);
	h2_trigsumTDC_shid->Fill(idblkSH, BBtrigsumTDCelemID[ihit]);
      }
    }
  }
  
  // h2_bbcal_hcal_corr->SetMinimum(-0.1);
  // h2_bbcal_hcal_corr->SetStats(0);
  // h2_bbcal_hcal_corr->Draw("colz");

  padptr->cd(4);
  //gStyle->SetOptStat(1100);
  h_coin->Draw();
  double maxcount = h_coin->GetBinContent(h_coin->GetMaximumBin());
  // draw lines to indicate the safe region of coincidence peak
  TLine lb, lf;
  lb.SetLineColor(2); lb.SetLineWidth(2); lb.SetLineStyle(9); lb.DrawLine(460,0,460,maxcount);
  lf.SetLineColor(2); lf.SetLineWidth(2); lf.SetLineStyle(9); lf.DrawLine(580,0,580,maxcount);

  padptr->cd(5);
  gStyle->SetOptStat(0);
  // h2_trigsum_shrow->SetMinimum(-0.1);
  // h2_trigsum_shrow->Draw("colz");
  h2_trigsum_shid->SetMinimum(-0.1);
  h2_trigsum_shid->Draw("colz");

  padptr->cd(6);
  gStyle->SetOptStat(0);
  h2_trigsumTDC_shid->SetMinimum(-0.1);
  h2_trigsumTDC_shid->Draw("colz");
}
