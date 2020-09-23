//#define gmn_dig_tree_cxx
#include "gmn_dig_tree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

gmn_dig_tree::gmn_dig_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("digitized/simdigtest_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("digitized/simdigtest_0.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gmn_dig_tree::~gmn_dig_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gmn_dig_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gmn_dig_tree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gmn_dig_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Earm_BBGEM_hit_plane = 0;
   Earm_BBGEM_hit_strip = 0;
   Earm_BBGEM_hit_x = 0;
   Earm_BBGEM_hit_y = 0;
   Earm_BBGEM_hit_z = 0;
   Earm_BBGEM_hit_polx = 0;
   Earm_BBGEM_hit_poly = 0;
   Earm_BBGEM_hit_polz = 0;
   Earm_BBGEM_hit_t = 0;
   Earm_BBGEM_hit_trms = 0;
   Earm_BBGEM_hit_tmin = 0;
   Earm_BBGEM_hit_tmax = 0;
   Earm_BBGEM_hit_tx = 0;
   Earm_BBGEM_hit_ty = 0;
   Earm_BBGEM_hit_xin = 0;
   Earm_BBGEM_hit_yin = 0;
   Earm_BBGEM_hit_zin = 0;
   Earm_BBGEM_hit_xout = 0;
   Earm_BBGEM_hit_yout = 0;
   Earm_BBGEM_hit_zout = 0;
   Earm_BBGEM_hit_txp = 0;
   Earm_BBGEM_hit_typ = 0;
   Earm_BBGEM_hit_xg = 0;
   Earm_BBGEM_hit_yg = 0;
   Earm_BBGEM_hit_zg = 0;
   Earm_BBGEM_hit_trid = 0;
   Earm_BBGEM_hit_mid = 0;
   Earm_BBGEM_hit_pid = 0;
   Earm_BBGEM_hit_vx = 0;
   Earm_BBGEM_hit_vy = 0;
   Earm_BBGEM_hit_vz = 0;
   Earm_BBGEM_hit_p = 0;
   Earm_BBGEM_hit_edep = 0;
   Earm_BBGEM_hit_beta = 0;
   Earm_BBGEM_hit_otridx = 0;
   Earm_BBGEM_hit_ptridx = 0;
   Earm_BBGEM_hit_sdtridx = 0;
   Earm_BBGEM_Track_TID = 0;
   Earm_BBGEM_Track_PID = 0;
   Earm_BBGEM_Track_MID = 0;
   Earm_BBGEM_Track_NumHits = 0;
   Earm_BBGEM_Track_NumPlanes = 0;
   Earm_BBGEM_Track_NDF = 0;
   Earm_BBGEM_Track_Chi2fit = 0;
   Earm_BBGEM_Track_Chi2true = 0;
   Earm_BBGEM_Track_X = 0;
   Earm_BBGEM_Track_Y = 0;
   Earm_BBGEM_Track_Xp = 0;
   Earm_BBGEM_Track_Yp = 0;
   Earm_BBGEM_Track_T = 0;
   Earm_BBGEM_Track_P = 0;
   Earm_BBGEM_Track_Sx = 0;
   Earm_BBGEM_Track_Sy = 0;
   Earm_BBGEM_Track_Sz = 0;
   Earm_BBGEM_Track_Xfit = 0;
   Earm_BBGEM_Track_Yfit = 0;
   Earm_BBGEM_Track_Xpfit = 0;
   Earm_BBGEM_Track_Ypfit = 0;
   Earm_BBGEM_Track_otridx = 0;
   Earm_BBGEM_Track_ptridx = 0;
   Earm_BBGEM_Track_sdtridx = 0;
   Earm_BBHodoScint_hit_row = 0;
   Earm_BBHodoScint_hit_col = 0;
   Earm_BBHodoScint_hit_cell = 0;
   Earm_BBHodoScint_hit_plane = 0;
   Earm_BBHodoScint_hit_wire = 0;
   Earm_BBHodoScint_hit_xcell = 0;
   Earm_BBHodoScint_hit_ycell = 0;
   Earm_BBHodoScint_hit_zcell = 0;
   Earm_BBHodoScint_hit_xcellg = 0;
   Earm_BBHodoScint_hit_ycellg = 0;
   Earm_BBHodoScint_hit_zcellg = 0;
   Earm_BBHodoScint_hit_xhit = 0;
   Earm_BBHodoScint_hit_yhit = 0;
   Earm_BBHodoScint_hit_zhit = 0;
   Earm_BBHodoScint_hit_xhitg = 0;
   Earm_BBHodoScint_hit_yhitg = 0;
   Earm_BBHodoScint_hit_zhitg = 0;
   Earm_BBHodoScint_hit_sumedep = 0;
   Earm_BBHodoScint_hit_tavg = 0;
   Earm_BBHodoScint_hit_trms = 0;
   Earm_BBHodoScint_hit_tmin = 0;
   Earm_BBHodoScint_hit_tmax = 0;
   Earm_BBHodoScint_hit_otridx = 0;
   Earm_BBHodoScint_hit_ptridx = 0;
   Earm_BBHodoScint_hit_sdtridx = 0;
   Earm_BBPS_hit_PMT = 0;
   Earm_BBPS_hit_row = 0;
   Earm_BBPS_hit_col = 0;
   Earm_BBPS_hit_plane = 0;
   Earm_BBPS_hit_xcell = 0;
   Earm_BBPS_hit_ycell = 0;
   Earm_BBPS_hit_zcell = 0;
   Earm_BBPS_hit_xgcell = 0;
   Earm_BBPS_hit_ygcell = 0;
   Earm_BBPS_hit_zgcell = 0;
   Earm_BBPS_hit_NumPhotoelectrons = 0;
   Earm_BBPS_hit_Time_avg = 0;
   Earm_BBPS_hit_Time_rms = 0;
   Earm_BBPS_hit_Time_min = 0;
   Earm_BBPS_hit_Time_max = 0;
   Earm_BBPS_hit_otridx = 0;
   Earm_BBPS_hit_ptridx = 0;
   Earm_BBPS_hit_sdtridx = 0;
   Earm_BBPSTF1_hit_row = 0;
   Earm_BBPSTF1_hit_col = 0;
   Earm_BBPSTF1_hit_cell = 0;
   Earm_BBPSTF1_hit_plane = 0;
   Earm_BBPSTF1_hit_wire = 0;
   Earm_BBPSTF1_hit_xcell = 0;
   Earm_BBPSTF1_hit_ycell = 0;
   Earm_BBPSTF1_hit_zcell = 0;
   Earm_BBPSTF1_hit_xcellg = 0;
   Earm_BBPSTF1_hit_ycellg = 0;
   Earm_BBPSTF1_hit_zcellg = 0;
   Earm_BBPSTF1_hit_xhit = 0;
   Earm_BBPSTF1_hit_yhit = 0;
   Earm_BBPSTF1_hit_zhit = 0;
   Earm_BBPSTF1_hit_xhitg = 0;
   Earm_BBPSTF1_hit_yhitg = 0;
   Earm_BBPSTF1_hit_zhitg = 0;
   Earm_BBPSTF1_hit_sumedep = 0;
   Earm_BBPSTF1_hit_tavg = 0;
   Earm_BBPSTF1_hit_trms = 0;
   Earm_BBPSTF1_hit_tmin = 0;
   Earm_BBPSTF1_hit_tmax = 0;
   Earm_BBPSTF1_hit_otridx = 0;
   Earm_BBPSTF1_hit_ptridx = 0;
   Earm_BBPSTF1_hit_sdtridx = 0;
   Earm_BBSH_hit_PMT = 0;
   Earm_BBSH_hit_row = 0;
   Earm_BBSH_hit_col = 0;
   Earm_BBSH_hit_plane = 0;
   Earm_BBSH_hit_xcell = 0;
   Earm_BBSH_hit_ycell = 0;
   Earm_BBSH_hit_zcell = 0;
   Earm_BBSH_hit_xgcell = 0;
   Earm_BBSH_hit_ygcell = 0;
   Earm_BBSH_hit_zgcell = 0;
   Earm_BBSH_hit_NumPhotoelectrons = 0;
   Earm_BBSH_hit_Time_avg = 0;
   Earm_BBSH_hit_Time_rms = 0;
   Earm_BBSH_hit_Time_min = 0;
   Earm_BBSH_hit_Time_max = 0;
   Earm_BBSH_hit_otridx = 0;
   Earm_BBSH_hit_ptridx = 0;
   Earm_BBSH_hit_sdtridx = 0;
   Earm_BBSHTF1_hit_row = 0;
   Earm_BBSHTF1_hit_col = 0;
   Earm_BBSHTF1_hit_cell = 0;
   Earm_BBSHTF1_hit_plane = 0;
   Earm_BBSHTF1_hit_wire = 0;
   Earm_BBSHTF1_hit_xcell = 0;
   Earm_BBSHTF1_hit_ycell = 0;
   Earm_BBSHTF1_hit_zcell = 0;
   Earm_BBSHTF1_hit_xcellg = 0;
   Earm_BBSHTF1_hit_ycellg = 0;
   Earm_BBSHTF1_hit_zcellg = 0;
   Earm_BBSHTF1_hit_xhit = 0;
   Earm_BBSHTF1_hit_yhit = 0;
   Earm_BBSHTF1_hit_zhit = 0;
   Earm_BBSHTF1_hit_xhitg = 0;
   Earm_BBSHTF1_hit_yhitg = 0;
   Earm_BBSHTF1_hit_zhitg = 0;
   Earm_BBSHTF1_hit_sumedep = 0;
   Earm_BBSHTF1_hit_tavg = 0;
   Earm_BBSHTF1_hit_trms = 0;
   Earm_BBSHTF1_hit_tmin = 0;
   Earm_BBSHTF1_hit_tmax = 0;
   Earm_BBSHTF1_hit_otridx = 0;
   Earm_BBSHTF1_hit_ptridx = 0;
   Earm_BBSHTF1_hit_sdtridx = 0;
   Earm_GRINCH_hit_PMT = 0;
   Earm_GRINCH_hit_row = 0;
   Earm_GRINCH_hit_col = 0;
   Earm_GRINCH_hit_xpmt = 0;
   Earm_GRINCH_hit_ypmt = 0;
   Earm_GRINCH_hit_zpmt = 0;
   Earm_GRINCH_hit_xgpmt = 0;
   Earm_GRINCH_hit_ygpmt = 0;
   Earm_GRINCH_hit_zgpmt = 0;
   Earm_GRINCH_hit_NumPhotoelectrons = 0;
   Earm_GRINCH_hit_Time_avg = 0;
   Earm_GRINCH_hit_Time_rms = 0;
   Earm_GRINCH_hit_Time_min = 0;
   Earm_GRINCH_hit_Time_max = 0;
   Earm_GRINCH_hit_mTrackNo = 0;
   Earm_GRINCH_hit_xhit = 0;
   Earm_GRINCH_hit_yhit = 0;
   Earm_GRINCH_hit_zhit = 0;
   Earm_GRINCH_hit_pxhit = 0;
   Earm_GRINCH_hit_pyhit = 0;
   Earm_GRINCH_hit_pzhit = 0;
   Earm_GRINCH_hit_pvx = 0;
   Earm_GRINCH_hit_pvy = 0;
   Earm_GRINCH_hit_pvz = 0;
   Earm_GRINCH_hit_ppx = 0;
   Earm_GRINCH_hit_ppy = 0;
   Earm_GRINCH_hit_ppz = 0;
   Earm_GRINCH_hit_volume_flag = 0;
   Earm_GRINCH_hit_otridx = 0;
   Earm_GRINCH_hit_ptridx = 0;
   Earm_GRINCH_hit_sdtridx = 0;
   Harm_HCal_hit_PMT = 0;
   Harm_HCal_hit_row = 0;
   Harm_HCal_hit_col = 0;
   Harm_HCal_hit_plane = 0;
   Harm_HCal_hit_xcell = 0;
   Harm_HCal_hit_ycell = 0;
   Harm_HCal_hit_zcell = 0;
   Harm_HCal_hit_xgcell = 0;
   Harm_HCal_hit_ygcell = 0;
   Harm_HCal_hit_zgcell = 0;
   Harm_HCal_hit_NumPhotoelectrons = 0;
   Harm_HCal_hit_Time_avg = 0;
   Harm_HCal_hit_Time_rms = 0;
   Harm_HCal_hit_Time_min = 0;
   Harm_HCal_hit_Time_max = 0;
   Harm_HCal_hit_otridx = 0;
   Harm_HCal_hit_ptridx = 0;
   Harm_HCal_hit_sdtridx = 0;
   Harm_HCalScint_hit_row = 0;
   Harm_HCalScint_hit_col = 0;
   Harm_HCalScint_hit_cell = 0;
   Harm_HCalScint_hit_plane = 0;
   Harm_HCalScint_hit_wire = 0;
   Harm_HCalScint_hit_xcell = 0;
   Harm_HCalScint_hit_ycell = 0;
   Harm_HCalScint_hit_zcell = 0;
   Harm_HCalScint_hit_xcellg = 0;
   Harm_HCalScint_hit_ycellg = 0;
   Harm_HCalScint_hit_zcellg = 0;
   Harm_HCalScint_hit_xhit = 0;
   Harm_HCalScint_hit_yhit = 0;
   Harm_HCalScint_hit_zhit = 0;
   Harm_HCalScint_hit_xhitg = 0;
   Harm_HCalScint_hit_yhitg = 0;
   Harm_HCalScint_hit_zhitg = 0;
   Harm_HCalScint_hit_sumedep = 0;
   Harm_HCalScint_hit_tavg = 0;
   Harm_HCalScint_hit_trms = 0;
   Harm_HCalScint_hit_tmin = 0;
   Harm_HCalScint_hit_tmax = 0;
   Harm_HCalScint_hit_otridx = 0;
   Harm_HCalScint_hit_ptridx = 0;
   Harm_HCalScint_hit_sdtridx = 0;
   OTrack_TID = 0;
   OTrack_MID = 0;
   OTrack_PID = 0;
   OTrack_posx = 0;
   OTrack_posy = 0;
   OTrack_posz = 0;
   OTrack_momx = 0;
   OTrack_momy = 0;
   OTrack_momz = 0;
   OTrack_polx = 0;
   OTrack_poly = 0;
   OTrack_polz = 0;
   OTrack_Etot = 0;
   OTrack_T = 0;
   PTrack_TID = 0;
   PTrack_PID = 0;
   PTrack_posx = 0;
   PTrack_posy = 0;
   PTrack_posz = 0;
   PTrack_momx = 0;
   PTrack_momy = 0;
   PTrack_momz = 0;
   PTrack_polx = 0;
   PTrack_poly = 0;
   PTrack_polz = 0;
   PTrack_Etot = 0;
   PTrack_T = 0;
   SDTrack_TID = 0;
   SDTrack_MID = 0;
   SDTrack_PID = 0;
   SDTrack_posx = 0;
   SDTrack_posy = 0;
   SDTrack_posz = 0;
   SDTrack_momx = 0;
   SDTrack_momy = 0;
   SDTrack_momz = 0;
   SDTrack_polx = 0;
   SDTrack_poly = 0;
   SDTrack_polz = 0;
   SDTrack_Etot = 0;
   SDTrack_T = 0;
   Earm_BBPS_dighit_chan = 0;
   Earm_BBPS_dighit_adc = 0;
   Earm_BBSH_dighit_chan = 0;
   Earm_BBSH_dighit_adc = 0;
   Earm_BBHodo_dighit_chan = 0;
   Earm_BBHodo_dighit_adc = 0;
   Earm_BBHodo_dighit_tdc_l = 0;
   Earm_BBHodo_dighit_tdc_t = 0;
   Earm_GRINCH_dighit_chan = 0;
   Earm_GRINCH_dighit_adc = 0;
   Earm_GRINCH_dighit_tdc_l = 0;
   Earm_GRINCH_dighit_tdc_t = 0;
   Earm_BBGEM_1x_dighit_strip = 0;
   Earm_BBGEM_1x_dighit_adc_0 = 0;
   Earm_BBGEM_1x_dighit_adc_1 = 0;
   Earm_BBGEM_1x_dighit_adc_2 = 0;
   Earm_BBGEM_1x_dighit_adc_3 = 0;
   Earm_BBGEM_1x_dighit_adc_4 = 0;
   Earm_BBGEM_1x_dighit_adc_5 = 0;
   Earm_BBGEM_1y_dighit_strip = 0;
   Earm_BBGEM_1y_dighit_adc_0 = 0;
   Earm_BBGEM_1y_dighit_adc_1 = 0;
   Earm_BBGEM_1y_dighit_adc_2 = 0;
   Earm_BBGEM_1y_dighit_adc_3 = 0;
   Earm_BBGEM_1y_dighit_adc_4 = 0;
   Earm_BBGEM_1y_dighit_adc_5 = 0;
   Earm_BBGEM_2x_dighit_strip = 0;
   Earm_BBGEM_2x_dighit_adc_0 = 0;
   Earm_BBGEM_2x_dighit_adc_1 = 0;
   Earm_BBGEM_2x_dighit_adc_2 = 0;
   Earm_BBGEM_2x_dighit_adc_3 = 0;
   Earm_BBGEM_2x_dighit_adc_4 = 0;
   Earm_BBGEM_2x_dighit_adc_5 = 0;
   Earm_BBGEM_2y_dighit_strip = 0;
   Earm_BBGEM_2y_dighit_adc_0 = 0;
   Earm_BBGEM_2y_dighit_adc_1 = 0;
   Earm_BBGEM_2y_dighit_adc_2 = 0;
   Earm_BBGEM_2y_dighit_adc_3 = 0;
   Earm_BBGEM_2y_dighit_adc_4 = 0;
   Earm_BBGEM_2y_dighit_adc_5 = 0;
   Earm_BBGEM_3x_dighit_strip = 0;
   Earm_BBGEM_3x_dighit_adc_0 = 0;
   Earm_BBGEM_3x_dighit_adc_1 = 0;
   Earm_BBGEM_3x_dighit_adc_2 = 0;
   Earm_BBGEM_3x_dighit_adc_3 = 0;
   Earm_BBGEM_3x_dighit_adc_4 = 0;
   Earm_BBGEM_3x_dighit_adc_5 = 0;
   Earm_BBGEM_3y_dighit_strip = 0;
   Earm_BBGEM_3y_dighit_adc_0 = 0;
   Earm_BBGEM_3y_dighit_adc_1 = 0;
   Earm_BBGEM_3y_dighit_adc_2 = 0;
   Earm_BBGEM_3y_dighit_adc_3 = 0;
   Earm_BBGEM_3y_dighit_adc_4 = 0;
   Earm_BBGEM_3y_dighit_adc_5 = 0;
   Earm_BBGEM_4x_dighit_strip = 0;
   Earm_BBGEM_4x_dighit_adc_0 = 0;
   Earm_BBGEM_4x_dighit_adc_1 = 0;
   Earm_BBGEM_4x_dighit_adc_2 = 0;
   Earm_BBGEM_4x_dighit_adc_3 = 0;
   Earm_BBGEM_4x_dighit_adc_4 = 0;
   Earm_BBGEM_4x_dighit_adc_5 = 0;
   Earm_BBGEM_4y_dighit_strip = 0;
   Earm_BBGEM_4y_dighit_adc_0 = 0;
   Earm_BBGEM_4y_dighit_adc_1 = 0;
   Earm_BBGEM_4y_dighit_adc_2 = 0;
   Earm_BBGEM_4y_dighit_adc_3 = 0;
   Earm_BBGEM_4y_dighit_adc_4 = 0;
   Earm_BBGEM_4y_dighit_adc_5 = 0;
   Earm_BBGEM_5x_dighit_strip = 0;
   Earm_BBGEM_5x_dighit_adc_0 = 0;
   Earm_BBGEM_5x_dighit_adc_1 = 0;
   Earm_BBGEM_5x_dighit_adc_2 = 0;
   Earm_BBGEM_5x_dighit_adc_3 = 0;
   Earm_BBGEM_5x_dighit_adc_4 = 0;
   Earm_BBGEM_5x_dighit_adc_5 = 0;
   Earm_BBGEM_5y_dighit_strip = 0;
   Earm_BBGEM_5y_dighit_adc_0 = 0;
   Earm_BBGEM_5y_dighit_adc_1 = 0;
   Earm_BBGEM_5y_dighit_adc_2 = 0;
   Earm_BBGEM_5y_dighit_adc_3 = 0;
   Earm_BBGEM_5y_dighit_adc_4 = 0;
   Earm_BBGEM_5y_dighit_adc_5 = 0;
   Harm_HCal_dighit_chan = 0;
   Harm_HCal_dighit_adc_0 = 0;
   Harm_HCal_dighit_adc_1 = 0;
   Harm_HCal_dighit_adc_2 = 0;
   Harm_HCal_dighit_adc_3 = 0;
   Harm_HCal_dighit_adc_4 = 0;
   Harm_HCal_dighit_adc_5 = 0;
   Harm_HCal_dighit_adc_6 = 0;
   Harm_HCal_dighit_adc_7 = 0;
   Harm_HCal_dighit_adc_8 = 0;
   Harm_HCal_dighit_adc_9 = 0;
   Harm_HCal_dighit_adc_10 = 0;
   Harm_HCal_dighit_adc_11 = 0;
   Harm_HCal_dighit_adc_12 = 0;
   Harm_HCal_dighit_adc_13 = 0;
   Harm_HCal_dighit_adc_14 = 0;
   Harm_HCal_dighit_adc_15 = 0;
   Harm_HCal_dighit_adc_16 = 0;
   Harm_HCal_dighit_adc_17 = 0;
   Harm_HCal_dighit_adc_18 = 0;
   Harm_HCal_dighit_adc_19 = 0;
   Harm_HCal_dighit_tdc = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("Earm.BBGEM.hit.nhits", &Earm_BBGEM_hit_nhits, &b_Earm_BBGEM_hit_nhits);
   fChain->SetBranchAddress("Earm.BBGEM.hit.plane", &Earm_BBGEM_hit_plane, &b_Earm_BBGEM_hit_plane);
   fChain->SetBranchAddress("Earm.BBGEM.hit.strip", &Earm_BBGEM_hit_strip, &b_Earm_BBGEM_hit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.hit.x", &Earm_BBGEM_hit_x, &b_Earm_BBGEM_hit_x);
   fChain->SetBranchAddress("Earm.BBGEM.hit.y", &Earm_BBGEM_hit_y, &b_Earm_BBGEM_hit_y);
   fChain->SetBranchAddress("Earm.BBGEM.hit.z", &Earm_BBGEM_hit_z, &b_Earm_BBGEM_hit_z);
   fChain->SetBranchAddress("Earm.BBGEM.hit.polx", &Earm_BBGEM_hit_polx, &b_Earm_BBGEM_hit_polx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.poly", &Earm_BBGEM_hit_poly, &b_Earm_BBGEM_hit_poly);
   fChain->SetBranchAddress("Earm.BBGEM.hit.polz", &Earm_BBGEM_hit_polz, &b_Earm_BBGEM_hit_polz);
   fChain->SetBranchAddress("Earm.BBGEM.hit.t", &Earm_BBGEM_hit_t, &b_Earm_BBGEM_hit_t);
   fChain->SetBranchAddress("Earm.BBGEM.hit.trms", &Earm_BBGEM_hit_trms, &b_Earm_BBGEM_hit_trms);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tmin", &Earm_BBGEM_hit_tmin, &b_Earm_BBGEM_hit_tmin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tmax", &Earm_BBGEM_hit_tmax, &b_Earm_BBGEM_hit_tmax);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tx", &Earm_BBGEM_hit_tx, &b_Earm_BBGEM_hit_tx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.ty", &Earm_BBGEM_hit_ty, &b_Earm_BBGEM_hit_ty);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xin", &Earm_BBGEM_hit_xin, &b_Earm_BBGEM_hit_xin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yin", &Earm_BBGEM_hit_yin, &b_Earm_BBGEM_hit_yin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zin", &Earm_BBGEM_hit_zin, &b_Earm_BBGEM_hit_zin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xout", &Earm_BBGEM_hit_xout, &b_Earm_BBGEM_hit_xout);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yout", &Earm_BBGEM_hit_yout, &b_Earm_BBGEM_hit_yout);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zout", &Earm_BBGEM_hit_zout, &b_Earm_BBGEM_hit_zout);
   fChain->SetBranchAddress("Earm.BBGEM.hit.txp", &Earm_BBGEM_hit_txp, &b_Earm_BBGEM_hit_txp);
   fChain->SetBranchAddress("Earm.BBGEM.hit.typ", &Earm_BBGEM_hit_typ, &b_Earm_BBGEM_hit_typ);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xg", &Earm_BBGEM_hit_xg, &b_Earm_BBGEM_hit_xg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yg", &Earm_BBGEM_hit_yg, &b_Earm_BBGEM_hit_yg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zg", &Earm_BBGEM_hit_zg, &b_Earm_BBGEM_hit_zg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.trid", &Earm_BBGEM_hit_trid, &b_Earm_BBGEM_hit_trid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.mid", &Earm_BBGEM_hit_mid, &b_Earm_BBGEM_hit_mid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.pid", &Earm_BBGEM_hit_pid, &b_Earm_BBGEM_hit_pid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vx", &Earm_BBGEM_hit_vx, &b_Earm_BBGEM_hit_vx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vy", &Earm_BBGEM_hit_vy, &b_Earm_BBGEM_hit_vy);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vz", &Earm_BBGEM_hit_vz, &b_Earm_BBGEM_hit_vz);
   fChain->SetBranchAddress("Earm.BBGEM.hit.p", &Earm_BBGEM_hit_p, &b_Earm_BBGEM_hit_p);
   fChain->SetBranchAddress("Earm.BBGEM.hit.edep", &Earm_BBGEM_hit_edep, &b_Earm_BBGEM_hit_edep);
   fChain->SetBranchAddress("Earm.BBGEM.hit.beta", &Earm_BBGEM_hit_beta, &b_Earm_BBGEM_hit_beta);
   fChain->SetBranchAddress("Earm.BBGEM.hit.otridx", &Earm_BBGEM_hit_otridx, &b_Earm_BBGEM_hit_otridx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.ptridx", &Earm_BBGEM_hit_ptridx, &b_Earm_BBGEM_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.sdtridx", &Earm_BBGEM_hit_sdtridx, &b_Earm_BBGEM_hit_sdtridx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.ntracks", &Earm_BBGEM_Track_ntracks, &b_Earm_BBGEM_Track_ntracks);
   fChain->SetBranchAddress("Earm.BBGEM.Track.TID", &Earm_BBGEM_Track_TID, &b_Earm_BBGEM_Track_TID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.PID", &Earm_BBGEM_Track_PID, &b_Earm_BBGEM_Track_PID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.MID", &Earm_BBGEM_Track_MID, &b_Earm_BBGEM_Track_MID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NumHits", &Earm_BBGEM_Track_NumHits, &b_Earm_BBGEM_Track_NumHits);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NumPlanes", &Earm_BBGEM_Track_NumPlanes, &b_Earm_BBGEM_Track_NumPlanes);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NDF", &Earm_BBGEM_Track_NDF, &b_Earm_BBGEM_Track_NDF);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2fit", &Earm_BBGEM_Track_Chi2fit, &b_Earm_BBGEM_Track_Chi2fit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2true", &Earm_BBGEM_Track_Chi2true, &b_Earm_BBGEM_Track_Chi2true);
   fChain->SetBranchAddress("Earm.BBGEM.Track.X", &Earm_BBGEM_Track_X, &b_Earm_BBGEM_Track_X);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Y", &Earm_BBGEM_Track_Y, &b_Earm_BBGEM_Track_Y);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xp", &Earm_BBGEM_Track_Xp, &b_Earm_BBGEM_Track_Xp);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Yp", &Earm_BBGEM_Track_Yp, &b_Earm_BBGEM_Track_Yp);
   fChain->SetBranchAddress("Earm.BBGEM.Track.T", &Earm_BBGEM_Track_T, &b_Earm_BBGEM_Track_T);
   fChain->SetBranchAddress("Earm.BBGEM.Track.P", &Earm_BBGEM_Track_P, &b_Earm_BBGEM_Track_P);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sx", &Earm_BBGEM_Track_Sx, &b_Earm_BBGEM_Track_Sx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sy", &Earm_BBGEM_Track_Sy, &b_Earm_BBGEM_Track_Sy);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sz", &Earm_BBGEM_Track_Sz, &b_Earm_BBGEM_Track_Sz);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xfit", &Earm_BBGEM_Track_Xfit, &b_Earm_BBGEM_Track_Xfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Yfit", &Earm_BBGEM_Track_Yfit, &b_Earm_BBGEM_Track_Yfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xpfit", &Earm_BBGEM_Track_Xpfit, &b_Earm_BBGEM_Track_Xpfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Ypfit", &Earm_BBGEM_Track_Ypfit, &b_Earm_BBGEM_Track_Ypfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.otridx", &Earm_BBGEM_Track_otridx, &b_Earm_BBGEM_Track_otridx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.ptridx", &Earm_BBGEM_Track_ptridx, &b_Earm_BBGEM_Track_ptridx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.sdtridx", &Earm_BBGEM_Track_sdtridx, &b_Earm_BBGEM_Track_sdtridx);
   fChain->SetBranchAddress("Earm.BBHodoScint.det.esum", &Earm_BBHodoScint_det_esum, &b_Earm_BBHodoScint_det_esum);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.nhits", &Earm_BBHodoScint_hit_nhits, &b_Earm_BBHodoScint_hit_nhits);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.row", &Earm_BBHodoScint_hit_row, &b_Earm_BBHodoScint_hit_row);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.col", &Earm_BBHodoScint_hit_col, &b_Earm_BBHodoScint_hit_col);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.cell", &Earm_BBHodoScint_hit_cell, &b_Earm_BBHodoScint_hit_cell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.plane", &Earm_BBHodoScint_hit_plane, &b_Earm_BBHodoScint_hit_plane);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.wire", &Earm_BBHodoScint_hit_wire, &b_Earm_BBHodoScint_hit_wire);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xcell", &Earm_BBHodoScint_hit_xcell, &b_Earm_BBHodoScint_hit_xcell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.ycell", &Earm_BBHodoScint_hit_ycell, &b_Earm_BBHodoScint_hit_ycell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zcell", &Earm_BBHodoScint_hit_zcell, &b_Earm_BBHodoScint_hit_zcell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xcellg", &Earm_BBHodoScint_hit_xcellg, &b_Earm_BBHodoScint_hit_xcellg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.ycellg", &Earm_BBHodoScint_hit_ycellg, &b_Earm_BBHodoScint_hit_ycellg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zcellg", &Earm_BBHodoScint_hit_zcellg, &b_Earm_BBHodoScint_hit_zcellg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xhit", &Earm_BBHodoScint_hit_xhit, &b_Earm_BBHodoScint_hit_xhit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.yhit", &Earm_BBHodoScint_hit_yhit, &b_Earm_BBHodoScint_hit_yhit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zhit", &Earm_BBHodoScint_hit_zhit, &b_Earm_BBHodoScint_hit_zhit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xhitg", &Earm_BBHodoScint_hit_xhitg, &b_Earm_BBHodoScint_hit_xhitg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.yhitg", &Earm_BBHodoScint_hit_yhitg, &b_Earm_BBHodoScint_hit_yhitg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zhitg", &Earm_BBHodoScint_hit_zhitg, &b_Earm_BBHodoScint_hit_zhitg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.sumedep", &Earm_BBHodoScint_hit_sumedep, &b_Earm_BBHodoScint_hit_sumedep);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.tavg", &Earm_BBHodoScint_hit_tavg, &b_Earm_BBHodoScint_hit_tavg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.trms", &Earm_BBHodoScint_hit_trms, &b_Earm_BBHodoScint_hit_trms);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.tmin", &Earm_BBHodoScint_hit_tmin, &b_Earm_BBHodoScint_hit_tmin);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.tmax", &Earm_BBHodoScint_hit_tmax, &b_Earm_BBHodoScint_hit_tmax);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.otridx", &Earm_BBHodoScint_hit_otridx, &b_Earm_BBHodoScint_hit_otridx);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.ptridx", &Earm_BBHodoScint_hit_ptridx, &b_Earm_BBHodoScint_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.sdtridx", &Earm_BBHodoScint_hit_sdtridx, &b_Earm_BBHodoScint_hit_sdtridx);
   fChain->SetBranchAddress("Earm.BBPS.hit.nhits", &Earm_BBPS_hit_nhits, &b_Earm_BBPS_hit_nhits);
   fChain->SetBranchAddress("Earm.BBPS.hit.PMT", &Earm_BBPS_hit_PMT, &b_Earm_BBPS_hit_PMT);
   fChain->SetBranchAddress("Earm.BBPS.hit.row", &Earm_BBPS_hit_row, &b_Earm_BBPS_hit_row);
   fChain->SetBranchAddress("Earm.BBPS.hit.col", &Earm_BBPS_hit_col, &b_Earm_BBPS_hit_col);
   fChain->SetBranchAddress("Earm.BBPS.hit.plane", &Earm_BBPS_hit_plane, &b_Earm_BBPS_hit_plane);
   fChain->SetBranchAddress("Earm.BBPS.hit.xcell", &Earm_BBPS_hit_xcell, &b_Earm_BBPS_hit_xcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.ycell", &Earm_BBPS_hit_ycell, &b_Earm_BBPS_hit_ycell);
   fChain->SetBranchAddress("Earm.BBPS.hit.zcell", &Earm_BBPS_hit_zcell, &b_Earm_BBPS_hit_zcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.xgcell", &Earm_BBPS_hit_xgcell, &b_Earm_BBPS_hit_xgcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.ygcell", &Earm_BBPS_hit_ygcell, &b_Earm_BBPS_hit_ygcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.zgcell", &Earm_BBPS_hit_zgcell, &b_Earm_BBPS_hit_zgcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.NumPhotoelectrons", &Earm_BBPS_hit_NumPhotoelectrons, &b_Earm_BBPS_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_avg", &Earm_BBPS_hit_Time_avg, &b_Earm_BBPS_hit_Time_avg);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_rms", &Earm_BBPS_hit_Time_rms, &b_Earm_BBPS_hit_Time_rms);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_min", &Earm_BBPS_hit_Time_min, &b_Earm_BBPS_hit_Time_min);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_max", &Earm_BBPS_hit_Time_max, &b_Earm_BBPS_hit_Time_max);
   fChain->SetBranchAddress("Earm.BBPS.hit.otridx", &Earm_BBPS_hit_otridx, &b_Earm_BBPS_hit_otridx);
   fChain->SetBranchAddress("Earm.BBPS.hit.ptridx", &Earm_BBPS_hit_ptridx, &b_Earm_BBPS_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBPS.hit.sdtridx", &Earm_BBPS_hit_sdtridx, &b_Earm_BBPS_hit_sdtridx);
   fChain->SetBranchAddress("Earm.BBPSTF1.det.esum", &Earm_BBPSTF1_det_esum, &b_Earm_BBPSTF1_det_esum);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.nhits", &Earm_BBPSTF1_hit_nhits, &b_Earm_BBPSTF1_hit_nhits);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.row", &Earm_BBPSTF1_hit_row, &b_Earm_BBPSTF1_hit_row);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.col", &Earm_BBPSTF1_hit_col, &b_Earm_BBPSTF1_hit_col);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.cell", &Earm_BBPSTF1_hit_cell, &b_Earm_BBPSTF1_hit_cell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.plane", &Earm_BBPSTF1_hit_plane, &b_Earm_BBPSTF1_hit_plane);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.wire", &Earm_BBPSTF1_hit_wire, &b_Earm_BBPSTF1_hit_wire);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xcell", &Earm_BBPSTF1_hit_xcell, &b_Earm_BBPSTF1_hit_xcell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.ycell", &Earm_BBPSTF1_hit_ycell, &b_Earm_BBPSTF1_hit_ycell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zcell", &Earm_BBPSTF1_hit_zcell, &b_Earm_BBPSTF1_hit_zcell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xcellg", &Earm_BBPSTF1_hit_xcellg, &b_Earm_BBPSTF1_hit_xcellg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.ycellg", &Earm_BBPSTF1_hit_ycellg, &b_Earm_BBPSTF1_hit_ycellg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zcellg", &Earm_BBPSTF1_hit_zcellg, &b_Earm_BBPSTF1_hit_zcellg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xhit", &Earm_BBPSTF1_hit_xhit, &b_Earm_BBPSTF1_hit_xhit);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.yhit", &Earm_BBPSTF1_hit_yhit, &b_Earm_BBPSTF1_hit_yhit);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zhit", &Earm_BBPSTF1_hit_zhit, &b_Earm_BBPSTF1_hit_zhit);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xhitg", &Earm_BBPSTF1_hit_xhitg, &b_Earm_BBPSTF1_hit_xhitg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.yhitg", &Earm_BBPSTF1_hit_yhitg, &b_Earm_BBPSTF1_hit_yhitg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zhitg", &Earm_BBPSTF1_hit_zhitg, &b_Earm_BBPSTF1_hit_zhitg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.sumedep", &Earm_BBPSTF1_hit_sumedep, &b_Earm_BBPSTF1_hit_sumedep);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.tavg", &Earm_BBPSTF1_hit_tavg, &b_Earm_BBPSTF1_hit_tavg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.trms", &Earm_BBPSTF1_hit_trms, &b_Earm_BBPSTF1_hit_trms);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.tmin", &Earm_BBPSTF1_hit_tmin, &b_Earm_BBPSTF1_hit_tmin);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.tmax", &Earm_BBPSTF1_hit_tmax, &b_Earm_BBPSTF1_hit_tmax);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.otridx", &Earm_BBPSTF1_hit_otridx, &b_Earm_BBPSTF1_hit_otridx);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.ptridx", &Earm_BBPSTF1_hit_ptridx, &b_Earm_BBPSTF1_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.sdtridx", &Earm_BBPSTF1_hit_sdtridx, &b_Earm_BBPSTF1_hit_sdtridx);
   fChain->SetBranchAddress("Earm.BBSH.hit.nhits", &Earm_BBSH_hit_nhits, &b_Earm_BBSH_hit_nhits);
   fChain->SetBranchAddress("Earm.BBSH.hit.PMT", &Earm_BBSH_hit_PMT, &b_Earm_BBSH_hit_PMT);
   fChain->SetBranchAddress("Earm.BBSH.hit.row", &Earm_BBSH_hit_row, &b_Earm_BBSH_hit_row);
   fChain->SetBranchAddress("Earm.BBSH.hit.col", &Earm_BBSH_hit_col, &b_Earm_BBSH_hit_col);
   fChain->SetBranchAddress("Earm.BBSH.hit.plane", &Earm_BBSH_hit_plane, &b_Earm_BBSH_hit_plane);
   fChain->SetBranchAddress("Earm.BBSH.hit.xcell", &Earm_BBSH_hit_xcell, &b_Earm_BBSH_hit_xcell);
   fChain->SetBranchAddress("Earm.BBSH.hit.ycell", &Earm_BBSH_hit_ycell, &b_Earm_BBSH_hit_ycell);
   fChain->SetBranchAddress("Earm.BBSH.hit.zcell", &Earm_BBSH_hit_zcell, &b_Earm_BBSH_hit_zcell);
   fChain->SetBranchAddress("Earm.BBSH.hit.xgcell", &Earm_BBSH_hit_xgcell, &b_Earm_BBSH_hit_xgcell);
   fChain->SetBranchAddress("Earm.BBSH.hit.ygcell", &Earm_BBSH_hit_ygcell, &b_Earm_BBSH_hit_ygcell);
   fChain->SetBranchAddress("Earm.BBSH.hit.zgcell", &Earm_BBSH_hit_zgcell, &b_Earm_BBSH_hit_zgcell);
   fChain->SetBranchAddress("Earm.BBSH.hit.NumPhotoelectrons", &Earm_BBSH_hit_NumPhotoelectrons, &b_Earm_BBSH_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.BBSH.hit.Time_avg", &Earm_BBSH_hit_Time_avg, &b_Earm_BBSH_hit_Time_avg);
   fChain->SetBranchAddress("Earm.BBSH.hit.Time_rms", &Earm_BBSH_hit_Time_rms, &b_Earm_BBSH_hit_Time_rms);
   fChain->SetBranchAddress("Earm.BBSH.hit.Time_min", &Earm_BBSH_hit_Time_min, &b_Earm_BBSH_hit_Time_min);
   fChain->SetBranchAddress("Earm.BBSH.hit.Time_max", &Earm_BBSH_hit_Time_max, &b_Earm_BBSH_hit_Time_max);
   fChain->SetBranchAddress("Earm.BBSH.hit.otridx", &Earm_BBSH_hit_otridx, &b_Earm_BBSH_hit_otridx);
   fChain->SetBranchAddress("Earm.BBSH.hit.ptridx", &Earm_BBSH_hit_ptridx, &b_Earm_BBSH_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBSH.hit.sdtridx", &Earm_BBSH_hit_sdtridx, &b_Earm_BBSH_hit_sdtridx);
   fChain->SetBranchAddress("Earm.BBSHTF1.det.esum", &Earm_BBSHTF1_det_esum, &b_Earm_BBSHTF1_det_esum);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.nhits", &Earm_BBSHTF1_hit_nhits, &b_Earm_BBSHTF1_hit_nhits);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.row", &Earm_BBSHTF1_hit_row, &b_Earm_BBSHTF1_hit_row);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.col", &Earm_BBSHTF1_hit_col, &b_Earm_BBSHTF1_hit_col);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.cell", &Earm_BBSHTF1_hit_cell, &b_Earm_BBSHTF1_hit_cell);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.plane", &Earm_BBSHTF1_hit_plane, &b_Earm_BBSHTF1_hit_plane);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.wire", &Earm_BBSHTF1_hit_wire, &b_Earm_BBSHTF1_hit_wire);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.xcell", &Earm_BBSHTF1_hit_xcell, &b_Earm_BBSHTF1_hit_xcell);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.ycell", &Earm_BBSHTF1_hit_ycell, &b_Earm_BBSHTF1_hit_ycell);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.zcell", &Earm_BBSHTF1_hit_zcell, &b_Earm_BBSHTF1_hit_zcell);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.xcellg", &Earm_BBSHTF1_hit_xcellg, &b_Earm_BBSHTF1_hit_xcellg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.ycellg", &Earm_BBSHTF1_hit_ycellg, &b_Earm_BBSHTF1_hit_ycellg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.zcellg", &Earm_BBSHTF1_hit_zcellg, &b_Earm_BBSHTF1_hit_zcellg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.xhit", &Earm_BBSHTF1_hit_xhit, &b_Earm_BBSHTF1_hit_xhit);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.yhit", &Earm_BBSHTF1_hit_yhit, &b_Earm_BBSHTF1_hit_yhit);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.zhit", &Earm_BBSHTF1_hit_zhit, &b_Earm_BBSHTF1_hit_zhit);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.xhitg", &Earm_BBSHTF1_hit_xhitg, &b_Earm_BBSHTF1_hit_xhitg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.yhitg", &Earm_BBSHTF1_hit_yhitg, &b_Earm_BBSHTF1_hit_yhitg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.zhitg", &Earm_BBSHTF1_hit_zhitg, &b_Earm_BBSHTF1_hit_zhitg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.sumedep", &Earm_BBSHTF1_hit_sumedep, &b_Earm_BBSHTF1_hit_sumedep);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.tavg", &Earm_BBSHTF1_hit_tavg, &b_Earm_BBSHTF1_hit_tavg);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.trms", &Earm_BBSHTF1_hit_trms, &b_Earm_BBSHTF1_hit_trms);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.tmin", &Earm_BBSHTF1_hit_tmin, &b_Earm_BBSHTF1_hit_tmin);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.tmax", &Earm_BBSHTF1_hit_tmax, &b_Earm_BBSHTF1_hit_tmax);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.otridx", &Earm_BBSHTF1_hit_otridx, &b_Earm_BBSHTF1_hit_otridx);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.ptridx", &Earm_BBSHTF1_hit_ptridx, &b_Earm_BBSHTF1_hit_ptridx);
   fChain->SetBranchAddress("Earm.BBSHTF1.hit.sdtridx", &Earm_BBSHTF1_hit_sdtridx, &b_Earm_BBSHTF1_hit_sdtridx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.nhits", &Earm_GRINCH_hit_nhits, &b_Earm_GRINCH_hit_nhits);
   fChain->SetBranchAddress("Earm.GRINCH.hit.PMT", &Earm_GRINCH_hit_PMT, &b_Earm_GRINCH_hit_PMT);
   fChain->SetBranchAddress("Earm.GRINCH.hit.row", &Earm_GRINCH_hit_row, &b_Earm_GRINCH_hit_row);
   fChain->SetBranchAddress("Earm.GRINCH.hit.col", &Earm_GRINCH_hit_col, &b_Earm_GRINCH_hit_col);
   fChain->SetBranchAddress("Earm.GRINCH.hit.xpmt", &Earm_GRINCH_hit_xpmt, &b_Earm_GRINCH_hit_xpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ypmt", &Earm_GRINCH_hit_ypmt, &b_Earm_GRINCH_hit_ypmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.zpmt", &Earm_GRINCH_hit_zpmt, &b_Earm_GRINCH_hit_zpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.xgpmt", &Earm_GRINCH_hit_xgpmt, &b_Earm_GRINCH_hit_xgpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ygpmt", &Earm_GRINCH_hit_ygpmt, &b_Earm_GRINCH_hit_ygpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.zgpmt", &Earm_GRINCH_hit_zgpmt, &b_Earm_GRINCH_hit_zgpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.NumPhotoelectrons", &Earm_GRINCH_hit_NumPhotoelectrons, &b_Earm_GRINCH_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_avg", &Earm_GRINCH_hit_Time_avg, &b_Earm_GRINCH_hit_Time_avg);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_rms", &Earm_GRINCH_hit_Time_rms, &b_Earm_GRINCH_hit_Time_rms);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_min", &Earm_GRINCH_hit_Time_min, &b_Earm_GRINCH_hit_Time_min);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_max", &Earm_GRINCH_hit_Time_max, &b_Earm_GRINCH_hit_Time_max);
   fChain->SetBranchAddress("Earm.GRINCH.hit.mTrackNo", &Earm_GRINCH_hit_mTrackNo, &b_Earm_GRINCH_hit_mTrackNo);
   fChain->SetBranchAddress("Earm.GRINCH.hit.xhit", &Earm_GRINCH_hit_xhit, &b_Earm_GRINCH_hit_xhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.yhit", &Earm_GRINCH_hit_yhit, &b_Earm_GRINCH_hit_yhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.zhit", &Earm_GRINCH_hit_zhit, &b_Earm_GRINCH_hit_zhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pxhit", &Earm_GRINCH_hit_pxhit, &b_Earm_GRINCH_hit_pxhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pyhit", &Earm_GRINCH_hit_pyhit, &b_Earm_GRINCH_hit_pyhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pzhit", &Earm_GRINCH_hit_pzhit, &b_Earm_GRINCH_hit_pzhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pvx", &Earm_GRINCH_hit_pvx, &b_Earm_GRINCH_hit_pvx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pvy", &Earm_GRINCH_hit_pvy, &b_Earm_GRINCH_hit_pvy);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pvz", &Earm_GRINCH_hit_pvz, &b_Earm_GRINCH_hit_pvz);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ppx", &Earm_GRINCH_hit_ppx, &b_Earm_GRINCH_hit_ppx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ppy", &Earm_GRINCH_hit_ppy, &b_Earm_GRINCH_hit_ppy);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ppz", &Earm_GRINCH_hit_ppz, &b_Earm_GRINCH_hit_ppz);
   fChain->SetBranchAddress("Earm.GRINCH.hit.volume_flag", &Earm_GRINCH_hit_volume_flag, &b_Earm_GRINCH_hit_volume_flag);
   fChain->SetBranchAddress("Earm.GRINCH.hit.otridx", &Earm_GRINCH_hit_otridx, &b_Earm_GRINCH_hit_otridx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ptridx", &Earm_GRINCH_hit_ptridx, &b_Earm_GRINCH_hit_ptridx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.sdtridx", &Earm_GRINCH_hit_sdtridx, &b_Earm_GRINCH_hit_sdtridx);
   fChain->SetBranchAddress("Harm.HCal.hit.nhits", &Harm_HCal_hit_nhits, &b_Harm_HCal_hit_nhits);
   fChain->SetBranchAddress("Harm.HCal.hit.PMT", &Harm_HCal_hit_PMT, &b_Harm_HCal_hit_PMT);
   fChain->SetBranchAddress("Harm.HCal.hit.row", &Harm_HCal_hit_row, &b_Harm_HCal_hit_row);
   fChain->SetBranchAddress("Harm.HCal.hit.col", &Harm_HCal_hit_col, &b_Harm_HCal_hit_col);
   fChain->SetBranchAddress("Harm.HCal.hit.plane", &Harm_HCal_hit_plane, &b_Harm_HCal_hit_plane);
   fChain->SetBranchAddress("Harm.HCal.hit.xcell", &Harm_HCal_hit_xcell, &b_Harm_HCal_hit_xcell);
   fChain->SetBranchAddress("Harm.HCal.hit.ycell", &Harm_HCal_hit_ycell, &b_Harm_HCal_hit_ycell);
   fChain->SetBranchAddress("Harm.HCal.hit.zcell", &Harm_HCal_hit_zcell, &b_Harm_HCal_hit_zcell);
   fChain->SetBranchAddress("Harm.HCal.hit.xgcell", &Harm_HCal_hit_xgcell, &b_Harm_HCal_hit_xgcell);
   fChain->SetBranchAddress("Harm.HCal.hit.ygcell", &Harm_HCal_hit_ygcell, &b_Harm_HCal_hit_ygcell);
   fChain->SetBranchAddress("Harm.HCal.hit.zgcell", &Harm_HCal_hit_zgcell, &b_Harm_HCal_hit_zgcell);
   fChain->SetBranchAddress("Harm.HCal.hit.NumPhotoelectrons", &Harm_HCal_hit_NumPhotoelectrons, &b_Harm_HCal_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_avg", &Harm_HCal_hit_Time_avg, &b_Harm_HCal_hit_Time_avg);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_rms", &Harm_HCal_hit_Time_rms, &b_Harm_HCal_hit_Time_rms);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_min", &Harm_HCal_hit_Time_min, &b_Harm_HCal_hit_Time_min);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_max", &Harm_HCal_hit_Time_max, &b_Harm_HCal_hit_Time_max);
   fChain->SetBranchAddress("Harm.HCal.hit.otridx", &Harm_HCal_hit_otridx, &b_Harm_HCal_hit_otridx);
   fChain->SetBranchAddress("Harm.HCal.hit.ptridx", &Harm_HCal_hit_ptridx, &b_Harm_HCal_hit_ptridx);
   fChain->SetBranchAddress("Harm.HCal.hit.sdtridx", &Harm_HCal_hit_sdtridx, &b_Harm_HCal_hit_sdtridx);
   fChain->SetBranchAddress("Harm.HCalScint.det.esum", &Harm_HCalScint_det_esum, &b_Harm_HCalScint_det_esum);
   fChain->SetBranchAddress("Harm.HCalScint.hit.nhits", &Harm_HCalScint_hit_nhits, &b_Harm_HCalScint_hit_nhits);
   fChain->SetBranchAddress("Harm.HCalScint.hit.row", &Harm_HCalScint_hit_row, &b_Harm_HCalScint_hit_row);
   fChain->SetBranchAddress("Harm.HCalScint.hit.col", &Harm_HCalScint_hit_col, &b_Harm_HCalScint_hit_col);
   fChain->SetBranchAddress("Harm.HCalScint.hit.cell", &Harm_HCalScint_hit_cell, &b_Harm_HCalScint_hit_cell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.plane", &Harm_HCalScint_hit_plane, &b_Harm_HCalScint_hit_plane);
   fChain->SetBranchAddress("Harm.HCalScint.hit.wire", &Harm_HCalScint_hit_wire, &b_Harm_HCalScint_hit_wire);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xcell", &Harm_HCalScint_hit_xcell, &b_Harm_HCalScint_hit_xcell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ycell", &Harm_HCalScint_hit_ycell, &b_Harm_HCalScint_hit_ycell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zcell", &Harm_HCalScint_hit_zcell, &b_Harm_HCalScint_hit_zcell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xcellg", &Harm_HCalScint_hit_xcellg, &b_Harm_HCalScint_hit_xcellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ycellg", &Harm_HCalScint_hit_ycellg, &b_Harm_HCalScint_hit_ycellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zcellg", &Harm_HCalScint_hit_zcellg, &b_Harm_HCalScint_hit_zcellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xhit", &Harm_HCalScint_hit_xhit, &b_Harm_HCalScint_hit_xhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.yhit", &Harm_HCalScint_hit_yhit, &b_Harm_HCalScint_hit_yhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zhit", &Harm_HCalScint_hit_zhit, &b_Harm_HCalScint_hit_zhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xhitg", &Harm_HCalScint_hit_xhitg, &b_Harm_HCalScint_hit_xhitg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.yhitg", &Harm_HCalScint_hit_yhitg, &b_Harm_HCalScint_hit_yhitg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zhitg", &Harm_HCalScint_hit_zhitg, &b_Harm_HCalScint_hit_zhitg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.sumedep", &Harm_HCalScint_hit_sumedep, &b_Harm_HCalScint_hit_sumedep);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tavg", &Harm_HCalScint_hit_tavg, &b_Harm_HCalScint_hit_tavg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.trms", &Harm_HCalScint_hit_trms, &b_Harm_HCalScint_hit_trms);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tmin", &Harm_HCalScint_hit_tmin, &b_Harm_HCalScint_hit_tmin);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tmax", &Harm_HCalScint_hit_tmax, &b_Harm_HCalScint_hit_tmax);
   fChain->SetBranchAddress("Harm.HCalScint.hit.otridx", &Harm_HCalScint_hit_otridx, &b_Harm_HCalScint_hit_otridx);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ptridx", &Harm_HCalScint_hit_ptridx, &b_Harm_HCalScint_hit_ptridx);
   fChain->SetBranchAddress("Harm.HCalScint.hit.sdtridx", &Harm_HCalScint_hit_sdtridx, &b_Harm_HCalScint_hit_sdtridx);
   fChain->SetBranchAddress("OTrack.ntracks", &OTrack_ntracks, &b_OTrack_ntracks);
   fChain->SetBranchAddress("OTrack.TID", &OTrack_TID, &b_OTrack_TID);
   fChain->SetBranchAddress("OTrack.MID", &OTrack_MID, &b_OTrack_MID);
   fChain->SetBranchAddress("OTrack.PID", &OTrack_PID, &b_OTrack_PID);
   fChain->SetBranchAddress("OTrack.posx", &OTrack_posx, &b_OTrack_posx);
   fChain->SetBranchAddress("OTrack.posy", &OTrack_posy, &b_OTrack_posy);
   fChain->SetBranchAddress("OTrack.posz", &OTrack_posz, &b_OTrack_posz);
   fChain->SetBranchAddress("OTrack.momx", &OTrack_momx, &b_OTrack_momx);
   fChain->SetBranchAddress("OTrack.momy", &OTrack_momy, &b_OTrack_momy);
   fChain->SetBranchAddress("OTrack.momz", &OTrack_momz, &b_OTrack_momz);
   fChain->SetBranchAddress("OTrack.polx", &OTrack_polx, &b_OTrack_polx);
   fChain->SetBranchAddress("OTrack.poly", &OTrack_poly, &b_OTrack_poly);
   fChain->SetBranchAddress("OTrack.polz", &OTrack_polz, &b_OTrack_polz);
   fChain->SetBranchAddress("OTrack.Etot", &OTrack_Etot, &b_OTrack_Etot);
   fChain->SetBranchAddress("OTrack.T", &OTrack_T, &b_OTrack_T);
   fChain->SetBranchAddress("PTrack.ntracks", &PTrack_ntracks, &b_PTrack_ntracks);
   fChain->SetBranchAddress("PTrack.TID", &PTrack_TID, &b_PTrack_TID);
   fChain->SetBranchAddress("PTrack.PID", &PTrack_PID, &b_PTrack_PID);
   fChain->SetBranchAddress("PTrack.posx", &PTrack_posx, &b_PTrack_posx);
   fChain->SetBranchAddress("PTrack.posy", &PTrack_posy, &b_PTrack_posy);
   fChain->SetBranchAddress("PTrack.posz", &PTrack_posz, &b_PTrack_posz);
   fChain->SetBranchAddress("PTrack.momx", &PTrack_momx, &b_PTrack_momx);
   fChain->SetBranchAddress("PTrack.momy", &PTrack_momy, &b_PTrack_momy);
   fChain->SetBranchAddress("PTrack.momz", &PTrack_momz, &b_PTrack_momz);
   fChain->SetBranchAddress("PTrack.polx", &PTrack_polx, &b_PTrack_polx);
   fChain->SetBranchAddress("PTrack.poly", &PTrack_poly, &b_PTrack_poly);
   fChain->SetBranchAddress("PTrack.polz", &PTrack_polz, &b_PTrack_polz);
   fChain->SetBranchAddress("PTrack.Etot", &PTrack_Etot, &b_PTrack_Etot);
   fChain->SetBranchAddress("PTrack.T", &PTrack_T, &b_PTrack_T);
   fChain->SetBranchAddress("SDTrack.ntracks", &SDTrack_ntracks, &b_SDTrack_ntracks);
   fChain->SetBranchAddress("SDTrack.TID", &SDTrack_TID, &b_SDTrack_TID);
   fChain->SetBranchAddress("SDTrack.MID", &SDTrack_MID, &b_SDTrack_MID);
   fChain->SetBranchAddress("SDTrack.PID", &SDTrack_PID, &b_SDTrack_PID);
   fChain->SetBranchAddress("SDTrack.posx", &SDTrack_posx, &b_SDTrack_posx);
   fChain->SetBranchAddress("SDTrack.posy", &SDTrack_posy, &b_SDTrack_posy);
   fChain->SetBranchAddress("SDTrack.posz", &SDTrack_posz, &b_SDTrack_posz);
   fChain->SetBranchAddress("SDTrack.momx", &SDTrack_momx, &b_SDTrack_momx);
   fChain->SetBranchAddress("SDTrack.momy", &SDTrack_momy, &b_SDTrack_momy);
   fChain->SetBranchAddress("SDTrack.momz", &SDTrack_momz, &b_SDTrack_momz);
   fChain->SetBranchAddress("SDTrack.polx", &SDTrack_polx, &b_SDTrack_polx);
   fChain->SetBranchAddress("SDTrack.poly", &SDTrack_poly, &b_SDTrack_poly);
   fChain->SetBranchAddress("SDTrack.polz", &SDTrack_polz, &b_SDTrack_polz);
   fChain->SetBranchAddress("SDTrack.Etot", &SDTrack_Etot, &b_SDTrack_Etot);
   fChain->SetBranchAddress("SDTrack.T", &SDTrack_T, &b_SDTrack_T);
   fChain->SetBranchAddress("Earm.BBPS.dighit.nchan", &Earm_BBPS_dighit_nchan, &b_Earm_BBPS_dighit_nchan);
   fChain->SetBranchAddress("Earm.BBPS.dighit.chan", &Earm_BBPS_dighit_chan, &b_Earm_BBPS_dighit_chan);
   fChain->SetBranchAddress("Earm.BBPS.dighit.adc", &Earm_BBPS_dighit_adc, &b_Earm_BBPS_dighit_adc);
   fChain->SetBranchAddress("Earm.BBSH.dighit.nchan", &Earm_BBSH_dighit_nchan, &b_Earm_BBSH_dighit_nchan);
   fChain->SetBranchAddress("Earm.BBSH.dighit.chan", &Earm_BBSH_dighit_chan, &b_Earm_BBSH_dighit_chan);
   fChain->SetBranchAddress("Earm.BBSH.dighit.adc", &Earm_BBSH_dighit_adc, &b_Earm_BBSH_dighit_adc);
   fChain->SetBranchAddress("Earm.BBHodo.dighit.nchan", &Earm_BBHodo_dighit_nchan, &b_Earm_BBHodo_dighit_nchan);
   fChain->SetBranchAddress("Earm.BBHodo.dighit.chan", &Earm_BBHodo_dighit_chan, &b_Earm_BBHodo_dighit_chan);
   fChain->SetBranchAddress("Earm.BBHodo.dighit.adc", &Earm_BBHodo_dighit_adc, &b_Earm_BBHodo_dighit_adc);
   fChain->SetBranchAddress("Earm.BBHodo.dighit.tdc_l", &Earm_BBHodo_dighit_tdc_l, &b_Earm_BBHodo_dighit_tdc_l);
   fChain->SetBranchAddress("Earm.BBHodo.dighit.tdc_t", &Earm_BBHodo_dighit_tdc_t, &b_Earm_BBHodo_dighit_tdc_t);
   fChain->SetBranchAddress("Earm.GRINCH.dighit.nchan", &Earm_GRINCH_dighit_nchan, &b_Earm_GRINCH_dighit_nchan);
   fChain->SetBranchAddress("Earm.GRINCH.dighit.chan", &Earm_GRINCH_dighit_chan, &b_Earm_GRINCH_dighit_chan);
   fChain->SetBranchAddress("Earm.GRINCH.dighit.adc", &Earm_GRINCH_dighit_adc, &b_Earm_GRINCH_dighit_adc);
   fChain->SetBranchAddress("Earm.GRINCH.dighit.tdc_l", &Earm_GRINCH_dighit_tdc_l, &b_Earm_GRINCH_dighit_tdc_l);
   fChain->SetBranchAddress("Earm.GRINCH.dighit.tdc_t", &Earm_GRINCH_dighit_tdc_t, &b_Earm_GRINCH_dighit_tdc_t);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.nstrips", &Earm_BBGEM_1x_dighit_nstrips, &b_Earm_BBGEM_1x_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.strip", &Earm_BBGEM_1x_dighit_strip, &b_Earm_BBGEM_1x_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.adc_0", &Earm_BBGEM_1x_dighit_adc_0, &b_Earm_BBGEM_1x_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.adc_1", &Earm_BBGEM_1x_dighit_adc_1, &b_Earm_BBGEM_1x_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.adc_2", &Earm_BBGEM_1x_dighit_adc_2, &b_Earm_BBGEM_1x_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.adc_3", &Earm_BBGEM_1x_dighit_adc_3, &b_Earm_BBGEM_1x_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.adc_4", &Earm_BBGEM_1x_dighit_adc_4, &b_Earm_BBGEM_1x_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.1x.dighit.adc_5", &Earm_BBGEM_1x_dighit_adc_5, &b_Earm_BBGEM_1x_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.nstrips", &Earm_BBGEM_1y_dighit_nstrips, &b_Earm_BBGEM_1y_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.strip", &Earm_BBGEM_1y_dighit_strip, &b_Earm_BBGEM_1y_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.adc_0", &Earm_BBGEM_1y_dighit_adc_0, &b_Earm_BBGEM_1y_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.adc_1", &Earm_BBGEM_1y_dighit_adc_1, &b_Earm_BBGEM_1y_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.adc_2", &Earm_BBGEM_1y_dighit_adc_2, &b_Earm_BBGEM_1y_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.adc_3", &Earm_BBGEM_1y_dighit_adc_3, &b_Earm_BBGEM_1y_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.adc_4", &Earm_BBGEM_1y_dighit_adc_4, &b_Earm_BBGEM_1y_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.1y.dighit.adc_5", &Earm_BBGEM_1y_dighit_adc_5, &b_Earm_BBGEM_1y_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.nstrips", &Earm_BBGEM_2x_dighit_nstrips, &b_Earm_BBGEM_2x_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.strip", &Earm_BBGEM_2x_dighit_strip, &b_Earm_BBGEM_2x_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.adc_0", &Earm_BBGEM_2x_dighit_adc_0, &b_Earm_BBGEM_2x_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.adc_1", &Earm_BBGEM_2x_dighit_adc_1, &b_Earm_BBGEM_2x_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.adc_2", &Earm_BBGEM_2x_dighit_adc_2, &b_Earm_BBGEM_2x_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.adc_3", &Earm_BBGEM_2x_dighit_adc_3, &b_Earm_BBGEM_2x_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.adc_4", &Earm_BBGEM_2x_dighit_adc_4, &b_Earm_BBGEM_2x_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.2x.dighit.adc_5", &Earm_BBGEM_2x_dighit_adc_5, &b_Earm_BBGEM_2x_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.nstrips", &Earm_BBGEM_2y_dighit_nstrips, &b_Earm_BBGEM_2y_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.strip", &Earm_BBGEM_2y_dighit_strip, &b_Earm_BBGEM_2y_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.adc_0", &Earm_BBGEM_2y_dighit_adc_0, &b_Earm_BBGEM_2y_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.adc_1", &Earm_BBGEM_2y_dighit_adc_1, &b_Earm_BBGEM_2y_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.adc_2", &Earm_BBGEM_2y_dighit_adc_2, &b_Earm_BBGEM_2y_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.adc_3", &Earm_BBGEM_2y_dighit_adc_3, &b_Earm_BBGEM_2y_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.adc_4", &Earm_BBGEM_2y_dighit_adc_4, &b_Earm_BBGEM_2y_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.2y.dighit.adc_5", &Earm_BBGEM_2y_dighit_adc_5, &b_Earm_BBGEM_2y_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.nstrips", &Earm_BBGEM_3x_dighit_nstrips, &b_Earm_BBGEM_3x_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.strip", &Earm_BBGEM_3x_dighit_strip, &b_Earm_BBGEM_3x_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.adc_0", &Earm_BBGEM_3x_dighit_adc_0, &b_Earm_BBGEM_3x_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.adc_1", &Earm_BBGEM_3x_dighit_adc_1, &b_Earm_BBGEM_3x_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.adc_2", &Earm_BBGEM_3x_dighit_adc_2, &b_Earm_BBGEM_3x_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.adc_3", &Earm_BBGEM_3x_dighit_adc_3, &b_Earm_BBGEM_3x_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.adc_4", &Earm_BBGEM_3x_dighit_adc_4, &b_Earm_BBGEM_3x_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.3x.dighit.adc_5", &Earm_BBGEM_3x_dighit_adc_5, &b_Earm_BBGEM_3x_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.nstrips", &Earm_BBGEM_3y_dighit_nstrips, &b_Earm_BBGEM_3y_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.strip", &Earm_BBGEM_3y_dighit_strip, &b_Earm_BBGEM_3y_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.adc_0", &Earm_BBGEM_3y_dighit_adc_0, &b_Earm_BBGEM_3y_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.adc_1", &Earm_BBGEM_3y_dighit_adc_1, &b_Earm_BBGEM_3y_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.adc_2", &Earm_BBGEM_3y_dighit_adc_2, &b_Earm_BBGEM_3y_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.adc_3", &Earm_BBGEM_3y_dighit_adc_3, &b_Earm_BBGEM_3y_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.adc_4", &Earm_BBGEM_3y_dighit_adc_4, &b_Earm_BBGEM_3y_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.3y.dighit.adc_5", &Earm_BBGEM_3y_dighit_adc_5, &b_Earm_BBGEM_3y_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.nstrips", &Earm_BBGEM_4x_dighit_nstrips, &b_Earm_BBGEM_4x_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.strip", &Earm_BBGEM_4x_dighit_strip, &b_Earm_BBGEM_4x_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.adc_0", &Earm_BBGEM_4x_dighit_adc_0, &b_Earm_BBGEM_4x_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.adc_1", &Earm_BBGEM_4x_dighit_adc_1, &b_Earm_BBGEM_4x_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.adc_2", &Earm_BBGEM_4x_dighit_adc_2, &b_Earm_BBGEM_4x_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.adc_3", &Earm_BBGEM_4x_dighit_adc_3, &b_Earm_BBGEM_4x_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.adc_4", &Earm_BBGEM_4x_dighit_adc_4, &b_Earm_BBGEM_4x_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.4x.dighit.adc_5", &Earm_BBGEM_4x_dighit_adc_5, &b_Earm_BBGEM_4x_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.nstrips", &Earm_BBGEM_4y_dighit_nstrips, &b_Earm_BBGEM_4y_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.strip", &Earm_BBGEM_4y_dighit_strip, &b_Earm_BBGEM_4y_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.adc_0", &Earm_BBGEM_4y_dighit_adc_0, &b_Earm_BBGEM_4y_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.adc_1", &Earm_BBGEM_4y_dighit_adc_1, &b_Earm_BBGEM_4y_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.adc_2", &Earm_BBGEM_4y_dighit_adc_2, &b_Earm_BBGEM_4y_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.adc_3", &Earm_BBGEM_4y_dighit_adc_3, &b_Earm_BBGEM_4y_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.adc_4", &Earm_BBGEM_4y_dighit_adc_4, &b_Earm_BBGEM_4y_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.4y.dighit.adc_5", &Earm_BBGEM_4y_dighit_adc_5, &b_Earm_BBGEM_4y_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.nstrips", &Earm_BBGEM_5x_dighit_nstrips, &b_Earm_BBGEM_5x_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.strip", &Earm_BBGEM_5x_dighit_strip, &b_Earm_BBGEM_5x_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.adc_0", &Earm_BBGEM_5x_dighit_adc_0, &b_Earm_BBGEM_5x_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.adc_1", &Earm_BBGEM_5x_dighit_adc_1, &b_Earm_BBGEM_5x_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.adc_2", &Earm_BBGEM_5x_dighit_adc_2, &b_Earm_BBGEM_5x_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.adc_3", &Earm_BBGEM_5x_dighit_adc_3, &b_Earm_BBGEM_5x_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.adc_4", &Earm_BBGEM_5x_dighit_adc_4, &b_Earm_BBGEM_5x_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.5x.dighit.adc_5", &Earm_BBGEM_5x_dighit_adc_5, &b_Earm_BBGEM_5x_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.nstrips", &Earm_BBGEM_5y_dighit_nstrips, &b_Earm_BBGEM_5y_dighit_nstrips);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.strip", &Earm_BBGEM_5y_dighit_strip, &b_Earm_BBGEM_5y_dighit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.adc_0", &Earm_BBGEM_5y_dighit_adc_0, &b_Earm_BBGEM_5y_dighit_adc_0);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.adc_1", &Earm_BBGEM_5y_dighit_adc_1, &b_Earm_BBGEM_5y_dighit_adc_1);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.adc_2", &Earm_BBGEM_5y_dighit_adc_2, &b_Earm_BBGEM_5y_dighit_adc_2);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.adc_3", &Earm_BBGEM_5y_dighit_adc_3, &b_Earm_BBGEM_5y_dighit_adc_3);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.adc_4", &Earm_BBGEM_5y_dighit_adc_4, &b_Earm_BBGEM_5y_dighit_adc_4);
   fChain->SetBranchAddress("Earm.BBGEM.5y.dighit.adc_5", &Earm_BBGEM_5y_dighit_adc_5, &b_Earm_BBGEM_5y_dighit_adc_5);
   fChain->SetBranchAddress("Earm.BBGEM.dighit.nchan", &Earm_BBGEM_dighit_nchan, &b_Earm_BBGEM_dighit_nchan);
   fChain->SetBranchAddress("Harm.HCal.dighit.chan", &Harm_HCal_dighit_chan, &b_Harm_HCal_dighit_chan);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_0", &Harm_HCal_dighit_adc_0, &b_Harm_HCal_dighit_adc_0);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_1", &Harm_HCal_dighit_adc_1, &b_Harm_HCal_dighit_adc_1);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_2", &Harm_HCal_dighit_adc_2, &b_Harm_HCal_dighit_adc_2);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_3", &Harm_HCal_dighit_adc_3, &b_Harm_HCal_dighit_adc_3);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_4", &Harm_HCal_dighit_adc_4, &b_Harm_HCal_dighit_adc_4);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_5", &Harm_HCal_dighit_adc_5, &b_Harm_HCal_dighit_adc_5);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_6", &Harm_HCal_dighit_adc_6, &b_Harm_HCal_dighit_adc_6);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_7", &Harm_HCal_dighit_adc_7, &b_Harm_HCal_dighit_adc_7);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_8", &Harm_HCal_dighit_adc_8, &b_Harm_HCal_dighit_adc_8);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_9", &Harm_HCal_dighit_adc_9, &b_Harm_HCal_dighit_adc_9);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_10", &Harm_HCal_dighit_adc_10, &b_Harm_HCal_dighit_adc_10);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_11", &Harm_HCal_dighit_adc_11, &b_Harm_HCal_dighit_adc_11);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_12", &Harm_HCal_dighit_adc_12, &b_Harm_HCal_dighit_adc_12);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_13", &Harm_HCal_dighit_adc_13, &b_Harm_HCal_dighit_adc_13);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_14", &Harm_HCal_dighit_adc_14, &b_Harm_HCal_dighit_adc_14);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_15", &Harm_HCal_dighit_adc_15, &b_Harm_HCal_dighit_adc_15);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_16", &Harm_HCal_dighit_adc_16, &b_Harm_HCal_dighit_adc_16);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_17", &Harm_HCal_dighit_adc_17, &b_Harm_HCal_dighit_adc_17);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_18", &Harm_HCal_dighit_adc_18, &b_Harm_HCal_dighit_adc_18);
   fChain->SetBranchAddress("Harm.HCal.dighit.adc_19", &Harm_HCal_dighit_adc_19, &b_Harm_HCal_dighit_adc_19);
   fChain->SetBranchAddress("Harm.HCal.dighit.tdc", &Harm_HCal_dighit_tdc, &b_Harm_HCal_dighit_tdc);
   Notify();
}

Bool_t gmn_dig_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gmn_dig_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gmn_dig_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
  return 1;
}

void gmn_dig_tree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L gmn_dig_tree.C
//      root> gmn_dig_tree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   }
}
