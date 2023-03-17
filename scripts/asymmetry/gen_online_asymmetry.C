


void gen_online_asymmetry(){

  int high_run_num = 5000;

  int good_files = 0;
  int first_run = 0;
  int last_run = 0;
  TString first_run_time;
  TString last_run_time;

  TString root_dir = "/chafs2/work1/sbs/Rootfiles";

  TChain *T = new TChain("T");

    
  for(int irun = high_run_num; irun > 3000; irun--){

    TString root_name = root_dir + Form("/gen_replayed_%i_200k_events.root",irun);
   
    if(!gSystem->AccessPathName(root_name)){
      TFile *file = new TFile(root_name,"read");
      
      if(file->IsZombie()) continue;


      TH1D *h = (TH1D*)file->Get("hdidhitx_bb_gem_m0");

      if(h->GetEntries() == 0) continue;

      T->Add(root_name);

      if(good_files == 0){
	last_run = irun;
	last_run_time = Form("%i-%i-%i",file->GetCreationDate().GetYear(),file->GetCreationDate().GetMonth(),file->GetCreationDate().GetDay());

      }

      first_run = irun;
      first_run_time = Form("%i-%i-%i",file->GetCreationDate().GetYear(),file->GetCreationDate().GetMonth(),file->GetCreationDate().GetDay());
      good_files++;
      
    }
    
    if(good_files == 30) break;
  }
  

  
  TH1D *hpim = new TH1D("hpim","",3,-1.5,1.5);
  TH1D *hpi0 = new TH1D("hpi0","",3,-1.5,1.5);

  T->Draw("scalhel.hel>>hpim","bb.tr.n==1&&abs(bb.tr.vz[0])<0.27&&sbs.hcal.e>0.025&&bb.ps.e<0.2&&bb.ps.e>0.02","goff");
  T->Draw("scalhel.hel>>hpi0","bb.tr.n==0&&sbs.hcal.e>0.025&&bb.ps.e<0.2&&bb.ps.e>0.02","goff");
  

  int Ym = hpim->GetBinContent(1) + hpi0->GetBinContent(1); 
  int Yp = hpim->GetBinContent(3) + hpi0->GetBinContent(3); 
  

  double A = (Yp - Ym)*1.0/(Yp + Ym);
  double A_err = sqrt((1 - A*A)/(Yp + Ym));

  A *= 100;
  A_err *= 100;

  

  cout<<"\n\nOldest file in this replay is run "<<first_run<<" taken on "<<first_run_time<<endl;
  cout<<"Newest file in this replay is run "<<last_run<<" taken on "<<last_run_time<<endl;
  cout<<"Total "<<T->GetEntries()<<" Events Replayed\n\n"<<endl;
  cout<<"Asymmetry = "<<A<<"% +/- "<<A_err<<" %\n"<<endl;
  
}

