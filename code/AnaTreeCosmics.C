//
//     Program to read and analyse the TREE created 
//     by AliCosmics.cxx from ESD for 
//     cosmic muons in the central barrel  (AB)   
//     First Version : 20/Jan/2016 
//     Last Version : 
//

void AnaTreeCosmics(const char* filename){

  // Open the root file in which there is the TREE
//  TFile *cosm  = new TFile("/home/bruno/CosmicGrid/cosmic_run_2015/lhc15a/LHC15a_Bminus.root","READ");
  TFile *cosm  = new TFile(filename,"READ");
//     Open root file in which save the plots
  TFile *saveplot = new TFile("histo_lhc15a_bminus.root","RECREATE","histo_lhc15a_bm");


  //    Create specific output files
//  fileout = fopen("lhc15a_bminus.out","at");
//  fileint = fopen("intev_lhc15a_bminus.out","at");

  //       TREE variables
  TTree *tcosm = (TTree*)cosm->Get("fDatree");
  Int_t nESDentries = (Int_t)tcosm->GetEntries();
  Int_t nRun,nChu,nEve;
  UInt_t timestp;
  Float_t magfield; 
  Int_t nTracks;
  Float_t fracnMuTrk,meanDist;
  Int_t flagintev;
  Int_t nMuons,triggerFlag;
  Int_t uptracmu[500], dwtracmu[500];
  Float_t upthetaCosmic[500],upphiCosmic[500],upphiCosmiccor[500];
  Float_t uptheta[500],upphi[500];
  Float_t dwthetaCosmic[500],dwphiCosmic[500],dwphiCosmiccor[500];
  Float_t dwtheta[500],dwphi[500];
  Float_t upSign[500],dwSign[500];
  Float_t upxv[500],upyv[500],upzv[500],dwxv[500],dwyv[500],dwzv[500];
  Float_t uptpcSignal[500],dwtpcSignal[500];
  Float_t Pup[500],Pdw[500];
  Float_t Sig1Ptup[500],Sig1Ptdw[500];
  Float_t uppx[500],uppy[500],uppz[500],dwpx[500],dwpy[500],dwpz[500];
  Float_t Pmu[500],PmuMed[500],Pcov[500],Pcovupdw[500],Pcovdwup[500],Xinv[500],PRes[500];
  Int_t uptpcNcls[500],dwtpcNcls[500];
  Float_t uptpcChi2[500],dwtpcChi2[500];
  Float_t corrTimeTOFup[500],corrTimeTOFdw[500],ToT_TOFup[500];
  Float_t ToT_TOFdw[500],DzTOFup[500],DzTOFdw[500],DxTOFup[500],DxTOFdw[500];
  Float_t Xinup[500],Yinup[500],Zinup[500],Xoutup[500],Youtup[500],Zoutup[500];
  Float_t Xindw[500],Yindw[500],Zindw[500],Xoutdw[500],Youtdw[500],Zoutdw[500];
  TObjString *fFileName; 
  TObjString *fFileNameChunk;
  Int_t *fEventInFile;
  Int_t hitsACO[60],mcnACO,flagV0;
  Float_t upxdir[500],upydir[500],upzdir[500],uptetadir[500],upphidir[500];
  Float_t dwxdir[500],dwydir[500],dwzdir[500],dwtetadir[500],dwphidir[500];

  //        Read the TREE variables (90 elements)
  tcosm->SetBranchAddress("nRun",&nRun);
  tcosm->SetBranchAddress("nChu",&nChu);
  tcosm->SetBranchAddress("nEve",&nEve);
  tcosm->SetBranchAddress("timestp",&timestp);
  tcosm->SetBranchAddress("magfield",&magfield);
  tcosm->SetBranchAddress("nTrk",&nTracks);
  tcosm->SetBranchAddress("frTrkMu",&fracnMuTrk);
  tcosm->SetBranchAddress("mDist",&meanDist);
  tcosm->SetBranchAddress("fintev",&flagintev);
  tcosm->SetBranchAddress("nMuons",&nMuons);
  tcosm->SetBranchAddress("trigFlag",&triggerFlag);
  tcosm->SetBranchAddress("uptracmu",&uptracmu); 
  tcosm->SetBranchAddress("dwtracmu",&dwtracmu); 
  tcosm->SetBranchAddress("upthetaCos",upthetaCosmic);
  tcosm->SetBranchAddress("upphiCos",upphiCosmic);
  tcosm->SetBranchAddress("upphiCoscor",upphiCosmiccor);
  tcosm->SetBranchAddress("uptheta",uptheta);
  tcosm->SetBranchAddress("upphi",upphi);
  tcosm->SetBranchAddress("dwthetaCos",dwthetaCosmic);
  tcosm->SetBranchAddress("dwphiCos",dwphiCosmic);
  tcosm->SetBranchAddress("dwphiCoscor",dwphiCosmiccor);
  tcosm->SetBranchAddress("dwtheta",dwtheta);
  tcosm->SetBranchAddress("dwphi",dwphi);  
  tcosm->SetBranchAddress("upSign",upSign);
  tcosm->SetBranchAddress("dwSign",dwSign);  
  tcosm->SetBranchAddress("upxv",upxv);
  tcosm->SetBranchAddress("upyv",upyv);
  tcosm->SetBranchAddress("upzv",upzv);
  tcosm->SetBranchAddress("dwxv",dwxv);
  tcosm->SetBranchAddress("dwyv",dwyv);
  tcosm->SetBranchAddress("dwzv",dwzv);
  tcosm->SetBranchAddress("updEdx",uptpcSignal);
  tcosm->SetBranchAddress("dwdEdx",dwtpcSignal);  
  tcosm->SetBranchAddress("Pup",Pup);
  tcosm->SetBranchAddress("Pdw",Pdw);
  tcosm->SetBranchAddress("Sig1Ptup",Sig1Ptup);
  tcosm->SetBranchAddress("Sig1Ptdw",Sig1Ptdw);
  tcosm->SetBranchAddress("uppx",uppx);
  tcosm->SetBranchAddress("uppy",uppy);
  tcosm->SetBranchAddress("uppz",uppz);
  tcosm->SetBranchAddress("dwpx",dwpx);
  tcosm->SetBranchAddress("dwpy",dwpy);
  tcosm->SetBranchAddress("dwpz",dwpz);  
  tcosm->SetBranchAddress("Pw",Pmu);
  tcosm->SetBranchAddress("PMed",PmuMed);
  tcosm->SetBranchAddress("Pcov",Pcov);
  tcosm->SetBranchAddress("Pcovupdw",Pcovupdw);
  tcosm->SetBranchAddress("Pcovdwup",Pcovdwup);
  tcosm->SetBranchAddress("Pull",Xinv);
  tcosm->SetBranchAddress("PRes",PRes);
  tcosm->SetBranchAddress("NclsUp",uptpcNcls);
  tcosm->SetBranchAddress("NclsDw",dwtpcNcls);
  tcosm->SetBranchAddress("Chi2Up",uptpcChi2);
  tcosm->SetBranchAddress("Chi2Dw",dwtpcChi2);
  tcosm->SetBranchAddress("corrTimeTOFup",corrTimeTOFup);
  tcosm->SetBranchAddress("corrTimeTOFdw",corrTimeTOFdw);
  tcosm->SetBranchAddress("ToT_TOFup",ToT_TOFup);
  tcosm->SetBranchAddress("ToT_TOFdw",ToT_TOFdw);
  tcosm->SetBranchAddress("DzTOFup",DzTOFup);
  tcosm->SetBranchAddress("DzTOFdw",DzTOFdw);
  tcosm->SetBranchAddress("DxTOFup",DxTOFup);
  tcosm->SetBranchAddress("DxTOFdw",DxTOFdw);
  tcosm->SetBranchAddress("Xinup",Xinup);
  tcosm->SetBranchAddress("Yinup",Yinup);
  tcosm->SetBranchAddress("Zinup",Zinup);
  tcosm->SetBranchAddress("Xoutup",Xoutup);
  tcosm->SetBranchAddress("Youtup",Youtup);
  tcosm->SetBranchAddress("Zoutup",Zoutup);
  tcosm->SetBranchAddress("Xindw",Xindw);
  tcosm->SetBranchAddress("Yindw",Yindw);
  tcosm->SetBranchAddress("Zindw",Zindw);
  tcosm->SetBranchAddress("Xoutdw",Xoutdw);
  tcosm->SetBranchAddress("Youtdw",Youtdw);
  tcosm->SetBranchAddress("Zoutdw",Zoutdw);
  tcosm->SetBranchAddress("FileName",&fFileName); 
  tcosm->SetBranchAddress("FileNameChunk",&fFileNameChunk); 
  tcosm->SetBranchAddress("EventInFile",&fEventInFile); 
  tcosm->SetBranchAddress("hitsACO",hitsACO);
  tcosm->SetBranchAddress("mcnACO",&mcnACO);
  tcosm->SetBranchAddress("flagV0",&flagV0);
  tcosm->SetBranchAddress("upxdir",upxdir);
  tcosm->SetBranchAddress("upydir",upydir);
  tcosm->SetBranchAddress("upzdir",upzdir);
  tcosm->SetBranchAddress("uptetadir",uptetadir);
  tcosm->SetBranchAddress("upphidir",upphidir);  
  tcosm->SetBranchAddress("dwxdir",dwxdir);
  tcosm->SetBranchAddress("dwydir",dwydir);
  tcosm->SetBranchAddress("dwzdir",dwzdir);
  tcosm->SetBranchAddress("dwtetadir",dwtetadir);
  tcosm->SetBranchAddress("dwphidir",dwphidir);


  // Definition other variables
  Float_t  tof_mu; 

                //  All the plots done with Nmu>0           
 
  //    Plots to count the different triggers
  TH1F *trigger = new TH1F("trigger"," Trigger Flag ",1000,0,1001); 

  //    Plots for the ratio between "interaction events" and EAS events (Nmu>0)
  TH1F *flaginteas = new TH1F("flaginteas"," Flag events: EAS=0  Interaction=10",21,-5,15); 

  //  Plots Muon Multiplicity Distribution for different triggers
  TH1F *mumult = new TH1F("mumult"," Muon Multiplicity Distribution (all triggers) ",300,-0.5,299.5);
  TH1F *mumultOB3 = new TH1F("mumultOB3"," Muon Multiplicity Distribution (only TOF OB3) ",300,-0.5,299.5);
  TH1F *mumultAMU = new TH1F("mumultAMU"," Muon Multiplicity Distribution (only ACORDE) ",300,-0.5,299.5);

  //  Plots spatial distribution of the muons at Y=0 plane
  TH2F *hspatialup = new TH2F("hspatialup","Spatial Dist. of muons at Y=0 plane (up)",200,-300,300,200,-300,300);
  TH2F *hspatialdw = new TH2F("hspatialdw","Spatial Dist. of muons at Y=0 plane (dw)",200,-300,300,200,-300,300);
  TH1F *yzeroup = new TH1F("yzeroup"," Fluctuation plane Y=0 (up)",200,-100,100); 
  TH1F *yzerodw = new TH1F("yzerodw"," Fluctuation plane Y=0 (dw)",200,-100,100); 
  TH1F *deltay = new TH1F("deltay"," DealtaY (Yup-Ydw) ",100,-10,10); 


  //  Plots angular distribution of the muons (up and down)
  TH2F *tetaphiup= new TH2F("tetaphiup","Theta vs Phi (up)",100,0,360,90,0,90);
  TH2F *tetaphidw= new TH2F("tetaphidw","Theta vs Phi (dw)",100,0,360,90,0,90);

  //  Plots momentum distribution of the muons and charge
  TH1F *momup = new TH1F("momup","Momentum Distribution (up) ",400,0,400); 
  TH1F *momdw = new TH1F("momdw","Momentum Distribution (dw) ",400,0,400);
  TH1F *momcov = new TH1F("momcov","Momentum Distribution (Pcov=matched) ",400,0,400);
  TH1F *mucharge = new TH1F("mucharge","Charge of the muons (=5 wrong) ",20,-10,10); 

  //            Study the DeltaTime of each muon with the TOF information
  TH1F *tofmu = new TH1F("tofmu"," Time of Flight mu (dw-up) [ns]",400,-200,200); 


  //     ---- Loop over events ----
  cout << " Total Number of events in ESD =  " << nESDentries << endl;
 
  for(Int_t iEntry=0;iEntry<nESDentries;iEntry++){
//  for(Int_t iEntry=0;iEntry<10000;iEntry++){

//    cout << "Actual Event Number =  " << iEntry << endl;

    tcosm->GetEntry(iEntry);

    if(nMuons>0){ 

//      cout << "Actual Event Number (Nmuons>0) =  " << iEntry << endl;

      trigger->Fill(triggerFlag);
  
      flaginteas->Fill(flagintev*10);

      mumult->Fill(nMuons);
      if(triggerFlag==512)mumultOB3->Fill(nMuons);
      if(triggerFlag==4)mumultAMU->Fill(nMuons);

      for(Int_t j=0; j<nMuons ; j++){
        if(upyv[j]!=-999&&dwyv[j]!=-999){
          hspatialup->Fill(upzv[j],upxv[j]);
          yzeroup->Fill(upyv[j]);
          hspatialdw->Fill(dwzv[j],dwxv[j]);
          yzerodw->Fill(dwyv[j]);
          deltay->Fill(upyv[j]-dwyv[j]);
        }  // END if(upyv[j]!=-999&&dwyv[j]!=-999) 
      
        tetaphiup->Fill(upphidir[j],uptetadir[j]);
        tetaphidw->Fill(dwphidir[j],dwtetadir[j]);
     
        if(magfield!=0){
          if(Pup[j]>0)momup->Fill(Pup[j]); 
          if(Pdw[j]>0)momdw->Fill(Pdw[j]);
          if(Pcov[j]>0){
            momcov->Fill(Pcov[j]);
            if(upSign[j]==-1&&dwSign[j]==1) {      // mu+
              mucharge->Fill(1);
            }
            if(upSign[j]==1&&dwSign[j]==-1) {      // mu-
              mucharge->Fill(-1);
            }  
            if((upSign[j]==1&&dwSign[j]==1)||(upSign[j]==-1&&dwSign[j]==-1)){
              mucharge->Fill(5);    // wrong charge undefined
            }
          }   // END if(Pcov[j]>0)
        } // END if(magfield>0)

        if(corrTimeTOFdw[j]!=-999&&corrTimeTOFdw[j]!=-10&&corrTimeTOFup[j]!=-999&&corrTimeTOFup[j]!=-10){
//        if(corrTimeTOFdw[j]>0&&corrTimeTOFup[j]>0){
          tof_mu= (corrTimeTOFdw[j]-corrTimeTOFup[j])/1000; // time of flight in ns
          tofmu->Fill(tof_mu);
        }   

      }  //  END for(Int_t j=0; j<nMuons ; j++) 

    }  //  END if(nMuons>0)

  }  // END for(Int_t iEntry=0;iEntry<nESDentries;iEntry++)

  //    -----  Draw the plots  ----

  //     Trigger
  TCanvas *t1 = new TCanvas("Trigger Flag All Events ","",200,10,700,500);
  t1->cd();
  trigger->GetXaxis()->SetTitle("Trigger Flag");
  trigger->GetXaxis()->SetLabelSize(0.03);
  trigger->GetXaxis()->SetTickLength(0.02);
  trigger->GetYaxis()->SetTitle("Number of events");
  trigger->GetYaxis()->SetLabelSize(0.03);
  trigger->GetYaxis()->SetTickLength(0.02);
  trigger->SetMarkerStyle(21);
  trigger->SetMarkerColor(kBlue);
  trigger->SetLineColor(kBlue);
  trigger->Draw();

  //   Count interaction events and EAS events
  TCanvas *ie1 = new TCanvas("Flag interaction=10 and eas=0 events","",200,10,700,500);
  ie1->cd();
  flaginteas->GetXaxis()->SetTitle("Flag");
  flaginteas->GetXaxis()->SetLabelSize(0.03);
  flaginteas->GetXaxis()->SetTickLength(0.02);
  flaginteas->GetYaxis()->SetTitle("Number of events");
  flaginteas->GetYaxis()->SetLabelSize(0.03);
  flaginteas->GetYaxis()->SetTickLength(0.02);
  flaginteas->SetMarkerStyle(21);
  flaginteas->SetMarkerColor(kBlue);
  flaginteas->SetLineColor(kBlue);
  flaginteas->Draw();

  //     Muon Multiplicity Distribution (MMD)
  TCanvas *m1 = new TCanvas("Muon Multiplicity Distribution ","",200,10,700,500);
  m1->Divide(1,3);
  m1->cd(1);
  gPad->SetLogy();
  mumult->GetXaxis()->SetTitle("Number of muons");
  mumult->GetXaxis()->SetLabelSize(0.03);
  mumult->GetXaxis()->SetTickLength(0.02);
  mumult->GetYaxis()->SetTitle("Number of events");
  mumult->GetYaxis()->SetLabelSize(0.03);
  mumult->GetYaxis()->SetTickLength(0.02);
  mumult->SetMarkerStyle(21);
  mumult->SetMarkerColor(kRed);
  mumult->SetLineColor(kRed);
  mumult->Draw("E");
  m1->cd(2);
  gPad->SetLogy();  
  mumultOB3->GetXaxis()->SetTitle("Number of muons");
  mumultOB3->GetXaxis()->SetLabelSize(0.03);
  mumultOB3->GetXaxis()->SetTickLength(0.02);
  mumultOB3->GetYaxis()->SetTitle("Number of events");
  mumultOB3->GetYaxis()->SetLabelSize(0.03);
  mumultOB3->GetYaxis()->SetTickLength(0.02);
  mumultOB3->SetMarkerStyle(22);
  mumultOB3->SetMarkerColor(kBlue);
  mumultOB3->SetLineColor(kBlue);
  mumultOB3->Draw("E");
  m1->cd(3);
  gPad->SetLogy();  
  mumultAMU->GetXaxis()->SetTitle("Number of muons");
  mumultAMU->GetXaxis()->SetLabelSize(0.03);
  mumultAMU->GetXaxis()->SetTickLength(0.02);
  mumultAMU->GetYaxis()->SetTitle("Number of events");
  mumultAMU->GetYaxis()->SetLabelSize(0.03);
  mumultAMU->GetYaxis()->SetTickLength(0.02);
  mumultAMU->SetMarkerStyle(23);
  mumultAMU->SetMarkerColor(kGreen);
  mumultAMU->SetLineColor(kGreen);
  mumultAMU->Draw("E");

  //   Spatial Distribution at plane Y=0 of the muons
  TCanvas *sp1= new TCanvas("sp1","Spatial Distribution mu ",200,10,700,500);
  sp1->Divide(1,2);
  sp1->cd(1);  
  hspatialup->SetXTitle("Z [cm]");
  hspatialup->SetYTitle("X [cm]");
  hspatialup->SetMarkerStyle(1);
  hspatialup->SetMarkerColor(1);
  hspatialup->Draw();
  sp1->cd(2);  
  hspatialdw->SetXTitle("Z [cm]");
  hspatialdw->SetYTitle("X [cm]");
  hspatialdw->SetMarkerStyle(1);
  hspatialdw->SetMarkerColor(1);
  hspatialdw->Draw();

  TCanvas *sp2= new TCanvas("sp2","Projection Spatial Distribution mu ",200,10,700,500);
  sp2->Divide(2,2);
  sp2->cd(1);
  hspatialup->SetLineColor(1);
  hspatialup->ProjectionX()->Draw();
  hspatialup->SetLineColor(2);
  hspatialup->ProjectionY()->Draw("same");
  sp2->cd(2);
  hspatialdw->SetLineColor(1);
  hspatialdw->ProjectionX()->Draw();
  hspatialdw->SetLineColor(2);
  hspatialdw->ProjectionY()->Draw("same");
  sp2->cd(3);
  yzeroup->GetXaxis()->SetTitle("Y [cm]");
  yzeroup->GetYaxis()->SetTitle("Counts");
  yzeroup->SetLineColor(1);
  yzeroup->Draw();
  yzerodw->SetLineColor(2);
  yzerodw->Draw("same");
  sp2->cd(4);
  deltay->GetXaxis()->SetTitle("DeltaY [cm]");
  deltay->GetYaxis()->SetTitle("Counts");
  deltay->SetLineColor(4);
  deltay->Draw(); 


  //     Angular distribution : theta vs phi
  TCanvas *ang1= new TCanvas("ang1","Theta vs Phi ",200,10,700,500);
  ang1->Divide(1,2);
  ang1->cd(1);  
  tetaphiup->SetXTitle("Phi [deg]");
  tetaphiup->SetYTitle("Theta [deg]");
  tetaphiup->SetMarkerStyle(1);
  tetaphiup->SetMarkerColor(2);
  tetaphiup->Draw("COLZ");
  ang1->cd(2);  
  tetaphidw->SetXTitle("Phi [deg]");
  tetaphidw->SetYTitle("Theta [deg]");
  tetaphidw->SetMarkerStyle(1);
  tetaphidw->SetMarkerColor(2);
  tetaphidw->Draw("COLZ");

  TCanvas *ang2= new TCanvas("ang2","Projection Theta vs Phi ",200,10,700,500);
  ang2->Divide(2,2);
  ang2->cd(1);
  tetaphiup->ProjectionX()->Draw();
  ang2->cd(2);
  tetaphiup->ProjectionY()->Draw();
  ang2->cd(3);
  tetaphidw->ProjectionX()->Draw();
  ang2->cd(4);
  tetaphidw->ProjectionY()->Draw();

    //       Momentum distribution and charge 
    TCanvas *mom1= new TCanvas("mom1","Momentum distribution and charge of muons",200,10,700,500);
    mom1->Divide(2,2);
    mom1->cd(1); 
    gPad->SetLogy();   
    momup->GetXaxis()->SetTitle("Pup [GeV/c]");
    momup->GetXaxis()->SetLabelSize(0.03);
    momup->GetXaxis()->SetTickLength(0.02);
    momup->GetYaxis()->SetTitle("Counts");
    momup->GetYaxis()->SetLabelSize(0.03);
    momup->GetYaxis()->SetTickLength(0.02);
    momup->SetMarkerStyle(21);
    momup->SetMarkerColor(kBlue);
    momup->SetLineColor(kBlue);
    momup->Draw();
    mom1->cd(2);
    gPad->SetLogy();  
    momdw->GetXaxis()->SetTitle("Pdw [GeV/c]");
    momdw->GetXaxis()->SetLabelSize(0.03);
    momdw->GetXaxis()->SetTickLength(0.02);
    momdw->GetYaxis()->SetTitle("Counts");
    momdw->GetYaxis()->SetLabelSize(0.03);
    momdw->GetYaxis()->SetTickLength(0.02);
    momdw->SetMarkerStyle(21);
    momdw->SetMarkerColor(kBlue);
    momdw->SetLineColor(kBlue);
    momdw->Draw();
    mom1->cd(3);
    gPad->SetLogy();  
    momcov->GetXaxis()->SetTitle("Pcov [GeV/c]");
    momcov->GetXaxis()->SetLabelSize(0.03);
    momcov->GetXaxis()->SetTickLength(0.02);
    momcov->GetYaxis()->SetTitle("Counts");
    momcov->GetYaxis()->SetLabelSize(0.03);
    momcov->GetYaxis()->SetTickLength(0.02);
    momcov->SetMarkerStyle(21);
    momcov->SetMarkerColor(kBlue);
    momcov->SetLineColor(kBlue);
    momcov->Draw();
    mom1->cd(4);
    mucharge->GetXaxis()->SetTitle("charge [GeV/c]");
    mucharge->GetXaxis()->SetLabelSize(0.03);
    mucharge->GetXaxis()->SetTickLength(0.02);
    mucharge->GetYaxis()->SetTitle("Counts");
    mucharge->GetYaxis()->SetLabelSize(0.03);
    mucharge->GetYaxis()->SetTickLength(0.02);
    mucharge->SetMarkerStyle(21);
    mucharge->SetMarkerColor(kBlue);
    mucharge->SetLineColor(kBlue);
    mucharge->Draw();

    //        Plot Time of Flight of muons and other TOF information
    TCanvas *mutof = new TCanvas("Time of Flight of muons","",200,10,700,500);
    mutof->cd(1);
    tofmu->GetXaxis()->SetTitle("Time of Flight of muons");
    tofmu->GetXaxis()->SetLabelSize(0.03);
    tofmu->GetXaxis()->SetTickLength(0.02);
    tofmu->GetYaxis()->SetTitle("Num. of muons");
    tofmu->GetXaxis()->SetTitle("tof [ns]");
    tofmu->GetYaxis()->SetLabelSize(0.03);
    tofmu->GetYaxis()->SetTickLength(0.02);
    tofmu->Draw();

  saveplot->Write();

} // END void AnaTreeCosmics()
