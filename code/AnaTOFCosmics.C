//
//     Program to read and analyse the TREE created 
//     by AliCosmics.cxx from ESD for 
//     cosmic muons in the central barrel  (AB)   
//     First Version : 20/Jan/2016 
//     Last Version : 
//

void AnaTOFCosmics(const char* filename){

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
  Float_t uptrackLength[500], dwtrackLength[500], upfTexp[500], dwfTexp[500];

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
  tcosm->SetBranchAddress("uptrackLength",uptrackLength);
  tcosm->SetBranchAddress("dwtrackLength",dwtrackLength);
  tcosm->SetBranchAddress("upfTexp",upfTexp);
  tcosm->SetBranchAddress("dwfTexp",dwfTexp);


  // Definition other variables
  Float_t  tof_mu;
  Int_t nSingle = 0;

                //  All the plots done with Nmu>0           
 
  //            Study the DeltaTime of each muon with the TOF information
  TH1F *tofmu = new TH1F("tofmu"," Time of Flight mu (dw-up) [ns]",400,-200,200); 
  TH1F *intlen = new TH1F("intlen","Integrated length for single 2-track muons", 100, 0., 5000.);
  TH1F *toflen = new TH1F("toflen","TOF time length",1000,-200000.,200000.);

  //     ---- Loop over events ----
  cout << " Total Number of events in ESD =  " << nESDentries << endl;
 
  for(Int_t iEntry=0;iEntry<nESDentries;iEntry++){
    tcosm->GetEntry(iEntry);

    if(nMuons == 1 && nTracks == 2){ 
      nSingle++;

      Float_t intLen = uptrackLength[0] + dwtrackLength[0];
      intlen->Fill(intLen);

      Float_t tofLen = corrTimeTOFdw[0] - corrTimeTOFup[0];
      toflen->Fill(tofLen);

    }  //  END if(nMuons==1)

  }  // END for(Int_t iEntry=0;iEntry<nESDentries;iEntry++)

  //    -----  Draw the plots  ----

  printf("Number of single 2-track muons: %d\n",nSingle);

  TCanvas *c1 = new TCanvas();
  c1->cd();
  intlen->Draw();

  TCanvas *c2 = new TCanvas();
  c2->cd();
  toflen->Draw();

  saveplot->Write();

} // END void AnaTreeCosmics()
