//////////////////////////////////////////////////////////////////////////////////////////////////////
//           Analysis Macro for Cosmics Data
//           Matching of the tracks to reconstruct the muons (up and down track)
//           Single track muons
//           Count the total number of muon per event
//           Write the measured variables of each muon in a TREE (root file)
//           Original version : 17/Jan/2011
//            M.Subieta INFN TORINO
//            B.Alessandro INFN TORINO
//		Update: May 30th 2012
//		Adding:
//--> branches: chunk path, # of event, ACORDE hits, ACORDE MCN and VZERO BG rejection
//	From: Mario Rodriguez Cahuantiz (FCFM-BUAP), Puebla-MX <mrodrigu@mail.cern.ch>
//		===================================
//		Update: Oct 15th 2014
//		Modified:
//			--> fixed the double tree filling: fDatree->Fill()
//              Update: Oct 20th 2015 : Teta e Phi with direction cosines
//              Add 5 variables up and 5 variables down at the tree
//              Update: Jan 27th 2016 : Info in the tree only when Nmu>0
//              Reduced the variables in the TREE
//              Update: 12th February 2016
//              Added TOF variables in the TREE -> Pad hits for each track
//              Eliminated TRD and others in TriggerFlag and added ASL=16 
//
//              Update: Jul 13th 2017
//              Small corrections on vector vecmatch keeping into account
//              the vecmatch that has to be reset at -10 see: 13/jul/2017 
//   
//		
//		Update: October 29th 2017
//		removing obsolete getter for ESD event. Adding GetESDEvent() (Mario Rodríguez)
//
//        6/Aug/2018 :
//        Added   covariant charge: utrack->Charge() in the tree
//        Delete some variables in the tree (see  //%%2018)
//
/////////////////////////////////////////////////////////////////////////////////////////////////////
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include <sstream>
#include <string>
#include "AliStack.h"
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TFile.h"
#include "TArrayI.h"
#include "TVectorD.h"
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TLatex.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisManager.h"

#include "AliESDVertex.h"
#include "AliESDtrack.h"

#include "AliESDEvent.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliTracker.h"
#include "AliCosmics.h"

#include <AliRunLoader.h>
#include <AliCluster.h>
//#include <AliCosmicTracker.h>
#include <AliESDCosmicTrack.h>

#include "AliMagF.h"
//#include "AliTPCExBExact.h"

#include <stdio.h>
#include "AliESDVZERO.h"
#include "AliESDACORDE.h"
ClassImp(AliCosmics)

//________________________________________________________________________
AliCosmics::AliCosmics(const char *name, const Bool_t isMC) 
: AliAnalysisTask(name,""),
  fESD(0),
  fListHist(0),
  fDatree(0),
  fTrackThetaCosmicOB0(0),
  fTrackThetaCosmicOB1(0),
  fTrackThetaCosmicOB3(0),
  fTrackPhiCosmic(0),
  fTriggerMask(0),
  fnTrksSPD(0),
  TotalTrksSPD(0),
  fnTrksTOFOB0(0),
  TotalTrksTOFOB0(0),
  fnTrksTOFOB1(0),
  TotalTrksTOFOB1(0),
  fnTrksTOFOB3(0),
  TotalTrksTOFOB3(0),
  fnTrksACORDEAMU(0),
  TotalTrksACORDEAMU(0),
  fnTrksTOFOB0Muons(0),
  fnTrksTOFOB1Muons(0),
  fnTrksTOFOB3Muons(0),
  fnTrksAMUMuons(0),
  fnTrksASLMuons(0),
  fnTrksSCOMuons(0),
  fmeanDistMu3(0),
  fmeanDistMu4(0),
  fdistHist(0),
  Debtrackused(0),
  Debtrackusedmatch(0),
  Debtrackusedsing(0),
  fmumult(0), 
  fphiDist(0), 
  fthetaDist(0),
  fupdwMuDist(0),
  fmatchZX(0), 
  fsingleZX(0), 
  fsingleupZX(0), 
  fsingledwZX(0),
  triggerMask(0),
  triggerFlag(0),
  nEvent(0),
  nRun(0),
  nTracks(0),
  nMuons(0),
  nMuonsMC(0),
  ntpar(0),
  maxntpar(0),
  ntmaxpar(0),
  ntup(0),
  ntdw(0),
  flagintev(0),
  itag(0), ntmatch(0), ntsingle(0),
  magfield(0),
  timestp(0),
  meanDist(0),fracnMuTrk(0),rmsDist(0),rdist1t2(0),rdistmatch(0),cost1t2(0),rDist(0),
  //%%pxz(0),E(0),P(0),Pt1(0),Pt2(0),Et1(0),Et2(0),
  E(0),P(0),Pt1(0),Pt2(0),Et1(0),Et2(0),
  fFileName(""), 
  fFileNameChunk(""),
  eventNumberInFile(0),
  fEventInFile(0),
  mcnACO(0),
  flagV0(0),
  mcData(0)
 // nMuonsMC(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(0, TTree::Class());
  DefineOutput(1, TList::Class());

  currchunk[0] = '\0';
  prevchunk[0] = '\0';

  mcData = isMC;
}
//________________________________________________________________________
AliCosmics::~AliCosmics(){
	// Destructor
	if (fListHist){
	 delete fListHist;
	 fListHist = 0;
	}
}
/*
//________________________________________________________________________
Bool_t AliCosmics::GetESDEvent()
{
        fESD = (AliESDEvent*)InputEvent();
        if (!fESD)
        {
                Printf("ERROR: ESD event not available ... check your inputs .. !!!!");
                return kFALSE;
        }
        return kTRUE;
}

//________________________________________________________________________
*/
void AliCosmics::ConnectInputData(Option_t *) 
{
  
  //================================================================================================ 
  // 			         Connect ESD or AOD here
  //================================================================================================ 

  TTree* tree = dynamic_cast<TTree*> (GetInputData(0));
  if (!tree) {
    Printf("ERROR: Could not read chain from input slot 0");
  } else {

    AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*> (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    esdH->SetNeedField();

    if (!esdH) {
      Printf("ERROR: Could not get ESDInputHandler");
    } else
 
      // Getting the pointer of ESD event from input Handler 
 
 //!     fESD   = dynamic_cast<AliESDEvent*>(AliESDEvent*)esdH->GetEvent();
 //!        fESD = esdH->GetEvent();
 //      fESD = (AliESDEvent*)InputEvent();
      fESD   = (AliESDEvent*)esdH->GetEvent();

  }
}

//________________________________________________________________________
void AliCosmics::CreateOutputObjects()
{

 //============================================================================================ 
  // 				Creating output files
  //============================================================================================ 

  filehmmeout = fopen ("highnumtracks.out","at");  
  intereve = fopen ("interaction_events.out","at");
  fprintf(intereve,"N.Int.|N.run|N.ev|Timest|N.tr|N.mu|Fracmu|MDist|File \n");  
 
  //================================================================================================= 
  // 				Create histograms in TList
  //================================================================================================= 

  fListHist = new TList();
  
  // Debug plots
  Debtrackused = new TH1F("Freq.track used","Frequency each track; Track Number; N.times used",702,-1,701); 
  Debtrackusedmatch = new TH1F("Freq.track used MATCHED","Frequency each track; Track Number; N.times used",702,-1,701); 
  Debtrackusedsing = new TH1F("Freq.track used SINGLE","Frequency each track; Track Number; N.times used",702,-1,701);   

//   Other plots 
 
  fnTrksSPD = new TH1F("nTrksSPD","Track Multiplicity SCO; No. of Tracks; No. of Events",200,-0.5,199.5);
  fListHist->Add(fnTrksSPD);
 
  fnTrksTOFOB0 = new TH1F("nTrksTOFOB0","Track Multiplicity OB0; No. of Tracks; No. of Events",200,-0.5,199.5);
  fListHist->Add(fnTrksTOFOB0);
  
  fnTrksTOFOB1 = new TH1F("nTrksTOFOB1","Track Multiplicity OB1; No. of Tracks; No. of Events",200,-0.5,199.5);
  fListHist->Add(fnTrksTOFOB1);

  fnTrksTOFOB3 = new TH1F("nTrksTOFOB3","Track Multiplicity OB3; No. of Tracks; No. of Events",200,-0.5,199.5);
  fListHist->Add(fnTrksTOFOB3);
  
  fnTrksACORDEAMU = new TH1F("nTrksACORDEAMU","Track Multiplicity AMU; No. of Tracks; No. of Events",200,-0.5,199.5);
  fListHist->Add(fnTrksACORDEAMU);

  fnTrksTOFOB0Muons = new TH1F("nTrksTOFOB0Muons","Muon Multiplicity OB0; No. of Muons; No. of Events",500,0.5,500.5);
  fListHist->Add(fnTrksTOFOB0Muons);
 
  fnTrksTOFOB1Muons = new TH1F("nTrksTOFOB1Muons","Muon Multiplicity OB1; No. of Muons; No. of Events",500,0.5,500.5);
  fListHist->Add(fnTrksTOFOB1Muons);

  fnTrksTOFOB3Muons = new TH1F("nTrksTOFOB3Muons","Muon Multiplicity OB3; No. of Muons; No. of Events",500,0.5,500.5);
  fListHist->Add(fnTrksTOFOB3Muons);
 
  fnTrksSCOMuons = new TH1F("nTrksSCOMuons","Muon Multiplicity SCO; No. of Muons; No. of Events",500,0.5,500.5);
  fListHist->Add(fnTrksSCOMuons);

  fnTrksAMUMuons = new TH1F("nTrksAMUMuons","Muon Multiplicity AMU; No. of Muons; No. of Events",500,0.5,500.5);
  fListHist->Add(fnTrksAMUMuons);
 
  fnTrksASLMuons = new TH1F("nTrksASLMuons","Number of Muons Triggered by ACORDE; Muon Multiplicity; No. of Events",500,0.5,500.5);
  fListHist->Add(fnTrksASLMuons);
 
  fTrackThetaCosmicOB0 = new TH1F("TrackThetaCosmicOB0","Muon Zenithal Angle Distribution; Theta(degrees)",90,0,90);
  fListHist->Add(fTrackThetaCosmicOB0);
 
  fTrackThetaCosmicOB1 = new TH1F("TrackThetaCosmicOB1","Muon Zenith Angle Distribution; Theta(degrees)",90,0,90);
  fListHist->Add(fTrackThetaCosmicOB1);

  fTrackThetaCosmicOB3 = new TH1F("TrackThetaCosmicOB3","Muon Zenith Angle Distribution; Theta(degrees)",90,0,90);
  fListHist->Add(fTrackThetaCosmicOB3);
 
  fTrackPhiCosmic = new TH1F("TrackPhiCosmic","Muon Azimuthal Angle Distribution; Phi(degrees)",360,0,360);
  fListHist->Add(fTrackPhiCosmic);

  fdistHist = new TH1F("distHist","Distance Parallel Tracks Respect to t1; distance (cm); Parallel tracks",500,0,500);
  fListHist->Add(fdistHist);

  fTriggerMask = new TH1F("TriggerMask","Trigger Mask; TriggerMask number; Number of Events",30,1,30);
  fListHist->Add(fTriggerMask);

  const Int_t nbins = 500;
  Double_t xmin = 0.1;
  Double_t xmax = 30;
  Double_t logxmin = TMath::Log10(xmin);
  Double_t logxmax = TMath::Log10(xmax);
  Double_t binwidth = (logxmax-logxmin)/nbins;
  Double_t xbins[nbins+1];
  xbins[0] = xmin;
  for (Int_t i=1;i<=nbins;i++){
    xbins[i]=xmin+TMath::Power(10,logxmin+i*binwidth);
  }

//       Study of the interaction events i comparison with multimuons

  fmeanDistMu3 = new TH2F("meanDistMu3"," Interaction Events; Mean_{Distance}(cm); N_{Tracks}/N_{mu}",500,0,500,nbins,xbins);
  fmeanDistMu4 = new TH2F("meanDistMu4","Multimuon Events; Mean_{Distance}(cm); N_{Tracks}/N_{mu}",500,0,500,nbins,xbins);
  fListHist->Add(fmeanDistMu3);   //Intercation events
  fListHist->Add(fmeanDistMu4);   // Multimuons
  
  // Plots prepared by Katherin (Thesis) 
  
  fmumult = new TH1F("mumult","Muon Mult. Dist.",300,-0.5,299.5);   // "mumult"
  fListHist->Add(fmumult);
  
  fthetaDist = new TH1F("ThetaDist","",50,0.5,51.5);  //30,-0.5,59.5   " Cosmic Tracks Zenith Angle Dist. (Deg) "
  fListHist->Add(fthetaDist);
 
  fphiDist = new TH1F("PhiDist","",90,0,360);  //100,100.5,300.5 " Cosmic Tracks Azimuth Angle Dist. (Deg) "
  fListHist->Add(fphiDist);
  
  fupdwMuDist = new TH1F("updwMuDist"," Distance between matched tracks up and down (cm) ",35, -0.05, 3.45);  //35, 0, 3.5
  fListHist->Add(fupdwMuDist);
  
  fmatchZX = new TH2F("matchZX","Spatial dist. of matched muons at y = 0 plane; Z[cm]; X[cm]",200,-300,300,200,-300,300);
  fListHist->Add(fmatchZX);
  
  fsingleZX = new TH2F("singleZX","Spatial dist. of single muons at y = 0 plane; Z[cm]; X[cm]",200,-300,300,200,-300,300);
  fListHist->Add(fsingleZX);
  
  fsingleupZX = new TH2F("singleupZX","Spatial dist. of single muons UP at y = 0 plane; Z[cm]; X[cm]",200,-300,300,200,-300,300);
  fListHist->Add(fsingleupZX);
  
  fsingledwZX = new TH2F("singledwZX","Spatial dist. of single muons DW at y = 0 plane; Z[cm]; X[cm]",200,-300,300,200,-300,300);
  fListHist->Add(fsingledwZX);
  
  

  //============================================================================================ 
  // 				Creating the data tree
  //============================================================================================ 

  fDatree = new TTree("fDatree","dataCosmicsTree");
 
  fDatree->Branch("nRun",&nRun,"nRun/I");
  fDatree->Branch("chnk",chnk,"chnk[50]/C");
  fDatree->Branch("nEve",&nEvent,"nEve/I");
  
  fDatree->Branch("timestp",&timestp,"timestp/i");
  fDatree->Branch("magfield",&magfield,"magfield/F");
  
  fDatree->Branch("nTrk",&nTracks,"nTrk/I");
  fDatree->Branch("frTrkMu",&fracnMuTrk,"frTrkMu/F");
  fDatree->Branch("mDist",&meanDist,"mDist/F");
  fDatree->Branch("fintev",&flagintev,"fintev/I");
  
  fDatree->Branch("nMuons",&nMuons,"nMuons/I");
  fDatree->Branch("trigFlag",&triggerFlag,"trigFlag/I");
  
  // uptracmu[rmu] ==> is the track number of the muon rmu for up
  //  ex.  Supposing we have reconstructed 10 muons , some are matched, 
  //  others are single track (or up or down), 
  //  uptracmu[5]=16   the track 16 is the track up of the muon n.5
  // if dwtracmu[5]=-10 means the muon n. 5 is not matched with only
  // track up, if dwtracmu[5]=22 means that track 16 and 22 are matched 
  // and form the muon n.5. The same things for the other variables. 
  fDatree->Branch("uptracmu",&uptracmu,"uptracmu[nMuons]/I"); 
  fDatree->Branch("dwtracmu",&dwtracmu,"dwtracmu[nMuons]/I"); 
  
  fDatree->Branch("upthetaCos",upthetaCosmic,"upthetaCos[nMuons]/F");
  fDatree->Branch("upphiCos",upphiCosmic,"upphiCos[nMuons]/F");
  //%%fDatree->Branch("upphiCoscor",upphiCosmiccor,"upphiCoscor[nMuons]/F");
  fDatree->Branch("uptheta",uptheta,"uptheta[nMuons]/F");
  fDatree->Branch("upphi",upphi,"upphi[nMuons]/F");
  fDatree->Branch("dwthetaCos",dwthetaCosmic,"dwthetaCos[nMuons]/F");
  fDatree->Branch("dwphiCos",dwphiCosmic,"dwphiCos[nMuons]/F");
  //%%fDatree->Branch("dwphiCoscor",dwphiCosmiccor,"dwphiCoscor[nMuons]/F");
  fDatree->Branch("dwtheta",dwtheta,"dwtheta[nMuons]/F");
  fDatree->Branch("dwphi",dwphi,"dwphi[nMuons]/F");
  
  fDatree->Branch("upSign",upSign,"upSign[nMuons]/F");
  fDatree->Branch("dwSign",dwSign,"dwSign[nMuons]/F");
  
  fDatree->Branch("upxv",upxv,"upxv[nMuons]/F");
  fDatree->Branch("upyv",upyv,"upyv[nMuons]/F");
  fDatree->Branch("upzv",upzv,"upzv[nMuons]/F");
  fDatree->Branch("dwxv",dwxv,"dwxv[nMuons]/F");
  fDatree->Branch("dwyv",dwyv,"dwyv[nMuons]/F");
  fDatree->Branch("dwzv",dwzv,"dwzv[nMuons]/F");
  
 //%%2018  fDatree->Branch("updEdx",uptpcSignal,"updEdx[nMuons]/F");
 //%%2018  fDatree->Branch("dwdEdx",dwtpcSignal,"dwdEdx[nMuons]/F");
  
  fDatree->Branch("Pup",Pup,"Pup[nMuons]/F");
  fDatree->Branch("Pdw",Pdw,"Pdw[nMuons]/F");
  fDatree->Branch("Sig1Ptup",Sig1Ptup,"Sig1Ptup[nMuons]/F");
  fDatree->Branch("Sig1Ptdw",Sig1Ptdw,"Sig1Ptdw[nMuons]/F");
  fDatree->Branch("uppx",uppx,"uppx[nMuons]/F");
  fDatree->Branch("uppy",uppy,"uppy[nMuons]/F");
  fDatree->Branch("uppz",uppz,"uppz[nMuons]/F");
  fDatree->Branch("dwpx",dwpx,"dwpx[nMuons]/F");
  fDatree->Branch("dwpy",dwpy,"dwpy[nMuons]/F");
  fDatree->Branch("dwpz",dwpz,"dwpz[nMuons]/F");
  
  //%%2018 fDatree->Branch("Pw",Pmu,"Pw[nMuons]/F");
  //%%fDatree->Branch("PMed",PmuMed,"PMed[nMuons]/F");
  fDatree->Branch("Pcov",Pcov,"Pcov[nMuons]/F");
  fDatree->Branch("Chargecov",Chargecov,"Chargecov[nMuons]/F");
  //%%fDatree->Branch("Pcovupdw",Pcovupdw,"Pcovupdw[nMuons]/F");
  //%%fDatree->Branch("Pcovdwup",Pcovdwup,"Pcovdwup[nMuons]/F");
  //%%2018  fDatree->Branch("Pull",Xinv,"Pull[nMuons]/F");
  //%%2018  fDatree->Branch("PRes",PRes,"PRes[nMuons]/F");
  
  fDatree->Branch("NclsUp",uptpcNcls,"NclsUp[nMuons]/I");
  fDatree->Branch("NclsDw",dwtpcNcls,"NclsDw[nMuons]/I");
  fDatree->Branch("Chi2Up",uptpcChi2,"Chi2Up[nMuons]/F");
  fDatree->Branch("Chi2Dw",dwtpcChi2,"Chi2Dw[nMuons]/F");
  
  fDatree->Branch("corrTimeTOFup",corrTimeTOFup,"corrTimeTOFup[nMuons]/F");
  fDatree->Branch("corrTimeTOFdw",corrTimeTOFdw,"corrTimeTOFdw[nMuons]/F");
  //%%fDatree->Branch("ToT_TOFup",ToT_TOFup,"ToT_TOFup[nMuons]/F");
  //%%fDatree->Branch("ToT_TOFdw",ToT_TOFdw,"ToT_TOFdw[nMuons]/F");
  //%%fDatree->Branch("DzTOFup",DzTOFup,"DzTOFup[nMuons]/F");
  //%%fDatree->Branch("DzTOFdw",DzTOFdw,"DzTOFdw[nMuons]/F");
  //%%fDatree->Branch("DxTOFup",DxTOFup,"DxTOFup[nMuons]/F");
  //%%fDatree->Branch("DxTOFdw",DxTOFdw,"DxTOFdw[nMuons]/F");
//MS-28/9/18
  fDatree->Branch("nMatchup",nMatchup,"nMatchup[nMuons]/b");
  fDatree->Branch("nMatchdw",nMatchdw,"nMatchdw[nMuons]/b");
  fDatree->Branch("indexTOFup",indexTOFup,"indexTOFup[nMuons]/I");
  fDatree->Branch("indexTOFdw",indexTOFdw,"indexTOFdw[nMuons]/I");
//MS-28/9/18
  
  fDatree->Branch("Xinup",Xinup,"Xinup[nMuons]/F");
  fDatree->Branch("Yinup",Yinup,"Yinup[nMuons]/F");
  fDatree->Branch("Zinup",Zinup,"Zinup[nMuons]/F");
  fDatree->Branch("Xoutup",Xoutup,"Xoutup[nMuons]/F");
  fDatree->Branch("Youtup",Youtup,"Youtup[nMuons]/F");
  fDatree->Branch("Zoutup",Zoutup,"Zoutup[nMuons]/F");
  fDatree->Branch("Xindw",Xindw,"Xindw[nMuons]/F");
  fDatree->Branch("Yindw",Yindw,"Yindw[nMuons]/F");
  fDatree->Branch("Zindw",Zindw,"Zindw[nMuons]/F");
  fDatree->Branch("Xoutdw",Xoutdw,"Xoutdw[nMuons]/F");
  fDatree->Branch("Youtdw",Youtdw,"Youtdw[nMuons]/F");
  fDatree->Branch("Zoutdw",Zoutdw,"Zoutdw[nMuons]/F");
  fDatree->Branch("FileName",&fFileName); 
  fDatree->Branch("FileNameChunk",&fFileNameChunk); 
  fDatree->Branch("EventInFile",&fEventInFile); 
  fDatree->Branch("hitsACO",hitsACO,"hitsACO[60]/I");
  fDatree->Branch("mcnACO",&mcnACO,"mcnACO/I");
  fDatree->Branch("flagV0",&flagV0,"flagV0/I");

  fDatree->Branch("upxdir",upxdir,"upxdir[nMuons]/F");
  fDatree->Branch("upydir",upydir,"upydir[nMuons]/F");
  fDatree->Branch("upzdir",upzdir,"upzdir[nMuons]/F");
  fDatree->Branch("uptetadir",uptetadir,"uptetadir[nMuons]/F");
  fDatree->Branch("upphidir",upphidir,"upphidir[nMuons]/F");  
  fDatree->Branch("dwxdir",dwxdir,"dwxdir[nMuons]/F");
  fDatree->Branch("dwydir",dwydir,"dwydir[nMuons]/F");
  fDatree->Branch("dwzdir",dwzdir,"dwzdir[nMuons]/F");
  fDatree->Branch("dwtetadir",dwtetadir,"dwtetadir[nMuons]/F");
  fDatree->Branch("dwphidir",dwphidir,"dwphidir[nMuons]/F");

  //! Monte Carlo branches
  if (mcData) {
    fDatree->Branch("nMuonsMC",&nMuonsMC,"nMuonsMC/I");
    fDatree->Branch("pMC",pMC,"pMC[nMuons]/F"); // p momentum 
    fDatree->Branch("pxMC",pxMC,"pxMC[nMuons]/F"); // px momentum 
    fDatree->Branch("pyMC",pyMC,"pyMC[nMuons]/F"); // py momentum 
    fDatree->Branch("pzMC",pzMC,"pzMC[nMuons]/F"); // pz momentum 
    fDatree->Branch("xMC",xMC,"xMC[nMuons]/F"); // x-muon 
    fDatree->Branch("yMC",yMC,"yMC[nMuons]/F"); // y-muon 
    fDatree->Branch("zMC",zMC,"zMC[nMuons]/F"); // z-muon 
    fDatree->Branch("tetaMC",tetaMC,"tetaMC[nMuons]/F"); // theta
    fDatree->Branch("energyMC",energyMC,"energyMC[nMuons]/F"); // energy
    fDatree->Branch("pdgCode",pdgCode,"pdgCode[nMuons]/F"); // pdg code	
  }
  // Post data here as request by AliAnalysisManager - M.S.
  PostData(0, fDatree);
  PostData(1, fListHist);
}

//________________________________________________________________________
void AliCosmics::Exec(Option_t *) 
{

  if (!fESD) {
    Printf("ERROR: fESD not available");
    return;
  }
/*
 if (!GetESDEvent())
        {
                PostData(0,fDatree);
                PostData(1,fListHist);
                return;
        }
*/
        //! pointer to analysis manager to access to the PID utils

        AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
        AliInputEventHandler *inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());

        fESD->InitMagneticField();


  	for (Int_t i=0; i<60; i++) hitsACO[i]=0;
	flagV0 = 0;
  
	// ACORDE information to the fDatree

	AliESDACORDE *acordeESD = fESD->GetACORDEData();
	Int_t acoMulti=0;
	mcnACO = 0;
	for(Int_t i=0; i<60; i++){
		hitsACO[i] = acordeESD->GetHitChannel(i);
		if (acordeESD->GetHitChannel(i)){
			acoMulti++;
		}
	}mcnACO = acoMulti;
	flagV0+=1;
	// VZero BG information to reject BG

	AliESDVZERO* vzeroring = fESD->GetVZEROData();

	if(AliCosmics::VZeroBG(vzeroring)) flagV0+=2;
	if(AliCosmics::VZeroBG(vzeroring) && acoMulti>0) flagV0+=4;


  Int_t nevdeb1 = -10;   // lower event number to be debugged 
  Int_t nevdeb2 = -10;   // upper event number to be debugged

  Int_t ntrkdeb = 400;   // debug events having more than ntrkdeb tracks 
//  Int_t ntrkdeb = 4;   // debug events having more than ntrkdeb tracks 
// Initialization for interaction events
  flagintev = 0;  // tag interaction events with flagintev=1 

  
  triggerMask  = fESD->GetTriggerMask();
  nEvent = fESD->GetEventNumberInFile();
  nRun   = fESD->GetRunNumber();
  nTracks = fESD->GetNumberOfTracks();
  magfield = fESD->GetMagneticField();
  timestp= fESD->GetTimeStamp();

  if (magfield == 0) minCutEnergy = 0.;  // M.S. 15 dec 2015
  
  TFile *curfile = AliAnalysisManager::GetAnalysisManager()->GetTree()->GetCurrentFile();
  eventNumberInFile++;

  TString fileName = curfile->GetName();
        stringstream ss;//create a stringstream
 //       ss << eventNumberInFile-1;//add number to the stream
        //printf("Event number: %s\n",ss);
        //cout << "SS: " << ss.str() << endl;
        //TString chunk = ss.str();


	fEventInFile = fESD->GetEventNumberInFile();
	ss << fEventInFile;
        TString chunk = ss.str();
        //TString sumaChunks = fileName+chunk;
        TString sumaChunks = fileName+chunk;
        fFileName.SetString(fileName.Data());
        fFileNameChunk.SetString(sumaChunks.Data());
//	cout << "fFileName: " << fileName.Data() << endl;
//	cout << "fFileNameChunk: " << sumaChunks.Data() << endl;
	cout << "Event number " << chunk.Data() << endl;

     
  nMuons = 0;
  nMuonsMC = 0;
 
  Float_t ptotup,ptotdw;
  ptotup = -999;
  ptotdw = -999;
       
  for (Int_t i = 0; i < 500; i++){
    xMC[i] = -9999.;
    yMC[i] = -9999.;
    zMC[i] = -9999.;

    pMC[i] = -9999.;
    pxMC[i] = -9999.;
    pyMC[i] = -9999.;
    pzMC[i] = -9999.;


    tetaMC[i] = -9999.;
    energyMC[i] = -9999.;
    pdgCode[i] = -9999.;

 //!   xMC[i] = yMC[i] = zMC[i] = pMC[i] = pxMC[i] = pyMC[i] = pzMC[i] = tetaMC[i] = energyMC[i] = pdgCode[i] = -999.;
    uptracmu[i] = -10;
    dwtracmu[i] = -10; 
    upthetaCosmic[i] = -10.;
    upphiCosmic[i] = -10.;
 //%%   upphiCosmiccor[i] = -10.;
    uptheta[i] = -10.;
    upphi[i] = -10.;
    dwthetaCosmic[i] = -10.;
    dwphiCosmic[i] = -10.;
 //%%   dwphiCosmiccor[i] = -10.;   
    dwtheta[i] = -10.;
    dwphi[i] = -10.;
    upSign[i] = -10.;
    dwSign[i] = -10.;
    upxv[i] = -999.;
    upyv[i] = -999.;
    upzv[i] = -999.;
    dwxv[i] = -999.;
    dwyv[i] = -999.;
    dwzv[i] = -999.;
    uppx[i] = -999.;
    uppy[i] = -999.;
    uppz[i] = -999.;
    dwpx[i] = -999.;
    dwpy[i] = -999.;
    dwpz[i] = -999.;
//%%2018    uptpcSignal[i] = -10.;	
//%%2018    dwtpcSignal[i] = -10.;	
    Pup[i] = -10.;
    Pdw[i] = -10.;
//%%2018    Pmu[i] = -10.;
//%%    PmuMed[i] = -10.;
    Pcov[i] = -10.;
    Chargecov[i] = -10.;
//%%    Pcovupdw[i] = -10.;
//%%    Pcovupdw[i] = -10.;
//%%2018    PRes[i] = -10.;
    Sig1Ptup[i] = -10.;
    Sig1Ptdw[i] = -10.;
//%%2018    Xinv[i] = -999.;
    SigUp2[i] = -10.;
    SigDw2[i] = -10.;
    uptpcNcls[i] = -10;
    dwtpcNcls[i] = -10; 
    uptpcChi2[i] = -10.;
    dwtpcChi2[i] = -10.;
//    upratioChi2cls[i] = -10.;
//    dwratioChi2cls[i] = -10.;
    corrTimeTOFup[i] = -999;
    corrTimeTOFdw[i] = -999;	
//%%    ToT_TOFup[i] = -999.;
//%%    ToT_TOFdw[i] = -999.;
//%%    DzTOFup[i] = -999.;
//%%    DzTOFdw[i] = -999.;
//%%    DxTOFup[i] = -999.;
//%%    DxTOFdw[i] = -999.;
//MS-28/9/18
    nMatchup[i] = -10;
    nMatchdw[i] = -10;
    indexTOFup[i] = -10.;
    indexTOFdw[i] = -10.;
//MS-28/9/18
    Xinup[i] = -999;
    Yinup[i] = -999;
    Zinup[i] = -999;
    Xoutup[i] = -999;
    Youtup[i] = -999;
    Zoutup[i] = -999;
    Xindw[i] = -999;
    Yindw[i] = -999;
    Zindw[i] = -999;
    Xoutdw[i] = -999;
    Youtdw[i] = -999;
    Zoutdw[i] = -999;
    
    MuDistupdw[i] = -10;
    Xpos[i] = -999;
    Zpos[i] = -999;
    
    upxdir[i] = -999.;
    upydir[i] = -999.;
    upzdir[i] = -999.;
    dwxdir[i] = -999.;
    dwydir[i] = -999.;
    dwzdir[i] = -999.;
    uptetadir[i] = -999.;
    upphidir[i] = -999.;
    dwtetadir[i] = -999.;
    dwphidir[i] = -999.;
    
    
  }
  
  //============================================================================================ 
  //			   Printing to output files & Filling Histograms 
  //============================================================================================ 
  
  strcpy(currchunk,curfile->GetName());
/*
  if (strcmp(currchunk,prevchunk) != 0) {
    if (neverFill)
      fprintf(stderr,"No filled events for chunk %s\n",prevchunk);
      
    printf("===> currchunk is %s\n",currchunk);

    // Extract the chunk substring from the whole file name - M.S. 01 feb 2017
    Int_t startChnk, endChnk;
    for (Int_t i=42; i<strlen(currchunk); i++)
      if (currchunk[i] == '/') {
	startChnk = i+1;
	break;
      }
    for (Int_t i=startChnk; i<strlen(currchunk); i++)
      if (currchunk[i] == '/') {
	endChnk = i-1;
	break;
      }
    for (Int_t i=0; i<(endChnk-startChnk+1); i++)
      chnk[i] = currchunk[startChnk+i];
    chnk[endChnk-startChnk+1] = '\0';
    printf("===> chnk is %s\n",chnk);

    strcpy(prevchunk, currchunk);
    neverFill = kTRUE;
  }
*/
  
  TString nameTrigger = fESD->GetFiredTriggerClasses();
  Char_t Triggername[1000];

  sprintf(Triggername,"%s",nameTrigger.Data());
  
 // triggerFlag  for 2009-2010 :  OB3=1  OCP=2  AMU=4  SCO=8
 // triggerFlag for 2011 data : OB0=1  OB1=2  AMU=4  SCO=8  we add the TRD trigger, HWV=16
 // triggerFlag  important for 2015   AMU=4  and  OB3=512 
  
  triggerFlag = 0;

//  Fill number of Tracks in the event also if zero. This count the
//  number of triggers (entries of the plot)

  if ((nameTrigger.Contains("OB0"))){
    fnTrksTOFOB0->Fill(nTracks);
    triggerFlag=triggerFlag+1;  
  }

  if ((nameTrigger.Contains("OB1"))){
    fnTrksTOFOB1->Fill(nTracks); 
    triggerFlag=triggerFlag+2; 
  }
  
  if ((nameTrigger.Contains("AMU"))){
    fnTrksACORDEAMU->Fill(nTracks);
    triggerFlag=triggerFlag+4;  
  }	
  
  if ((nameTrigger.Contains("SCO"))){
    fnTrksSPD->Fill(nTracks); 
    triggerFlag=triggerFlag+8; 
  }

  if ((nameTrigger.Contains("ASL"))){
    fnTrksSPD->Fill(nTracks); 
    triggerFlag=triggerFlag+16; 
  }

  if ((nameTrigger.Contains("OB3"))){
    fnTrksTOFOB3->Fill(nTracks);
    triggerFlag=triggerFlag+512;
  }

  if ((nameTrigger.Contains("OBE"))){
    triggerFlag=triggerFlag+1024;
  }


if (nameTrigger.Contains("LSR")) return;


  //    Analyze only Events with these triggers : 
  
//  if ((nameTrigger.Contains("OB3"))||(nameTrigger.Contains("OBE"))||(nameTrigger.Contains("SH2"))||((nameTrigger.Contains("SCO1")))||(nameTrigger.Contains("EMC"))||(nameTrigger.Contains("HWU"))||(nameTrigger.Contains("TRD")) || (nameTrigger.Contains("OB1")) || (nameTrigger.Contains("OB0")) || (nameTrigger.Contains("AMU")) || (nameTrigger.Contains("SCO")) ||  (nameTrigger.Contains("ASL"))){
  //if ((nameTrigger.Contains("OB3"))||(nameTrigger.Contains("TRD"))||(nameTrigger.Contains("OB0"))||(nameTrigger.Contains("AMU"))||(nameTrigger.Contains("ASL"))){

//if ((nameTrigger.Contains("OB0"))||(nameTrigger.Contains("OB1"))||(nameTrigger.Contains("OB2"))||(nameTrigger.Contains("OB3"))||(nameTrigger.Contains("AMU"))||(nameTrigger.Contains("ASL")) ||  (nameTrigger.Contains("SCO")))   
{

    // write output event if it has more than ntrkdeb tracks
  if(nTracks >= ntrkdeb) {
     nevdeb1 = nEvent;
     nevdeb2 = nEvent;  
  }
          
  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  
      cout << endl;
      cout << " DEBUG EVENT : " << nEvent << "   that has " << nTracks << " TRACKS " << endl; 
      cout << " ===== RUN : " << nRun << " Chunk : " << chnk << " Event " << nEvent << " ===== " << endl;   
      cout << "Triggername :  " << Triggername << " Triggermask : " << triggerMask << "  TriggeFLag : " << triggerFlag << endl;
      cout << "Time Stamp :   " << timestp << "   Magnetic Field :  " << magfield << endl;
      cout << " Used CUTS for muons:  N.of TPC cluster matched tracks= " << cutnumcls << "   N.of TPC cluster single-track= " << cutsingncls << "   Energy cut in GeV= " << minCutEnergy << "   Parallelis parameter= " << fCutDirPar << "   Max Distance matched tracks in cm= " << minCutDist << endl;     
      cout << " Current File : " << curfile->GetName() << endl;
      cout << endl;
  }
    
    //      Analyze Events with at least one track in the TPC 
    if (nTracks > 0){
      
      //============================================================================================ 
      // 			         Loop on Tracks 				     	 
      //============================================================================================ 
      
      Int_t vecpar[nTracks], vecmatch[nTracks], vecup[nTracks], vecdw[nTracks], vecpar_u[nTracks], vecup_u[nTracks], vecdw_u[nTracks], tagtrk[nTracks];
      vector<Int_t> vecmu, trkmatch;
      Float_t distmatch[nTracks];
      Float_t distmatch_p[nTracks];
      Int_t freqtrack[nTracks];
      
      Double_t dirt1[3], dirt2[3]; 
      Int_t nTracksm05;
      Int_t oldt1,oldt2;
      
      // Initialize values for each vector
      // Explanation of the matching vectors and variables
      // ntmaxpar = track with the maximum number of parallel tracks
      // maxntpar = max number of parallel tracks ==> ntmaxpar has maxntpar parallel tracks
      // vecmatch[t1] = the track t1 matches with the track vecmatch[t1]
      // vecdist[t1] = distance between t1 and the matched track vecmatch[t1]
      // vecpar[i] = list of the tracks parallel to ntmaxpar
      // vecup[i] = list of the up tracks parallel to ntmaxpar
      // vecdw[i] = list of the down tracks parallel to ntmaxpar
      // vecpar_u, vecup_u, vecdw_u same meaning but provisional vectors
      // vecmu[i] = final list of all the parallel tracks used and tagged as good 
      // vecmu[0] = ntmaxpar , vecmu[1]=vecpar[0], ....
      // trkmatch[i] = list of all matched tracks to check and cut when already used
      
      for (Int_t it = 0; it < nTracks; it++){
	
	vecpar[it]   = -10;
	vecup[it]    = -10;
	vecdw[it]    = -10;
	vecmatch[it] = -10;
        distmatch[it] = -10;
	distmatch_p[it]  = -10;
	vecpar_u[it] = -10;
	vecup_u[it]  = -10;
	vecdw_u[it]  = -10;
	tagtrk[it]   = -10;
	freqtrack[it] = 0;
      }

      maxntpar =  0; //max num. parallel tracks
      ntmaxpar = -10 ; //num. of the track with more parallel tracks
      
//      if(nEvent>=nevdeb1&&nEvent<=nevdeb2){   
//        cout << " #*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#* " << endl;
//	cout << " nTracks > 0 --- Starting various cases ----"<< " nTracks = " << nTracks << " nEvent : " << nEvent << endl; 
//      }
      
      //  $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      // Explanation of the matching algorithm and the muon reconstruction
      // Particular case : event with 1 track. If the track passes the cuts
      // the track is considered a muon and Nmuons = 1.
      // Standard case more than 1 track    if (nTracks > 1).
      // Loop on track t1 and the Loop on track t2=t1+1.
      // Sign the eventual track t2 matching t1  vecmatch[t1]=t2.
      // List of the tracks parallel to t1 vecpar[i].
      // After all the loop on t2 and on t1 we have choosen the
      // tracks with the maximum number of parallel tracks = ntmaxpar
      // and we have filled the vector for this track (ntmaxpar).
      // =======================================================
      // If ntmaxpar >=0 (exist) ==> 
      // Calculus of ntmatch
      // vecmatch[ntmatch] if exist and
      // loop on vecpar[i] to find the eventual matching track
      // of each parallel track : vecmatch[vecpar[i]]
      // at the end we have ntmatch the number of match 
      // every two tracks matched ntmatch increases by 1
      // nMuons = maxntpar + 1 - ntmatch - itag   is the number of muons of the event 
      // ========================================================
      // If nMuons > 0   ====>
      // Sign multimuon events and Interaction events and use only multimuon events
      // For each muon write the variables in the tree
      
      // $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
      
      
      //  Start with the particular case of having only 1 Track in the TPC
     
      if (nTracks == 1){
	
	for (Int_t t1 = 0; t1 < nTracks; t1++){	
	  
	  AliESDtrack *tr1 = fESD->GetTrack(t1);
	  
	  if(tr1 == NULL || tr1->GetOuterParam() == 0 || tr1->GetInnerParam() == 0)
	    continue;
	  
	  if(tr1->GetTPCNcls() < cutsingncls)
	    continue;
	  
//	  if((tr1->GetTPCchi2()/tr1->GetTPCNcls())>=cutchi2track)
//	    continue;
	  
	  Float_t xvert = tr1->Xv(); 
//	  Float_t yvert = tr1->Yv(); 
	  Float_t zvert = tr1->Zv(); 
	  E             = tr1->E();
	  P             = tr1->P();

//         No this cut with B=OFF	  
	  if (E < minCutEnergy) continue; 
	  
	  if (TMath::Abs(zvert) > zcut || TMath::Abs(xvert) > xcut) 
	    continue;
	  
	  TVectorD pxpypz(3);
	  TVectorD xyzv(3);
	  TVectorD pxpypz1(3);
	  TVectorD xyzv1(3);
	  TVectorD xyzout(3);
          TVectorD xyzin(3);	 
 
	  
	  Int_t rmu = 0 ;    // Count the number of muons for the index vector
          Int_t rmureal = 1 ;  // Count real number of muons
	  fracnMuTrk = 1.; // In this case 1 track and 1 muon
          meanDist = -10.;  // Inizialization
	  
	  nMuons = rmureal;
	  
	  if (tr1->GetOuterParam()->GetAlpha() > 0)  // up track tr1
	    {
	      AliExternalTrackParam *uptr1 = (AliExternalTrackParam *)tr1;      
	      uptr1->PxPyPz(pxpypz.GetMatrixArray());
	      uptr1->XvYvZv(xyzv.GetMatrixArray());

              tr1->GetOuterXYZ(xyzout.GetMatrixArray());
	      tr1->GetInnerXYZ(xyzin.GetMatrixArray()); 
	      tr1->GetDirection(dirt1); 

	      uptracmu[rmu] = t1;   // track number 
	      upxv[rmu] = xyzv[0];
	      upyv[rmu] = xyzv[1];
	      upzv[rmu] = xyzv[2];
	      uptheta[rmu] = RadGrad*uptr1->Theta();
	      upphi[rmu] = RadGrad*uptr1->Phi();
	      uppx[rmu] = pxpypz[0];  
	      uppy[rmu] = pxpypz[1];
	      uppz[rmu] = pxpypz[2];  
	      Pup[rmu] = uptr1->P();
//%%2018	      Pmu[rmu] = Pup[rmu];
//%%	      PmuMed[rmu] = Pup[rmu];
//%%2018	      Xinv[rmu] = -666.;  
//%%	      ToT_TOFup[rmu] = tr1->GetTOFsignalToT();		     		     
	      corrTimeTOFup[rmu] = tr1->GetTOFsignal();	
//%%	      DzTOFup[rmu] = tr1->GetTOFsignalDz();		     
//%%	      DxTOFup[rmu] = tr1->GetTOFsignalDx();		     
//MS-28/9/18
	      nMatchup[rmu] = tr1->GetTOFclusterN();
	      indexTOFup[rmu] = tr1->GetTOFCalChannel();
//MS-28/9/18
	      Xoutup[rmu] = xyzout[0];	
	      Youtup[rmu] = xyzout[1];	
	      Zoutup[rmu] = xyzout[2];	
	      Xinup[rmu] = xyzin[0];	
	      Yinup[rmu] = xyzin[1];	
	      Zinup[rmu] = xyzin[2];	
	      upxdir[rmu] = dirt1[0];
	      upydir[rmu] = dirt1[1];
	      upzdir[rmu] = dirt1[2];

//    Calculate the teta and phi of the up track with directions 
              uptetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(upydir[rmu]));
              if(upydir[rmu]>=0){   
              upphidir[rmu] = RadGrad*TMath::ACos(upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(upxdir[rmu]<0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(upxdir[rmu]>0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }
              if(upydir[rmu]<0){   
              upphidir[rmu] = RadGrad*TMath::ACos(-upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(-upxdir[rmu]<0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(-upxdir[rmu]>0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }


	      if (Pup[rmu]>0) upthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(uppy[rmu])/Pup[rmu]);
	      if (uppz[rmu]>0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu]);
	      if (uppz[rmu]<0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu])+180;
	      
//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(uppx[rmu]*uppx[rmu]+uppz[rmu]*uppz[rmu]);

//%%
/*
              if(uppz[rmu]>0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= RadGrad*TMath::ACos(uppx[rmu]/pxz);      

              if(uppz[rmu]>0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(uppz[rmu]/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(uppx[rmu])/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(uppz[rmu])/pxz);  
*/
//%%              
	      
	      upSign[rmu] = uptr1->GetSign();
//%%2018	      uptpcSignal[rmu] = tr1->GetTPCsignal();
	      uptpcNcls[rmu] = tr1->GetTPCNcls();
	      uptpcChi2[rmu] = tr1->GetTPCchi2();
	      //		upratioChi2cls[rmu] = uptpcChi2[rmu]/uptpcNcls[rmu];
	      Sig1Ptup[rmu] = tr1->GetSigma1Pt2();
	      // Variable used but not to record in the tree
	      Ptupinv[rmu] = uptr1->OneOverPt();  
	      Pup2[rmu] = uppx[rmu]*uppx[rmu]+uppy[rmu]*uppy[rmu]+uppz[rmu]*uppz[rmu];
	      
	      SigUp2[rmu] = (uptr1->GetSigmaY2()*uppy[rmu]*uppy[rmu]/Pup2[rmu])+(uptr1->GetSigmaZ2()*uppz[rmu]*uppz[rmu]/Pup2[rmu]);  
	      	      
	      ntmaxpar = 0;

//  Fill number of Muons (nTracks=nMuons=1 in this case) in the event with 
//  at least one muon. This count the number of trigger with at least one 
//  muon in TPC (entries of the plot)
    	      
	      if ((nameTrigger.Contains("OB0"))){
	        fnTrksTOFOB0Muons->Fill(nTracks);
	      }		
	      if ((nameTrigger.Contains("OB1"))){
		fnTrksTOFOB1Muons->Fill(nTracks);
	      }
	      if ((nameTrigger.Contains("OB3"))){
		fnTrksTOFOB3Muons->Fill(nTracks);
	      }
	      if ((nameTrigger.Contains("SCO"))){
	        fnTrksSCOMuons->Fill(nTracks);
	      }
	      if ((nameTrigger.Contains("AMU"))){
	        fnTrksAMUMuons->Fill(nTracks);
	      }

	    }  // END up track tr1

	  if (tr1->GetOuterParam()->GetAlpha() < 0)  // down track tr1
	    {
	      AliExternalTrackParam *dwtr1 = (AliExternalTrackParam *)tr1;      
	      dwtr1->PxPyPz(pxpypz1.GetMatrixArray());
	      dwtr1->XvYvZv(xyzv1.GetMatrixArray());

	      tr1->GetOuterXYZ(xyzout.GetMatrixArray());
              tr1->GetInnerXYZ(xyzin.GetMatrixArray());	
  	      tr1->GetDirection(dirt1); 

	      dwtracmu[rmu] = t1;   // track number 
	      dwxv[rmu] = xyzv1[0];
	      dwyv[rmu] = xyzv1[1];
	      dwzv[rmu] = xyzv1[2];
	      dwtheta[rmu] = RadGrad*dwtr1->Theta();
	      dwphi[rmu] = RadGrad*dwtr1->Phi();
	      dwpx[rmu] = pxpypz1[0];  
	      dwpy[rmu] = pxpypz1[1];
	      dwpz[rmu] = pxpypz1[2];  
	      Pdw[rmu] = dwtr1->P();
//%%2018	      Pmu[rmu] = Pdw[rmu];
//%%	      PmuMed[rmu] = Pdw[rmu];
//%%2018	      Xinv[rmu] = -666;
//%%	      ToT_TOFdw[rmu] = tr1->GetTOFsignalToT();		     		     
	      corrTimeTOFdw[rmu] = tr1->GetTOFsignal();	
//%%	      DzTOFdw[rmu] = tr1->GetTOFsignalDz();		     
//%%	      DxTOFdw[rmu] = tr1->GetTOFsignalDx();		     
//MS-28/9/18
	      nMatchdw[rmu] = tr1->GetTOFclusterN();
	      indexTOFdw[rmu] = tr1->GetTOFCalChannel();
//MS-28/9/18
	      Xoutdw[rmu] = xyzout[0];	
	      Youtdw[rmu] = xyzout[1];	
	      Zoutdw[rmu] = xyzout[2];	
	      Xindw[rmu] = xyzin[0];	
	      Yindw[rmu] = xyzin[1];	
	      Zindw[rmu] = xyzin[2];	
	      dwxdir[rmu] = dirt1[0];
	      dwydir[rmu] = dirt1[1];
	      dwzdir[rmu] = dirt1[2];

//    Calculate the teta and phi of the down track with directions for down change the sign
//       of the vectors of the direction in the formulas (dwxdir,dwydir,dwzdir)
//     Sometimes the TPC reconstruction gives the wrong sign to dwydir ??? so we check
//     it both for dw and up track.
              dwtetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwydir[rmu]));
              if(dwydir[rmu]<=0){
              dwphidir[rmu] = RadGrad*TMath::ACos(-dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(-dwxdir[rmu]<0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(-dwxdir[rmu]>0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }
         //    if wrong sign given by TPC
              if(dwydir[rmu]>0){
              dwphidir[rmu] = RadGrad*TMath::ACos(dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(dwxdir[rmu]<0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(dwxdir[rmu]>0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }


	      if (Pdw[rmu]>0) dwthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwpy[rmu])/Pdw[rmu]);
	      if (dwpz[rmu]<0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu]);
	      if (dwpz[rmu]>0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu])+180;

//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(dwpx[rmu]*dwpx[rmu]+dwpz[rmu]*dwpz[rmu]);

//%%
/*               
              if(dwpz[rmu]<0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]<0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);      
*/
//%%

	      dwSign[rmu] = tr1->GetSign();
//%%2018	      dwtpcSignal[rmu] = tr1->GetTPCsignal();
	      dwtpcNcls[rmu] = tr1->GetTPCNcls();
	      dwtpcChi2[rmu] = tr1->GetTPCchi2();
	      //		dwratioChi2cls[rmu] = dwtpcChi2[rmu]/dwtpcNcls[rmu];
	      Sig1Ptdw[rmu]= dwtr1->GetSigma1Pt2();
	      // Variable used but not to record in the tree
	      Ptdwinv[rmu] = dwtr1->OneOverPt();  
	      Pdw2[rmu]        = dwpx[rmu]*dwpx[rmu]+dwpy[rmu]*dwpy[rmu]+dwpz[rmu]*dwpz[rmu];
	      SigDw2[rmu] = (dwtr1->GetSigmaY2()*dwpy[rmu]*dwpy[rmu]/Pdw2[rmu])+(dwtr1->GetSigmaZ2()*dwpz[rmu]*dwpz[rmu]/Pdw2[rmu]); 

	      ntmaxpar = 0;


//  Fill number of Muons (nTracks=nMuons=1 in this case) in the event with 
//  at least one muon. This count the number of trigger with at least one 
//  muon in TPC (entries of the plot)
	
	      if ((nameTrigger.Contains("OB0"))){
	        fnTrksTOFOB0Muons->Fill(nTracks);
	      }		
	      if ((nameTrigger.Contains("OB1"))){
	        fnTrksTOFOB1Muons->Fill(nTracks);
	      }
	      if ((nameTrigger.Contains("OB3"))){
	        fnTrksTOFOB3Muons->Fill(nTracks);
	      }
	      if ((nameTrigger.Contains("SCO"))){
	        fnTrksSCOMuons->Fill(nTracks);
	      }
	      if ((nameTrigger.Contains("AMU"))){
	        fnTrksAMUMuons->Fill(nTracks);
	      }

	    }  // END down track tr1

	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){
	    cout << " uptracmu[rmu] = " << uptracmu[rmu] << "   dwtracmu[rmu] =  " << dwtracmu[rmu] << endl; 
	    cout << " upSign[rmu] = " << upSign[rmu] << "   dwSign[rmu] =  " << dwSign[rmu] << endl;   
	    cout << " ---End  Particular case only 1 Track ----"<< " nTracks = " << nTracks << " nEvent : " << nEvent << endl; 
	  }



	}   // END Loop for (Int_t t1 = 0; t1 < nTracks; t1++) (only one track in this case)

	if(nMuons>0){                   // Fill the TREE if Nmu>0
	  fDatree->Fill();
	  neverFill = kFALSE;
	}
	
      }//---------------------->>> END if Track == 1 <<<------------------------- 


      //   =========  STARTING LOOP ON THE TRACKS FOR MATCHING TRACKS num.tracks > 1 =========
      //      Loop on tracks :  loop on t1   and  loop on t2=t1+1
      
      if (nTracks > 1){ 

	nTracksm05 = 0;
	itag = 0; // count the tracks to be deleted
 
	for (Int_t t1 = 0; t1 < nTracks; t1++){
 
	  AliESDtrack *itrack = fESD->GetTrack(t1);
		    
          rdistmatch = 10000;  //init check distance t1 t2 for matching	    
		    
	  if(itrack == NULL || itrack->GetOuterParam() == 0 || itrack->GetInnerParam() == 0)
	    continue;

	  if(itrack->GetTPCNcls() < cutnumcls)  continue;
	  
//	  if((itrack->GetTPCchi2()/itrack->GetTPCNcls())>=cutchi2track)
//	    continue;

     
	  Float_t t1xvert = itrack->Xv(); 
	  Float_t t1yvert = itrack->Yv(); 
	  Float_t t1zvert = itrack->Zv(); 
	  Et1     	  = itrack->E();
	  Pt1     	  = itrack->P();


//         No this cut with B=OFF          
	  if(Et1 < minCutEnergy) continue; 
  
	  if (TMath::Abs(t1zvert) > zcut || TMath::Abs(t1xvert) > xcut) 
	    continue;
      
	nTracksm05++;
  
	itrack->GetDirection(dirt1);
  
	ntpar = 0; // count the parallel tracks to t1
	ntup  = 0;
	ntdw  = 0;

//   Loop on t2 to find all parallel tracks at t1

	for (Int_t t2 = t1 + 1; t2 < nTracks; t2++){

	  AliESDtrack *jtrack = fESD->GetTrack(t2);

	  if (jtrack == NULL || jtrack->GetOuterParam() == 0 || jtrack->GetInnerParam() == 0)	
	    continue;
    
	  if (jtrack->GetTPCNcls() < cutnumcls) 
	    continue;
	    
//	  if((jtrack->GetTPCchi2()/jtrack->GetTPCNcls())>=cutchi2track)
//	    continue;

   
	  Float_t t2xvert = jtrack->Xv(); 
	  Float_t t2yvert = jtrack->Yv(); 
	  Float_t t2zvert = jtrack->Zv(); 
	  Et2     	  = jtrack->E();
	  Pt2     	  = jtrack->P();

//         No this cut with B=OFF   
	  if(Et2 < minCutEnergy) continue; 
    
	  if (TMath::Abs(t2zvert) > zcut || TMath::Abs(t2xvert) > xcut) 
	    continue;

	  jtrack->GetDirection(dirt2);

	  // Parallelism --> Compute cost1t2   dot product for direction t1 // t2

	  cost1t2 = (dirt1[0]*dirt2[0]+dirt1[1]*dirt2[1]+dirt1[2]*dirt2[2]); 

	  if (TMath::Abs(cost1t2) < fCutDirPar) continue;   
 
	  //Check if t1 and t2 match (distance between them) 
	  // Compute rdist1t2 the distance between t1 and t2
  
	  rdist1t2 = sqrt((t1xvert-t2xvert)*(t1xvert-t2xvert) + (t1yvert-t2yvert)*(t1yvert-t2yvert) + (t1zvert-t2zvert)*(t1zvert-t2zvert));

//        MATCHING : Find and Tag the track t2 that matches t1 (vecmatch[t1]=t2)

// track t2 has passed all the cuts and is parallel to t1 ==> decide if t2 matches with t1 or not.
// t2 matches with t1 if t2 is the closest track to t1 and is different (up or down) from t1.
// Sign in tagtrk[itag] the matched track close to another matched track but more distant to match.
// Example track 1 match with track 0 , then we find that also track 2 matches with track 0,
// but track 2 is more distant than track 1. So the matching is 0=>1 and we sign in vector
// tagtrk[itag]=2 that has to be eliminated from parallel track vector vecmu. 

          oldt1 = -10;  // initialization
          oldt2 = -10;  // initialization
	  if(rdist1t2<minCutDist){  // First check if the distance between t1 and t2 is < minCutDist
	    Float_t t1alpha=itrack->GetOuterParam()->GetAlpha() ;
	    Float_t t2alpha=jtrack->GetOuterParam()->GetAlpha() ;
	    if((t1alpha>0&&t2alpha<0)||(t1alpha<0&&t2alpha>0)){ // t1=up and t2=down or viceversa
	      if(rdist1t2<rdistmatch){  // Check if the distance between t1 and t2 is the lowest
	         distmatch_p[t1] = rdist1t2;
	         distmatch_p[t2] = rdist1t2;
	         oldt1 = vecmatch[t2];
	         oldt2 = vecmatch[t1];
	         
	         // Check that new tracks t1 and t2 were not already used for a previous match
	         if(oldt1>=0){    // t2 already matched with old t1
	           if(distmatch_p[t1]<distmatch[oldt1]){
	             // new match t1-->t2 good and closer, cancel old match and tag old t1
	 // correction 13/jul/2017 tagtrk[itag] = vecmatch[oldt1];  // tag the old t1 matched track
                     tagtrk[itag] = oldt1;  // tag the old t1 matched track
	             itag++;
	           // 13/jul/2017  vecmatch[vecmatch[oldt1]] = -10; // delete old t2 match
	             vecmatch[oldt1] = -10;  // delete old t1 match
	             vecmatch[t1] = t2;  // new closer match
	             vecmatch[t2] = t1;  // new closer match
	           }else{
	             // new match t1-->t2 more distant than previous not taken and tag t1
	             tagtrk[itag] = t1;  // tag the new t1 
	             itag++;
	           }// END if(distmatch_p[t1]<distmatch[oldt1])
	         } // END  if(oldt1>0)
	         
	         if(oldt2>=0){    // t1 already matched with old t2
	           if(distmatch_p[t2]<distmatch[oldt2]){
	             // new match t1-->t2 good and closer, cancel old match and tag old t2
	       // correction 13/jul/2017      tagtrk[itag] = vecmatch[oldt2];  // tag the old t2 matched track
	             tagtrk[itag] = oldt2;  // tag the old t2 matched track
	             itag++;
	            // 13/jul/2017 vecmatch[vecmatch[oldt2]] = -10; // delete old t1 match
	             vecmatch[oldt2] = -10;  // delete old t2 match
	             vecmatch[t1] = t2;  // new closer match
	             vecmatch[t2] = t1;  // new closer match
	           }else{
	             // new macth t1-->t2 more distant than previous not taken and tag t2
	             tagtrk[itag] = t2;  // tag the new t2 
	             itag++;
	           }// END if(distmatch_p[t2]<distmatch[oldt2])
	         } // END  if(oldt2>0)
	         //   END Check new tracks t1 and t2 used or not
	         
	         if(oldt1<0&&oldt2<0){
	           vecmatch[t1] = t2;
	           vecmatch[t2] = t1;
	           distmatch_p[t1] = rdist1t2;
	           distmatch_p[t2] = rdist1t2;
	         }
	         
	      }  // END  if(rdist1t2<rdistmatch)
            } // END  if((t1alpha>0&&t2alpha<0)||(t1alpha<0&&t2alpha>0))
	  }else{       // if (rdist1t2<minCutDist)
	    if (jtrack->GetTPCNcls() < cutsingncls) {
	      continue; // cut in number of clusters for unmatched t2
	    }            
          }  // END if(rdist1t2<minCutDist)
          
       //   END MATCHING : Terminated the part of the matching 

	         
   //  Save the track t2 that is parallel to t1 in the vectors : vecpar, vecup, vecdw          

 	  vecpar_u[ntpar] = t2; // track t2 parallel to track t1
	  ntpar++; // increment the number of parallel tracks to t1
         
	  if (jtrack->GetOuterParam()->GetAlpha() > 0){
	    vecup_u[ntup] = t2; // uptrack t2
	    ntup++;
	  } 
	  
	  if (jtrack->GetOuterParam()->GetAlpha() < 0){
	    vecdw_u[ntdw] = t2; // downtrack t2
	    ntdw++;
	  }
  
	}// END loop t2   for (Int_t t2 = t1 + 1; t2 < nTracks; t2++) 
 
	// Store the information of the track with max. number of parallel tracks (Main Track)      

	if (ntpar>maxntpar){
 
	  maxntpar = ntpar ;  // ntpar = number of parallel tracks to t1  
	  ntmaxpar = t1 ; // ntmaxpar=track with max. num. of parallel tracks (Main Track)

	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  
	      cout << " Actual Main track (ntmaxpar) : " << ntmaxpar << "     Num. parallel tracks (maxntpar) : " << maxntpar << endl; 
	      cout << " Num. par. tracks up : " << ntup << "       Num. par. tracks dw  : " << ntdw << endl;
	  } 

	  for(Int_t iw = 0; iw < maxntpar; iw++){	
	    vecpar[iw]=vecpar_u[iw] ; // store in vecpar the num. of all the tracks parallel to ntmaxpar
	  }

	  for(Int_t iup = 0; iup < ntup; iup++){
	    vecup[iup]=vecup_u[iup] ; // store in vecup the num. of all the tracks up par. to ntmaxpar 
	  }
 	
	  for(Int_t idw = 0; idw < ntdw; idw++){
	    vecdw[idw]=vecdw_u[idw] ; // store in vecdw the num. of all the tracks down par. to ntmaxpar 
//           if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  
//	      cout << " idw : " << idw << " vecdw[idw] : " << vecdw[idw] << endl; 
//	    }
	  }
   
	}//END if ntpar>maxntpar
  
      }// END loop t1    (Int_t t1 = 0; t1 < nTracks; t1++)

      if(nEvent>=nevdeb1&&nEvent<=nevdeb2){   
	cout << "++++++++++++++++ Ending loop t1 ----------- " << endl;  
      }

//  We have terminated the double loop over all the tracks  [t1 --> (t2 --> End) End]
//  We have the following information :
//  The Main Track : ntmaxpar    track with maximum number of parallel tracks
//  A list of all the parallel tracks : vecpar[i]  i=0,maxntpar-1  
//  For each track t1 the matched track t2 (if it exists) : vecmatch[t1]=t2 and distance distmatch[t1] between them


//======= LOOP on parallel tracks :  ntmaxpar ==> itrack    track with max. number of parallel tracks ====

      if (ntmaxpar >= 0){
     
        if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  
     	  cout << endl; 
	  cout << " $$$$$$$$$ FOUND PARALLEL TRACKS ==> SUMMARY EVENT : " << nEvent  << " $$$$$$$$$ " << endl; 
	  cout << " Main Track (ntmaxpar) " << ntmaxpar << "    Num. parallel tracks (maxntpar) : " << maxntpar << endl;
        }
 		
	fdistHist->Reset();
	  
	AliESDtrack *itrack = fESD->GetTrack(ntmaxpar); //itrack=track with max. num. of parallel tracks

        
     //   ==========   Calculus of meanDist : multimu or interaction =============
     //  This variable is useful to decide if the event is a multimuon or interaction 

 	Float_t t1xvert  = itrack->Xv(); 
	Float_t t1yvert  = itrack->Yv(); 
	Float_t t1zvert  = itrack->Zv(); 
	itrack->GetDirection(dirt1);
  
	if (vecpar[0] >= 0 && maxntpar > 0){
		
	  for (Int_t i = 0; i < maxntpar; i++){
	    // vecpar[i] contains the list of all the tracks parallel to itrack
	    AliESDtrack *ttrack = fESD->GetTrack(vecpar[i]);

	    Float_t t2xvert = ttrack->Xv(); 
	    Float_t t2yvert = ttrack->Yv(); 
	    Float_t t2zvert = ttrack->Zv(); 
	    ttrack->GetDirection(dirt2);
  
  // Distance between Main Track and every parallel track  		
	    rDist = sqrt((t1xvert-t2xvert)*(t1xvert-t2xvert) + (t1yvert-t2yvert)*(t1yvert-t2yvert) + (t1zvert-t2zvert)*(t1zvert-t2zvert));
   
//	    cost1t2 = (dirt1[0]*dirt2[0]+dirt1[1]*dirt2[1]+dirt1[2]*dirt2[2]);
 
	    fdistHist->Fill(rDist);

	  }  // end loop   for (Int_t i = 0; i < maxntpar; i++)

	  meanDist = fdistHist->GetMean();	
	  rmsDist = fdistHist->GetRMS();	

	}   // end loop (vecpar[0] >= 0 && maxntpar > 0)	  
 
      // =================  End calculus meanDist ====================        	    


  // ==== Count the Number of Matched Tracks  ===
		
	ntmatch = 0;
	vecmu.push_back(ntmaxpar);
	// Main track ntmaxpar : vecmatch[ntmaxpar] is the track that matches with ntmaxpar
	if (vecmatch[ntmaxpar] > 0){
	  ntmatch++;
	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){ 
	    cout << endl; 
	    cout << "**** $$$$  N.match : " << ntmatch << " Main Track   " << ntmaxpar << "  matches with  " << vecmatch[ntmaxpar] << endl; 
	  }
	}
//                 --------------------------    13/jul/2017    -------------------
//      Check vecmatch that belongs to parallel tracks  13/jul/2017
//      If vecmatch[imat] does not belong to any parallel track (vecpar) set vecmtach[imat]=-10
//      If vecmatch[imat] belongs to a track that has to be deleted (tagtrk) set vecmtach[imat]=-10
      for(Int_t imat=0; imat<nTracks; imat++){
        if(vecmatch[imat]>=0){
          Bool_t found=kFALSE; 
          for(Int_t jmat=0; jmat<maxntpar; jmat++){
            if(vecmatch[imat]==vecpar[jmat])found=kTRUE;
          }
          if(found){
            for(Int_t idel=0; idel<itag ; idel++){
              if(vecmatch[imat]==tagtrk[idel]) vecmatch[imat]=-10;
            }  
          } 
          else {
            vecmatch[imat]=-10;
          }
        }
      }    	
//           ---------------- end    13/jul/2017       ----------------------------

      // === LOOP ON PARALLEL TRACKS (jpar) ; main track is *itrack ===	
      
	for(Int_t im = 0; im < maxntpar; im++){
	  // jpar = vecpar[im] contains the list of all the tracks parallel to ntmaxpar
	  Int_t jpar = vecpar[im];
      //   * Put in the vector vecmu all the good parallel tracks	  		
	  Int_t jdel=0; 
	  for(Int_t idel=0; idel<itag ; idel++){
	    if(jpar==tagtrk[idel]) jdel=1;
	  }
	  if(jdel==0){  // good parallel track
	    vecmu.push_back(vecpar[im]); 
	   //  * Count the matched tracks and write the matching
	  //  vecmatch[jpar] is the track that matches with jpar
	    if(vecmatch[jpar] > 0 && vecmatch[jpar] > jpar){
	      ntmatch++;  	      	        
	      if(nEvent>=nevdeb1&&nEvent<=nevdeb2)
		{  
	        cout << " N.match : " << ntmatch << " Track       " << jpar << " matches with " << vecmatch[jpar] << endl; 
	      }
	    }  // END if(vecmatch[jpar] > 0 && vecmatch[jpar] > jpar)	    	    
	  } // END if(jdel==0)			      	      	  		
	}  // END loop   for(Int_t im = 0; im < maxntpar; im++)
	
	if(nEvent>=nevdeb1&&nEvent<=nevdeb2)
	{ 
	  cout << " Number of matching :  " << ntmatch << endl;
	  cout << endl;
	  cout << "  List of all the parallel tracks and matching tracks " << endl;
	  for(UInt_t l = 0; l < vecmu.size(); l++){
            Int_t jvmu=vecmu[l]; 
	    cout << l << "   vecmu[l] " << vecmu[l] << "    matched track= " << vecmatch[jvmu] << endl ;
	  }
	  cout << endl;
	  cout << "  List of all the deleted tracks " << endl;	  
	  for(Int_t idel=0; idel<itag ; idel++){
	    cout << idel << "   tagtrk[idel] " << tagtrk[idel]  << endl;	     
	  }
	}      
	      
          nMuons = vecmu.size() - ntmatch ;

        //!if(nEvent>=nevdeb1&&nEvent<=nevdeb2)
	{  
            cout << endl;
	    cout << " RECONSTRUCTED NUMBER OF MUONS : nMuons = vecmu.size() - ntmatch  =  " << vecmu.size() << " - " << ntmatch << " = " << nMuons << endl;
	    cout << endl; 
        // Write in a special file only the numerical values for study with event display	
           fprintf(filehmmeout," ============================================ \n");
           fprintf(filehmmeout,"N.run|N.ev|Timest|N.Aco|N.tr|N.par.tr|Fracmu|MDist|N.matchmu|N.mu|File \n");     
	   fprintf(filehmmeout,"%6d %6d %u %2d %4d %3d %4.2f %5.2f %3d %3d %s \n",nRun,nEvent,timestp,mcnACO,nTracks,maxntpar,fracnMuTrk,meanDist,ntmatch,nMuons,curfile->GetName()); 	    	    
        }

        
 // -------- Tag the interaction events (mu inter. in iron) flagintev=1 ------- 
	Float_t nTrk = nTracks;
	fracnMuTrk = -10.;
	flagintev = 0; 
 
	if(nMuons>0) fracnMuTrk = nTrk/nMuons;        
		
 	if (nMuons > 3 && meanDist > 50 && meanDist < 120 && fracnMuTrk > 2.2){ 
	  nInter++;
          flagintev=1;
	  fmeanDistMu3->Fill(meanDist,fracnMuTrk);
	}else if (nMuons > 3 && meanDist > 0 && meanDist < 50){
	  nInter++;
	  flagintev=1;
	  fmeanDistMu3->Fill(meanDist,fracnMuTrk);	
	}
//       --------- End interaction events -------------


	// Counts the multi-muon events (multimuev) they have at least 4 muons and are not flagged as interaction events
	if (nMuons > 3 && flagintev==0){ 
	  multimuev++;
	  fmeanDistMu4->Fill(meanDist,fracnMuTrk);
	}
	
      }//END if ntmaxpar >= 0

     //      -------------------- END OF MATCHING ALGORITHM --------------------

         
      // =========   START TO WRITE THE VARIABLE IN THE TREE FOR EACH MUON FOUND =========	
 

	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){
	    if (nMuons == 0) {
	    cout << " No Muons Found " << "  " << "Event = " << nEvent << endl;
	    } 
	    if (flagintev==1) {   
	    cout << " @@@ INTERACTION EVENT @@@  " << "    meanDist = " << meanDist << " frac = " << fracnMuTrk << " nTracks = " << nTracks << "   Num. of muons reconst. = "   << nMuons << endl;  
	    }
       	  }
 
/*     Insert after the TREE for Interaction Events      	 
//   -------- Interaction Events :  TREE not written, only some variables written in a file -------- 
       if(flagintev==1){
         fprintf(intereve," ============================================ \n");
         fprintf(intereve,"%3d %6d %6d %u %4d %3d %4.2f %5.2f %s\n",nInter,nRun,nEvent,timestp,nTracks,nMuons,fracnMuTrk,meanDist,curfile->GetName());
          //    Fill the TREE for Interaction Events 
          fDatree->Fill();   // Fill the TREE  
	  neverFill = kFALSE;
       } // END   if(flagintev==1)        
//   ---------------END  TREE for Interaction Events --------------
*/        


//      if (nMuons > 0 && flagintev==0){
      if(nMuons > 0){

//  Fill number of Muons in the event with at least one muon. This count the
//  number of trigger with at least one muon in TPC (entries of the plot)
//  excluding the Interaction Events
 
	if (nameTrigger.Contains("OB0")) fnTrksTOFOB0Muons->Fill(nMuons);		
	if (nameTrigger.Contains("OB1")) fnTrksTOFOB1Muons->Fill(nMuons);
	if (nameTrigger.Contains("OB3")) fnTrksTOFOB3Muons->Fill(nMuons);
	if (nameTrigger.Contains("SCO")) fnTrksSCOMuons->Fill(nMuons);
	if (nameTrigger.Contains("AMU")) fnTrksAMUMuons->Fill(nMuons);
	if (nameTrigger.Contains("ASL")) fnTrksASLMuons->Fill(nMuons);
  
	//=====  Getting information of each muon and FILL THE TREE  =======
    
	if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  
	  cout << " ====  Getting information of each muon  =====" << endl;
	}

	TVectorD pxpypz(3);
	TVectorD xyzv(3);
	TVectorD pxpypz1(3);
	TVectorD xyzv1(3);
        TVectorD xyzout(3);
        TVectorD xyzin(3);

	Int_t rmu = 0 ;    // Count the number of muons for the index vector
	Int_t rmureal = 0 ;  // Count real number of muons
        Int_t rmumatch = 0;  // Count num. muon matched
        Int_t rmusingtrk = 0; // Count num. muon single track
        Int_t rmuplus  = 0; // Count num. mu+
        Int_t rmuminus = 0; // Count num mu-
        Int_t rmuplusminus = 0; // Count the number of muon in which the charged is ambigous
        Int_t srmuplus  = 0; // Count num. mu+ single track
        Int_t srmuminus = 0; // Count num mu- single track
	

	for (UInt_t j = 0; j < vecmu.size(); j++){ // Loop on tracks in vecmu (jtr1)

	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  //use this to check
            cout << "Loop on vecmu: " << j << "   vecmu.size()= " << vecmu.size() << endl;
          } 
	
	  Int_t jtr1 = vecmu[j];	
	  Int_t jmatch = vecmatch[jtr1];

	if(nEvent>=nevdeb1&&nEvent<=nevdeb2)
	{  //use this to check
            cout << "Track jtr1= " << vecmu[j] << "  Track matched jmatch= " << vecmatch[jtr1] << endl;
          }
	  
	  if(jtr1>jmatch&&jmatch>=0) continue;  // already matched jmatch-->jtr1

//        Sign who is track up and track down : Initialization
          Int_t uptrk=-10;
          Int_t dwtrk=-10;

	  AliESDtrack *tr1 = fESD->GetTrack(jtr1);

	  if(jmatch>=0){ // Muon with 2 matched tracks (jtr1 and jmatch) pointer tr1 and tr2

	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  //use this to check
            cout << "Matched:   Track jtr1= " << vecmu[j] << "  Track matched jmatch= " << vecmatch[jtr1] << endl;
          } 

	    //         MATCHED TRACKS : jtr1 is a track with a matched track = jmatch   	     

	    AliESDtrack *tr2 = fESD->GetTrack(jmatch);

	    rmureal++;       // Number of muons
	    rmu=rmureal-1;       // Variable for the loop
	    rmumatch++;  // Number of muons matched

	    trkmatch.push_back(jmatch);
		 	
	    if (tr1->GetOuterParam()->GetAlpha() > 0)  // up track tr1
	      {	   
	        uptrk=jtr1;      // track jtr1 is up
		AliExternalTrackParam *uptr1 = (AliExternalTrackParam *)tr1;		
		uptr1->PxPyPz(pxpypz.GetMatrixArray());
		uptr1->XvYvZv(xyzv.GetMatrixArray());
		
		tr1->GetOuterXYZ(xyzout.GetMatrixArray());
                tr1->GetInnerXYZ(xyzin.GetMatrixArray());
  	        tr1->GetDirection(dirt1); 

		uptracmu[rmu] = jtr1;   // track number 
		upxv[rmu] = xyzv[0];
		upyv[rmu] = xyzv[1];
		upzv[rmu] = xyzv[2];
		uptheta[rmu] = RadGrad*uptr1->Theta();
		upphi[rmu] = RadGrad*uptr1->Phi();
		uppx[rmu] = pxpypz[0];  
		uppy[rmu] = pxpypz[1];
		uppz[rmu] = pxpypz[2];  
		Pup[rmu] = uptr1->P();
//%%	        ToT_TOFup[rmu] = tr1->GetTOFsignalToT();		  
		corrTimeTOFup[rmu] = tr1->GetTOFsignal();
//%%	        DzTOFup[rmu] = tr1->GetTOFsignalDz();		     
//%%	        DxTOFup[rmu] = tr1->GetTOFsignalDx();		     
//MS-28/9/18
		nMatchup[rmu] = tr1->GetTOFclusterN();
		indexTOFup[rmu] = tr1->GetTOFCalChannel();
//MS-28/9/18
	        Xoutup[rmu] = xyzout[0];	
	        Youtup[rmu] = xyzout[1];	
	        Zoutup[rmu] = xyzout[2];	
	        Xinup[rmu] = xyzin[0];	
	        Yinup[rmu] = xyzin[1];	
	        Zinup[rmu] = xyzin[2];	
                upxdir[rmu] = dirt1[0];
	        upydir[rmu] = dirt1[1];
	        upzdir[rmu] = dirt1[2];

//    Calculate the teta and phi of the up track with directions 
              uptetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(upydir[rmu]));
              if(upydir[rmu]>=0){   
              upphidir[rmu] = RadGrad*TMath::ACos(upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(upxdir[rmu]<0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(upxdir[rmu]>0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }
              if(upydir[rmu]<0){   
              upphidir[rmu] = RadGrad*TMath::ACos(-upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(-upxdir[rmu]<0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(-upxdir[rmu]>0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }


		if (Pup[rmu]>0)  upthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(uppy[rmu])/Pup[rmu]);
		if (uppz[rmu]>0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu]);
		if (uppz[rmu]<0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu])+180;
		
		//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(uppx[rmu]*uppx[rmu]+uppz[rmu]*uppz[rmu]);

//%%
/*  
              if(uppz[rmu]>0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= RadGrad*TMath::ACos(uppx[rmu]/pxz);      

              if(uppz[rmu]>0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(uppz[rmu]/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(uppx[rmu])/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(uppz[rmu])/pxz);  
*/
//%%
            
		upSign[rmu] = uptr1->GetSign();
//%%2018		uptpcSignal[rmu] = tr1->GetTPCsignal();
		uptpcNcls[rmu] = tr1->GetTPCNcls();
		uptpcChi2[rmu] = tr1->GetTPCchi2();
//		upratioChi2cls[rmu] = uptpcChi2[rmu]/uptpcNcls[rmu];		
		Sig1Ptup[rmu] = uptr1->GetSigma1Pt2();
		// Variable used but not to record in the tree
		Ptupinv[rmu] = uptr1->OneOverPt();  
		Pup2[rmu] = uppx[rmu]*uppx[rmu]+uppy[rmu]*uppy[rmu]+uppz[rmu]*uppz[rmu];
		SigUp2[rmu] = (uptr1->GetSigmaY2()*uppy[rmu]*uppy[rmu]/Pup2[rmu])+(uptr1->GetSigmaZ2()*uppz[rmu]*uppz[rmu]/Pup2[rmu]);  

	      }  // END up track tr1

	    if (tr1->GetOuterParam()->GetAlpha() < 0)  // down track tr1
	      {
	        dwtrk=jtr1;      // track jtr1 is down
		AliExternalTrackParam *dwtr1 = (AliExternalTrackParam *)tr1;
		dwtr1->PxPyPz(pxpypz1.GetMatrixArray());
		dwtr1->XvYvZv(xyzv1.GetMatrixArray());

		tr1->GetOuterXYZ(xyzout.GetMatrixArray());
                tr1->GetInnerXYZ(xyzin.GetMatrixArray());	
  	        tr1->GetDirection(dirt1); 

		dwtracmu[rmu] = jtr1;   // track number 
		dwxv[rmu] = xyzv1[0];
		dwyv[rmu] = xyzv1[1];
		dwzv[rmu] = xyzv1[2];
		dwtheta[rmu] = RadGrad*dwtr1->Theta();
		dwphi[rmu] = RadGrad*dwtr1->Phi();
		dwpx[rmu] = pxpypz1[0];  
		dwpy[rmu] = pxpypz1[1];
		dwpz[rmu] = pxpypz1[2];  
		Pdw[rmu] = dwtr1->P();
//%%	        ToT_TOFdw[rmu] = tr1->GetTOFsignalToT();				     
	        corrTimeTOFdw[rmu] = tr1->GetTOFsignal();
//%%		DzTOFdw[rmu] = tr1->GetTOFsignalDz();		     
//%%	        DxTOFdw[rmu] = tr1->GetTOFsignalDx();		     
//MS-28/9/18
		nMatchdw[rmu] = tr1->GetTOFclusterN();
		indexTOFdw[rmu] = tr1->GetTOFCalChannel();
//MS-28/9/18
	        Xoutdw[rmu] = xyzout[0];	
	        Youtdw[rmu] = xyzout[1];	
	        Zoutdw[rmu] = xyzout[2];	
	        Xindw[rmu] = xyzin[0];	
	        Yindw[rmu] = xyzin[1];	
	        Zindw[rmu] = xyzin[2];	
                dwxdir[rmu] = dirt1[0];
	        dwydir[rmu] = dirt1[1];
	        dwzdir[rmu] = dirt1[2];

//    Calculate the teta and phi of the down track with directions for down change the sign
//       of the vectors of the direction in the formulas (dwxdir,dwydir,dwzdir)
//     Sometimes the TPC reconstruction gives the wrong sign to dwydir ??? so we check
//     it both for dw and up track.
              dwtetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwydir[rmu]));
              if(dwydir[rmu]<=0){
              dwphidir[rmu] = RadGrad*TMath::ACos(-dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(-dwxdir[rmu]<0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(-dwxdir[rmu]>0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }
         //    if wrong sign given by TPC
              if(dwydir[rmu]>0){
              dwphidir[rmu] = RadGrad*TMath::ACos(dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(dwxdir[rmu]<0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(dwxdir[rmu]>0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }



		if (Pdw[rmu]>0) dwthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwpy[rmu])/Pdw[rmu]);
		if (dwpz[rmu]<0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu]);
		if (dwpz[rmu]>0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu])+180;

//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(dwpx[rmu]*dwpx[rmu]+dwpz[rmu]*dwpz[rmu]);

//%%
/* 
              if(dwpz[rmu]<0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]<0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);      
*/
//%%           

		dwSign[rmu] = dwtr1->GetSign();
//%%2018		dwtpcSignal[rmu] = tr1->GetTPCsignal();
		dwtpcNcls[rmu] = tr1->GetTPCNcls();
		dwtpcChi2[rmu] = tr1->GetTPCchi2();
		//		  dwratioChi2cls[rmu] = dwtpcChi2[rmu]/dwtpcNcls[rmu];
		Sig1Ptdw[rmu] = tr1->GetSigma1Pt2();
		// Variable used but not to record in the tree
		Ptdwinv[rmu] = dwtr1->OneOverPt();  
		Pdw2[rmu] = dwpx[rmu]*dwpx[rmu]+dwpy[rmu]*dwpy[rmu]+dwpz[rmu]*dwpz[rmu];
		SigDw2[rmu] = (dwtr1->GetSigmaY2()*dwpy[rmu]*dwpy[rmu]/Pdw2[rmu])+(dwtr1->GetSigmaZ2()*dwpz[rmu]*dwpz[rmu]/Pdw2[rmu]);
 
	      }  // END down track tr1
	    
	    if (tr2->GetOuterParam()->GetAlpha() > 0)  // up track tr2
	      { 
	        uptrk=jmatch;      // track jmatch is up
 		AliExternalTrackParam *uptr2 = (AliExternalTrackParam *)tr2; 			
		uptr2->PxPyPz(pxpypz.GetMatrixArray());
		uptr2->XvYvZv(xyzv.GetMatrixArray());

		tr2->GetOuterXYZ(xyzout.GetMatrixArray());
                tr2->GetInnerXYZ(xyzin.GetMatrixArray());
  	        tr2->GetDirection(dirt2); 

		uptracmu[rmu] = jmatch;   // track number 
		upxv[rmu] = xyzv[0];
		upyv[rmu] = xyzv[1];
		upzv[rmu] = xyzv[2];
		uptheta[rmu] = RadGrad*uptr2->Theta();
		upphi[rmu] = RadGrad*uptr2->Phi();
		uppx[rmu] = pxpypz[0];  
		uppy[rmu] = pxpypz[1];
		uppz[rmu] = pxpypz[2];  
		Pup[rmu] = uptr2->P();
//%%	        ToT_TOFup[rmu] = tr2->GetTOFsignalToT();			               
	        corrTimeTOFup[rmu] = tr2->GetTOFsignal();
//%%		DzTOFup[rmu] = tr2->GetTOFsignalDz();		     
//%%	        DxTOFup[rmu] = tr2->GetTOFsignalDx();		     
//MS-28/9/18
		nMatchup[rmu] = tr2->GetTOFclusterN();
		indexTOFup[rmu] = tr2->GetTOFCalChannel();
//MS-28/9/18
	        Xoutup[rmu] = xyzout[0];	
	        Youtup[rmu] = xyzout[1];	
	        Zoutup[rmu] = xyzout[2];	
	        Xinup[rmu] = xyzin[0];	
	        Yinup[rmu] = xyzin[1];	
	        Zinup[rmu] = xyzin[2];	
                upxdir[rmu] = dirt2[0];
	        upydir[rmu] = dirt2[1];
	        upzdir[rmu] = dirt2[2];

//    Calculate the teta and phi of the up track with directions 
              uptetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(upydir[rmu]));
              if(upydir[rmu]>=0){   
              upphidir[rmu] = RadGrad*TMath::ACos(upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(upxdir[rmu]<0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(upxdir[rmu]>0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }
              if(upydir[rmu]<0){   
              upphidir[rmu] = RadGrad*TMath::ACos(-upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(-upxdir[rmu]<0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(-upxdir[rmu]>0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }



		if (Pup[rmu]>0)  upthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(uppy[rmu])/Pup[rmu]);
		if (uppz[rmu]>0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu]);
		if (uppz[rmu]<0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu])+180;


//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(uppx[rmu]*uppx[rmu]+uppz[rmu]*uppz[rmu]);

//%%
/* 
              if(uppz[rmu]>0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= RadGrad*TMath::ACos(uppx[rmu]/pxz);      

              if(uppz[rmu]>0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(uppz[rmu]/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(uppx[rmu])/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(uppz[rmu])/pxz);  
*/
//%%             

		upSign[rmu] = uptr2->GetSign();
//%%2018		uptpcSignal[rmu] = tr2->GetTPCsignal();
		uptpcNcls[rmu] = tr2->GetTPCNcls();
		uptpcChi2[rmu] = tr2->GetTPCchi2();
		//		  upratioChi2cls[rmu] = uptpcChi2[rmu]/uptpcNcls[rmu];
		  
		Sig1Ptup[rmu] = tr2->GetSigma1Pt2();
		// Variable used but not to record in the tree
		Ptupinv[rmu] = uptr2->OneOverPt();  
		Pup2[rmu] = uppx[rmu]*uppx[rmu]+uppy[rmu]*uppy[rmu]+uppz[rmu]*uppz[rmu];
		SigUp2[rmu] = (uptr2->GetSigmaY2()*uppy[rmu]*uppy[rmu]/Pup2[rmu])+(uptr2->GetSigmaZ2()*uppz[rmu]*uppz[rmu]/Pup2[rmu]);  

	      }  // END up track tr2

	    if (tr2->GetOuterParam()->GetAlpha() < 0)  // down track tr2
	      {
	        dwtrk=jmatch;      // track jmatch is down
		AliExternalTrackParam *dwtr2 = (AliExternalTrackParam *)tr2;		
		dwtr2->PxPyPz(pxpypz1.GetMatrixArray());
		dwtr2->XvYvZv(xyzv1.GetMatrixArray());

		tr2->GetOuterXYZ(xyzout.GetMatrixArray());
                tr2->GetInnerXYZ(xyzin.GetMatrixArray());
  	        tr2->GetDirection(dirt2); 
	
		dwtracmu[rmu] = jmatch;   // track number 
		dwxv[rmu] = xyzv1[0];
		dwyv[rmu] = xyzv1[1];
		dwzv[rmu] = xyzv1[2];
		dwtheta[rmu] = RadGrad*dwtr2->Theta();
		dwphi[rmu] = RadGrad*dwtr2->Phi();
		dwpx[rmu] = pxpypz1[0];  
		dwpy[rmu] = pxpypz1[1];
		dwpz[rmu] = pxpypz1[2];  
		dwxv[rmu] = xyzv1[0];
		dwyv[rmu] = xyzv1[1];
		dwzv[rmu] = xyzv1[2];
		Pdw[rmu] = dwtr2->P();
//%%	        ToT_TOFdw[rmu] = tr2->GetTOFsignalToT();		  
	        corrTimeTOFdw[rmu] = tr2->GetTOFsignal();
//%%		DzTOFdw[rmu] = tr2->GetTOFsignalDz();		     
//%%	        DxTOFdw[rmu] = tr2->GetTOFsignalDx();		     
//MS-28/9/18
		nMatchdw[rmu] = tr2->GetTOFclusterN();
		indexTOFdw[rmu] = tr2->GetTOFCalChannel();
//MS-28/9/18
	        Xoutdw[rmu] = xyzout[0];	
	        Youtdw[rmu] = xyzout[1];	
	        Zoutdw[rmu] = xyzout[2];	
	        Xindw[rmu] = xyzin[0];	
	        Yindw[rmu] = xyzin[1];	
	        Zindw[rmu] = xyzin[2];	
                dwxdir[rmu] = dirt2[0];
	        dwydir[rmu] = dirt2[1];
	        dwzdir[rmu] = dirt2[2];

//    Calculate the teta and phi of the down track with directions for down change the sign
//       of the vectors of the direction in the formulas (dwxdir,dwydir,dwzdir)
//     Sometimes the TPC reconstruction gives the wrong sign to dwydir ??? so we check
//     it both for dw and up track.
              dwtetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwydir[rmu]));
              if(dwydir[rmu]<=0){
              dwphidir[rmu] = RadGrad*TMath::ACos(-dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(-dwxdir[rmu]<0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(-dwxdir[rmu]>0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }
         //    if wrong sign given by TPC
              if(dwydir[rmu]>0){
              dwphidir[rmu] = RadGrad*TMath::ACos(dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(dwxdir[rmu]<0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(dwxdir[rmu]>0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }



		if (Pdw[rmu]>0) dwthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwpy[rmu])/Pdw[rmu]);
		if (dwpz[rmu]<0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu]);
		if (dwpz[rmu]>0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu])+180;

//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(dwpx[rmu]*dwpx[rmu]+dwpz[rmu]*dwpz[rmu]);

//%%
/*  
              if(dwpz[rmu]<0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]<0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);    
*/
//%%

		dwSign[rmu] = dwtr2->GetSign();
//%%2018		dwtpcSignal[rmu] = tr2->GetTPCsignal();
		dwtpcNcls[rmu] = tr2->GetTPCNcls();
		dwtpcChi2[rmu] = tr2->GetTPCchi2();
		//		  dwratioChi2cls[rmu] =  dwtpcChi2[rmu]/dwtpcNcls[rmu];
		  
		Sig1Ptdw[rmu] = dwtr2->GetSigma1Pt2();
		// Variable used but not to record in the tree
		Ptdwinv[rmu] = dwtr2->OneOverPt();  
		Pdw2[rmu] = dwpx[rmu]*dwpx[rmu]+dwpy[rmu]*dwpy[rmu]+dwpz[rmu]*dwpz[rmu];
		SigDw2[rmu] = (dwtr2->GetSigmaY2()*dwpy[rmu]*dwpy[rmu]/Pdw2[rmu])+(dwtr2->GetSigmaZ2()*dwpz[rmu]*dwpz[rmu]/Pdw2[rmu]);  

	      }  // END down track tr2

             // Count the number of mu+ and mu- or bad charge (ambigous)
 	    if(upSign[rmu]==-1&&dwSign[rmu]==+1) rmuplus++;            
 	    if(upSign[rmu]==+1&&dwSign[rmu]==-1) rmuminus++; 
 	    if(upSign[rmu]==-1&&dwSign[rmu]==-1) rmuplusminus++; 
 	    if(upSign[rmu]==+1&&dwSign[rmu]==+1) rmuplusminus++;  	    

	    // Calculus Pcov for matching tracks

	    // Ask to have one track UP and the other DOWN both found

	    if(upSign[rmu]!=-10.&&dwSign[rmu]!=-10.) { 	      
	      //calculus Pcov   
          // substitute  tr1 with uptrk     tr2 with dwtrk
          AliESDtrack *uptrkp = fESD->GetTrack(uptrk);
          AliESDtrack *dwtrkp = fESD->GetTrack(dwtrk);             
	      AliExternalTrackParam *utrack= new AliExternalTrackParam(*uptrkp);
	      AliExternalTrackParam *par1R= new AliExternalTrackParam(*dwtrkp);
	      
	      // START Ruben instructions 
	      par1R->Invert();
	      if(par1R->Propagate(uptrkp->GetAlpha(),uptrkp->GetX(),AliTracker::GetBz()))      {	            	      
            AliTrackerBase::UpdateTrack(*utrack,*par1R);
	        Double_t *param=(Double_t*) utrack->GetParameter();
	        Pcov[rmu] =TMath::Sqrt(1.+ param[3]*param[3])/TMath::Abs(param[4]);	     
//%%2018 	        PRes[rmu] = (Ptupinv[rmu]-Ptdwinv[rmu])/(0.5*(Ptupinv[rmu]+Ptdwinv[rmu]));	      
	//!      cout << "charge covariant " << utrack->Charge() << endl;
	        Chargecov[rmu] = -utrack->Charge();
 	      } // END if(par1R->Propagate(tr1->GetAlpha(),tr1->GetX(),AliTracker::GetBz()))
          //  END calculus Pcov and Charge covariant  
          //  END Ruben instructions              

              
         //  Calculus distance between track up and down (Katherin)
          if( uptracmu[rmu] >= 0 && dwtracmu[rmu] >= 0 ){
	        AliESDtrack *upMuTrack = fESD->GetTrack(uptracmu[rmu]);
            AliESDtrack *dwMuTrack = fESD->GetTrack(dwtracmu[rmu]);
		
  		    Float_t tXup = upMuTrack->Xv(); 
		    Float_t tYup = upMuTrack->Yv(); 
		    Float_t tZup = upMuTrack->Zv();
		
		    Float_t tXdw = dwMuTrack->Xv(); 
		    Float_t tYdw = dwMuTrack->Yv(); 
		    Float_t tZdw = dwMuTrack->Zv();
		
		    MuDistupdw[rmu] = sqrt((tXup-tXdw)*(tXup-tXdw) + (tYup-tYdw)*(tYup-tYdw) + (tZup-tZdw)*(tZup-tZdw));

	      }    // END if( uptracmu[rmu] >= 0 && dwtracmu[rmu] >= 0 )

/*   //%%2018               
             
	      if (SigUp2[rmu]>0 && SigDw2[rmu]>0){ 
		Xinv[rmu] = (Ptupinv[rmu]-Ptdwinv[rmu])/sqrt(Sig1Ptup[rmu]+Sig1Ptdw[rmu]);         
		Pmu[rmu] = (SigUp2[rmu]/(SigUp2[rmu]+SigDw2[rmu]))*Pup[rmu] + (SigDw2[rmu]/(SigUp2[rmu]+SigDw2[rmu]))*Pdw[rmu];
	      } else {
		badmatch++;
	      }
*/

	      if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  //use this to check
//             Study directions of muons for event with 1 muon (check).
//              if(nMuons==1) {    
	        Debtrackused->Fill(jtr1); 
	        Debtrackused->Fill(vecmatch[jtr1]);
	        Debtrackusedmatch->Fill(jtr1); 
	        Debtrackusedmatch->Fill(vecmatch[jtr1]);
	        freqtrack[jtr1]++;
	        freqtrack[(vecmatch[jtr1])]++;
	        if(uppx[rmu]!=-999&&uppy[rmu]!=-999&&uppz[rmu]!=-999){
	        ptotup = TMath::Sqrt(uppx[rmu]*uppx[rmu]+uppz[rmu]*uppz[rmu]+uppy[rmu]*uppy[rmu]);   
	        } else {
	        ptotup = -999;
	        }
	        if(dwpx[rmu]!=-999&&dwpy[rmu]!=-999&&dwpz[rmu]!=-999){	        
	        ptotdw = TMath::Sqrt(dwpx[rmu]*dwpx[rmu]+dwpz[rmu]*dwpz[rmu]+dwpy[rmu]*dwpy[rmu]);   
	        } else {
	        ptotdw = -999;
	        }
	        	        	        
	        cout << endl;     	    	
	        cout << " ===== Matched Muons ===== Total Muons : " << rmureal  <<endl;	
               	cout << " uptrack : " << uptrk << " dwtrack : " << dwtrk << endl;		        
		cout << " N. Matched Mu: " << rmumatch <<  " upxv: " << upxv[rmu] <<  " upyv: " << upyv[rmu] << " upzv: " << upzv[rmu] << " dwxv: "  << dwxv[rmu] << " dwyv: " << dwyv[rmu] << " dwzv: " << dwzv[rmu] <<  " uppx: " << uppx[rmu] <<  " uppy: " << uppy[rmu] << " uppz: " << uppz[rmu] << " Pup: " << Pup[rmu] << " dwpx: "  << dwpx[rmu] << " dwpy: " << dwpy[rmu] << " dwpz: " << dwpz[rmu] << " Pdw:" << Pdw[rmu] << " upthetaCosmic "  << upthetaCosmic[rmu] << " upphiCosmic "  << upphiCosmic[rmu]  << "  dwthetaCosmic "  << dwthetaCosmic[rmu] << " dwphiCosmic "  << dwphiCosmic[rmu]  << " Ptotup, Ptotdw : " << ptotup << " , " << ptotdw << "  Pcov:  "  << Pcov[rmu] <<  "  Up Charge: " << upSign[rmu] << "  Dw Charge: " << dwSign[rmu] << "  Covariant Charge: " << Chargecov[rmu] <<  "   Ev: " << nEvent << endl;
		cout << " Dist. up-dw : " << MuDistupdw[rmu] << endl;
		cout << " Index rmu : " << rmu << endl;
		cout << endl; 
                cout <<   " &&&&&&&&&& Study the direction of the muons &&&&&&&& " << endl;
                cout << " upthetacosmic, upphicosmic = " << upthetaCosmic[rmu] << "  " << upphiCosmic[rmu] << endl;
                cout << " dwthetacosmic, dwphicosmic = " << dwthetaCosmic[rmu] << "  " << dwphiCosmic[rmu] << endl;
                cout << " upthetadir, upphidir = " << uptetadir[rmu] << "  " << upphidir[rmu] << endl;
                cout << " dwthetadir, dwphidir = " << dwtetadir[rmu] << "  " << dwphidir[rmu] << endl;
                cout << " upxv, upyv, upzv = " << upxv[rmu] << "  " << upyv[rmu] << "  " << upzv[rmu] << endl;
                cout << " dwxv, dwyv, dwzv = " << dwxv[rmu] << "  " << dwyv[rmu] << "  " << dwzv[rmu] << endl;
                cout << " upxdir, upydir, upzdir = " << upxdir[rmu] << "  " << upydir[rmu] << "  " << upzdir[rmu] << endl;
                cout << " dwxdir, dwydir, dwzdir = " << dwxdir[rmu] << "  " << dwydir[rmu] << "  " << dwzdir[rmu] << endl;
                cout << " xoutup, youtup, zoutup = " << Xoutup[rmu] << "  " << Youtup[rmu] << "  " << Zoutup[rmu] << endl;
                cout << " xinup, yinup, zinup = " << Xinup[rmu] << "  " << Yinup[rmu] << "  " << Zinup[rmu] << endl;
                cout << " xoutdw, youtdw, zoutdw = " << Xoutdw[rmu] << "  " << Youtdw[rmu] << "  " << Zoutdw[rmu] << endl;
                cout << " xindw, yindw, zindw = " << Xindw[rmu] << "  " << Yindw[rmu] << "  " << Zindw[rmu] << endl;
                cout <<   " &&&&&&&&&& END Study the direction of the muons &&&&&&&& " << endl;
                cout << endl << endl;   


    // Write in a special file only the numerical values for further fast analysis        
    fprintf(filehmmeout,"%3d %3d %3d %7.3f %7.3f %7.3f \n",rmureal,uptrk,dwtrk,ptotup,ptotdw,Pcov[rmu]);
		//============================================================
	      
	      // Filling Kath histograms (TList) for Thesis  (matched muon tracks)

	      fthetaDist->Fill(dwthetaCosmic[rmu]);   //angular (theta) distribution
	      fphiDist->Fill(dwphiCosmic[rmu]);       //angular (phi) distribution
	      fmatchZX->Fill(upzv[rmu],upxv[rmu]);    //spacial (ZX) distribution
	      fupdwMuDist->Fill(MuDistupdw[rmu]);     //dist. between matched tracks Up and Down
	      
	      ncountK++;
	      
	      Xpos[rmu] = upxv[rmu];
	      Zpos[rmu] = upzv[rmu];
	      
	      cout << endl;
	      cout << " Zpos, Xpos : " << Zpos[rmu] << " , " << Xpos[rmu] << endl;
	      cout << endl;
	      
	      } // END if(nEvent>=nevdeb1&&nEvent<=nevdeb2)
	      
	    }  // End   if(upSign[rmu]!=-10.&&dwSign[rmu]!=-10.)	     

	    
	  }else if( jmatch < 0 ){    // END if( jmatch >= 0 )
	
 // Check if the track jtr1 (the track of the loop) was taken before as matched track
	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  //use this to check
            cout << "Not matched:   Track jtr1= " << vecmu[j] << "  Track matched jmatch= " << vecmatch[jtr1] << endl;
          }

	    Int_t jtr1ana=0;  // Flag if track has been already used

	    for(Int_t jj=0; jj<trkmatch.size() ; jj++)  //Loop to check jtr1
	      {

		if(jtr1 == trkmatch[jj]) {
		  jtr1ana = 1 ;
		 break;  // if the track has been used before as matched track
		}
        
	      } // END Loop to check jtr1

	    if(jtr1ana == 1) continue ; // track already used as a matched track


	  if(nEvent>=nevdeb1&&nEvent<=nevdeb2){  //use this to check
            cout << "Not matched taken:  Track jtr1= " << vecmu[j] << "  Track matched jmatch= " << vecmatch[jtr1] << endl;
          }


     //      =============    MUON SINGLE TRACK  ====================
    // jtr1 is a MUON SINGLE TRACK  (not matched)  fill the variable for the TREE 

	    rmureal++;
	    rmu=rmureal-1;
	    rmusingtrk++;
	       
	    if (tr1->GetOuterParam()->GetAlpha() > 0)  // Single track is up 
	      {
		AliExternalTrackParam *uptr1 = (AliExternalTrackParam *)tr1;
		uptr1->PxPyPz(pxpypz.GetMatrixArray());
		uptr1->XvYvZv(xyzv.GetMatrixArray());
	
		tr1->GetOuterXYZ(xyzout.GetMatrixArray());
                tr1->GetInnerXYZ(xyzin.GetMatrixArray());
	        tr1->GetDirection(dirt1); 

		uptracmu[rmu] = jtr1;   // track number 
		upxv[rmu] = xyzv[0];
		upyv[rmu] = xyzv[1];
		upzv[rmu] = xyzv[2];
		uptheta[rmu] = RadGrad*uptr1->Theta();
		upphi[rmu] = RadGrad*uptr1->Phi();
		uppx[rmu] = pxpypz[0];  
		uppy[rmu] = pxpypz[1];
		uppz[rmu] = pxpypz[2];  
		Pup[rmu] = uptr1->P();
//%%2018		Pmu[rmu] = Pup[rmu];
//%%		PmuMed[rmu] = Pup[rmu];
//%%2018		Xinv[rmu] = -666.;  
//%%	        ToT_TOFup[rmu] = tr1->GetTOFsignalToT();
	        corrTimeTOFup[rmu] = tr1->GetTOFsignal();
//%%		DzTOFup[rmu] = tr1->GetTOFsignalDz();		     
//%%	        DxTOFup[rmu] = tr1->GetTOFsignalDx();		     
//MS-28/9/18
		nMatchup[rmu] = tr1->GetTOFclusterN();
		indexTOFup[rmu] = tr1->GetTOFCalChannel();
//MS-28/9/18
	        Xoutup[rmu] = xyzout[0];	
	        Youtup[rmu] = xyzout[1];	
	        Zoutup[rmu] = xyzout[2];	
	        Xinup[rmu] = xyzin[0];	
	        Yinup[rmu] = xyzin[1];	
	        Zinup[rmu] = xyzin[2];	
                upxdir[rmu] = dirt1[0];
	        upydir[rmu] = dirt1[1];
	        upzdir[rmu] = dirt1[2];

//    Calculate the teta and phi of the up track with directions 
              uptetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(upydir[rmu]));
              if(upydir[rmu]>=0){   
              upphidir[rmu] = RadGrad*TMath::ACos(upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(upxdir[rmu]<0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(upxdir[rmu]>0&&upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }
              if(upydir[rmu]<0){   
              upphidir[rmu] = RadGrad*TMath::ACos(-upxdir[rmu]/TMath::Sqrt(1-upydir[rmu]*upydir[rmu]));
              if(-upxdir[rmu]<0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu]; 
              if(-upxdir[rmu]>0&&-upzdir[rmu]<0)upphidir[rmu]=360-upphidir[rmu];
              }


		if (Pup[rmu]>0)  upthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(uppy[rmu])/Pup[rmu]);
		if (uppz[rmu]>0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu]);
		if (uppz[rmu]<0) upphiCosmic[rmu] = 90-RadGrad*TMath::ATan(uppx[rmu]/uppz[rmu])+180;

//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(uppx[rmu]*uppx[rmu]+uppz[rmu]*uppz[rmu]);

//%%
/*
              if(uppz[rmu]>0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= RadGrad*TMath::ACos(uppx[rmu]/pxz);      

              if(uppz[rmu]>0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(uppz[rmu]/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]<0)
              upphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(uppx[rmu])/pxz);      

              if(uppz[rmu]<0&&uppx[rmu]>0)
              upphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(uppz[rmu])/pxz);  
*/
//%%

		upSign[rmu] = uptr1->GetSign();
//%%2018		uptpcSignal[rmu] = tr1->GetTPCsignal();
		uptpcNcls[rmu] = tr1->GetTPCNcls();
		uptpcChi2[rmu] = tr1->GetTPCchi2();
		//		  upratioChi2cls[rmu] = uptpcChi2[rmu]/uptpcNcls[rmu];
		  
		Sig1Ptup[rmu] = tr1->GetSigma1Pt2();
		// Variable used but not to record in the tree
		Ptupinv[rmu] = uptr1->OneOverPt();  
		Pup2[rmu] = uppx[rmu]*uppx[rmu]+uppy[rmu]*uppy[rmu]+uppz[rmu]*uppz[rmu];
		SigUp2[rmu] = (uptr1->GetSigmaY2()*uppy[rmu]*uppy[rmu]/Pup2[rmu])+(uptr1->GetSigmaZ2()*uppz[rmu]*uppz[rmu]/Pup2[rmu]);  

	      }  // END Single track is up

	    if (tr1->GetOuterParam()->GetAlpha() < 0)  // Single track is down
	      { 
		AliExternalTrackParam *dwtr1 = (AliExternalTrackParam *)tr1;
		dwtr1->PxPyPz(pxpypz1.GetMatrixArray());
		dwtr1->XvYvZv(xyzv1.GetMatrixArray());
	
		tr1->GetOuterXYZ(xyzout.GetMatrixArray());
                tr1->GetInnerXYZ(xyzin.GetMatrixArray());
	        tr1->GetDirection(dirt1);
 	
		dwtracmu[rmu] = jtr1;   // track number 
		dwxv[rmu] = xyzv1[0];
		dwyv[rmu] = xyzv1[1];
		dwzv[rmu] = xyzv1[2];
		dwtheta[rmu] = RadGrad*dwtr1->Theta();
		dwphi[rmu] = RadGrad*dwtr1->Phi();
		dwpx[rmu] = pxpypz1[0];  
		dwpy[rmu] = pxpypz1[1];
		dwpz[rmu] = pxpypz1[2];  
		Pdw[rmu] = dwtr1->P();
//%%	        ToT_TOFdw[rmu] = tr1->GetTOFsignalToT();		 
	        corrTimeTOFdw[rmu] = tr1->GetTOFsignal();
//%%		DzTOFdw[rmu] = tr1->GetTOFsignalDz();		     
//%%	        DxTOFdw[rmu] = tr1->GetTOFsignalDx();		     
//MS-28/9/18
		nMatchdw[rmu] = tr1->GetTOFclusterN();
		indexTOFdw[rmu] = tr1->GetTOFCalChannel();
//MS-28/9/18
	        Xoutdw[rmu] = xyzout[0];	
	        Youtdw[rmu] = xyzout[1];	
	        Zoutdw[rmu] = xyzout[2];	
	        Xindw[rmu] = xyzin[0];	
	        Yindw[rmu] = xyzin[1];	
	        Zindw[rmu] = xyzin[2];	
                dwxdir[rmu] = dirt1[0];  // dirt1 instead dirt2 : Correction 13/jul/2017
	        dwydir[rmu] = dirt1[1];  // idem
	        dwzdir[rmu] = dirt1[2];  // idem

//    Calculate the teta and phi of the down track with directions for down change the sign
//       of the vectors of the direction in the formulas (dwxdir,dwydir,dwzdir)
//     Sometimes the TPC reconstruction gives the wrong sign to dwydir ??? so we check
//     it both for dw and up track.
              dwtetadir[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwydir[rmu]));
              if(dwydir[rmu]<=0){
              dwphidir[rmu] = RadGrad*TMath::ACos(-dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(-dwxdir[rmu]<0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(-dwxdir[rmu]>0&&-dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }
         //    if wrong sign given by TPC
              if(dwydir[rmu]>0){
              dwphidir[rmu] = RadGrad*TMath::ACos(dwxdir[rmu]/TMath::Sqrt(1-dwydir[rmu]*dwydir[rmu]));
              if(dwxdir[rmu]<0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu]; 
              if(dwxdir[rmu]>0&&dwzdir[rmu]<0)dwphidir[rmu]=360-dwphidir[rmu];
              }



		if (Pdw[rmu]>0) dwthetaCosmic[rmu] = RadGrad*TMath::ACos(TMath::Abs(dwpy[rmu])/Pdw[rmu]);
		if (dwpz[rmu]<0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu]);
		if (dwpz[rmu]>0)dwphiCosmic[rmu] = 90-RadGrad*TMath::ATan(dwpx[rmu]/dwpz[rmu])+180;


//  Correction of the phi cosmic angle
//%%              pxz = TMath::Sqrt(dwpx[rmu]*dwpx[rmu]+dwpz[rmu]*dwpz[rmu]);


//%%
/* 
              if(dwpz[rmu]<0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]<0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 90+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]>0)
              dwphiCosmiccor[rmu]= 180+RadGrad*TMath::ACos(TMath::Abs(dwpx[rmu])/pxz);      

              if(dwpz[rmu]>0&&dwpx[rmu]<0)
              dwphiCosmiccor[rmu]= 270+RadGrad*TMath::ACos(TMath::Abs(dwpz[rmu])/pxz);    
*/
//%%

//%%2018		Pmu[rmu] = Pdw[rmu];
//%%		PmuMed[rmu] = Pdw[rmu];
//%%2018		Xinv[rmu] = -666.;  
		dwSign[rmu] = tr1->GetSign();
//%%2018		dwtpcSignal[rmu] = tr1->GetTPCsignal();
		dwtpcNcls[rmu] = tr1->GetTPCNcls();
		dwtpcChi2[rmu] = tr1->GetTPCchi2();
		//		  dwratioChi2cls[rmu] = dwtpcChi2[rmu]/dwtpcNcls[rmu];
		  
		Sig1Ptdw[rmu]= dwtr1->GetSigma1Pt2();
		// Variable used but not to record in the tree
		Ptdwinv[rmu] = dwtr1->OneOverPt();  
		Pdw2[rmu]        = dwpx[rmu]*dwpx[rmu]+dwpy[rmu]*dwpy[rmu]+dwpz[rmu]*dwpz[rmu];
		SigDw2[rmu] = (dwtr1->GetSigmaY2()*dwpy[rmu]*dwpy[rmu]/Pdw2[rmu])+(dwtr1->GetSigmaZ2()*dwpz[rmu]*dwpz[rmu]/Pdw2[rmu]);  

	      }  // END Single track is down

             // Count the number of mu+ and mu- or bad charge (ambigous)
 	    if(upSign[rmu]==-1&&dwSign[rmu]==-10) srmuplus++;   
 	    if(upSign[rmu]==-10&&dwSign[rmu]==+1) srmuplus++;   	             
  	    if(upSign[rmu]==+1&&dwSign[rmu]==-10) srmuminus++;   
 	    if(upSign[rmu]==-10&&dwSign[rmu]==-1) srmuminus++;   	             

	      if(nEvent>=nevdeb1&&nEvent<=nevdeb2){   
//             Study directions of muons for event with 1 muon (check).
//              if(nMuons==1) {  

	        Debtrackused->Fill(jtr1); 
	        Debtrackusedsing->Fill(jtr1); 
	        freqtrack[jtr1]++;
	        if(uppx[rmu]!=-999&&uppy[rmu]!=-999&&uppz[rmu]!=-999){
	        ptotup = TMath::Sqrt(uppx[rmu]*uppx[rmu]+uppz[rmu]*uppz[rmu]+uppy[rmu]*uppy[rmu]);   
	        } else {
	        ptotup = -999;
	        }
	        if(dwpx[rmu]!=-999&&dwpy[rmu]!=-999&&dwpz[rmu]!=-999){	        
	        ptotdw = TMath::Sqrt(dwpx[rmu]*dwpx[rmu]+dwpz[rmu]*dwpz[rmu]+dwpy[rmu]*dwpy[rmu]);   
	        } else {
	        ptotdw = -999;
	        }	        
	        	        	        	        	      
	        cout << endl;     	    	
	        cout << " $$$$$ NOT Matched Muons $$$$$ Total muons : " << rmureal <<endl;	
               	cout << " jtr1 : " << jtr1 << " jmatch : " << vecmatch[jtr1] << endl;		       	   
		cout << " Num. Single track  Mu: " << rmusingtrk <<  " upxv: " << upxv[rmu] <<  " upyv: " << upyv[rmu] << " upzv: " << upzv[rmu] << " dwxv: "  << dwxv[rmu] << " dwyv: " << dwyv[rmu] << " dwzv: " << dwzv[rmu] <<  " uppx: " << uppx[rmu] <<  " uppy: " << uppy[rmu] << " uppz: " << uppz[rmu] << " dwpx: "  << dwpx[rmu] << " dwpy: " << dwpy[rmu] << " dwpz: " << dwpz[rmu] << " upthetaCosmic "  << upthetaCosmic[rmu] << " upphiCosmic "  << upphiCosmic[rmu] << " dwthetaCosmic "  << dwthetaCosmic[rmu] << " dwphiCosmic "  << dwphiCosmic[rmu] << "  Up Charge: " << upSign[rmu] << "  Dw Charge: " << dwSign[rmu] << " Ev: " << nEvent << endl;
		cout << endl;

		cout << endl; 
                cout <<   " &&&&&& SINGLE TRACK: Study the direction of the muons &&&&&&&& " << endl;
                cout << " upthetacosmic, upphicosmic = " << upthetaCosmic[rmu] << "  " << upphiCosmic[rmu] << endl;
                cout << " dwthetacosmic, dwphicosmic = " << dwthetaCosmic[rmu] << "  " << dwphiCosmic[rmu] << endl;
                cout << " upthetadir, upphidir = " << uptetadir[rmu] << "  " << upphidir[rmu] << endl;
                cout << " dwthetadir, dwphidir = " << dwtetadir[rmu] << "  " << dwphidir[rmu] << endl;
                cout << " upxv, upyv, upzv = " << upxv[rmu] << "  " << upyv[rmu] << "  " << upzv[rmu] << endl;
                cout << " dwxv, dwyv, dwzv = " << dwxv[rmu] << "  " << dwyv[rmu] << "  " << dwzv[rmu] << endl;
                cout << " upxdir, upydir, upzdir = " << upxdir[rmu] << "  " << upydir[rmu] << "  " << upzdir[rmu] << endl;
                cout << " dwxdir, dwydir, dwzdir = " << dwxdir[rmu] << "  " << dwydir[rmu] << "  " << dwzdir[rmu] << endl;
                cout << " xoutup, youtup, zoutup = " << Xoutup[rmu] << "  " << Youtup[rmu] << "  " << Zoutup[rmu] << endl;
                cout << " xinup, yinup, zinup = " << Xinup[rmu] << "  " << Yinup[rmu] << "  " << Zinup[rmu] << endl;
                cout << " xoutdw, youtdw, zoutdw = " << Xoutdw[rmu] << "  " << Youtdw[rmu] << "  " << Zoutdw[rmu] << endl;
                cout << " xindw, yindw, zindw = " << Xindw[rmu] << "  " << Yindw[rmu] << "  " << Zindw[rmu] << endl;
                cout <<   " &&&&&&&&&& END Study the direction of the muons &&&&&&&& " << endl;
                cout << endl << endl;   


		
		// Write in a special file only the numerical values for further fast analysis
		fprintf(filehmmeout,"%d %d %d %f %f %f \n",rmureal,jtr1,jmatch,ptotup,ptotdw,Pcov[rmu]);
		//=====================================================================
		
		// Filling Kath histograms (TList) for Thesis  (single muon tracks)
		
		if(upSign[rmu]==-1||upSign[rmu]==1){
		  fthetaDist->Fill(upthetaCosmic[rmu]);
		  fphiDist->Fill(upphiCosmic[rmu]);
		  fsingleZX->Fill(upzv[rmu],upxv[rmu]);
		  fsingleupZX->Fill(upzv[rmu],upxv[rmu]);
		  
		  Xpos[rmu] = upxv[rmu];
		  Zpos[rmu] = upzv[rmu];
		}
		if(dwSign[rmu]==-1||dwSign[rmu]==1){
		  fthetaDist->Fill(dwthetaCosmic[rmu]);
		  fphiDist->Fill(dwphiCosmic[rmu]);
		  fsingleZX->Fill(dwzv[rmu],dwxv[rmu]);
		  fsingledwZX->Fill(dwzv[rmu],dwxv[rmu]);
		  
		  Xpos[rmu] = dwxv[rmu];
		  Zpos[rmu] = dwzv[rmu];
		}
		
	      cout << endl;
	      cout << " Zpos, Xpos : " << Zpos[rmu] << " , " << Xpos[rmu] << endl;
	      cout << endl;

		
		
	    } // END (nEvent>=nevdeb1&&nEvent<=nevdeb2)
	    
	  } // END if( jmatch < 0 )
	  
	}  // END Loop on vector vecmu       for (Int_t j = 0; j < vecmu.size(); j++)
	
	
	          
		if(rmureal != nMuons) cout << " ERROR ATTENTION !!!! " << "rmureal = " << rmureal << "  nMuons = " << nMuons <<  endl;

	// Check of the TREE variables 

	if(nEvent>=nevdeb1&&nEvent<=nevdeb2)
	{ 
	  
	  cout << endl;   
	  cout << " nRun "  << nRun  <<  " nEvent " << nEvent <<  " nTracks " << nTracks << endl;
	  cout << " fracnMuTrk " << fracnMuTrk << " meanDist " << meanDist << " nMuons " << nMuons <<  " triggerFlag  " << triggerFlag << endl;
	  cout << " Num. muons : " << nMuons << "    Num. matched muons : " << rmumatch << "   Num. single track muons : " << rmusingtrk << endl;  
	  cout << " Num. matched mu+ : " << rmuplus << "    Num. matched mu- : " << rmuminus << "   Num. single track mu+ : " << srmuplus << "   Num. single track mu- :  " << srmuminus << "    Num. ambigous charge mu :  " << rmuplusminus << endl;  
//	  cout << " Xmin, Zmin : " << xvmin << " , " << zvmin << endl;
//	  cout << " Xmax, Zmax : " << xvmax << " , " << zvmax << endl;
          cout << "$$$$$$$$$ END SUMMARY EVENT : " << nEvent  << " $$$$$$$$$$$$ " << endl; 
  	  cout << " ***************************************" << endl;
  	  cout << endl;
  	  cout << "   Debug Frequency Track Used  " << endl;
  	    Int_t ikk=0;
  	    for (Int_t it = 0; it < nTracks; it++){
  	      if(freqtrack[it]>0){
  	        ikk++;
  	        cout << ikk << ")    Track Num. " << it << "  : " << freqtrack[it] << endl;
  	        if(freqtrack[it]>1){
  	        cout << ikk << ")    Track Num. " << it << "  : " << freqtrack[it] << "   USED MORE THAN 1 TIME " << endl;
  	        }  	        
  	      }  
  	    }
  	  // Kath calculations for Thesis  (extrem position of muons in ZX plane)  
  	  cout << " " <<  endl; 
  	  cout << " kath --- Num. muons : " << nMuons << "   event : " << nEvent << endl; 
	  for (Int_t j = 0; j < nMuons; j++){ 
//	    cout << " Xpos, Zpos : " << j << " " << Xpos[j] << " , " << Zpos[j] << endl;
	    if (Xpos[j] <= xvmin && Xpos[j] != -999){ 
	      xvmin = Xpos[j];
	    }
	    if (Zpos[j] <= zvmin  && Zpos[j] != -999){
	      zvmin = Zpos[j];
	    }
	    if (Xpos[j] > xvmax  && Xpos[j] != -999){ 
	      xvmax = Xpos[j];
	    }
	    if (Zpos[j] > zvmax  && Zpos[j] != -999){
	      zvmax = Zpos[j];
	    }
	  }
	 cout << " " <<  endl; 
         cout << " Xmin, Zmin : " << xvmin << " , " << zvmin << endl;
	 cout << " Xmax, Zmax : " << xvmax << " , " << zvmax << endl;
	 cout << " " <<  endl; 

  	    
	}  // END if(nEvent>=nevdeb1&&nEvent<=nevdeb2)


	
	/*============================================================
	 MC loop analysis
	=============================================================*/
	if (mcData) {
	  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
	  if (!eventHandler) {
		Printf("ERROR: Could not retrieve MC event handler");
		return;
	  }
	
	  AliMCEvent* mcEvent = eventHandler->MCEvent();
	  if (!mcEvent) {
		Printf("ERROR: Could not retrieve MC event");
		return;
	  }
	  Int_t muonesMC=0;
	  AliStack *mcStack = mcEvent->Stack();
	  if (!mcStack) return;
	  for (Int_t iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); iTracks++) 
	    {
		AliVParticle* track = mcEvent->GetTrack(iTracks);
		if ((!track)) {
			Printf("ERROR: Could not receive track %d (mc loop)", iTracks);
			continue;
		}
		if ((track->PdgCode()==13)||(track->PdgCode()==-13))
		{
			muonesMC++;
		}	
	    }
	  nMuonsMC = muonesMC;
	  cout << "******** No. of MUONS from MC: " << muonesMC << endl;
	  cout << "******** No. of MUONS from AnaCosmic (in TPC): " << nMuons << endl;
	  cout << "******* No. of trackas: " << fESD->GetNumberOfTracks() << endl;
/*

	for (UInt_t jtrack=0; jtrack < vecmu.size(); jtrack++)
	{
	    Int_t jtr1 = vecmu[jtrack];	
	    Int_t jmatch = vecmatch[jtr1];
            cout << "*** $$$$ Track jtr1= " << vecmu[jtrack] << "  Track matched jmatch= " << vecmatch[jtr1] << endl;
	    AliESDtrack *esdTrack = (AliESDtrack*)fESD->GetTrack(vecmu[jtrack]);
	    Int_t labelTrack = TMath::Abs(esdTrack->GetLabel());
	    cout << "**** ---- track label: " << labelTrack << " ---- track ID: " << esdTrack->GetID() <<endl;
	    cout << "**** ---- uptracmu: " << uptracmu[jtr1] << " dwtracmu: "<< dwtracmu[jtr1]<<endl;
	    if(jtr1>jmatch&&jmatch>=0) continue;  // already matched jmatch-->jtr1
	}
*/
	  AliESDtrack *UpESDTrack=0x0;
	  AliESDtrack *DwESDTrack=0x0;
	  AliESDtrack *ESDTrack=0x0;
	  for (Int_t imuon=0; imuon<nMuons; imuon++)
	    {
	//!	cout << "**** ---- imuon: " << imuon << " pCov: " << Pcov[imuon] << endl;
	//!        cout << "**** ---- uptracmu: " << uptracmu[imuon] << " dwtracmu: "<< dwtracmu[imuon]<<endl;
		if (uptracmu[imuon]==-10 && dwtracmu[imuon]>=0) //! single track (dwtrack)
		{
			ESDTrack = fESD->GetTrack(dwtracmu[imuon]);	
		}
		if (uptracmu[imuon]>=0 && dwtracmu[imuon]==-10) //! single track (uptrack)
		{
			ESDTrack = fESD->GetTrack(uptracmu[imuon]);	
		}
		if (uptracmu[imuon]>=0 && dwtracmu[imuon]>=0) //! matched tracks (up and dw track)
		{
			ESDTrack = fESD->GetTrack(uptracmu[imuon]);	
		}
		Int_t labelTrack = TMath::Abs(ESDTrack->GetLabel());
		//!Int_t labelTrack = TMath::Abs(ESDTrack->GetID());
//!		cout << "**** ---- label track: " << labelTrack <<  " ---- label ID: " << ESDTrack->GetID() <<endl;
		AliMCParticle* mcParticle = (AliMCParticle*)mcEvent->GetTrack(labelTrack);
		TParticle *particle = mcStack->Particle(labelTrack);
//!		cout << "**** ---- pCovMC: " << mcParticle->P() << " --- pCovMC-stack: " << particle->P() << " --- pdg: " << particle->GetPdgCode()<<  endl;
	      	Int_t lp = mcStack->Particle(labelTrack)->GetUniqueID();
//!		cout << "**** ---- unique ID track: " << lp << endl;
		pMC[imuon] = particle->P();
		pxMC[imuon] = particle->Px();
		pyMC[imuon] = particle->Py();
		pzMC[imuon] = particle->Pz();
		xMC[imuon] = particle->Vx();
		yMC[imuon] = particle->Vy();
		zMC[imuon] = particle->Vz();
		tetaMC[imuon] = particle->Theta();
		energyMC[imuon] = particle->Energy();
		pdgCode[imuon] = particle->GetPdgCode();
	    }


	//! loop over tracks
	  cout << "**** $$$$ Summary of event $$$$ ****" << endl;
	  for (Int_t imuon=0; imuon<nMuons; imuon++)
	    {	
		cout << " --- imuon: " << imuon << " pCovRec: " << Pcov[imuon] <<  " pCovMC: "<< pMC[imuon] << " Pdg: " << pdgCode[imuon] << " charge cov: " << Chargecov[imuon]<< endl;
	    }

	} // if (mcData)

        fDatree->Fill();   // Fill the TREE
	neverFill = kFALSE;
	
	nRecMuons = nRecMuons + nMuons;
//	cout << "===== Kath Rec nMuons = " << nMuons <<  endl;
//	cout << "===== Kath # of fmatchZX = " << ncountK <<  endl;
	
//	 cout << " Xmin, Zmin : " << xvmin << " , " << zvmin << endl;
//	 cout << " Xmax, Zmax : " << xvmax << " , " << zvmax << endl;



      } // END if (nMuons > 0 && flagintev==0)

    } //---------------------->>> END if nTracks > 1 <<<------------------------- 
    
//    cout << "===== Kath # of Reconstructed (accumulated) Muons = " << nRecMuons <<  endl;
    
    if (nMuons > 0)fmumult->Fill(nMuons);

  }    //END if nTracks > 0
  
  fTriggerMask->Fill(triggerMask);
//fDatree->Fill();
} // END Triggers OB1 OB0 OB3 SCO AMU  

  //============================================================================================ 
  // 				Post output data.
  //============================================================================================ 

//PostData(0, fDatree);
//PostData(1, fListHist);

}      

//______________________________________________________________________________
AliExternalTrackParam *AliCosmics::MakeTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // Make a atrack using the kalman update of track0 and track1
  //
  AliExternalTrackParam *par1R= new AliExternalTrackParam(*track1);
  	      // Ruben instructions
	      par1R->Invert();
	      if(par1R->Propagate(track0->GetAlpha(),track0->GetX(),AliTracker::GetBz())){
             } else {
             }
  return par1R;
}

AliExternalTrackParam *AliCosmics::MakeCombinedTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // Make combined track
  //
  //
  AliExternalTrackParam * par1T = MakeTrack(track0,track1);
  AliExternalTrackParam * par0U = new AliExternalTrackParam(*track0);
  //
  UpdateTrack(*par0U,*par1T);
  delete par1T;
  return par0U;
}

//
//________________________________________________________________________
Double_t AliCosmics::GetPcov(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1){
  //
  // Make combined track
  //
  //
  AliExternalTrackParam * par1T = MakeTrack(track0,track1);
  AliExternalTrackParam * par0U = new AliExternalTrackParam(*track0);
  //
  Double_t pCovMI;
  pCovMI = GetPCovMI(*par0U,*par1T);
  delete par1T;
  return pCovMI;
}


//______________________________________________________________________________
Double_t AliCosmics::GetPCovMI(AliExternalTrackParam &track1, const AliExternalTrackParam &track2){
  //
  // Update track 1 with track 2
  //
  //
  //
  TMatrixD vecXk(5,1);    // X vector
  TMatrixD covXk(5,5);    // X covariance 
  TMatrixD matHk(5,5);    // vector to mesurement
  TMatrixD measR(5,5);    // measurement error 
  TMatrixD vecZk(5,1);    // measurement
  //
  TMatrixD vecYk(5,1);    // Innovation or measurement residual
  TMatrixD matHkT(5,5);
  TMatrixD matSk(5,5);    // Innovation (or residual) covariance
  TMatrixD matKk(5,5);    // Optimal Kalman gain
  TMatrixD mat1(5,5);     // update covariance matrix
  TMatrixD covXk2(5,5);   // 
  TMatrixD covOut(5,5);
  //

  Double_t param[4];
  Double_t pcovMI;


  Double_t *param1=(Double_t*) track1.GetParameter();
  Double_t *covar1=(Double_t*) track1.GetCovariance();
  Double_t *param2=(Double_t*) track2.GetParameter();
  Double_t *covar2=(Double_t*) track2.GetCovariance();
  //
  // copy data to the matrix
  for (Int_t ipar=0; ipar<5; ipar++){
    for (Int_t jpar=0; jpar<5; jpar++){
      covXk(ipar,jpar) = covar1[track1.GetIndex(ipar, jpar)];
      measR(ipar,jpar) = covar2[track2.GetIndex(ipar, jpar)];
      matHk(ipar,jpar)=0;
      mat1(ipar,jpar)=0;
    }
    vecXk(ipar,0) = param1[ipar];
    vecZk(ipar,0) = param2[ipar];
    matHk(ipar,ipar)=1;
    mat1(ipar,ipar)=0;
  }
  //
  //
  //
  //
  //
  vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
  vecXk += matKk*vecYk;                      //  updated vector 
  covXk2 = (mat1-(matKk*matHk));
  covOut =  covXk2*covXk; 
  //
  //
  //
  // copy from matrix to parameters
  if (0) {
    vecXk.Print();
    vecZk.Print();
    //
    measR.Print();
    covXk.Print();
    covOut.Print();
    //
    track1.Print();
    track2.Print();
  }

  for (Int_t ipar=0; ipar<5; ipar++){
    param1[ipar]= vecXk(ipar,0) ;
    param[ipar]=param1[ipar];

    for (Int_t jpar=0; jpar<5; jpar++){
      covar1[track1.GetIndex(ipar, jpar)]=covOut(ipar,jpar);
    }
  }

  pcovMI =TMath::Sqrt(1.+ param[3]*param[3])/TMath::Abs(param[4]);
  
  return pcovMI;

}

//______________________________________________________________________________
void AliCosmics::UpdateTrack(AliExternalTrackParam &track1, const AliExternalTrackParam &track2){
  //
  // Update track 1 with track 2
  //
  //
  //
  TMatrixD vecXk(5,1);    // X vector
  TMatrixD covXk(5,5);    // X covariance 
  TMatrixD matHk(5,5);    // vector to mesurement
  TMatrixD measR(5,5);    // measurement error 
  TMatrixD vecZk(5,1);    // measurement
  //
  TMatrixD vecYk(5,1);    // Innovation or measurement residual
  TMatrixD matHkT(5,5);
  TMatrixD matSk(5,5);    // Innovation (or residual) covariance
  TMatrixD matKk(5,5);    // Optimal Kalman gain
  TMatrixD mat1(5,5);     // update covariance matrix
  TMatrixD covXk2(5,5);   // 
  TMatrixD covOut(5,5);
  //
  Double_t *param1=(Double_t*) track1.GetParameter();
  Double_t *covar1=(Double_t*) track1.GetCovariance();
  Double_t *param2=(Double_t*) track2.GetParameter();
  Double_t *covar2=(Double_t*) track2.GetCovariance();
  //
  // copy data to the matrix
  for (Int_t ipar=0; ipar<5; ipar++){
    for (Int_t jpar=0; jpar<5; jpar++){
      covXk(ipar,jpar) = covar1[track1.GetIndex(ipar, jpar)];
      measR(ipar,jpar) = covar2[track2.GetIndex(ipar, jpar)];
      matHk(ipar,jpar)=0;
      mat1(ipar,jpar)=0;
    }
    vecXk(ipar,0) = param1[ipar];
    vecZk(ipar,0) = param2[ipar];
    matHk(ipar,ipar)=1;
    mat1(ipar,ipar)=0;
  }
  //
  //
  //
  //
  //
  vecYk = vecZk-matHk*vecXk;                 // Innovation or measurement residual
  matHkT=matHk.T(); matHk.T();
  matSk = (matHk*(covXk*matHkT))+measR;      // Innovation (or residual) covariance
  matSk.Invert();
  matKk = (covXk*matHkT)*matSk;              //  Optimal Kalman gain
  vecXk += matKk*vecYk;                      //  updated vector 
  covXk2 = (mat1-(matKk*matHk));
  covOut =  covXk2*covXk; 
  //
  //
  //
  // copy from matrix to parameters
  if (0) {
    vecXk.Print();
    vecZk.Print();
    //
    measR.Print();
    covXk.Print();
    covOut.Print();
    //
    track1.Print();
    track2.Print();
  }

  for (Int_t ipar=0; ipar<5; ipar++){
    param1[ipar]= vecXk(ipar,0) ;
    for (Int_t jpar=0; jpar<5; jpar++){
      covar1[track1.GetIndex(ipar, jpar)]=covOut(ipar,jpar);
    }
  }
}

//________________________________________________________________________
// Metodo para poner el style del macro de figuras de ALICE
void AliCosmics::SetStyle() {                      
    cout << "Setting style!" << endl;
    
    gStyle->Reset("Plain");
    //  gStyle->SetOptTitle(0);
    //  gStyle->SetOptStat(0);
    gStyle->SetCanvasColor(10);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetFrameFillColor(kWhite);
    gStyle->SetPadColor(10);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetHistLineWidth(1);
    gStyle->SetHistLineColor(kRed);
    gStyle->SetFuncWidth(2);
    gStyle->SetFuncColor(kGreen);
    gStyle->SetLineWidth(2);
    gStyle->SetLabelSize(0.045,"xyz");
    gStyle->SetLabelOffset(0.01,"y");
    gStyle->SetLabelOffset(0.01,"x");
    gStyle->SetLabelColor(kBlack,"xyz");
    gStyle->SetTitleSize(0.05,"xyz");
    gStyle->SetTitleOffset(1.25,"y");
    gStyle->SetTitleOffset(1.2,"x");
    gStyle->SetTitleFillColor(kWhite);
    gStyle->SetTextSizePixels(26);
    gStyle->SetTextFont(42);
  }
//________________________________________________________________________
void AliCosmics::Terminate(Option_t *) 
{

       //   Close the open files
  // In AliEn, Terminate does not know filehmmeout/intereve, so check to avoid segfaults - M.S.
  if (filehmmeout) fclose(filehmmeout);
  if (intereve) fclose(intereve);
  
//  TString fileCanvasOut = "r152599-ch020.12.root";
  TString CanvasName;
  
   // Draw some results in Canvases
  
  //=========================================================================================
  TCanvas *mumult = new TCanvas("Muon Multiplicity Distribution","",200,10,700,500);
  
  CanvasName = "canvas-mumult.root";
  mumult->cd();
  
  gPad->SetLogy();
    
  // In AliEn, Terminate does not know fmumult, so check and reaload to avoid segfaults - M.S.
  if (!fmumult) {
    fListHist = (TList*)GetOutputData(1);
    fmumult = (TH1F*)(fListHist->FindObject("mumult"));
  }

  fmumult->GetXaxis()->SetTitle("# reco #mu");
  fmumult->GetXaxis()->SetLabelSize(0.03);
  fmumult->GetXaxis()->SetTickLength(0.02);
  fmumult->GetYaxis()->SetTitle("Counts");
  fmumult->GetYaxis()->SetLabelSize(0.03);
  fmumult->GetYaxis()->SetTickLength(0.02);
  fmumult->GetYaxis()->SetTitleOffset(0.8);
  fmumult->SetMarkerStyle(20);		//20
  fmumult->SetMarkerSize(1);		//1.0
  fmumult->SetMarkerColor(2);
  fmumult->SetLineColor(2);
  fmumult->SetLineWidth(2);
  fmumult->Draw("E");			//"E"

  gPad->Update();
  TPaveStats *st_mumult = (TPaveStats*)fmumult->GetListOfFunctions()->FindObject("stats");
  st_mumult->SetX1NDC(0.69);
  st_mumult->SetY1NDC(0.72);
  st_mumult->SetX2NDC(0.89);
  st_mumult->SetY2NDC(0.88);
  gPad->Modified();
  gPad->Update();
  
  mumult->SaveAs(CanvasName);
  
  //=========================================================================================
  
  SetStyle();
  
  TCanvas *angDist = new TCanvas("Angular Distributions (Zenith & Azimuth)","",200,10,700,500);
  
  CanvasName = "canvas-AngDist.root";
  
  angDist->Divide(2,1);
 
  angDist->cd(1);
  
  // In AliEn, Terminate does not know fthetaDist, so check and reaload to avoid segfaults - M.S.
  if (!fthetaDist)
    fthetaDist = (TH1F*)(fListHist->FindObject("ThetaDist"));
  
  fthetaDist->SetTitle("");
  fthetaDist->GetXaxis()->SetTitle("Zenith Angle #theta (Deg)");
//  fthetaDist->GetXaxis()->SetRange(xmintheta,xmaxtheta);
//  fthetaDist->GetXaxis()->SetLimits(xmintheta,xmaxtheta);   //18,34
// fthetaDist->SetAxisRange(xmintheta,xmaxtheta, "X");
  fthetaDist->GetXaxis()->SetLabelSize(0.03);
//  fthetaDist->GetXaxis()->SetTickLength(0.02);
  fthetaDist->GetYaxis()->SetTitle("Counts");
//  fthetaDist->GetYaxis()->SetLabelSize(0.03);
//  fthetaDist->GetYaxis()->SetTickLength(0.02);
  fthetaDist->GetYaxis()->SetTitleOffset(1.4);
  fthetaDist->SetMarkerStyle(20);		//20
  fthetaDist->SetMarkerSize(1.0);		//1.0
  fthetaDist->SetMarkerColor(2);
  fthetaDist->SetLineColor(2);
  fthetaDist->SetLineWidth(2);
//  fthetaDist->Draw("p");			//"E" "p"
  
  TLatex * text_theta = new TLatex(0.20,.84,"ALICE");  //0.47,.83
  text_theta->SetNDC();
  text_theta->SetTextFont(42);
//  text_theta->Draw();
    
 angDist->cd(2);
  
  // In AliEn, Terminate does not know fphiDist, so check and reaload to avoid segfaults - M.S.
  if (!fphiDist)
    fphiDist = (TH1F*)(fListHist->FindObject("PhiDist"));
  
  fphiDist->SetTitle("");
  fphiDist->GetXaxis()->SetTitle("Azimuth Angle #phi (Deg)");
  fphiDist->GetXaxis()->SetLabelSize(0.03);
//  fphiDist->GetXaxis()->SetTickLength(0.02);
  fphiDist->GetYaxis()->SetTitle("Counts");
//  fphiDist->GetYaxis()->SetLabelSize(0.03);
//  fphiDist->GetYaxis()->SetTickLength(0.02);
  fphiDist->GetYaxis()->SetTitleOffset(1.4);
  fphiDist->SetMarkerStyle(0);
  fphiDist->SetMarkerSize(20);
  fphiDist->SetMarkerColor(2);
  fphiDist->SetLineColor(2);
  fphiDist->SetLineWidth(2);
//  fphiDist->Draw("p");
  
  TLatex * text_phi = new TLatex(0.20,.84,"ALICE");  //0.47,.83
  text_phi->SetNDC();
  text_phi->SetTextFont(42);  
//  text_phi->Draw();
    
//  angDist->SaveAs(CanvasName);  
   
 //=========================================================================================
  SetStyle();
  
  TCanvas *SpatialZX = new TCanvas("Muon Spatial Distribution at plane Y=0","",200,10,700,500);
  
  CanvasName = "canvas-ZXmuon.root";
  SpatialZX->cd();
  
//  gStyle->SetOptTitle(kFALSE);
  
  // In AliEn, Terminate does not know fmatchZX/fsingleZX, so check and reaload to avoid segfaults - M.S.
  if (!fmatchZX)
    fmatchZX = (TH2F*)(fListHist->FindObject("matchZX"));
  
 // fmatchZX->SetStats(kFALSE);
  fmatchZX->SetTitle("");
  fmatchZX->GetXaxis()->SetTitle("Z [cm]");
  fmatchZX->GetXaxis()->SetLabelSize(0.03);
  fmatchZX->GetXaxis()->SetTickLength(0.02);
  fmatchZX->GetYaxis()->SetTitle("X [cm]");
  fmatchZX->GetYaxis()->SetLabelSize(0.03);
  fmatchZX->GetYaxis()->SetTickLength(0.02);
  fmatchZX->GetYaxis()->SetTitleOffset(0.8);
  fmatchZX->SetMarkerStyle(20);	
  fmatchZX->SetMarkerSize(1.2);
  fmatchZX->SetMarkerColor(9);   //9-blue
//  fmatchZX->Draw();

  if (!fsingleZX)
    fsingleZX = (TH2F*)(fListHist->FindObject("singleZX"));
  
  fsingleZX->SetStats(kTRUE);
  fsingleZX->SetTitle("");
  fsingleZX->SetMarkerStyle(22);		//20
  fsingleZX->SetMarkerSize(1.3);		//1.0
  fsingleZX->SetMarkerColor(2);  //2-red
//  fsingleZX->Draw("sames");
  
  TLatex * text_ZX = new TLatex(0.48,.82,"ALICE");
  text_ZX->SetNDC();
  text_ZX->SetTextFont(42);
  text_ZX->SetTextSize(0.056);
  text_ZX->SetLineWidth(2);
//  text_ZX->Draw();
  
  TLegend *leg_ZX = new TLegend(0.74,0.76,0.86,0.87);
  leg_ZX->AddEntry(fmatchZX,"Matched tracks","p");
  leg_ZX->AddEntry(fsingleZX,"Single tracks","p");
  leg_ZX->SetTextSize(0.025);
  leg_ZX->SetLineWidth(0);
  leg_ZX->SetLineColor(0);
//  leg_ZX->Draw();
  
//  SpatialZX->SaveAs(CanvasName);  
  
  //=========================================================================================
  TCanvas *SpatialZXupdw = new TCanvas("Muon Spatial Distribution at plane Y=0","",200,10,700,500);
  
  CanvasName = "canvas-ZXupdwmuon.root";
  SpatialZXupdw->cd();
  
//  gStyle->SetOptTitle(kFALSE);
  
 // fmatchZX->SetStats(kFALSE);
  fmatchZX->SetTitle("");
  fmatchZX->GetXaxis()->SetTitle("Z [cm]");
  fmatchZX->GetXaxis()->SetLabelSize(0.03);
  fmatchZX->GetXaxis()->SetTickLength(0.02);
  fmatchZX->GetYaxis()->SetTitle("X [cm]");
  fmatchZX->GetYaxis()->SetLabelSize(0.03);
  fmatchZX->GetYaxis()->SetTickLength(0.02);
  fmatchZX->GetYaxis()->SetTitleOffset(0.8);
  fmatchZX->SetMarkerStyle(20);	
  fmatchZX->SetMarkerSize(1.2);
  fmatchZX->SetMarkerColor(9);   //9-blue
//  fmatchZX->Draw();
  
  // In AliEn, Terminate does not know fsingleupZX/fsingledwZX, so check and reaload to avoid segfaults - M.S.
  if (!fsingleupZX)
    fsingleupZX = (TH2F*)(fListHist->FindObject("singleupZX"));
  if (!fsingledwZX)
    fsingledwZX = (TH2F*)(fListHist->FindObject("singledwZX"));

//  fsingleupZX->SetStats(kTRUE);
  fsingleupZX->SetTitle("");
  fsingleupZX->SetMarkerStyle(22);
  fsingleupZX->SetMarkerSize(1.3);
  fsingleupZX->SetMarkerColor(6);  //2-red  6-hot pink 3-bright green
//  fsingleupZX->Draw("sames");
  
  fsingledwZX->SetTitle("");
  fsingledwZX->SetMarkerStyle(23);
  fsingledwZX->SetMarkerSize(1.3);
  fsingledwZX->SetMarkerColor(3);  //2-red  6-hot pink 3-bright green
//  fsingledwZX->Draw("sames");
      
  TLatex * text_ZXupdw = new TLatex(0.47,.82,"ALICE");
  text_ZXupdw->SetNDC();
  text_ZXupdw->SetTextFont(42);  //132
//  text_ZXupdw->Draw();
  
  TLegend *leg_ZXupdw = new TLegend(0.77,0.77,0.88,0.88);
  leg_ZXupdw->AddEntry(fmatchZX,"Matched Muons","p");
  leg_ZXupdw->AddEntry(fsingleupZX,"Single UP Muons","p");
  leg_ZXupdw->AddEntry(fsingledwZX,"Single DOWN Muons","p");
//  leg_ZXupdw->Draw();
  
//  SpatialZXupdw->SaveAs(CanvasName);  
  //=========================================================================================
  
  TCanvas *mudist = new TCanvas("Distance between matched tracks up and down","",200,10,700,500);
  
  CanvasName = "canvas-mudist.root";
  mudist->cd();
  
//  gPad->SetLogy();

  // In AliEn, Terminate does not know fupdwMuDist, so check and reaload to avoid segfaults - M.S.
  if (!fupdwMuDist)
    fupdwMuDist = (TH1F*)(fListHist->FindObject("updwMuDist"));
    
  fupdwMuDist->GetXaxis()->SetTitle("Track distance (cm)");
  fupdwMuDist->GetXaxis()->SetLabelSize(0.03);
  fupdwMuDist->GetXaxis()->SetTickLength(0.02);
  fupdwMuDist->GetYaxis()->SetTitle("Counts");
  fupdwMuDist->GetYaxis()->SetLabelSize(0.03);
  fupdwMuDist->GetYaxis()->SetTickLength(0.02);
  fupdwMuDist->GetYaxis()->SetTitleOffset(0.8);
//  fupdwMuDist->SetMarkerStyle(20);
//  fupdwMuDist->SetMarkerSize(1);
//  fupdwMuDist->SetMarkerColor(2);
  fupdwMuDist->SetLineColor(2);
  fupdwMuDist->SetLineWidth(2);
//  fupdwMuDist->Draw();			//"E"
  
  
  //TLatex * text_mudist = new TLatex(0.47,.82,"ALICE");
  TLatex * text_mudist = new TLatex(0.273628,0.8296089,"Distance between track up and track down in matched muons");
  text_mudist->SetNDC();
  text_mudist->SetTextFont(42);
  text_mudist->SetTextSize(0.03072626);
  text_mudist->SetLineWidth(2);
//  text_mudist->Draw();
  
//  mudist->SaveAs(CanvasName);
  
  //=========================================================================================
  
  
  

  // Draw result for the Trigger Rate and Rate of the events with at least 1 muon
  // In AliEn, Terminate does not know fnTrksTOFOB0, so check and reaload to avoid segfaults - M.S.
  if (!fnTrksTOFOB0){
    fnTrksTOFOB0 = (TH1F*)(fListHist->FindObject("nTrksTOFOB0"));
    fnTrksTOFOB1 = (TH1F*)(fListHist->FindObject("nTrksTOFOB1"));
    fnTrksTOFOB3 = (TH1F*)(fListHist->FindObject("nTrksTOFOB3"));
    fnTrksACORDEAMU = (TH1F*)(fListHist->FindObject("nTrksACORDEAMU"));
    fnTrksTOFOB0Muons = (TH1F*)(fListHist->FindObject("nTrksTOFOB0Muons"));
    fnTrksTOFOB1Muons = (TH1F*)(fListHist->FindObject("nTrksTOFOB1Muons"));
    fnTrksTOFOB3Muons = (TH1F*)(fListHist->FindObject("nTrksTOFOB3Muons"));
    fnTrksAMUMuons = (TH1F*)(fListHist->FindObject("nTrksAMUMuons"));
  }

  TCanvas *c0 = new TCanvas("Num. of Triggers","Trigger",20,20,510,510);
  TCanvas *c1 = new TCanvas("Muon Multiplicity","Muons",20,20,510,510);
//  TCanvas *c2 = new TCanvas("Freq. Track","Freq. Track used",20,20,510,510);
//  TCanvas *c3 = new TCanvas("Freq. Track Matched","Freq. Track used Matched",20,20,510,510);
//  TCanvas *c4 = new TCanvas("Freq. Track Sing.","Freq. Track used Single",20,20,510,510);
  
  c0->Divide(2,2);  
  c0->cd(1);
  gPad->SetLogy();
  fnTrksTOFOB0->Draw();
  c0->cd(2);
  gPad->SetLogy();
  fnTrksTOFOB1->Draw();
  c0->cd(3);
  gPad->SetLogy(); 
  fnTrksTOFOB3->Draw();
  c0->cd(4);
  gPad->SetLogy(); 
  fnTrksACORDEAMU->Draw();
  
  c0->SaveAs("TriggerNumber.root");

  c1->Divide(2,2);  
  c1->cd(1);
  gPad->SetLogy();
  fnTrksTOFOB0Muons->Draw();
  c1->cd(2);
  gPad->SetLogy();
  fnTrksTOFOB1Muons->Draw();
  c1->cd(3);
  gPad->SetLogy(); 
  fnTrksTOFOB3Muons->Draw();  
  c1->cd(4);
  gPad->SetLogy(); 
  fnTrksAMUMuons->Draw();

  c1->SaveAs("MuonNumber.root");
   
  // In AliEn, Terminate does not know fmeanDistMu3, so check and reaload to avoid segfaults - M.S.
  if (!fmeanDistMu3){
    fmeanDistMu3 = (TH2F*)(fListHist->FindObject("meanDistMu3"));
    fmeanDistMu4 = (TH2F*)(fListHist->FindObject("meanDistMu4"));
  }

TCanvas *c5 = new TCanvas("Interaction Events","Study Interac. Events",20,20,510,510);       
  c5->cd();
  fmeanDistMu3->GetXaxis()->SetTitle("Average_{Dist} (cm)");
  fmeanDistMu3->GetYaxis()->SetTitle("N_{Tracks}/N_{mu}");
  fmeanDistMu3->SetMarkerStyle(25);
  fmeanDistMu3->SetMarkerColor(kBlue);
  fmeanDistMu3->SetMarkerSize(0.4);
  fmeanDistMu3->Draw();
  fmeanDistMu4->SetMarkerStyle(25);
  fmeanDistMu4->SetMarkerColor(kRed);
  fmeanDistMu4->SetMarkerSize(0.4);
  fmeanDistMu4->Draw("same");
  TLegend *Mu = new TLegend(0.55,0.65,0.76,0.82);
  Mu->AddEntry(fmeanDistMu3,"Interaction events","p");
  Mu->AddEntry(fmeanDistMu4,"Multi-Muon events","p");
  Mu->Draw();
 
  c5->SaveAs("MuInt_MultiMu_nrun.root");
  
}


//**************************************************************************
Bool_t AliCosmics::VZeroBG(AliESDVZERO* vzeroring) const
{
	AliESDVZERO *vzero = vzeroring;
	for(Int_t i=0;i<32;i++){
		if (vzero->BGTriggerV0A(i) > 0) return kTRUE;
		if (vzero->BGTriggerV0C(i) > 0) return kTRUE;

	} // end loop over cells scintillators
	return kFALSE;
}
