//              Last update: May 30th 2012
//              Adding:
//                      --> branches: chunk path, # of event, ACORDE hits, ACORDE MCN and VZERO BG rejection
//              From: Mario Rodriguez Cahuantiz (FCFM-BUAP), Puebla-MX <mrodrigu@mail.cern.ch>

#ifndef AliCosmics_cxx 
#define AliCosmics_cxx


class TH1F;
class TH2F;
class AliESDEvent;
class TList;
class TTree;
class AliESDVZERO;

//#include <vector>;
#include "AliAnalysisTask.h"
#include "TObject.h"
#include "TTreeStream.h"


class AliCosmics : public AliAnalysisTask {
 public:
  AliCosmics(const char *name = "AliCosmics", const Bool_t isMC = kFALSE);
  virtual ~AliCosmics();
  
  virtual void   ConnectInputData(Option_t *); //! Connects input data to class analysis
  virtual void   CreateOutputObjects();        //! Creates output object (cosmicTree)
  virtual void   Exec(Option_t *option);       //! Execution class
  virtual void   Terminate(Option_t *);        //! Terminate class
  virtual Bool_t VZeroBG(AliESDVZERO* vzeroring)const; //! VZER0 info
  AliExternalTrackParam *MakeTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1);
  AliExternalTrackParam *MakeCombinedTrack(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1);
  Double_t GetPcov(const AliExternalTrackParam *track0, const AliExternalTrackParam *track1);

  void UpdateTrack(AliExternalTrackParam &track0, const AliExternalTrackParam &track1);
  Double_t GetPCovMI(AliExternalTrackParam &track0, const AliExternalTrackParam &track1);
  
  void SetStyle();

 private:
  AliESDEvent *fESD;    //ESD object
  TList       *fListHist;
  TTree       *fDatree; 
  TH1F        *fTrackThetaCosmicOB0;
  TH1F        *fTrackThetaCosmicOB1;
  TH1F        *fTrackThetaCosmicOB3;
  TH1F        *fTrackPhiCosmic;
  TH1F        *fTriggerMask; 
  TH1F        *fnTrksSPD; 
  TH1F        *TotalTrksSPD; 
  TH1F        *fnTrksTOFOB0; 
  TH1F        *TotalTrksTOFOB0; 
  TH1F        *fnTrksTOFOB1; 
  TH1F        *TotalTrksTOFOB1; 
  TH1F        *fnTrksTOFOB3; 
  TH1F        *TotalTrksTOFOB3; 
  TH1F        *fnTrksACORDEAMU; 
  TH1F        *TotalTrksACORDEAMU; 
  TH1F        *fnTrksTOFOB0Muons; 
  TH1F        *fnTrksTOFOB1Muons; 
  TH1F        *fnTrksTOFOB3Muons; 
  TH1F        *fnTrksAMUMuons; 
  TH1F        *fnTrksASLMuons; 
  TH1F        *fnTrksSCOMuons; 
  TH2F        *fmeanDistMu3; 
  TH2F        *fmeanDistMu4; 
  TH1F        *fdistHist;
  TH1F        *Debtrackused;
  TH1F        *Debtrackusedmatch;
  TH1F        *Debtrackusedsing;
  
  // ======== Thesis ========= //
  
  TH1F        *fmumult;
  
  TH1F        *fphiDist;           
  TH1F        *fthetaDist;
  
  TH1F        *fupdwMuDist;
  
  TH2F        *fmatchZX;
  TH2F        *fsingleZX;
  TH2F        *fsingleupZX;
  TH2F        *fsingledwZX;
  
  // ======== Thesis end========= //
  

  ULong64_t triggerMask;

  Int_t triggerFlag;

  Int_t nEvent,nRun,nTracks,nMuons,nMuonsMC,ntpar,maxntpar,ntmaxpar,ntup,ntdw;
  Int_t flagintev;
  Int_t itag, ntmatch, ntsingle;
  Float_t magfield;
  UInt_t timestp;
  
  Int_t uptpcNcls[500],dwtpcNcls[500];

  Int_t uptracmu[500],dwtracmu[500];

  //%%2018 Float_t PRes[500];

  Float_t uptpcChi2[500],dwtpcChi2[500];
  
  // Float_t  upratioChi2cls[500],dwratioChi2cls[500];
  
  //%%2018 Float_t upSign[500],dwSign[500],uptpcSignal[500],dwtpcSignal[500];
  Float_t upSign[500],dwSign[500];
   
  Float_t meanDist,fracnMuTrk,rmsDist,rdist1t2,rdistmatch,cost1t2,rDist; 

  Float_t uptheta[500],upphi[500],upthetaCosmic[500],upphiCosmic[500];

  Float_t dwtheta[500],dwphi[500],dwthetaCosmic[500],dwphiCosmic[500];
  
  //%% Float_t upphiCosmiccor[500], dwphiCosmiccor[500], pxz; 

  //%%2018 Float_t Xinv[500],E,P,Pt1,Pt2,Et1,Et2;
  Float_t E,P,Pt1,Pt2,Et1,Et2;

  Float_t Pup[500],Ptupinv[500],upxv[500],upyv[500],upzv[500],uppx[500],uppy[500],uppz[500],Pup2[500],SigUp2[500],Sig1Ptup[500];
 
  Float_t Pdw[500],Ptdwinv[500],dwxv[500],dwyv[500],dwzv[500],dwpx[500],dwpy[500],dwpz[500],Pdw2[500],SigDw2[500],Sig1Ptdw[500];
 
  //%%2018 Float_t Pmu[500]; //momentum mu calculated with sigma Pweight
  
  //%% Float_t Pcov[500], Pcovupdw[500], Pcovdwup[500]; //covariant momentum 
  Float_t Pcov[500], Chargecov[500];
 
  //%% Float_t PmuMed[500]; //Mean momentum mu 

  //%% Float_t uncTimeTOFup[500],uncTimeTOFdw[500],corrTimeTOFup[500],corrTimeTOFdw[500],ToT_TOFup[500],ToT_TOFdw[500],DzTOFup[500],DzTOFdw[500],DxTOFup[500],DxTOFdw[500];
  Float_t uncTimeTOFup[500],uncTimeTOFdw[500],corrTimeTOFup[500],corrTimeTOFdw[500];

  Float_t Xinup[500],Yinup[500],Zinup[500],Xindw[500],Yindw[500],Zindw[500],Xoutup[500],Youtup[500],Zoutup[500],Xoutdw[500],Youtdw[500],Zoutdw[500]; 	
 
  Float_t upxdir[500],upydir[500],upzdir[500];
  Float_t dwxdir[500],dwydir[500],dwzdir[500];  
  Float_t uptetadir[500],upphidir[500];
  Float_t dwtetadir[500],dwphidir[500];  
  Char_t nMatchup[500],nMatchdw[500];
  Int_t indexTOFup[500],indexTOFdw[500];

  TObjString fFileName; 
  TObjString fFileNameChunk;

  Int_t eventNumberInFile;
  Int_t fEventInFile; 

  Int_t mcnACO;
  Int_t hitsACO[60];
  Int_t flagV0;

  Bool_t mcData;
  
  ClassDef(AliCosmics, 1); // example of analysis
};

  //  Initialization of some variables

Int_t nInter    = 0;   // number of interaction events
char chnk[50];         // chunk name
Int_t multimuev = 0;   // number of multimuon events (Nmu>3)
Int_t badmatch = 0;    // number of bad matching
Float_t fCutDirPar = 0.990; //cut direction cos(theta) = 8 Deg
Float_t minCutEnergy = 0.5; // cut energy tracks (Gev)
Float_t zcut = 250;// cut in dimension of z for TPC
Float_t xcut = 250;// cut in dimension of x for TPC
Int_t cutnumcls = 50; // cut number of clusters
Int_t cutsingncls = 50; // cut number of clusters for single track
Float_t cutchi2track = 4;  //cut on Chi2 of the track
Float_t minCutDist = 6; // cut of the distance between t1 and t2 (cm) -- Usamos 3cm --
Float_t RadGrad = 57.297469;  // factor to trasform radiant in degree
FILE *filehmmeout;    // file info of high number of track events
FILE *intereve;   // file info of interaction events

char currchunk[200], prevchunk[200];
Bool_t neverFill=kFALSE;

// ======== Thesis ========= //     (Add by Kath)
Int_t ncountK = 0;                                      // Count the number of anything
Int_t nRecMuons = 0;                                      // Count the number of reconstructed muons
Float_t MuDistupdw[500];                                  //Distance between track up and down in matched muons
Float_t xvmin = 0;                                       //Xmin, Xmax, Zmin, Zmax area projection of the TPC
Float_t xvmax = 0;
Float_t zvmin = 0;
Float_t zvmax = 0;

Float_t Xpos[500], Zpos[500];    


// ======== Thesis end========= //

//! *** Adding MC information ***

Float_t pMC[500],pxMC[500],pyMC[500],pzMC[500],tetaMC[500],energyMC[500],pdgCode[500];
Float_t xMC[500],yMC[500],zMC[500];
#endif
