void countChnks(const char* filename){
  TFile *cosm  = new TFile(filename,"READ");
  TTree *tcosm = (TTree*)cosm->Get("fDatree");
  UInt_t timestp, startTime, lastTime;
  Int_t nRun,nEve;
  char chnkName[50];

  Int_t nESDentries = (Int_t)tcosm->GetEntries();

  tcosm->SetBranchAddress("nRun",&nRun);
  tcosm->SetBranchAddress("nEve",&nEve);
  tcosm->SetBranchAddress("chnk",chnkName);
  tcosm->SetBranchAddress("timestp",&timestp);

  cout << "Total Number of events in ESD =  " << nESDentries << endl << endl;
 
  Int_t totEv=0, currRun=0, totChnks;
  char currChnk[50];
  currChnk[0] = '\0';
  for(Int_t iEntry=0;iEntry<nESDentries;iEntry++){
    tcosm->GetEntry(iEntry);

    if (currRun == 0) {
      currRun = nRun;
      startTime = timestp;
       lastTime = timestp;
      totChnks = 0;
    }

    if (currRun != nRun) {
      printf("Run %d : Events %d Chunks %d Time %u\n",currRun, totEv, totChnks,
	     lastTime-startTime);
      currRun = nRun;
      startTime = timestp;
       lastTime = timestp;
      totChnks = 0;
      totEv = 0;
    }

    if (currChnk[0] == '\0') {
      strcpy(currChnk,chnkName);
//      printf("Chunk: %s\n",chnkName);
      totChnks++;
    }

    if (strcmp(currChnk,chnkName) != 0) {
//      printf("\nChunk: %s\n",chnkName);
      strcpy(currChnk,chnkName);
      totChnks++;
    }

    if (timestp < startTime) startTime = timestp;
    if (timestp >  lastTime)  lastTime = timestp;
    totEv++;
  }

  printf("Run %d : Events %d Chunks %d Time %u\n",currRun, totEv, totChnks,
	 lastTime-startTime);
}
