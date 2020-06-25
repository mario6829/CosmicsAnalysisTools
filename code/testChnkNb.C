void testChnkNb(const char* filename, const Int_t runnum){
  TFile *cosm  = new TFile(filename,"READ");
  TTree *tcosm = (TTree*)cosm->Get("fDatree");
  Int_t nRun,nEve;
  char chnkName[50];

  Int_t nESDentries = (Int_t)tcosm->GetEntries();

  tcosm->SetBranchAddress("nRun",&nRun);
  tcosm->SetBranchAddress("nEve",&nEve);
  tcosm->SetBranchAddress("chnk",chnkName);

  cout << "Total Number of events in ESD =  " << nESDentries << endl << endl;
 
  Int_t totEv=0, nChunks=0;
  char currChnk[50];
  currChnk[0] = '\0';
  for(Int_t iEntry=0;iEntry<nESDentries;iEntry++){
    tcosm->GetEntry(iEntry);

    if (nRun == runnum) {
      if (currChnk[0] == '\0') {
	strcpy(currChnk,chnkName);
	printf("Chunk: %s\n",chnkName);
	nChunks++;
      }

      if (strcmp(currChnk,chnkName) != 0) {
	printf("Run %d Events %d\n",nRun,totEv);
	printf("\nChunk: %s\n",chnkName);
	strcpy(currChnk,chnkName);
	totEv = 0;
	nChunks++;
      }
    }

    totEv++;
  }
  if (nRun == runnum)
    printf("Run %d Events %d\n",nRun,totEv);

    printf("Summary for Run %d : %d chunks found\n",runnum,nChunks);
}
