#include <time.h>

void checkRun(const char* filename, const Int_t run){
  TFile *cosm  = new TFile(filename,"READ");
  TTree *tcosm = (TTree*)cosm->Get("fDatree");
//  UInt_t timestp, startTime, lastTime;
  UInt_t timestp;
  time_t startTime, lastTime;
  Int_t nRun,nEve;
  char chnkName[50];
  struct tm lt;
  char res[32];

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
      startTime = 2000000000;
       lastTime = 0;
      totChnks = 0;
    }

    if (currChnk[0] == '\0' && nRun == run) {
//      printf("Chunk: %s\n",chnkName);
      strcpy(currChnk,chnkName);
      totChnks++;
    }

    if (strcmp(currChnk,chnkName) != 0 && nRun == run) {
//      printf("\nChunk: %s\n",chnkName);
      strcpy(currChnk,chnkName);
      totChnks++;
    }

    if (nRun == run) {
      if (timestp < startTime) startTime = timestp;
      if (timestp >  lastTime)  lastTime = timestp;
      totEv++;
    }
  }

  printf("Run %d : Events %d Chunks %d Time %u (first %u last %u)\n",
	 run, totEv, totChnks, lastTime-startTime, startTime, lastTime);
//  localtime_r(&startTime, &lt);
  lt = localtime(&startTime);
  strftime(res, sizeof(res), "%a %b %d %Y, %H:%M:%S", &lt);
  printf("Start time: %s\n",res);
//  localtime_r(&lastTime, &lt);
  lt = localtime(&lastTime);
  strftime(res, sizeof(res), "%a %b %d %Y, %H:%M:%S", &lt);
  printf("Last time: %s\n",res);

}
