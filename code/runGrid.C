void runGrid(Int_t nRun, char *todo="F", char* period="LHC15a"){
// load libs needed to run from root 
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase.so");
  gSystem->Load("libESD.so");
  gSystem->Load("libAOD.so");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libOADB");
  gSystem->Load("libANALYSISalice");
  gSystem->Load("libXMLParser.so");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");

  // load private class for cosmic analysis

  gROOT->LoadMacro("AliCosmics.cxx++g");

  // Create the analysis manager
  AliAnalysisManager *mgr = new AliAnalysisManager("testAnalysis");

  // Create and configure the alien handler plugin
  gROOT->LoadMacro("CreateAlienHandler.C");
  AliAnalysisGrid *alienHandler = CreateAlienHandler(nRun, todo[0], period);
  if (!alienHandler) return;

  // Connect plug-in to the analysis manager
  mgr->SetGridHandler(alienHandler);
  char outputFileName[120];
  sprintf(outputFileName,"cosmicAnaRun_%d.root",nRun);

  // input handler
  AliVEventHandler* esdH = new AliESDInputHandler();
//  ((AliESDInputHandler*)esdH->GetInputEventHandler())->SetNeedField();
  mgr->SetInputEventHandler(esdH);
  
  // ==> Cosmic analysis task
  AliAnalysisTask *task = new AliCosmics("TaskAnaCosmicTO");
  mgr->AddTask((AliAnalysisTaskSE *) task);

  AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
  AliAnalysisDataContainer *coutput = mgr->CreateContainer("fDatree",TTree::Class(),AliAnalysisManager::kOutputContainer,outputFileName);
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer("fListHist",TList::Class(),AliAnalysisManager::kOutputContainer,outputFileName);

  // Connect input/output
  mgr->ConnectInput(task, 0, cinput);
  mgr->ConnectOutput(task, 0, coutput);
  mgr->ConnectOutput(task, 1, coutput1);

  // Enable debug printouts
  mgr->SetDebugLevel(0);

  if (!mgr->InitAnalysis())
    return;

  mgr->PrintStatus();

  // Start analysis in grid.
  mgr->StartAnalysis("grid");

  // Save the link run number <-> masterjob Id
  TString jobnum = ((AliAnalysisAlien*)alienHandler)->GetGridJobIDs();
  FILE* jobtab = fopen("jobs.tmp","a");
  fprintf(jobtab,"%d  %s\n",nRun,jobnum.Data());
  fclose(jobtab);
};
