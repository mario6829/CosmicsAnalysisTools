AliAnalysisGrid* CreateAlienHandler(Int_t nRun,
				    char todo='F', char* period="LHC15a")
{
// Check if user has a valid token, otherwise make one. This has limitations.
// One can always follow the standard procedure of calling alien-token-init then
//   source /tmp/gclient_env_$UID in the current shell.
//   if (!AliAnalysisGrid::CreateToken()) return NULL;
   AliAnalysisAlien *plugin = new AliAnalysisAlien();
   plugin->SetOverwriteMode();

// Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
   switch (todo) {
     case 'S':
       plugin->SetRunMode("submit");
       break;
 //    case 'F':
 //  plugin->SetRunMode("test");
 //    break;
     case 'T':
       plugin->SetRunMode("terminate");
       break;
     case 'F':
     default:
       plugin->SetRunMode("full");
       break;
   }

// Set versions of used packages
   plugin->SetAPIVersion("V1.1x");
//   plugin->SetROOTVersion("v5-30-06-1");
//   plugin->SetROOTVersion("v5-34-08-6");  //v5-34-08-7 v5-34-02-1   
//   plugin->SetROOTVersion("v5-34-08-7");  //v5-34-08-7 v5-34-02-1   
//   plugin->SetROOTVersion("v5-34-30-1");
//   plugin->SetROOTVersion("v5-34-30-alice-7");
//   plugin->SetROOTVersion("v5-34-30-alice10-1");
   plugin->SetROOTVersion("v5-34-30-alice10-12");
//   plugin->SetAliROOTVersion("vAN-20150118");  //v5-05-Rev-28 v5-04-01-AN
//   plugin->SetAliROOTVersion("v5-06-09");
//   plugin->SetAliROOTVersion("v5-06-08");
//   plugin->SetAliROOTVersion("v5-05-Rev-28");  //v5-05-Rev-28 v5-04-01-AN
//   plugin->SetAliROOTVersion("v5-03-09-AN");
//   plugin->SetAliROOTVersion("v5-06-16");
//   plugin->SetAliROOTVersion("v5-06-19");
//   plugin->SetAliROOTVersion("v5-06-33");
//   plugin->SetAliROOTVersion("v5-06-39");
//   plugin->SetAliROOTVersion("v5-07-18-2");
//   plugin->SetAliROOTVersion("v5-08-00-1");
//   plugin->SetAliROOTVersion("v5-09-24-1");
   plugin->SetAliROOTVersion("v5-09-45-1");

// Declare input data to be processed.
// Method 1: Create automatically XML collections using alien 'find' command.
// Define production directory LFN
   plugin->SetCheckCopy(kFALSE);
   char dataDir[120];
   sprintf(dataDir,"/alice/data/2018/%s",period);
   plugin->SetGridDataDir(dataDir); // /alice/data/2015/LHC15a
// Set data search pattern
//   plugin->SetDataPattern("*ESDs.root");
//   plugin->SetDataPattern("*pass1/*ESDs.root");
   plugin->SetDataPattern("*cosmics_pass1/*ESDs.root");
//   plugin->SetDataPattern("*calo_pass1/*ESDs.root");
//   plugin->SetDataPattern("*cosmics_pass2/*ESDs.root");
//   plugin->SetDataPattern("*calo_pass2/*ESDs.root");
//   plugin->SetDataPattern("*cosmics/*ESDs.root");
   plugin->SetRunPrefix("000");   // real data
// ...then add run numbers to be considered
 //  plugin->AddRunNumber(125020);
   plugin->AddRunNumber(nRun);  // real data
//   plugin->SetOutputSingleFolder("output");
//   plugin->SetOutputToRunNo();

// Method 2: Declare existing data files (raw collections, xml collections, root file)
// If no path mentioned data is supposed to be in the work directory (see SetGridWorkingDir())
// XML collections added via this method can be combined with the first method if
// the content is compatible (using or not tags)
//   plugin->AddDataFile("tag.xml");
//   plugin->AddDataFile("/alice/data/2008/LHC08c/000057657/raw/Run57657.Merged.RAW.tag.root");
// Define alien work directory where all files will be copied. Relative to alien $HOME.

   char txt[120];
   char txt1[120];
   char txt2[120];
   char txt3[120];
   char txt4[120];
   sprintf(txt,"anacosmic/%s/%d",period,nRun);
   sprintf(txt1,"anaCosmic_%d.C",nRun);
   sprintf(txt2,"anaCosmic_%d.sh",nRun);
   sprintf(txt3,"anaCosmic_%d.jdl",nRun);
   sprintf(txt4,"cosmicAnaRun_%d.root",nRun);

   plugin->SetGridWorkingDir(txt);
 
// Declare alien output directory. Relative to working directory.
   plugin->SetGridOutputDir("output"); // In this case will be $HOME/work/output
// Declare the analysis source files names separated by blancs. To be compiled runtime
// using ACLiC on the worker nodes.
   plugin->SetAnalysisSource("AliCosmics.cxx");

// Declare all libraries (other than the default ones for the framework. These will be
// loaded by the generated analysis macro. Add all extra files (task .cxx/.h) here.
   plugin->SetAdditionalLibs("AliCosmics.h AliCosmics.cxx");

// Declare the output file names separated by blancs.
// (can be like: file.root or file.root@ALICE::Niham::File)
  //  plugin->SetDefaultOutputs();
   plugin->SetDefaultOutputs(kFALSE);  // explicitly call to manually set output files
   plugin->SetOutputFiles(txt4);  // cosmicAnaRun_<runnum>.root

// Keep the log files (otherwise they get deleted by the validation script)
   plugin->SetKeepLogs(kTRUE);

// Optionally define the files to be archived.
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@ALICE::NIHAM::File root_archive.zip:*.root@ALICE::NIHAM::File");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@disk=1 root_archive.zip:cosmicAnaRun*.root@disk=1 events_archive.zip:*.out@disk=1");
   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@disk=1 events_archive.zip:*.out@disk=2");
//   plugin->SetOutputArchive("log_archive.zip:stdout,stderr@disk=1 root_archive.zip:cosmicAnaRun*.root@ALICE::Bari::SE events_archive.zip:*.out@disk=1");

// Optionally set a name for the generated analysis macro (default MyAnalysis.C)
   plugin->SetAnalysisMacro(txt1);

// Optionally set maximum number of input files/subjob (default 100, put 0 to ignore)
   plugin->SetSplitMaxInputFileNumber(50);
   //plugin->SetSplitMaxInputFileNumber(0);

// Optionally modify the executable name (default analysis.sh)
   plugin->SetExecutable(txt2);

// Optionally set number of failed jobs that will trigger killing waiting sub-jobs.
   plugin->SetMaxInitFailed(5);

// Optionally resubmit threshold.
   plugin->SetMasterResubmitThreshold(90);

// Optionally set time to live (default 30000 sec)
   plugin->SetTTL(86300);

// Optionally set input format (default xml-single)
   plugin->SetInputFormat("xml-single");

// Optionally modify the name of the generated JDL (default analysis.jdl)
   plugin->SetJDLName(txt3);

// Optionally modify job price (default 1)
   plugin->SetPrice(1);      

// Optionally modify split mode (default 'se')    
   plugin->SetSplitMode("se");

// Optionally do not wait for terminate (i.e. do no start aliensh)
   plugin->SetDropToShell(kFALSE);

   return plugin;
}
