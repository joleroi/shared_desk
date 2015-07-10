/*
authors Johannes King  king.j@mpi-hd.mpg.de
        Johannes.Veh@fau.de
some code is taken from Vincents Scripts in calibrationmakers/scripts

This is a preliminary run selection script for HESS2
 
    Input is a runlist HAP style

  
*/

// STL
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>


// ROOT
#include <TFile.h>
#include <TArrayD.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TFitResultPtr.h>

// HESS
#include <sash/HESSArray.hh>
#include <sash/DataSet.hh>
#include <hdcalibration/TelescopeParticipation.hh>
#include <sash/MakerChain.hh>
#include <sash/TreeWriter.hh>
#include <sashfile/HandlerC.hh>
#include <sashfile/Utils.hh>
#include <summary/ParticipationFractionMaker.hh>
#include <utilities/TextStyle.hh>

// for nice debug output
#define DEBUG 1 
// 0 no output 
// 1  output
#include <utilities/debugging.hh>



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// We dont have optimised cut atm so you should use this values for your cuts

// @JK kannst du die cuts comentieren? 
// gibt es einen grund das du double nutzt? wenn es geht immer Float da es performanter ist
// aufjedenfall sollten wir es einheitlich machen wenn du nix dagegen hast aender ich alles auf float
Float_t mean_cut = 2;
Float_t rms_cut = 0.5;
Float_t garbage_cut = 62;
Float_t mean_fit_cut = 0.6;
Float_t outliers_cut = 225;










//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TH1D* make_ppf_histogram(int runnumber);
TObjArray* get_list_of_histograms(std::vector<int> runlist);
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

int test_participation_fraction(std::string RunlistName){
  
  /*
    Input is a runlist HAP style
  */
  
  //  gStyle->SetOptStat(0);

  // read the run list from the file and store run numbers and telescope lists
  // in the runlist and tellist vectors
  std::vector<int> runlist;
  std::vector<std::string> tellist;
  int success = SashFile::ReadRunList( RunlistName, runlist, tellist );
  
  std::cout << "read " << success << " runs." << std::endl;
  
  if( !success )
    {
      std::cout << "cant read runlist file: " << RunlistName << std::endl;
    }

  TObjArray * histos = get_list_of_histograms(runlist);
  std::string outRunlist = RunlistName.substr(0,RunlistName.find_last_of("."))+"_passedPPFcuts.lis";
 
  ofstream saveFile;
  saveFile.open(outRunlist.c_str());

  TF1 *fit;
  Double_t mean, rms, garbage, mean_fit, sigma, outliers, bin1, bin2, overflow;
  TH1D * histo;
  TH1D * min_bin      = new TH1D("Garbage_Bin","",1000,0,1000);  //run has 94352 n_min = 468!
  TH2D * rms_vs_mean  = new TH2D("rms_vs_mean","",200,0,2,300,-1,2);
  TH1D * outliershist = new TH1D("outlierst","",2500,0,2500); 
  TH1D * rmshist      = new TH1D("rms","",200,0,2); 
  TH1D * meanhist     = new TH1D("mean","",300,-1,2); 
  
  for (int i=0;i<histos->GetEntriesFast();i++){
    histo = (TH1D*)histos->At(i);
    cout << histo->GetName();

    mean=histo->GetMean();
    meanhist->Fill(mean);
    rms=histo->GetRMS();
    rmshist->Fill(rms);
    rms_vs_mean->Fill(rms,mean);
    
    garbage=histo->GetBinContent(1);  //Get Content of Garbage Bin
    min_bin->Fill(garbage);

    overflow=histo->GetBinContent(histo->GetNbinsX()+1);  //Get Content of Garbage Bin
    if(overflow != 0){
      WARN_OUT   << "Overflow Bin in PPF distribution not empty -> Check Binning" << endl;
      return 1;
    }//Code not yet able to handle this
    
    if(mean > mean_cut) {
      cout << "\t Cut on Mean (>"<<mean_cut<<")" << endl;
      saveFile << runlist[i] << " " << 0 << " #failed cut on Mean(>"<< mean_cut<<")\n";
      continue;
    }
    if(rms > rms_cut) {
      cout << "\t Cut on RMS (>"<<rms_cut<<")" << endl;
      saveFile << runlist[i] << " " << 0 << " #failed cut on RMS (>"<<rms_cut<<")\n";
      continue;
    }
    histo->GetXaxis()->SetRangeUser(-.5,.5);
    histo->Fit("gaus","0Q"); //do not store graphics
    fit  = histo->GetFunction("gaus");
    mean_fit = fit->GetParameter(1);

    if(mean_fit < mean_fit_cut){
      cout << "\t Mean of fit: " << mean_fit << "(<"<<mean_fit_cut<<"), Mean: " << mean << endl;
      saveFile << runlist[i] << " " << 0 << " #failed cut mean_fit (<"<<mean_fit_cut<<")\n";
      continue;
    }
    if(garbage > garbage_cut){ // = 1/2 bus = 4 drawers
      cout << "\t Garbage bin: " << garbage << endl;
      saveFile << runlist[i] << " " << 0 << " #failed cut on missing pixel (>"<<garbage_cut<<")\n";
      continue;
     }
    sigma = fit->GetParameter(2);
    bin1=histo->GetXaxis()->FindBin(mean_fit-2*sigma);           
    bin2=histo->GetXaxis()->FindBin(mean_fit+2*sigma);          
    outliers=histo->Integral(0,bin1)+histo->Integral(bin2,2);
    outliershist->Fill(outliers);
    if(outliers > outliers_cut){
      cout << "\t Outliers    " << outliers << endl;
      saveFile << runlist[i] << " " << 0 << " #failed cut on outliers (>"<<outliers_cut<<")\n";
      continue;
    }
    // if(rms < 0.2){
    //   cout << "\t Cut on RMS (>0.2)" << endl;
    //   continue;
    // }
    cout << "\t ok" << endl;
    saveFile << runlist[i] << " " << tellist[i] << "\n";
  }

  saveFile.close();
  std::cout << "Your Runlist will be saved to: " << outRunlist<< std::endl;
  //Diagnostic plots
  TCanvas *canvas = new TCanvas("Default","Diagnostics Plots",1400,1000);
  canvas->Divide(3,2);

  canvas->cd(1);
  min_bin->GetXaxis()->SetTitle("Content of Garbage Bin");
  min_bin->GetYaxis()->SetTitle("# entries");
  min_bin->Draw();

  canvas->cd(2);
  rms_vs_mean->SetMarkerStyle(8);
  rms_vs_mean->GetXaxis()->SetTitle("RMS");
  rms_vs_mean->GetYaxis()->SetTitle("Mean");
  rms_vs_mean->Draw();
  
  canvas->cd(3);
  //  outliershist->Fit("gaus");
  outliershist->Draw();
  
  canvas->cd(4);
  //canvas2.SetLogy();
  //rmshist->Fit("gaus");
  rmshist->Draw();

  canvas->cd(5);
  meanhist->Draw();
  
  return 1;

}
void MakeParticipationFraction(Int_t nrun, Bool_t save=true) {

  std::string analysisname = "ParticipationFraction";

  int FirstEvent = 0;
  int EventLimit = -1;

  SashFile::HandlerC han("DST");
  std::string dstfilename = SashFile::MakeFileName(nrun,"DST",SashFile::GetDSTPath(),1,false);
  if ( !han.ConnectFileByName( dstfilename ) ) {
    std::cout << Utilities::TextStyle::Red() << "Can't open DST file [" << dstfilename << "]" << Utilities::TextStyle::Reset() << std::endl;
    return;
  }
  
  Sash::DataSet *ds_events = han.GetPrimaryDataSet();
  if (!ds_events) {
    std::cout << Utilities::TextStyle::Red() << "Can't find the DataSet [DST] in the file [" << dstfilename << "]" << Utilities::TextStyle::Reset() << std::endl;
    return;
  }
    
  std::ostringstream oss_outfilename;
  oss_outfilename << "run_" << nrun << "_" << analysisname << "_002.root"; // to have the same name as in HESSDATA.
  
  TFile *outfile =0;
  if (save) {
    outfile = new TFile(oss_outfilename.str().c_str(),"RECREATE");
    outfile->SetCompressionLevel(9);
    gROOT->cd();
  }
  
  std::cout << "Building maker chain for analysis type [" << analysisname << "]" << std::endl;
  
  // Construct the maker chain
  Sash::MakerChain *fChain = new Sash::MakerChain("Calibration Chain",true,false);
  
  //  Summary::ParticipationFractionMaker *spf = new Summary::ParticipationFractionMaker("^DST$",analysisname.c_str(),"Clean0407",5.,true);
  Summary::ParticipationFractionMaker *spf = new Summary::ParticipationFractionMaker("^DST$",analysisname.c_str(),"Extended0407",5.,true);
  
  fChain->UseMaker( spf );
  
  // A TreeWriter to Save Data :
  Sash::TreeWriter *fWriter = new Sash::TreeWriter("^run$","PedestalAnalysis",10000000u);
  fWriter->AddInputFolder( analysisname.c_str() );
  fWriter->SetOutDir(outfile);
  if (save) { fChain->UseMaker( fWriter ); }


  // Loop on Events :
  TStopwatch glob_watch;
  glob_watch.Start();
  int nread = ds_events->EventLoop(fChain,FirstEvent,EventLimit,true,true); // last true : send endrun and endanalysis
  int neventstot = (int)ds_events->GetEntries();
  std::cout << "nrun = " << nrun << " nread = " << nread << " neventstot = " << neventstot << std::endl;
  
  if (save) {
    fChain->Finish(Sash::Folder::GetFolder("Finish"));
    delete outfile;
    outfile=0;
    fWriter->SetOutDir(outfile);
  }

  glob_watch.Stop();

  std::cout << "Time to make ParticipationFraction = " << glob_watch.RealTime() << " s (" << glob_watch.CpuTime() << " s)" << std::endl;

  std::cout << "Deleting fChain" << std::endl;
  if (save) delete fChain;

}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


TObjArray* get_list_of_histograms(std::vector<int> runlist){


 int counter = 0;
 TObjArray * histos = new TObjArray;

  for( unsigned int i=0; i<runlist.size(); ++i ) {

      histos->Add(make_ppf_histogram(runlist[i]));
      counter++;
    }
  std::cout << counter << " files read" << std::endl;
  histos->SaveAs("ppf/ppf_histograms.root");

  return histos;
}


TH1D* make_ppf_histogram(int runnumber){

  //convert runnumber into ppt file name
  // run_94354_ParticipationFraction_002.root
  std::ostringstream fileName;
  fileName << "run_" << runnumber << "_ParticipationFraction_002.root";

  std::ostringstream histname;
  histname<<  "h_ppd_" << runnumber;
  //  std::cout << "histname: " << histname.str().c_str() << std::endl;
  //  TH1D* hist = new TH1D(histname.c_str(),"Pixel Participation Distribution",150,-1,2);
  TH1D* hist = new TH1D(histname.str().c_str(),"Pixel Participation Distribution",6000,-1,5);


  if(!SashFile::FileExists(fileName.str().c_str()))
    {
     MakeParticipationFraction(runnumber,true);
    }

  // this is a dirty hack to make sure we dont crash in case there is no DST file/PPF file
  // I dont like to connect to files that often it makes us slow
  // maybe we can implement a return value to MakeParticipationFraction ad check for that value
  // also It would be nice to have a mode in MPF that doesnt open the result canvas slow again ;)
  if(!SashFile::FileExists(fileName.str().c_str()))
    {
      WARN_OUT << "Could not find the DST of ppf file for run: " << runnumber
		   << endl;
      return hist;  
    }


  TFile* f = new TFile(fileName.str().c_str(),"READ"); 

  Sash::DataSet* data =  dynamic_cast<Sash::DataSet*>(f->Get("ParticipationFraction"));
  if( !data ) 
    return NULL;
  else{
    data->GetEntry(0);
    //std::cout << "Loading Participation Fraction Data Set" << std::endl;
    //data->Print();
  }

  Sash::DataSet* run =  dynamic_cast<Sash::DataSet*>(f->Get("run"));
  if( !run ) 
    return NULL;
  else{
    run->GetEntry(0);
    //std::cout << "Loading run Data Set" << std::endl;
    //run->Print();
  }


  Sash::HESSArray* hess = &Sash::HESSArray::GetHESSArray();  
  Sash::Telescope* tel = hess->GetTel(5);
  
  const HDCalibration::TelescopeParticipation *ppf = tel->Get<HDCalibration::TelescopeParticipation>();
 
  TArrayD* pix_values = ppf->Flatten("Fraction");
  //std::cout << "#entries in camera: " << pix_values->GetSize() << std::endl;
   
  //  int pos = fileName.find("run");
  // std::string runnumber = fileName.substr(pos+4,5);
 

  for (int i=0;i<pix_values->GetSize();i++){
    hist->Fill(pix_values->GetAt(i));
  }
  return hist;
  
}


