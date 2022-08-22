#include "TH1.h"
#include "TCanvas.h"
#include "TNtuple.h"
#include "TPad.h"
#include "TTree.h"
#include "TMath.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "TGraphErrors.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TROOT.h"
#include <vector>
#include "unistd.h"
#include "TSystem.h"
#include "TBranch.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TH2.h"
#include "TCut.h"
#include "Riostream.h"
#include "TNtupleD.h"
#include "TCut.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

#include "Fiducial.h"
#include "Acceptance.h"
#include "Hconstants.h"
#include "Hfunctions.h"


int main(){
  
  const int maxParticles=50;
   TFile *fout = new TFile("output/out_MF_scratch.root","RECREATE");
	gRandom->SetSeed(4357);

   TCanvas *canvas = new TCanvas("canvas","output plots", 700, 700);
   canvas->SetFillColor(0);
   canvas->SetBorderMode(0);
   canvas->SetBorderSize(0);
   canvas->SetFrameFillColor(0);
   canvas->SetFrameBorderMode(0);
   //canvas->SetLogy();
   std::string pdf_file_name = "output/out_MF_scratch.pdf";
   gROOT->SetBatch(true);
   gStyle->SetOptStat(0);
   
   TFile *f_3p = new TFile("../../../data/19jun20/skim_He3_p_Jun20.root");
   TFile *f_3n = new TFile("../../../data/19jun20/skim_He3_n_Mar5.root");
   TFile *f_4p = new TFile("../../../data/19jun20/skim_He4_p_Jun20.root");
   TFile *f_4n = new TFile("../../../data/19jun20/skim_He4_n_Mar5.root");
   TFile *f_5p = new TFile("../../../data/19jun20/skim_C12_p_Jun20.root");
   TFile *f_5n = new TFile("../../../data/19jun20/skim_C12_n_Mar5.root");
  
  
  
  //////////////////////////////////
  // Plots/////////////////////////
  /////////////////////////////////
  TH1F* h_emom[3][3]; 
  TH1F* h_etheta[3][3];
  TH1F* h_q2[3][3];
  TH1F* h_nu[3][3];
  TH1F* h_xb[3][3];
  TH1F* h_pmiss[3][3];
  TH2F* h_emVpm[3][3];
  TH1F* h_thetapq[3][3];
  TH1F* h_nmom[3][3];
  TH1F* h_ntheta[3][3];
  TH1D * h_emiss[3][3]; 
  TH1D * h_mmiss[3][3];


  //target: 3He, 4He, 12C
  for (int ii=0;ii<3;ii++){
    //lead: proton, smeared proton, neutron
    for(int jj=0;jj<3;jj++){
      h_emom[ii][jj]=new TH1F(Form("h_emom_%i_%i",ii,jj),";P_{e-} [GeV/c] (with lead cuts);",50,0,6); 
      h_etheta[ii][jj]=new TH1F(Form("h_etheta_%i_%i",ii,jj),";#theta_{e-} [deg] (with lead cuts);",50,0,45);
      h_q2[ii][jj]=new TH1F(Form("h_q2_%i_%i",ii,jj),";Q^{2} (with lead cuts);",50,0,6);
      h_nu[ii][jj]=new TH1F(Form("h_nu_%i_%i",ii,jj),";#nu [GeV] (with lead cuts);",50,0,3);
      h_xb[ii][jj]=new TH1F(Form("h_xb_%i_%i",ii,jj),";x_{B} (with lead cuts);",50,0,2);
      h_pmiss[ii][jj]=new TH1F(Form("h_pmiss_%i_%i",ii,jj),";P_{miss} [GeV/c];",50,0,1.5);
      h_emVpm[ii][jj]=new TH2F(Form("h_emVpm_%i_%i",ii,jj),";pmiss;emiss",50,0,0.8,50,-0.5,0.5);
      h_thetapq[ii][jj]=new TH1F(Form("h_thetapq_%i_%i",ii,jj),";#theta_{lead N/q} [deg];",50,0,25);
      h_nmom[ii][jj]=new TH1F(Form("h_nmom_%i_%i",ii,jj),";P_{lead N} [GeV/c];",50,0,3);
      h_ntheta[ii][jj]=new TH1F(Form("h_ntheta_%i_%i",ii,jj),";#theta_{lead N} [deg];",50,0,90);
      h_emiss[ii][jj]= new TH1D(Form("h_emiss_%i_%i",ii,jj),";missing energy",50,-1,2);
      h_mmiss[ii][jj]= new TH1D(Form("h_mmiss_%i_%i",ii,jj),";missing mass [GeV/c^{2}]",50,0,3);
    }
  }
  
 
  //////////////////////////////////
  // Set some target dep parameters
  /////////////////////////////////
  double tgt_min, tgt_max, tgt_M, tgt_R;
  int tgt_Z, tgt_A;

  ////////////////////////////////
  //Helium-3
  ////////////////////////////////
  TString target = "3He";
  Acceptance myMapHe3("He3", 4461, 2250, "p");

  tgt_min = vtx_min_He3; tgt_max = vtx_max_He3; tgt_Z = 2; tgt_A = 3; tgt_M = M3He; tgt_R = M2H;
  Fiducial fid_params(4461,2250,5996,"3He",true);

  //open the 3He proton file
  TTree * intree = (TTree*)f_3p->Get("T");
  int nParticles, type[maxParticles];
  double QSq, xB, mom_x[maxParticles], mom_y[maxParticles], mom_z[maxParticles], Nu, ec_x[maxParticles], ec_y[maxParticles], vtxZ[maxParticles], ec_u[maxParticles], ec_v[maxParticles], ec_w[maxParticles], sc_time[maxParticles], ec_time[maxParticles], sc_path[maxParticles], ec_path[maxParticles],ec_in[maxParticles], ec_out[maxParticles];
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  double theta_prev = 0;

    //get (e,e'p)
    for (int event =0 ; event < intree->GetEntries() ; event++)
    //    for (int event =0 ; event < 100000 ; event++)
      {
         if (event%10000==0)
        	cout << "Working on event " << event <<endl;

        intree->GetEvent(event);
        
        //if (nRun!=18020){continue;}
        
        if (QSq<Q2_min){continue;}

         // Get the electron and q vectors
        TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
        double Ebeam = Nu + mom_e.Mag();
        TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
        
        //mean field stuff:
        double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mp);

        // Loop over particles looking for a proton most likely to pass lead cuts
        int leadingID = -1;
        for (int part=1 ; part < nParticles ; part++)
  	{   
  	  // Only look at protons
  	  if (type[part] != 2212)
  	    continue;

  	  //vertex cut-target dep.
  	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
  	  
  	  // Get a 3-vector for this proton's momentum
  	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

  	  // Angle between p and q
  	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
  	  if (theta_pq>thetapq_MFcut){continue;}
  	  
  	  //first proton of interest in event
  	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

  	  /////////////////////////////////////////////////////////////////////
  	  //Here I selected the proton with angle closest to q-vector:	  
  	  //current p/q is larger and angle is smaller than previous
  	  if (theta_pq<theta_prev){leadingID=part;}
  	  /////////////////////////////////////////////////////////////////////	  
  	 
  	}

        //Now we have a leading proton!
        if (leadingID > 0)
  	{
  	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
  	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

  	 // if (mom_lead.Mag()>2.34){continue;}
 	  //if (mom_lead.Mag()>1.95){continue;}
  
  	  TVector3 mom_miss = mom_lead - q;
  	  double theta = mom_lead.Theta()*180./M_PI;
  	  double phi = mom_lead.Phi();
  	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());
  	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
  	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);

  	  //smear the proton momentum resolution
  	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
  	  TVector3 u1 = mom_lead.Unit();
  	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

  	  //recalculate the smeared factors
  	  TVector3 mom_missSmear = mom_smear - q;
  	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
  	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
  	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
  	  double p_over_q_sm = mom_smear.Mag()/q.Mag();
  	  
  	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
  	  phi *= 180./M_PI; 

  	  //apply proton fiducial cuts:

  	  //apply EC 10cm edge projection cut:
  	  //this doesn't work, but I don't know why:
  	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
  	  // Cutting top of EC
  	  if(theta>44.) continue;
  	  
  	  //passing lead proton cuts
  	  double p_over_q = mom_smear.Mag() / q.Mag();

  	  //check if in good acceptance
  	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
  	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}
  	  
		 double p_eff = myMapHe3.get_acc(mom_lead);

  	    
  	  //select the 3He proton sample
  	  	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut && mom_miss.Mag()<pmiss_MFnomCut && emiss<emiss_nomCut
  	  			  &&mom_lead.Mag()<1.95 && p_points_to_EC_fid(mom_lead,vtxZ[0],10.)&&fid_params.pFiducialCut(mom_lead)){
  	  	    h_emom[0][0]->Fill(mom_e.Mag()); 
  	  	    h_etheta[0][0]->Fill(mom_e.Theta()*180./M_PI);
  	  	    h_q2[0][0]->Fill(QSq);
  	  	    h_nu[0][0]->Fill(Nu);
  	  	    h_xb[0][0]->Fill(xB);
  	  	    h_pmiss[0][0]->Fill(mom_miss.Mag());
  	  	    h_emVpm[0][0]->Fill(mom_miss.Mag(),emiss);
  	  	    h_thetapq[0][0]->Fill(mom_lead.Angle(q)*180./M_PI);
  	  	    h_nmom[0][0]->Fill(mom_lead.Mag());
  	  	    h_ntheta[0][0]->Fill(theta);
  	  	    h_emiss[0][0]->Fill(emiss);
  	  	    h_mmiss[0][0]->Fill(mmiss);	    
  	  	  }

  	  	  //select the 3He smeared proton sample
  	  	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_smear.Angle(q)*180./M_PI<thetapq_MFcut && mom_missSmear.Mag()<pmiss_MFnomCut_n && emissSmear<emiss_nomCut_n
  	  			  &&mom_smear.Mag()<1.95 && p_points_to_EC_fid(mom_smear,vtxZ[0],10.)&&fid_params.pFiducialCut(mom_smear)){
  	  	    h_emom[0][1]->Fill(mom_e.Mag(),1./p_eff); 
  	  	    h_etheta[0][1]->Fill(mom_e.Theta()*180./M_PI,1./p_eff);
  	  	    h_q2[0][1]->Fill(QSq,1./p_eff);
  	  	    h_nu[0][1]->Fill(Nu,1./p_eff);
  	  	    h_xb[0][1]->Fill(xB,1./p_eff);
  	  	    h_pmiss[0][1]->Fill(mom_missSmear.Mag(),1./p_eff);
  	  	    h_emVpm[0][1]->Fill(mom_missSmear.Mag(),emissSmear,1./p_eff);
  	  	    h_thetapq[0][1]->Fill(mom_smear.Angle(q)*180./M_PI,1./p_eff);
  	  	    h_nmom[0][1]->Fill(mom_smear.Mag(),1./p_eff);
  	  	    h_ntheta[0][1]->Fill(theta,1./p_eff);
  	  	    h_emiss[0][1]->Fill(emissSmear,1./p_eff);
  	  	    h_mmiss[0][1]->Fill(mmissSmear,1./p_eff);
  	  	    
  	  	  }

  	}

  	}

  
 
  intree->Delete();
  
  //open the 3He neutron file
  intree = (TTree*)f_3n->Get("T");
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  for (int event =0 ; event < intree->GetEntries() ; event++)
     // for (int event =0 ; event < 100000 ; event++)
     {
      
       if (event%10000==0){cout << "Working on event " << event <<endl;}

       intree->GetEvent(event);
       //if (nRun!=18020){continue;}

       if (QSq<Q2_min){continue;}

        // Get the electron and q vectors
       TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
       double Ebeam = Nu + mom_e.Mag();
       TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
       
       //mean field stuff:
       double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mn);

       // Loop over particles looking for a proton most likely to pass lead cuts
       int leadingID = -1;
       for (int part=1 ; part < nParticles ; part++)
 	{   
 	  // Only look at neutrons
 	  if (type[part] != 2112)
 	    continue;

 	  //vertex cut-target dep.
 	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
 	  
 	  // Get a 3-vector for this proton's momentum
 	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

 	  // Angle between p and q
 	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
 	  if (theta_pq>thetapq_MFcut){continue;}
 	  
 	  //first proton of interest in event
 	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

 	  /////////////////////////////////////////////////////////////////////
 	  //Here I selected the proton with angle closest to q-vector:	  
 	  //current p/q is larger and angle is smaller than previous
 	  if (theta_pq<theta_prev){leadingID=part;}
 	  /////////////////////////////////////////////////////////////////////	  
 	  	  
 	}

       //Now we have a leading neutron!
       if (leadingID > 0)
 	{
 	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
 	  mom_lead = correctPnVec(mom_lead,ec_in[leadingID],ec_out[leadingID]);
 	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

 	  //calculate the nominal values
 	  TVector3 mom_miss = mom_lead - q;
 	  double theta = mom_lead.Theta()*180./M_PI;
 	  double phi = mom_lead.Phi();
 	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
 	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
 	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
 	  
 	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
 	  phi *= 180./M_PI; 
 	  
 	  //apply proton fiducial cuts:
 	  if (!fid_params.pFiducialCut(mom_lead)) {continue;}

 	  //apply EC 10cm edge projection cut:
 	  //this doesn't work, but I don't know why:
 	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
 	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
 	  // Cutting top of EC
 	  if(theta>44.) continue;
 	  if(mom_lead.Mag()>1.95) continue;
 	  
 	  //passing lead proton cuts
 	  double p_over_q = mom_lead.Mag() / q.Mag();

 	  //check if in good acceptance
 	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
 	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}


 	  //select the 3He neutron sample
 	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut 
 			  && mom_miss.Mag()<pmiss_MFnomCut_n && emiss<emiss_nomCut_n
 			  //&& mom_lead.Mag()>1.4 && mom_lead.Mag()<1.8
 			  ){
 		  
 		 	    h_emom[0][2]->Fill(mom_e.Mag()); 
 		 	    h_etheta[0][2]->Fill(mom_e.Theta()*180./M_PI);
 		 	    h_q2[0][2]->Fill(QSq);
 		 	    h_nu[0][2]->Fill(Nu);
 		 	    h_xb[0][2]->Fill(xB);
 		 	    h_pmiss[0][2]->Fill(mom_miss.Mag());
 		 	    h_emVpm[0][2]->Fill(mom_miss.Mag(),emiss);
 		 	    h_thetapq[0][2]->Fill(mom_lead.Angle(q)*180./M_PI);
 		 	    h_nmom[0][2]->Fill(mom_lead.Mag());
 		 	    h_ntheta[0][2]->Fill(theta);
 		 	    h_emiss[0][2]->Fill(emiss);
 		 	    h_mmiss[0][2]->Fill(mmiss);
 		 	    
 		 	  }	  
  
 	}//end leading
     }//end event
    
  
  ////////////////////////////////
  //Helium-4
  ////////////////////////////////
  intree->Delete();
	gRandom->SetSeed(4357);

  target = "4He";
  Acceptance myMapHe4("He4", 4461, 2250, "p");

  tgt_min = vtx_min_He4; tgt_max = vtx_max_He4; tgt_Z = 2; tgt_A = 4; tgt_M = M4He; tgt_R = M3H;
  Fiducial fid_params4He(4461,2250,5996,"4He",true);
  //open the 4He proton file
  intree = (TTree*)f_4p->Get("T");
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);

  //get (e,e'p)
   for (int event =0 ; event < intree->GetEntries() ; event++)
   //    for (int event =0 ; event < 100000 ; event++)
     {
        if (event%10000==0)
       	cout << "Working on event " << event <<endl;

       intree->GetEvent(event);
       
       //if (nRun!=18020){continue;}
       
       if (QSq<Q2_min){continue;}

        // Get the electron and q vectors
       TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
       double Ebeam = Nu + mom_e.Mag();
       TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
       
       //mean field stuff:
       double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mp);

       // Loop over particles looking for a proton most likely to pass lead cuts
       int leadingID = -1;
       for (int part=1 ; part < nParticles ; part++)
 	{   
 	  // Only look at protons
 	  if (type[part] != 2212)
 	    continue;

 	  //vertex cut-target dep.
 	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
 	  
 	  // Get a 3-vector for this proton's momentum
 	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

 	  // Angle between p and q
 	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
 	  if (theta_pq>thetapq_MFcut){continue;}
 	  
 	  //first proton of interest in event
 	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

 	  /////////////////////////////////////////////////////////////////////
 	  //Here I selected the proton with angle closest to q-vector:	  
 	  //current p/q is larger and angle is smaller than previous
 	  if (theta_pq<theta_prev){leadingID=part;}
 	  /////////////////////////////////////////////////////////////////////	  
 	 
 	}

       //Now we have a leading proton!
       if (leadingID > 0)
 	{
 	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
 	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

 	  //if (mom_lead.Mag()>2.34){continue;}
 	  //if (mom_lead.Mag()>1.95){continue;}
  
 	  TVector3 mom_miss = mom_lead - q;
 	  double theta = mom_lead.Theta()*180./M_PI;
 	  double phi = mom_lead.Phi();
 	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());
 	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
 	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);

 	  //smear the proton momentum resolution
 	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
 	  TVector3 u1 = mom_lead.Unit();
 	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

 	  //recalculate the smeared factors
 	  TVector3 mom_missSmear = mom_smear - q;
 	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
 	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
 	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
 	  double p_over_q_sm = mom_smear.Mag()/q.Mag();
 	  
 	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
 	  phi *= 180./M_PI; 

 	  //apply proton fiducial cuts:

 	  //apply EC 10cm edge projection cut:
 	  //this doesn't work, but I don't know why:
 	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
 	  // Cutting top of EC
 	  if(theta>44.) continue;
 	  
 	  //passing lead proton cuts
 	  double p_over_q = mom_smear.Mag() / q.Mag();

 	  //check if in good acceptance
 	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
 	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}
 	  
		 double p_eff = myMapHe4.get_acc(mom_lead);

 	  //select the 4He proton sample
 	  	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut && mom_miss.Mag()<pmiss_MFnomCut && emiss<emiss_nomCut
 	  			  &&mom_lead.Mag()<1.95 && p_points_to_EC_fid(mom_lead,vtxZ[0],10.)&&fid_params4He.pFiducialCut(mom_lead)){
 	  	    h_emom[1][0]->Fill(mom_e.Mag()); 
 	  	    h_etheta[1][0]->Fill(mom_e.Theta()*180./M_PI);
 	  	    h_q2[1][0]->Fill(QSq);
 	  	    h_nu[1][0]->Fill(Nu);
 	  	    h_xb[1][0]->Fill(xB);
 	  	    h_pmiss[1][0]->Fill(mom_miss.Mag());
 	  	    h_emVpm[1][0]->Fill(mom_miss.Mag(),emiss);
 	  	    h_thetapq[1][0]->Fill(mom_lead.Angle(q)*180./M_PI);
 	  	    h_nmom[1][0]->Fill(mom_lead.Mag());
 	  	    h_ntheta[1][0]->Fill(theta);
 	  	    h_emiss[1][0]->Fill(emiss);
 	  	    h_mmiss[1][0]->Fill(mmiss);	    
 	  	  }

 	  	  //select the 4He smeared proton sample
 	  	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_smear.Angle(q)*180./M_PI<thetapq_MFcut && mom_missSmear.Mag()<pmiss_MFnomCut_n && emissSmear<emiss_nomCut_n
 	  			  &&mom_smear.Mag()<1.95 && p_points_to_EC_fid(mom_smear,vtxZ[0],10.)&&fid_params4He.pFiducialCut(mom_smear)){
 	  	    h_emom[1][1]->Fill(mom_e.Mag(),1./p_eff); 
 	  	    h_etheta[1][1]->Fill(mom_e.Theta()*180./M_PI,1./p_eff);
 	  	    h_q2[1][1]->Fill(QSq,1./p_eff);
 	  	    h_nu[1][1]->Fill(Nu,1./p_eff);
 	  	    h_xb[1][1]->Fill(xB,1./p_eff);
 	  	    h_pmiss[1][1]->Fill(mom_missSmear.Mag(),1./p_eff);
 	  	    h_emVpm[1][1]->Fill(mom_missSmear.Mag(),emissSmear,1./p_eff);
 	  	    h_thetapq[1][1]->Fill(mom_smear.Angle(q)*180./M_PI,1./p_eff);
 	  	    h_nmom[1][1]->Fill(mom_smear.Mag(),1./p_eff);
 	  	    h_ntheta[1][1]->Fill(theta,1./p_eff);
 	  	    h_emiss[1][1]->Fill(emissSmear,1./p_eff);
 	  	    h_mmiss[1][1]->Fill(mmissSmear,1./p_eff);
 	  	    
 	  	  }

 	}

 	}

 
  
  intree->Delete();
  //open the 4He neutron file
  intree = (TTree*)f_4n->Get("T");
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  for (int event =0 ; event < intree->GetEntries() ; event++)
     // for (int event =0 ; event < 100000 ; event++)
     {
      
       if (event%10000==0){cout << "Working on event " << event <<endl;}

       intree->GetEvent(event);
       //if (nRun!=18020){continue;}

       if (QSq<Q2_min){continue;}

        // Get the electron and q vectors
       TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
       double Ebeam = Nu + mom_e.Mag();
       TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
       
       //mean field stuff:
       double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mn);

       // Loop over particles looking for a proton most likely to pass lead cuts
       int leadingID = -1;
       for (int part=1 ; part < nParticles ; part++)
 	{   
 	  // Only look at neutrons
 	  if (type[part] != 2112)
 	    continue;

 	  //vertex cut-target dep.
 	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
 	  
 	  // Get a 3-vector for this proton's momentum
 	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

 	  // Angle between p and q
 	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
 	  if (theta_pq>thetapq_MFcut){continue;}
 	  
 	  //first proton of interest in event
 	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

 	  /////////////////////////////////////////////////////////////////////
 	  //Here I selected the proton with angle closest to q-vector:	  
 	  //current p/q is larger and angle is smaller than previous
 	  if (theta_pq<theta_prev){leadingID=part;}
 	  /////////////////////////////////////////////////////////////////////	  
 	  	  
 	}

       //Now we have a leading neutron!
       if (leadingID > 0)
 	{
 	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
 	  mom_lead = correctPnVec(mom_lead,ec_in[leadingID],ec_out[leadingID]);
 	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

 	  //calculate the nominal values
 	  TVector3 mom_miss = mom_lead - q;
 	  double theta = mom_lead.Theta()*180./M_PI;
 	  double phi = mom_lead.Phi();
 	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
 	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
 	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
 	  
 	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
 	  phi *= 180./M_PI; 
 	  
 	  //apply proton fiducial cuts:
 	  if (!fid_params4He.pFiducialCut(mom_lead)) {continue;}

 	  //apply EC 10cm edge projection cut:
 	  //this doesn't work, but I don't know why:
 	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
 	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
 	  // Cutting top of EC
 	  if(theta>44.) continue;
 	  if(mom_lead.Mag()>1.95) continue;
 	  
 	  //passing lead proton cuts
 	  double p_over_q = mom_lead.Mag() / q.Mag();

 	  //check if in good acceptance
 	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
 	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}

 	  double n_eff = getNEff(mom_lead.Mag());


 	  //select the 4He neutron sample
 	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut 
 			  && mom_miss.Mag()<pmiss_MFnomCut_n && emiss<emiss_nomCut_n
 			  //&& mom_lead.Mag()>1.4 && mom_lead.Mag()<1.8
 			  ){
 		  
 		 	    h_emom[1][2]->Fill(mom_e.Mag(),1./n_eff); 
 		 	    h_etheta[1][2]->Fill(mom_e.Theta()*180./M_PI,1./n_eff);
 		 	    h_q2[1][2]->Fill(QSq,1./n_eff);
 		 	    h_nu[1][2]->Fill(Nu,1./n_eff);
 		 	    h_xb[1][2]->Fill(xB,1./n_eff);
 		 	    h_pmiss[1][2]->Fill(mom_miss.Mag(),1./n_eff);
 		 	    h_emVpm[1][2]->Fill(mom_miss.Mag(),emiss,1./n_eff);
 		 	    h_thetapq[1][2]->Fill(mom_lead.Angle(q)*180./M_PI,1./n_eff);
 		 	    h_nmom[1][2]->Fill(mom_lead.Mag(),1./n_eff);
 		 	    h_ntheta[1][2]->Fill(theta,1./n_eff);
 		 	    h_emiss[1][2]->Fill(emiss,1./n_eff);
 		 	    h_mmiss[1][2]->Fill(mmiss,1./n_eff);
 		 	    
 		 	  }	  
  
 	}//end leading
     }//end event

  
  ////////////////////////////////
  //Carbon-12
  ////////////////////////////////
  target = "12C";
  Acceptance myMapC("C12", 4461, 2250, "p");

  tgt_min = vtx_min_C; tgt_max = vtx_max_C; tgt_Z = 6; tgt_A = 12; tgt_M = M12C; tgt_R = M11B;
  Fiducial fid_paramsC(4461,2250,5996,"12C",true);
  intree->Delete();
	gRandom->SetSeed(4357);

  //open the 12C proton file
  intree = (TTree*)f_5p->Get("T");
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);

  //get (e,e'p)
   for (int event =0 ; event < intree->GetEntries() ; event++)
   //    for (int event =0 ; event < 100000 ; event++)
     {
        if (event%10000==0)
       	cout << "Working on event " << event <<endl;

       intree->GetEvent(event);
       
       //if (nRun!=18020){continue;}
       
       if (QSq<Q2_min){continue;}

        // Get the electron and q vectors
       TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
       double Ebeam = Nu + mom_e.Mag();
       TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
       
       //mean field stuff:
       double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mp);

       // Loop over particles looking for a proton most likely to pass lead cuts
       int leadingID = -1;
       for (int part=1 ; part < nParticles ; part++)
 	{   
 	  // Only look at protons
 	  if (type[part] != 2212)
 	    continue;

 	  //vertex cut-target dep.
 	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
 	  
 	  // Get a 3-vector for this proton's momentum
 	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

 	  // Angle between p and q
 	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
 	  if (theta_pq>thetapq_MFcut){continue;}
 	  
 	  //first proton of interest in event
 	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

 	  /////////////////////////////////////////////////////////////////////
 	  //Here I selected the proton with angle closest to q-vector:	  
 	  //current p/q is larger and angle is smaller than previous
 	  if (theta_pq<theta_prev){leadingID=part;}
 	  /////////////////////////////////////////////////////////////////////	  
 	 
 	}

       //Now we have a leading proton!
       if (leadingID > 0)
 	{
 	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
 	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

 	  //if (mom_lead.Mag()>2.34){continue;}
 	  if (mom_lead.Mag()>1.95){continue;}
 	  
 	  TVector3 mom_miss = mom_lead - q;
 	  double theta = mom_lead.Theta()*180./M_PI;
 	  double phi = mom_lead.Phi();
 	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());
 	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
 	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);

 	  //smear the proton momentum resolution
 	  double smearFactor = gRandom->Gaus(mom_lead.Mag(),f_delta_p(mom_lead.Mag()));
 	  TVector3 u1 = mom_lead.Unit();
 	  TVector3 mom_smear(smearFactor*u1.X(),smearFactor*u1.Y(),smearFactor*u1.Z());

 	  //recalculate the smeared factors
 	  TVector3 mom_missSmear = mom_smear - q;
 	  double epSmear = sqrt(pow(Mp,2.0)+mom_smear.Mag2());
 	  double emissSmear = fn_Emiss(mom_missSmear.Mag(), Nu, tgt_M, epSmear, Mp);
 	  double mmissSmear = fn_Mmiss(Nu, Mp, epSmear, q, mom_smear);
 	  double p_over_q_sm = mom_smear.Mag()/q.Mag();
 	  
 	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
 	  phi *= 180./M_PI; 

 	  //apply proton fiducial cuts:

 	  //apply EC 10cm edge projection cut:
 	  //this doesn't work, but I don't know why:
 	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
 	  // Cutting top of EC
 	  if(theta>44.) continue;
 	  
 	  //passing lead proton cuts
 	  double p_over_q = mom_smear.Mag() / q.Mag();

 	  //check if in good acceptance
 	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
 	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}
 	  
		 double p_eff = myMapC.get_acc(mom_lead);

 	    
 	  //select the 12C proton sample
 	  	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut && mom_miss.Mag()<pmiss_MFnomCut && emiss<emiss_nomCut
 	  			  &&mom_lead.Mag()<1.95 && p_points_to_EC_fid(mom_lead,vtxZ[0],10.)&&fid_params.pFiducialCut(mom_lead)){
 	  	    h_emom[2][0]->Fill(mom_e.Mag()); 
 	  	    h_etheta[2][0]->Fill(mom_e.Theta()*180./M_PI);
 	  	    h_q2[2][0]->Fill(QSq);
 	  	    h_nu[2][0]->Fill(Nu);
 	  	    h_xb[2][0]->Fill(xB);
 	  	    h_pmiss[2][0]->Fill(mom_miss.Mag());
 	  	    h_emVpm[2][0]->Fill(mom_miss.Mag(),emiss);
 	  	    h_thetapq[2][0]->Fill(mom_lead.Angle(q)*180./M_PI);
 	  	    h_nmom[2][0]->Fill(mom_lead.Mag());
 	  	    h_ntheta[2][0]->Fill(theta);
 	  	    h_emiss[2][0]->Fill(emiss);
 	  	    h_mmiss[2][0]->Fill(mmiss);	    
 	  	  }

 	  	  //select the 12C smeared proton sample
 	  	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_smear.Angle(q)*180./M_PI<thetapq_MFcut && mom_missSmear.Mag()<pmiss_MFnomCut_n && emissSmear<emiss_nomCut_n
 	  			  &&mom_smear.Mag()<1.95 && p_points_to_EC_fid(mom_smear,vtxZ[0],10.)&&fid_params.pFiducialCut(mom_smear)){
 	  	    h_emom[2][1]->Fill(mom_e.Mag(),1./p_eff); 
 	  	    h_etheta[2][1]->Fill(mom_e.Theta()*180./M_PI,1./p_eff);
 	  	    h_q2[2][1]->Fill(QSq,1./p_eff);
 	  	    h_nu[2][1]->Fill(Nu,1./p_eff);
 	  	    h_xb[2][1]->Fill(xB,1./p_eff);
 	  	    h_pmiss[2][1]->Fill(mom_missSmear.Mag(),1./p_eff);
 	  	    h_emVpm[2][1]->Fill(mom_missSmear.Mag(),emissSmear,1./p_eff);
 	  	    h_thetapq[2][1]->Fill(mom_smear.Angle(q)*180./M_PI,1./p_eff);
 	  	    h_nmom[2][1]->Fill(mom_smear.Mag(),1./p_eff);
 	  	    h_ntheta[2][1]->Fill(theta,1./p_eff);
 	  	    h_emiss[2][1]->Fill(emissSmear,1./p_eff);
 	  	    h_mmiss[2][1]->Fill(mmissSmear,1./p_eff);
 	  	    
 	  	  }

 	}

 	}

 
  
  intree->Delete();
  //open the 12C neutron file
  intree = (TTree*)f_5n->Get("T");
  intree->SetBranchAddress("nParticles",&nParticles);
  intree->SetBranchAddress("Part_type",type);
  intree->SetBranchAddress("Q2",&QSq);
  intree->SetBranchAddress("Xb",&xB);
  intree->SetBranchAddress("Nu",&Nu);
  intree->SetBranchAddress("mom_x",mom_x);
  intree->SetBranchAddress("mom_y",mom_y);
  intree->SetBranchAddress("mom_z",mom_z);
  intree->SetBranchAddress("ec_x",ec_x);
  intree->SetBranchAddress("ec_y",ec_y);
  intree->SetBranchAddress("ec_u",ec_u);
  intree->SetBranchAddress("ec_v",ec_v);
  intree->SetBranchAddress("ec_w",ec_w);
  intree->SetBranchAddress("ec_in",ec_in);
  intree->SetBranchAddress("ec_out",ec_out);
  intree->SetBranchAddress("vtx_z_cor",vtxZ);
  intree->SetBranchAddress("sc_time",sc_time);
  intree->SetBranchAddress("ec_time",ec_time);
  intree->SetBranchAddress("sc_path",sc_path);
  intree->SetBranchAddress("ec_path",ec_path);
  
  for (int event =0 ; event < intree->GetEntries() ; event++)
     // for (int event =0 ; event < 100000 ; event++)
     {
      
       if (event%10000==0){cout << "Working on event " << event <<endl;}

       intree->GetEvent(event);
       //if (nRun!=18020){continue;}

       if (QSq<Q2_min){continue;}

        // Get the electron and q vectors
       TVector3 mom_e(mom_x[0],mom_y[0],mom_z[0]);
       double Ebeam = Nu + mom_e.Mag();
       TVector3 q = TVector3(0.,0.,Ebeam) - mom_e;
       
       //mean field stuff:
       double y_sc = y_scale_func(q,Nu,tgt_M,tgt_R,Mn);

       // Loop over particles looking for a proton most likely to pass lead cuts
       int leadingID = -1;
       for (int part=1 ; part < nParticles ; part++)
 	{   
 	  // Only look at neutrons
 	  if (type[part] != 2112)
 	    continue;

 	  //vertex cut-target dep.
 	  if (vtxZ[part]<tgt_min || vtxZ[part]>tgt_max){continue;}
 	  
 	  // Get a 3-vector for this proton's momentum
 	  TVector3 this_mom(mom_x[part],mom_y[part],mom_z[part]);

 	  // Angle between p and q
 	  double theta_pq = this_mom.Angle(q)*180.0/M_PI;
 	  if (theta_pq>thetapq_MFcut){continue;}
 	  
 	  //first proton of interest in event
 	  if (leadingID==-1){theta_prev=theta_pq; leadingID=part;}

 	  /////////////////////////////////////////////////////////////////////
 	  //Here I selected the proton with angle closest to q-vector:	  
 	  //current p/q is larger and angle is smaller than previous
 	  if (theta_pq<theta_prev){leadingID=part;}
 	  /////////////////////////////////////////////////////////////////////	  
 	  	  
 	}

       //Now we have a leading neutron!
       if (leadingID > 0)
 	{
 	  TVector3 mom_lead(mom_x[leadingID],mom_y[leadingID],mom_z[leadingID]);
 	  mom_lead = correctPnVec(mom_lead,ec_in[leadingID],ec_out[leadingID]);
 	  //mom_lead = P_Coulomb_Corrected(mom_lead,tgt_Z,tgt_A);

 	  //calculate the nominal values
 	  TVector3 mom_miss = mom_lead - q;
 	  double theta = mom_lead.Theta()*180./M_PI;
 	  double phi = mom_lead.Phi();
 	  double ep = sqrt(pow(Mp,2.0)+mom_lead.Mag2());	  
 	  double emiss = fn_Emiss(mom_miss.Mag(), Nu, tgt_M, ep, Mp);
 	  double mmiss = fn_Mmiss(Nu, Mp, ep, q, mom_lead);
 	  
 	  if (phi < -M_PI/6.) phi+= 2.*M_PI;
 	  phi *= 180./M_PI; 
 	  
 	  //apply proton fiducial cuts:
 	  if (!fid_paramsC.pFiducialCut(mom_lead)) {continue;}

 	  //apply EC 10cm edge projection cut:
 	  //this doesn't work, but I don't know why:
 	  //if (!fid_params.CutUVW(pECUVW,10.0)) {continue;}
 	  if(!p_points_to_EC_fid(mom_lead,vtxZ[0],10.)) {continue;}
 	  // Cutting top of EC
 	  if(theta>44.) continue;
 	  if(mom_lead.Mag()>1.95) continue;
 	  
 	  //passing lead proton cuts
 	  double p_over_q = mom_lead.Mag() / q.Mag();

 	  //check if in good acceptance
 	  if (phi>100 && phi<150 && theta>25 && theta<35){continue;}
 	  if (phi>150 && phi<200 && theta>35 && theta<40){continue;}

 	  double n_eff = getNEff(mom_lead.Mag());

 	  
 	  //select the 12C neutron sample
 	  if (y_sc>y_cut_min && y_sc<y_cut_max && nu_min<Nu && Nu<nu_max && mom_lead.Angle(q)*180./M_PI<thetapq_MFcut 
 			  && mom_miss.Mag()<pmiss_MFnomCut_n && emiss<emiss_nomCut_n
 			  //&& mom_lead.Mag()>1.4 && mom_lead.Mag()<1.8
 			  ){
 		  
 		 	    h_emom[2][2]->Fill(mom_e.Mag(),1./n_eff); 
 		 	    h_etheta[2][2]->Fill(mom_e.Theta()*180./M_PI,1./n_eff);
 		 	    h_q2[2][2]->Fill(QSq,1./n_eff);
 		 	    h_nu[2][2]->Fill(Nu,1./n_eff);
 		 	    h_xb[2][2]->Fill(xB,1./n_eff);
 		 	    h_pmiss[2][2]->Fill(mom_miss.Mag(),1./n_eff);
 		 	    h_emVpm[2][2]->Fill(mom_miss.Mag(),emiss,1./n_eff);
 		 	    h_thetapq[2][2]->Fill(mom_lead.Angle(q)*180./M_PI,1./n_eff);
 		 	    h_nmom[2][2]->Fill(mom_lead.Mag(),1./n_eff);
 		 	    h_ntheta[2][2]->Fill(theta,1./n_eff);
 		 	    h_emiss[2][2]->Fill(emiss,1./n_eff);
 		 	    h_mmiss[2][2]->Fill(mmiss,1./n_eff);
 		 	    
 		 	  }	  
  
 	}//end leading
     }//end event
  
  //f->Close();
  int protons_3He, protons_4He, protons_12C, smeared_3He, smeared_4He, smeared_12C, neutrons_3He, neutrons_4He, neutrons_12C;
  protons_3He = h_emom[0][0]->GetEntries();
  neutrons_3He = h_emom[0][2]->GetEntries();
  smeared_3He =  h_emom[0][1]->GetEntries();
  protons_4He = h_emom[1][0]->GetEntries();
  neutrons_4He = h_emom[1][2]->GetEntries();
  smeared_4He =  h_emom[1][1]->GetEntries();
  protons_12C = h_emom[2][0]->GetEntries();
  neutrons_12C = h_emom[2][2]->GetEntries();
  smeared_12C =  h_emom[2][1]->GetEntries();
  
  canvas->Update();

  h_emVpm[0][0]->SetTitle("3He protons");
  h_emVpm[0][0]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[0][1]->SetTitle("3He smeared protons");
  h_emVpm[0][1]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[0][2]->SetTitle("3He neutrons");
  h_emVpm[0][2]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm[1][0]->SetTitle("4He protons");
  h_emVpm[1][0]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[1][1]->SetTitle("4He smeared protons");
  h_emVpm[1][1]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[1][2]->SetTitle("4He neutrons");
  h_emVpm[1][2]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emVpm[2][0]->SetTitle("12C protons");
  h_emVpm[2][0]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[2][1]->SetTitle("12C smeared protons");
  h_emVpm[2][1]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_emVpm[2][2]->SetTitle("12C neutrons");
  h_emVpm[2][2]->Draw("colz");
  canvas->Print( (pdf_file_name + "(").c_str());

  ///////////////////////////////////////
  //Compare neutrons and smeared protons
  ///////////////////////////////////////
  //Helium-3
  h_emom[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_emom[0][1]->SetLineColor(kRed);
  h_emom[0][2]->SetLineColor(kBlue);
  h_emom[0][1]->Scale(1./h_emom[0][1]->Integral());
  h_emom[0][2]->Scale(1./h_emom[0][2]->Integral());
  h_emom[0][1]->Draw("hist");
  h_emom[0][2]->Draw("histsame");
  double amax =  h_emom[0][1]->GetMaximum();
  double bmax =  h_emom[0][2]->GetMaximum();
  if (amax<bmax){h_emom[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_etheta[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_etheta[0][1]->SetLineColor(kRed);
  h_etheta[0][2]->SetLineColor(kBlue);
  h_etheta[0][1]->Scale(1./h_etheta[0][1]->Integral());
  h_etheta[0][2]->Scale(1./h_etheta[0][2]->Integral());
  h_etheta[0][1]->Draw("hist");
  h_etheta[0][2]->Draw("histsame");
  amax =  h_etheta[0][1]->GetMaximum();
  bmax =  h_etheta[0][2]->GetMaximum();
  if (amax<bmax){h_etheta[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_q2[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_q2[0][1]->SetLineColor(kRed);
  h_q2[0][2]->SetLineColor(kBlue);
  h_q2[0][1]->Scale(1./h_q2[0][1]->Integral());
  h_q2[0][2]->Scale(1./h_q2[0][2]->Integral());
  h_q2[0][1]->Draw("hist");
  h_q2[0][2]->Draw("histsame");
  amax =  h_q2[0][1]->GetMaximum();
  bmax =  h_q2[0][2]->GetMaximum();
  if (amax<bmax){h_q2[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_nu[0][1]->SetLineColor(kRed);
  h_nu[0][2]->SetLineColor(kBlue);
  h_nu[0][1]->Scale(1./h_nu[0][1]->Integral());
  h_nu[0][2]->Scale(1./h_nu[0][2]->Integral());
  h_nu[0][1]->Draw("hist");
  h_nu[0][2]->Draw("histsame");
  amax =  h_nu[0][1]->GetMaximum();
  bmax =  h_nu[0][2]->GetMaximum();
  if (amax<bmax){h_nu[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_xb[0][1]->SetLineColor(kRed);
  h_xb[0][2]->SetLineColor(kBlue);
  h_xb[0][1]->Scale(1./h_xb[0][1]->Integral());
  h_xb[0][2]->Scale(1./h_xb[0][2]->Integral());
  h_xb[0][1]->Draw("hist");
  h_xb[0][2]->Draw("histsame");
  amax =  h_xb[0][1]->GetMaximum();
  bmax =  h_xb[0][2]->GetMaximum();
  if (amax<bmax){h_xb[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_pmiss[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_pmiss[0][1]->SetLineColor(kRed);
  h_pmiss[0][2]->SetLineColor(kBlue);
  h_pmiss[0][1]->Scale(1./h_pmiss[0][1]->Integral());
  h_pmiss[0][2]->Scale(1./h_pmiss[0][2]->Integral());
  h_pmiss[0][1]->Draw("hist");
  h_pmiss[0][2]->Draw("histsame");
  amax =  h_pmiss[0][1]->GetMaximum();
  bmax =  h_pmiss[0][2]->GetMaximum();
  if (amax<bmax){h_pmiss[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_thetapq[0][1]->SetLineColor(kRed);
  h_thetapq[0][2]->SetLineColor(kBlue);
  h_thetapq[0][1]->Scale(1./h_thetapq[0][1]->Integral());
  h_thetapq[0][2]->Scale(1./h_thetapq[0][2]->Integral());
  h_thetapq[0][1]->Draw("hist");
  h_thetapq[0][2]->Draw("histsame");
  amax =  h_thetapq[0][1]->GetMaximum();
  bmax =  h_thetapq[0][2]->GetMaximum();
  if (amax<bmax){h_thetapq[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_nmom[0][1]->SetLineColor(kRed);
  h_nmom[0][2]->SetLineColor(kBlue);
  h_nmom[0][1]->Scale(1./h_nmom[0][1]->Integral());
  h_nmom[0][2]->Scale(1./h_nmom[0][2]->Integral());
  h_nmom[0][1]->Draw("hist");
  h_nmom[0][2]->Draw("histsame");
  amax =  h_nmom[0][1]->GetMaximum();
  bmax =  h_nmom[0][2]->GetMaximum();
  if (amax<bmax){h_nmom[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_ntheta[0][1]->SetLineColor(kRed);
  h_ntheta[0][2]->SetLineColor(kBlue);
  h_ntheta[0][1]->Scale(1./h_ntheta[0][1]->Integral());
  h_ntheta[0][2]->Scale(1./h_ntheta[0][2]->Integral());
  h_ntheta[0][1]->Draw("hist");
  h_ntheta[0][2]->Draw("histsame");
  amax =  h_ntheta[0][1]->GetMaximum();
  bmax =  h_ntheta[0][2]->GetMaximum();
  if (amax<bmax){h_ntheta[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emiss[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_emiss[0][1]->SetLineColor(kRed);
  h_emiss[0][2]->SetLineColor(kBlue);
  h_emiss[0][1]->Scale(1./h_emiss[0][1]->Integral());
  h_emiss[0][2]->Scale(1./h_emiss[0][2]->Integral());
  h_emiss[0][1]->Draw("hist");
  h_emiss[0][2]->Draw("histsame");
  amax =  h_emiss[0][1]->GetMaximum();
  bmax =  h_emiss[0][2]->GetMaximum();
  if (amax<bmax){h_emiss[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss[0][1]->SetTitle("3He (smeared p-r, neutrons-b)");
  h_mmiss[0][1]->SetLineColor(kRed);
  h_mmiss[0][2]->SetLineColor(kBlue);
  h_mmiss[0][1]->Scale(1./h_mmiss[0][1]->Integral());
  h_mmiss[0][2]->Scale(1./h_mmiss[0][2]->Integral());
  h_mmiss[0][1]->Draw("hist");
  h_mmiss[0][2]->Draw("histsame");
  amax =  h_mmiss[0][1]->GetMaximum();
  bmax =  h_mmiss[0][2]->GetMaximum();
  if (amax<bmax){h_mmiss[0][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  //Helium-4
  h_emom[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_emom[1][1]->SetLineColor(kRed);
  h_emom[1][2]->SetLineColor(kBlue);
  h_emom[1][1]->Scale(1./h_emom[1][1]->Integral());
  h_emom[1][2]->Scale(1./h_emom[1][2]->Integral());
  h_emom[1][1]->Draw("hist");
  h_emom[1][2]->Draw("histsame");
  amax =  h_emom[1][1]->GetMaximum();
  bmax =  h_emom[1][2]->GetMaximum();
  if (amax<bmax){h_emom[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_etheta[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_etheta[1][1]->SetLineColor(kRed);
  h_etheta[1][2]->SetLineColor(kBlue);
  h_etheta[1][1]->Scale(1./h_etheta[1][1]->Integral());
  h_etheta[1][2]->Scale(1./h_etheta[1][2]->Integral());
  h_etheta[1][1]->Draw("hist");
  h_etheta[1][2]->Draw("histsame");
  amax =  h_etheta[1][1]->GetMaximum();
  bmax =  h_etheta[1][2]->GetMaximum();
  if (amax<bmax){h_etheta[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_q2[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_q2[1][1]->SetLineColor(kRed);
  h_q2[1][2]->SetLineColor(kBlue);
  h_q2[1][1]->Scale(1./h_q2[1][1]->Integral());
  h_q2[1][2]->Scale(1./h_q2[1][2]->Integral());
  h_q2[1][1]->Draw("hist");
  h_q2[1][2]->Draw("histsame");
  amax =  h_q2[1][1]->GetMaximum();
  bmax =  h_q2[1][2]->GetMaximum();
  if (amax<bmax){h_q2[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_nu[1][1]->SetLineColor(kRed);
  h_nu[1][2]->SetLineColor(kBlue);
  h_nu[1][1]->Scale(1./h_nu[1][1]->Integral());
  h_nu[1][2]->Scale(1./h_nu[1][2]->Integral());
  h_nu[1][1]->Draw("hist");
  h_nu[1][2]->Draw("histsame");
  amax =  h_nu[1][1]->GetMaximum();
  bmax =  h_nu[1][2]->GetMaximum();
  if (amax<bmax){h_nu[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_xb[1][1]->SetLineColor(kRed);
  h_xb[1][2]->SetLineColor(kBlue);
  h_xb[1][1]->Scale(1./h_xb[1][1]->Integral());
  h_xb[1][2]->Scale(1./h_xb[1][2]->Integral());
  h_xb[1][1]->Draw("hist");
  h_xb[1][2]->Draw("histsame");
  amax =  h_xb[1][1]->GetMaximum();
  bmax =  h_xb[1][2]->GetMaximum();
  if (amax<bmax){h_xb[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_pmiss[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_pmiss[1][1]->SetLineColor(kRed);
  h_pmiss[1][2]->SetLineColor(kBlue);
  h_pmiss[1][1]->Scale(1./h_pmiss[1][1]->Integral());
  h_pmiss[1][2]->Scale(1./h_pmiss[1][2]->Integral());
  h_pmiss[1][1]->Draw("hist");
  h_pmiss[1][2]->Draw("histsame");
  amax =  h_pmiss[1][1]->GetMaximum();
  bmax =  h_pmiss[1][2]->GetMaximum();
  if (amax<bmax){h_pmiss[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_thetapq[1][1]->SetLineColor(kRed);
  h_thetapq[1][2]->SetLineColor(kBlue);
  h_thetapq[1][1]->Scale(1./h_thetapq[1][1]->Integral());
  h_thetapq[1][2]->Scale(1./h_thetapq[1][2]->Integral());
  h_thetapq[1][1]->Draw("hist");
  h_thetapq[1][2]->Draw("histsame");
  amax =  h_thetapq[1][1]->GetMaximum();
  bmax =  h_thetapq[1][2]->GetMaximum();
  if (amax<bmax){h_thetapq[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_nmom[1][1]->SetLineColor(kRed);
  h_nmom[1][2]->SetLineColor(kBlue);
  h_nmom[1][1]->Scale(1./h_nmom[1][1]->Integral());
  h_nmom[1][2]->Scale(1./h_nmom[1][2]->Integral());
  h_nmom[1][1]->Draw("hist");
  h_nmom[1][2]->Draw("histsame");
  amax =  h_nmom[1][1]->GetMaximum();
  bmax =  h_nmom[1][2]->GetMaximum();
  if (amax<bmax){h_nmom[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_ntheta[1][1]->SetLineColor(kRed);
  h_ntheta[1][2]->SetLineColor(kBlue);
  h_ntheta[1][1]->Scale(1./h_ntheta[1][1]->Integral());
  h_ntheta[1][2]->Scale(1./h_ntheta[1][2]->Integral());
  h_ntheta[1][1]->Draw("hist");
  h_ntheta[1][2]->Draw("histsame");
  amax =  h_ntheta[1][1]->GetMaximum();
  bmax =  h_ntheta[1][2]->GetMaximum();
  if (amax<bmax){h_ntheta[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emiss[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_emiss[1][1]->SetLineColor(kRed);
  h_emiss[1][2]->SetLineColor(kBlue);
  h_emiss[1][1]->Scale(1./h_emiss[1][1]->Integral());
  h_emiss[1][2]->Scale(1./h_emiss[1][2]->Integral());
  h_emiss[1][1]->Draw("hist");
  h_emiss[1][2]->Draw("histsame");
  amax =  h_emiss[1][1]->GetMaximum();
  bmax =  h_emiss[1][2]->GetMaximum();
  if (amax<bmax){h_emiss[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss[1][1]->SetTitle("4He (smeared p-r, neutrons-b)");
  h_mmiss[1][1]->SetLineColor(kRed);
  h_mmiss[1][2]->SetLineColor(kBlue);
  h_mmiss[1][1]->Scale(1./h_mmiss[1][1]->Integral());
  h_mmiss[1][2]->Scale(1./h_mmiss[1][2]->Integral());
  h_mmiss[1][1]->Draw("hist");
  h_mmiss[1][2]->Draw("histsame");
  amax =  h_mmiss[1][1]->GetMaximum();
  bmax =  h_mmiss[1][2]->GetMaximum();
  if (amax<bmax){h_mmiss[1][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  //Carbon-12
  h_emom[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_emom[2][1]->SetLineColor(kRed);
  h_emom[2][2]->SetLineColor(kBlue);
  h_emom[2][1]->Scale(1./h_emom[2][1]->Integral());
  h_emom[2][2]->Scale(1./h_emom[2][2]->Integral());
  h_emom[2][1]->Draw("hist");
  h_emom[2][2]->Draw("histsame");
  amax =  h_emom[2][1]->GetMaximum();
  bmax =  h_emom[2][2]->GetMaximum();
  if (amax<bmax){h_emom[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_etheta[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_etheta[2][1]->SetLineColor(kRed);
  h_etheta[2][2]->SetLineColor(kBlue);
  h_etheta[2][1]->Scale(1./h_etheta[2][1]->Integral());
  h_etheta[2][2]->Scale(1./h_etheta[2][2]->Integral());
  h_etheta[2][1]->Draw("hist");
  h_etheta[2][2]->Draw("histsame");
  amax =  h_etheta[2][1]->GetMaximum();
  bmax =  h_etheta[2][2]->GetMaximum();
  if (amax<bmax){h_etheta[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_q2[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_q2[2][1]->SetLineColor(kRed);
  h_q2[2][2]->SetLineColor(kBlue);
  h_q2[2][1]->Scale(1./h_q2[2][1]->Integral());
  h_q2[2][2]->Scale(1./h_q2[2][2]->Integral());
  h_q2[2][1]->Draw("hist");
  h_q2[2][2]->Draw("histsame");
  amax =  h_q2[2][1]->GetMaximum();
  bmax =  h_q2[2][2]->GetMaximum();
  if (amax<bmax){h_q2[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nu[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_nu[2][1]->SetLineColor(kRed);
  h_nu[2][2]->SetLineColor(kBlue);
  h_nu[2][1]->Scale(1./h_nu[2][1]->Integral());
  h_nu[2][2]->Scale(1./h_nu[2][2]->Integral());
  h_nu[2][1]->Draw("hist");
  h_nu[2][2]->Draw("histsame");
  amax =  h_nu[2][1]->GetMaximum();
  bmax =  h_nu[2][2]->GetMaximum();
  if (amax<bmax){h_nu[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_xb[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_xb[2][1]->SetLineColor(kRed);
  h_xb[2][2]->SetLineColor(kBlue);
  h_xb[2][1]->Scale(1./h_xb[2][1]->Integral());
  h_xb[2][2]->Scale(1./h_xb[2][2]->Integral());
  h_xb[2][1]->Draw("hist");
  h_xb[2][2]->Draw("histsame");
  amax =  h_xb[2][1]->GetMaximum();
  bmax =  h_xb[2][2]->GetMaximum();
  if (amax<bmax){h_xb[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  h_pmiss[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_pmiss[2][1]->SetLineColor(kRed);
  h_pmiss[2][2]->SetLineColor(kBlue);
  h_pmiss[2][1]->Scale(1./h_pmiss[2][1]->Integral());
  h_pmiss[2][2]->Scale(1./h_pmiss[2][2]->Integral());
  h_pmiss[2][1]->Draw("hist");
  h_pmiss[2][2]->Draw("histsame");
  amax =  h_pmiss[2][1]->GetMaximum();
  bmax =  h_pmiss[2][2]->GetMaximum();
  if (amax<bmax){h_pmiss[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_thetapq[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_thetapq[2][1]->SetLineColor(kRed);
  h_thetapq[2][2]->SetLineColor(kBlue);
  h_thetapq[2][1]->Scale(1./h_thetapq[2][1]->Integral());
  h_thetapq[2][2]->Scale(1./h_thetapq[2][2]->Integral());
  h_thetapq[2][1]->Draw("hist");
  h_thetapq[2][2]->Draw("histsame");
  amax =  h_thetapq[2][1]->GetMaximum();
  bmax =  h_thetapq[2][2]->GetMaximum();
  if (amax<bmax){h_thetapq[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_nmom[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_nmom[2][1]->SetLineColor(kRed);
  h_nmom[2][2]->SetLineColor(kBlue);
  h_nmom[2][1]->Scale(1./h_nmom[2][1]->Integral());
  h_nmom[2][2]->Scale(1./h_nmom[2][2]->Integral());
  h_nmom[2][1]->Draw("hist");
  h_nmom[2][2]->Draw("histsame");
  amax =  h_nmom[2][1]->GetMaximum();
  bmax =  h_nmom[2][2]->GetMaximum();
  if (amax<bmax){h_nmom[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_ntheta[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_ntheta[2][1]->SetLineColor(kRed);
  h_ntheta[2][2]->SetLineColor(kBlue);
  h_ntheta[2][1]->Scale(1./h_ntheta[2][1]->Integral());
  h_ntheta[2][2]->Scale(1./h_ntheta[2][2]->Integral());
  h_ntheta[2][1]->Draw("hist");
  h_ntheta[2][2]->Draw("histsame");
  amax =  h_ntheta[2][1]->GetMaximum();
  bmax =  h_ntheta[2][2]->GetMaximum();
  if (amax<bmax){h_ntheta[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_emiss[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_emiss[2][1]->SetLineColor(kRed);
  h_emiss[2][2]->SetLineColor(kBlue);
  h_emiss[2][1]->Scale(1./h_emiss[2][1]->Integral());
  h_emiss[2][2]->Scale(1./h_emiss[2][2]->Integral());
  h_emiss[2][1]->Draw("hist");
  h_emiss[2][2]->Draw("histsame");
  amax =  h_emiss[2][1]->GetMaximum();
  bmax =  h_emiss[2][2]->GetMaximum();
  if (amax<bmax){h_emiss[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());

  h_mmiss[2][1]->SetTitle("12C (smeared p-r, neutrons-b)");
  h_mmiss[2][1]->SetLineColor(kRed);
  h_mmiss[2][2]->SetLineColor(kBlue);
  h_mmiss[2][1]->Scale(1./h_mmiss[2][1]->Integral());
  h_mmiss[2][2]->Scale(1./h_mmiss[2][2]->Integral());
  h_mmiss[2][1]->Draw("hist");
  h_mmiss[2][2]->Draw("histsame");
  amax =  h_mmiss[2][1]->GetMaximum();
  bmax =  h_mmiss[2][2]->GetMaximum();
  if (amax<bmax){h_mmiss[2][1]->SetMaximum(bmax+0.03);}
  canvas->Print( (pdf_file_name + "(").c_str());
  
  //////////////////////////
  //Compare amongst targets
  //////////////////////////
  /*

  //h_emom[0][0]->SetTitle("protons (3He-r, 4He-b, 12C-p)");
  //h_emom[0][0]->SetLineColor(kRed);
  //h_emom[1][0]->SetLineColor(kBlue);
  //h_emom[2][0]->SetLineColor("kMagenta");
  //h_emom[0][0]->Scale(1./h_emom[0][0]->Integral());
  //h_emom[1][0]->Scale(1./h_emom[1][0]->Integral());
  //h_emom[2][0]->Scale(1./h_emom[2][0]->Integral());
  //h_emom[0][0]->Draw("hist");
  //h_emom[1][0]->Scale(1./h_emom[1][0]->Integral());
  //h_emom[2][0]->Scale(1./h_emom[2][0]->Integral());
  */
  canvas->Print( (pdf_file_name + ")").c_str());

  cout<<"3He protons: "<<protons_3He<<" 3He neutrons: "<<neutrons_3He<<" 3He smeared: "<<smeared_3He<<endl;
  cout<<"4He protons: "<<protons_4He<<" 4He neutrons: "<<neutrons_4He<<" 4He smeared: "<<smeared_4He<<endl;
  cout<<"12C protons: "<<protons_12C<<" 12C neutrons: "<<neutrons_12C<<" 12C smeared: "<<smeared_12C<<endl;

  

  fout->Write();
  fout->Close();

}
