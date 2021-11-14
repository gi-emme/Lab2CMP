//B --> K+pi

#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h" 
#include "TCanvas.h" 
#include "TFile.h"
#include "TLegend.h"
#include "TTree.h"

//Vector3D class with methods to extract angles ecc .. 
#include "Vector3D.h"

#include <iostream>
#include <vector>

using namespace std;

//this function returns the norm of the three-momentum of the products of the decay a -> b+c in the CoM system
double p_CoM(const double& m_a, const double& m_b, const double& m_c){ 
	return sqrt((m_a*m_a*m_a*m_a+(m_b*m_b-m_c*m_c)*(m_b*m_b-m_c*m_c)-2*m_a*m_a*(m_b*m_b+m_c*m_c))/(4*m_a*m_a));
}

int main(){

	int n_evts = 1e4;
	cout << "----Simulating " << n_evts << " B meson decays----" << endl;

	//masses
	double m_B = 5.279; //GeV
	double m_pi = 0.140; //GeV
	double m_K = 0.500; //GeV

	double resol = 0.03; //detector resolution

	TString rootfname("./output.root"); 
  	// Open TFile for output
  	// Overwite output file if it already exists (done with the RECREATE option)
  	TFile rfile(rootfname, "RECREATE"); 

  	// Open the ROOT file and make sure the file was opened successfully
  	if( !rfile.IsOpen() ) {
  	  cout << "problems creating root file. exiting... " << endl;
  	  exit(-1);
  	}
  	cout << "storing output in root file " << rootfname << endl;

  	// Create new histogram (this one for the invariant mass)
	int nbins1 = 100; 
	double xlo1 = m_B - 0.5;
  	double xhi1 = m_B + 0.5;
  	double binwidth1 = (xhi1-xlo1) / nbins1;

  	// Create new histogram (this one for the angle between the momenta)
	int nbins2 = 2000; 
	double xlo2 = 0;
  	double xhi2 = M_PI + 0.5;
  	double binwidth2 = (xhi2-xlo2) / nbins2;

  	// Create new histogram (this one for the measured invariant mass)
	int nbins3 = 100; 
	double xlo3 = m_B - 0.5;
  	double xhi3 = m_B + 0.5;
  	double binwidth3 = (xhi3-xlo3) / nbins3;


  	//true invariant mass histogram
  	TH1F h_inv_mass("true inv. mass", "distribution of invariant mass", nbins1, xlo1, xhi1);
  	//opening angle histogram
  	TH1F h_op_angle("\\phi (LAB)", "distribution of opening angles", nbins2, xlo2, xhi2);
  	//measured inv. mass hist
  	TH1F h_inv_mass_meas("meas. inv. mass", "distribution of measured invariant mass", nbins3, xlo3, xhi3);


	//B meson 4-momentum in the LAB frame
  	TLorentzVector p4_B; 
  	double p_B = 0.3; // GeV
  	// Flat metric, (- - - +) signature: m^2 = E^2 - p^2
  	p4_B.SetPxPyPzE(p_B, 0, 0, sqrt(p_B*p_B+m_B*m_B));

  	//Starting up a random generator
  	TRandom*  gen = new TRandom();
  	// ...exploiting the machine clock for the seed
  	gen->SetSeed(0);

  	double px_K, py_K, pz_K;
  	TLorentzVector p4CoM_K, p4CoM_pi;

  	Vector3D p3_K = Vector3D::Cartesian(0.,0.,0.);
  	Vector3D p3_pi = Vector3D::Cartesian(0.,0.,0.);

  	TLorentzVector p4Meas_K;
  	TLorentzVector p4Meas_pi;

  	double p_star = p_CoM(m_B, m_K, m_pi); //GeV
  	double sqrt_s = 0., sqrt_s_meas = 0.;
  	double phi = 0.; //opening angle between the products in the LAB
  	double p_K_0=0., p_pi_0=0.;
  	double p_K_meas=0., p_pi_meas=0.;

  	//loop over events 
  	for(int i = 0; i<n_evts; ++i){
  		//generating the three-momentum (of the K) in the CoM
  		gen->Sphere(px_K, py_K, pz_K, p_star); //random 3D vector with norm p_star -> stores the coordinates in px, py, pz

  		p3_K = Vector3D::Cartesian(px_K, py_K, pz_K);
  		p3_pi = (-1)*p3_K; //Conservation of momentum

  		//four-vectors of the products in the CoM
  		p4CoM_K.SetPxPyPzE(px_K, py_K, pz_K, sqrt(p3_K.scalarProduct(p3_K)+m_K*m_K));
  		p4CoM_pi.SetPxPyPzE(p3_pi.x(), p3_pi.y(), p3_pi.z(), sqrt(p3_pi.scalarProduct(p3_pi)+m_pi*m_pi));

  		//boosting to obtain 4-momenta in the LAB system
  		p4CoM_K.Boost(p4_B.BoostVector()); 
  		p4CoM_pi.Boost(p4_B.BoostVector()); 
  		//vector components:
  		TVector3 p_pi = p4CoM_pi.Vect();
  		TVector3 p_K = p4CoM_K.Vect();

  		//true norm of the momenta (in the LAB)
  		p_K_0 = sqrt(p_K.Dot(p_K));
  		p_pi_0 = sqrt(p_pi.Dot(p_pi));
  		//measured norm of the momenta: generated with gaussian model for the detector
  		p_K_meas = gen -> Gaus(p_K_0,resol*p_K_0);
  		p_pi_meas = gen -> Gaus(p_pi_0,resol*p_pi_0);

  		//new 4-momenta (LAB) with measured p:
  		p_pi *= (p_pi_meas/p_pi_0);
  		p_K *= (p_K_meas/p_K_0); 
  		p4Meas_pi.SetPxPyPzE(p_pi.X(), p_pi.Y(), p_pi.Z(), sqrt(m_pi*m_pi+p_pi_meas*p_pi_meas));
  		p4Meas_K.SetPxPyPzE(p_K.X(), p_K.Y(), p_K.Z(), sqrt(m_K*m_K+p_K_meas*p_K_meas));



  		TLorentzVector p4_tot = p4CoM_pi + p4CoM_K; //total 4-momentum in the LAB
  		TLorentzVector p4_tot_Meas = p4Meas_K + p4Meas_pi; //total measured 4-momentum

  		//invariant mass 
  		sqrt_s = sqrt(p4_tot*p4_tot);
  		//opening angle:
  		phi = p_pi.Angle(p_K);

  		//measured invariant mass
  		sqrt_s_meas = sqrt(p4_tot_Meas*p4_tot_Meas);


  		//fill histograms
  		h_inv_mass.Fill(sqrt_s);
  		h_op_angle.Fill(phi);
  		h_inv_mass_meas.Fill(sqrt_s_meas);
  	}

  	//plotting the results:

  	//nice features for overlayed histograms.. 
  	h_inv_mass.SetFillColor(kRed);
   	h_inv_mass_meas.SetFillColor(kBlue);

   	//creating canvas
  	TCanvas canv("canv", "canvas for plotting", 1280, 1024);

  	//true mass hist
  	h_inv_mass.GetXaxis()->SetTitle("\\sqrt{s} [GeV]");
  	h_inv_mass.Draw();
  	canv.SaveAs("./true-mass.pdf");

  	//drawing opening angles histogram
  	h_op_angle.GetXaxis()->SetTitle("opening angle \\phi [rad]");
  	h_op_angle.Draw();
  	canv.SaveAs("./opening-angle.pdf");

  	//drawing measured mass histogram
  	h_inv_mass_meas.GetXaxis()->SetTitle("\\sqrt{s} [GeV]");
  	h_inv_mass_meas.Draw();
  	canv.SaveAs("./measured-mass.pdf");
  	h_inv_mass.Draw("same");

  	//legend
  	TLegend *legend = new TLegend(0.1,0.9,0.4,0.8);
   	legend->SetHeader("","C"); 
   	legend->AddEntry(&h_inv_mass,"true \\sqrt{s}");
   	legend->AddEntry(&h_inv_mass_meas,"measured \\sqrt{s}");
   	legend->Draw();
  	canv.SaveAs("./invariant-mass.pdf");


  	//Now simulate the effect of different detectors with 1%, 5%, and 10% momentum resolution
  	double d_resol[3] = {0.01, 0.05, 0.10};

  	//storing data in a TTree
  	// Open a root file
  	TString Tfname("./data.root");
  	TFile orootfile(Tfname, "RECREATE"); 
  	if( !orootfile.IsOpen() ) {
  	  cout << "problems creating root file. exiting... " << endl;
  	  exit(-1);
  	}
  	cout << "storing output in root file " << Tfname << endl;


  	vector<TH1F> histo; //vector of histograms
  	TLorentzVector p4LAB_K;
	TLorentzVector p4LAB_pi;
	TVector3 p3LAB_K, p3LAB_pi;
	double p_K_LAB, p_pi_LAB, p_K_LAB_meas, p_pi_LAB_meas, theta_K, theta_pi, phi_K, phi_pi;

  	for(int i = 0; i<3 ; ++i){		//loop over det. resolutions
  		double resolution = d_resol[i];
  		TString histo_name("histo"+to_string(i)); 
		TH1F* h_mass = new TH1F(histo_name, "distribution of invariant mass", nbins1, 4., 6.6); 

		// Create a new TTree object
		TString tree_name("datatree"+to_string(i)); 
  		TTree* tree = new TTree(tree_name, "tree containing our data");

  		int nDau;
		double TpB;
		vector<double> Tnmass, Ttheta, Tphi, Tp;

  		//tree branches: Branch(branchname, address)
  		tree->Branch("B momentum", &TpB); 
  		tree->Branch("# daughters", &nDau); 
  		tree->Branch("Daughter masses", &Tnmass); 
  		tree->Branch("Daughter momenta (LAB)", &Tp);
  		tree->Branch("Theta angles (LAB)",&Ttheta);
  		tree->Branch("Phi angles (LAB)", &Tphi);

  		for(int j=0; j<n_evts ; ++j){		//loop over events

  			TpB = p_B;
  			nDau = 2;
  			Tnmass.push_back(m_K);
  			Tnmass.push_back(m_pi);

  			gen->Sphere(px_K, py_K, pz_K, p_star); //random 3D vector with norm p_star -> stores the coordinates in px, py, pz

  			p3_K = Vector3D::Cartesian(px_K, py_K, pz_K);
  			p3_pi = (-1)*p3_K; //Conservation of momentum

  			//four-vectors of the products in the CoM
  			p4CoM_K.SetPxPyPzE(px_K, py_K, pz_K, sqrt(p3_K.scalarProduct(p3_K)+m_K*m_K));
  			p4CoM_pi.SetPxPyPzE(p3_pi.x(), p3_pi.y(), p3_pi.z(), sqrt(p3_pi.scalarProduct(p3_pi)+m_pi*m_pi));

  			//boosting to obtain 4-momenta in the LAB system
  			p4CoM_K.Boost(p4_B.BoostVector()); 
  			p4CoM_pi.Boost(p4_B.BoostVector()); 

  			//3-momenta in the LAB
  			p3LAB_K = p4CoM_K.Vect();
  			p3LAB_pi = p4CoM_pi.Vect();
  			//theta (i.e. polar) angles:
  			theta_K = p3LAB_K.Theta();
  			theta_pi = p3LAB_pi.Theta();
  			Ttheta.push_back(theta_K);
  			Ttheta.push_back(theta_pi);
  			//phi (i.e. azimutal) angles:
  			theta_K = p3LAB_K.Phi();
  			theta_pi = p3LAB_pi.Phi();
  			Tphi.push_back(phi_K);
  			Tphi.push_back(phi_pi);

  			//p3LAB_K = p4LAB_K.Vect();
  			//p3LAB_pi = p4LAB_pi.Vect();
  			p_K_LAB = sqrt(p3LAB_K.Dot(p3LAB_K));
  			p_pi_LAB = sqrt(p3LAB_pi.Dot(p3LAB_pi));
  			Tp.push_back(p_K_LAB);
  			Tp.push_back(p_pi_LAB);
  			//measured norm:
  			p_K_LAB_meas = gen -> Gaus(p_K_LAB, p_K_LAB*resolution);
  			p_pi_LAB_meas = gen -> Gaus(p_pi_LAB, p_pi_LAB*resolution);

  			p3LAB_K *= (p_K_LAB_meas/p_K_LAB);
  			p3LAB_pi *= (p_pi_LAB_meas/p_pi_LAB);

  			//
  			p4LAB_K.SetPxPyPzE(p3LAB_K.X(), p3LAB_K.Y(), p3LAB_K.Z(), sqrt(m_K*m_K + p_K_LAB_meas*p_K_LAB_meas));
  			p4LAB_pi.SetPxPyPzE(p3LAB_pi.X(), p3LAB_pi.Y(), p3LAB_pi.Z(), sqrt(m_pi*m_pi + p_pi_LAB_meas*p_pi_LAB_meas));

  			TLorentzVector p4Meas = p4LAB_K + p4LAB_pi;

  			double inv_mass = sqrt(p4Meas*p4Meas);

  			h_mass->Fill(inv_mass);
  			tree->Fill();

  		}
  		//clearing vectors used for tree:
  		Tnmass.clear();
  		Tp.clear();
  		Ttheta.clear();
  		Tphi.clear();

		histo.push_back(*h_mass);
		tree->Write();

		tree->Print();
		cout << "\n";

  		delete h_mass;
  		delete tree;
  	}

  	histo[0].SetFillColor(kTeal);
  	histo[1].SetFillColor(kOrange);
  	histo[2].SetFillColor(kRed);

  	//true mass hist
  	histo[0].GetXaxis()->SetTitle("\\sqrt{s} [GeV]");
  	histo[0].GetYaxis()->SetTitle("# evts");
  	histo[0].Draw();
  	histo[1].Draw("same");
  	histo[2].Draw("same");

  	//legend
  	TLegend *legend2 = new TLegend(0.1,0.9,0.4,0.8);
   	legend2->SetHeader("","C"); 
   	legend2->AddEntry(&(histo[0]),"\\sigma = 1%");
   	legend2->AddEntry(&(histo[1]),"\\sigma = 5%");
   	legend2->AddEntry(&(histo[2]),"\\sigma = 10%");
   	legend2->Draw();
  	canv.SaveAs("./measured-mass-resolutions.pdf");


  	//Store histograms to file
  	h_inv_mass.Write();
  	h_op_angle.Write();
  	h_inv_mass_meas.Write();
  	histo[0].Write();
  	histo[1].Write();
  	histo[2].Write();
  	legend2->Write();


  	delete gen;
  	delete legend;
  	delete legend2;

  	//closing file
  	rfile.Close();
  	orootfile.Close();

	return 0;
}