//B --> K+pi

#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1F.h" 
#include "TCanvas.h" 
#include "TFile.h"

//Vector3D class with methods to extract angles ecc .. 
#include "Vector3D.h"

#include <iostream>
#include <cmath>

using std::cout;
using std::endl;

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
	double xlo1 = m_B - 0.00005;
  	double xhi1 = m_B + 0.00005;
  	double binwidth1 = (xhi1-xlo1) / nbins1;
  	cout << "# bins hist 1: " << nbins1 << "\t bin width: " << binwidth1 << endl;

  	// Create new histogram (this one for the angle between the momenta)
	int nbins2 = 100; 
	double xlo2 = -3*M_PI;
  	double xhi2 = 3*M_PI;
  	double binwidth2 = (xhi2-xlo2) / nbins2;
  	cout << "# bins hist 2: " << nbins2 << "\t bin width: " << binwidth2 << endl;

  	//invariant mass histogram
  	TH1F h_inv_mass("\\sqrt(s)", "distribution of invariant mass", nbins1, xlo1, xhi1);
  	//opening angle histogram
  	TH1F h_op_angle("\\phi", "distribution of opening angles", nbins2, xlo2, xhi2);


	//B meson 4-momentum in the LAB frame
  	TLorentzVector p4_B; 
  	double p_B = 0.300; // GeV
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

  	double p_star = p_CoM(m_B, m_K, m_pi); //GeV
  	double sqrt_s = 0.;
  	double phi = 0.; //opening angle between the products in the LAB

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

  		TLorentzVector p4_tot = p4CoM_pi + p4CoM_K; //total 4-momentum in the LAB

  		//invariant mass 
  		sqrt_s = sqrt(p4_tot*p4_tot);
  		//opening angle
  		phi = p3_pi.angle(p3_K); //LO DEVI FARE CON I VETTORI BOOSTATI!!!...

  		//fill histograms
  		h_inv_mass.Fill(sqrt_s);
  		h_op_angle.Fill(phi);
  	}


  	//plotting the results:
  	// * create canvas
  	TCanvas canv("canv", "canvas for plotting", 1280, 1024);
  	h_inv_mass.GetXaxis()->SetTitle("\\sqrt{s} [GeV]");
  	// * plot
  	h_inv_mass.Draw();
  	// * store to file in 2 formats
  	canv.SaveAs("./true-mass.pdf");

  	h_op_angle.GetXaxis()->SetTitle("opening angle \\phi [rad]");
  	h_op_angle.Draw();
  	canv.SaveAs("./opening-angle.pdf");



  	//Store histograms to file
  	h_inv_mass.Write();
  	h_op_angle.Write();

  	delete gen;

  	rfile.Close();

	return 0;
}