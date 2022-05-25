{
  TFile fileOutput1("Kaon_xi_Output_test.root","recreate");
  TFile file1("KL_flux_p_Pythia_24m.root");
  TFile file2("Maya_KLn_KPXi.root");
  gStyle->SetOptStat(0);
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// creating histograms
  // angular dependence histograms
  TH2D* h_kaon_plus=new TH2D("h_kaon_plus","Angular dependence of K^{+}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,200);                                    //k-plus
  TH2D* h_xi=new TH2D("h_xi","Angular dependence of #Xi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,15,200,0,180);                                                //xi
  TH2D* h_resc_xi=new TH2D("h_resc_xi","Angular dependence of rescattered #Xi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,180);                          //rescattered xi
  TH2D* h_resc_proton=new TH2D("h_resc_proton","Angular dependence of first rescattered proton by #Xi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,180);  //rescattered proton

  TH2D* h_lambda=new TH2D("h_lambda","Angular dependence of #Xi^{-}-decay produced #Lambda^{0}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,180);
  TH2D* h_pi=new TH2D("h_pi","Angular dependence of #Xi^{-}-decay produced #pi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,5,200,0,200);
  TH2D* h_lambda_proton=new TH2D("h_lambda_proton","Angular dependence of #Lambda^{0}-decay produced proton; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,180);
  TH2D* h_lambda_pion=new TH2D("h_lambda_pion","Angular dependence of #Lambda^{0}-decay produced #pi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,5,200,0,200);

  // momentum histograms
  TH1F* h_momentum_xi = new TH1F("h_momentum_xi","Momentum of #Xi^{-}; Momentum p [GeV/c]; Counts",200,0,10);
  TH1F* h_momentum_resc_xi = new TH1F("h_momentum_resc_xi","Momentum of rescattered #Xi^{-}; Momentum p [GeV/c]; Counts",200,0,10);
  TH1F* h24p = (TH1F*)file1.Get("h24p");
  TH1F* h_beamflux = new TH1F("h_beamflux","Momentum of beam w.r.t flux; Momentum p [GeV/c]; Events per second",12000,0,12);

  TH1F* h_kl_momentum = new TH1F("h_kl_momentum","Momentum of K-Long; Momentum p [GeV/c]; Counts",12000,0,12);
  TH1F* h_kl_momentumacceptance1 = new TH1F("h_kl_momentumacceptance1","Acceptance 1; Momentum p [GeV/c]; Counts",12000,0,12);
  TH1F* h_kl_momentumacceptance2 = new TH1F("h_kl_momentumacceptance2","Acceptance 2; Momentum p [GeV/c]; Counts",12000,0,12);
  TH1F* h_kl_momentumacceptance = new TH1F("h_kl_momentumacceptance","Acceptance of Momentum of K-Long; Momentum p [GeV/c]; Accepted Momentum",12000,0,12);

  // scattering events histograms
  TH1F* h_detected_xi = new TH1F("h_detected_xi","Number of detected #Xi^{-} scattering events per 100 days per GeV of Momentum; Momentum p [GeV/c]; Normalised Yield",12000,0,12);
  TH1F* h_scattered_xi = new TH1F("h_scattered_xi","Number of scattered #Xi^{-} scattering events per 100 days per GeV of Momentum; Momentum p [GeV/c]; Normalised Yield",12000,0,12);
  TH1F* hKPXi = (TH1F*)file2.Get("hKPXi");
  TH1F* hKPXi1 = new TH1F("hKPXi1", "K_L + p --> #Xi^{-}+ K^{+}",12000,0,12);
        // cross section scaling
        for (int j = 1; j < 12001; j++) {
          hKPXi1->SetBinContent(j,hKPXi->GetBinContent(j));
        }

  // path length plots
  TH2D* h_l_xy = new TH2D("h_l_xy","#phi of #Xi^{-} vs X-Y component of path length; #phi [radians]; length [cm]",200,-4,4,200,0,10);
  TH2D* h_l_z = new TH2D("h_l_z","#theta of #Xi^{-} vs Z component of path length; #theta [degrees]; length [cm]",200,0,50,200,0,50);
  TH2D* h_pathlength = new TH2D("h_pathlength","#theta of #Xi^{-} vs path length; #theta [degrees]; length [cm]",200,0,50,200,0,50);
  TH1F* h_filledpathlength = new TH1F("h_filledpathlength","Pathlength filled w.r.t Momentum; length [cm]; p [MeV/c]", 450,0,45);
  TF1* mydecay_plot = new TF1("mydecay_plot","TMath::Exp(-x/([0]*[1]*30*0.1639)); path length [cm]", 0, 100);  //0.1639 = mean lifetime of xi, ns
  TH2D* h_pathlengthdecay = new TH2D("h_pathlengthdecay","#theta of #Xi^{-} decay vs path length; #theta [degrees]; length [cm]",200,0,50,200,0,50);

  // constraint histograms
  TH2D* h_constrained_kaon_plus = new TH2D("h_constrained_kaon_plus","Constrained K^{+}",200,0,10,200,0,200);
  TH2D* h_constrained_resc_proton = new TH2D("h_constrained_resc_proton","Constrained proton", 200,0,10,200,0,180);
  TH2D* h_constrained_pi = new TH2D("h_constrained_pi","Constrained pion",200,0,5,200,0,200);
  TH2D* h_constrained_lambda_pion = new TH2D("h_constrained_lambda_pion","Constrained lambda pion",200,0,5,200,0,200);
  TH2D* h_constrained_lambda_proton = new TH2D("h_constrained_lambda_proton","Constrained lambda proton",200,0,10,200,0,180);

  // acceptance histograms
  TH2D* h_acceptance_kaon_plus = new TH2D("h_acceptance_kaon_plus","Acceptance of K^{+}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,200);
  TH2D* h_acceptance_resc_proton = new TH2D("h_acceptance_resc_proton","Acceptance of the rescattered proton; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,180);
  TH2D* h_acceptance_pi = new TH2D("h_acceptance_pi","Acceptance of #Xi^{-}-produced #pi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,5,200,0,200);
  TH2D* h_acceptance_lambda_pion = new TH2D("h_acceptance_lambda_pion","Acceptance of #Lambda^{0}-produced #pi^{-}; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,5,200,0,200);
  TH2D* h_acceptance_lambda_proton = new TH2D("h_acceptance_lambda_proton","Acceptance of #Lambda^{0}-produced proton; Momentum p [GeV/c]; Angle [#theta, degrees]",200,0,10,200,0,180);

  // rebinning
  h24p->Rebin(100);
  hKPXi1->Rebin(100);
  h_kl_momentum->Rebin(100);
  h_kl_momentumacceptance1->Rebin(100);
  h_kl_momentumacceptance2->Rebin(100);
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//// creating particles and terms

  // how many events to simulate and percentage completed
  Int_t nevents=250000;
  Int_t Percentage=nevents/100;

  // creating TLorentzVectors of particles
  TLorentzVector Beam, Target, Target_Proton; // beam and target
  TLorentzVector *xi, *kaon_plus; // first vertex particles
  TLorentzVector *resc_xi, *resc_proton; // second vertex particles

  TLorentzVector *lambda, *pi; // third vertex
  TLorentzVector *lambda_proton, *lambda_pion; // fourth vertex

  // setting TLorentzVectors for beam and target in GeV (Px,Py,Pz,M)
  Beam.SetXYZM(0,0,3.5,0.497); //Pz is 3.5, M is 0.497
  Target.SetXYZM(0,0,0,0.939); //M is 0.939
  Target_Proton.SetXYZM(0,0,0,0.938);

  // defining initial vertex and masses of particles
  /*
  mass of xi = 1.32171
  mass of k+ = 0.49368
  mass of proton = 0.93827
  mass of pi- = 0.13957
  mass of lambda = 1.115683
  */

  TLorentzVector V1 = Beam + Target;
  Double_t Masses_1[2] = {1.32171,0.49368};   // xi-, K+ (primary vertex)
  Double_t Masses_2[2] = {1.32171,0.93827};   // xi-, p (xi- interaction with proton)
  Double_t Masses_3[2] = {1.115683,0.13957};  // lambda, pi- (xi decay)
  Double_t Masses_4[2] = {0.93827,0.13957};   // proton, pi2- (lambda decay)

  // creating decay vertices
  TGenPhaseSpace Vertex_1, Vertex_2, Vertex_3, Vertex_4;

  // beam
  Double_t beamrandom;

  // target
  Double_t targetlength = 40; //in cm
  Double_t targetradius = 3;  //in cm

  // path length and components
  Double_t pathlength;
  Double_t l_xy;              //path length x-y component
  Double_t l_z;               //path length z component

  // angle components
  Double_t phi_v;
  Double_t phi_xi;
  Double_t theta_xi;
  Double_t alpha;             //= 180 - phi_v + phi_xi

  // vertex components
  Double_t x_v;
  Double_t y_v;
  Double_t z_v;
  Double_t r_v;               //radius from centre of target to vertex position
  Double_t z_exit;            //z vertex of exit position

  // decay components
  Double_t xi_beta;
  Double_t xi_gamma;

  Double_t decaypath;
  Double_t pathlengthdecay;   //variable dependent on decaypath

  TRandom3 *myrndm = new TRandom3();

  // luminosity and cross section components of beam momentum
  h_beamflux->Rebin(100);
  h_beamflux->Multiply(h24p,hKPXi1);
  h_beamflux->Scale((6.022E-9/2.014)*0.1624*targetlength);
  //               ((avogadro/mol.mass)*density of deuterium*length)
  cout<<"Total number of Xi produced over 100 days: "<< h_beamflux->Integral()*100*24*60*60 << endl;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//=========================================================================================================================================================================================================
//// looping over simulated events
//// start of loop
  for (Long64_t i=0;i<nevents;i++){
    if (i % Percentage == 0){
      fprintf (stderr, "%lld\r", i/Percentage); //shows percentage of simulated events completed
      fflush (stderr);
    }

    beamrandom = h24p->GetRandom();
    Beam.SetXYZM(0,0,beamrandom,0.497);
    h_kl_momentum->Fill(beamrandom);
    V1 = Beam + Target;

    /*
      creating the vertices for each branch of the interaction.
      Kl + n --> kaon_plus + Xi-                                   (vertex 1)
                    elastic scattering, --> resc_xi + resc_proton; (vertex 2)
                        xi decays --> Lambda + pi                  (vertex 3)
                            lambda decays --> p + pi               (vertex 4)
    */
    // first decay
    if (!Vertex_1.SetDecay(V1,2,Masses_1)) continue;
    //                    (total energy, no. particles, mass array)
    Double_t Phasespace_Weight_1 = Vertex_1.Generate();               // generating event and defining the phasespace weight for first vertex
    xi = Vertex_1.GetDecay(0);                                        // assigning the decay particles from array above
    kaon_plus = Vertex_1.GetDecay(1);

    // angular dependence of second vertex of scattered xi and proton
    TLorentzVector V2 = (TLorentzVector)*xi + Target_Proton;
    Vertex_2.SetDecay(V2,2,Masses_2);
    Double_t Phasespace_Weight_2 = Vertex_2.Generate();               // defining the phasespace for this scattering
    resc_xi = Vertex_2.GetDecay(0);
    resc_proton = Vertex_2.GetDecay(1);
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // fill momentum histogram
    h_momentum_xi->Fill(xi->Rho(),Phasespace_Weight_1*Phasespace_Weight_2);
    h_momentum_resc_xi->Fill(resc_xi->Rho(),Phasespace_Weight_1*Phasespace_Weight_2);
//-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // angular dependence of first decay
    TLorentzVector V3 = (TLorentzVector)*resc_xi;
    Vertex_3.SetDecay(V3,2,Masses_3);
    Double_t Phasespace_Weight_3 = Vertex_3.Generate();
    lambda = Vertex_3.GetDecay(0);
    pi = Vertex_3.GetDecay(1);

    // angular dependence of second decay
    TLorentzVector V4 = (TLorentzVector)*lambda;
    Vertex_4.SetDecay(V4,2,Masses_4);
    Double_t Phasespace_Weight_4 = Vertex_4.Generate();
    lambda_proton = Vertex_4.GetDecay(0);
    lambda_pion = Vertex_4.GetDecay(1);
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // path length calculations
    r_v = 5; //has to be smaller than the diameter of the target
    while(r_v>3){
      x_v = myrndm->Gaus(0,1); //in cm
      y_v = myrndm->Gaus(0,1); //in cm
      r_v=sqrt(x_v*x_v + y_v*y_v);
    }

    phi_v = TMath::ATan2(y_v,x_v);
    phi_xi = xi->Phi();
    theta_xi = xi->Theta();
    z_v = myrndm->Rndm()*targetlength;
    alpha = (TMath::ATan(1)*4) - phi_v + phi_xi;

    l_xy = (r_v*TMath::Cos(alpha)) + sqrt((r_v*r_v*(TMath::Cos(alpha))*(TMath::Cos(alpha))) + targetradius*targetradius - r_v*r_v);
    l_z = l_xy/TMath::Tan(theta_xi);
    z_exit = l_z + z_v;

    if (z_exit > targetlength){
      l_z = targetlength - z_v;
      l_xy = l_z*TMath::Tan(theta_xi);
      pathlength = l_z/TMath::Cos(theta_xi);
    }
    else {
      pathlength = l_xy/TMath::Sin(theta_xi);
    }

    // setting the decay constants for xi
    xi_beta = xi->Beta();
    xi_gamma = xi->Gamma();
    mydecay_plot->SetParameter(0, xi_beta);
    mydecay_plot->SetParameter(1, xi_gamma);
    decaypath = mydecay_plot->GetRandom();

    // plotting path length decay lengths
    if (decaypath < pathlength){
      pathlengthdecay = decaypath;
    }
    else {
      pathlengthdecay = pathlength;
    }
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // filling scattering and cross section histograms
    h_filledpathlength->Fill(pathlengthdecay);
    h_kl_momentumacceptance2->Fill(beamrandom,Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);
//---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    // constraints
    if (kaon_plus->Rho() > 0.2 && resc_proton->Rho() > 0.2 && pi->Rho() > 0.2 && lambda_pion->Rho() > 0.2 && lambda_proton->Rho() > 0.2){ //main if loop

      // if statement for each decay vertex particle
      if ((kaon_plus->Theta() > 3*TMath::DegToRad() && kaon_plus->Theta() < 15*TMath::DegToRad())
      || (kaon_plus->Theta() > 20*TMath::DegToRad() && kaon_plus->Theta() < 165*TMath::DegToRad())){              //1

        if ((resc_proton->Theta() > 3*TMath::DegToRad() && resc_proton->Theta() < 15*TMath::DegToRad())
        || (resc_proton->Theta() > 20*TMath::DegToRad() && resc_proton->Theta() < 165*TMath::DegToRad())){        //2

          if ((pi->Theta() > 3*TMath::DegToRad() && pi->Theta() < 15*TMath::DegToRad())
          || (pi->Theta() > 20*TMath::DegToRad() && pi->Theta() < 165*TMath::DegToRad())){                        //3

            if ((lambda_pion > 3*TMath::DegToRad() && lambda_pion->Theta() < 15*TMath::DegToRad())
            || (lambda_pion > 20*TMath::DegToRad() && lambda_pion->Theta() < 165*TMath::DegToRad())){             //4

              if ((lambda_proton > 3*TMath::DegToRad() && lambda_proton->Theta() < 15*TMath::DegToRad())
              || (lambda_proton > 20*TMath::DegToRad() && lambda_proton->Theta() < 165*TMath::DegToRad())){       //5

                h_constrained_kaon_plus->Fill(kaon_plus->Rho(),kaon_plus->Theta()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);
                h_constrained_resc_proton->Fill(resc_proton->Rho(),resc_proton->Theta()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);
                h_constrained_pi->Fill(pi->Rho(),pi->Theta()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);
                h_constrained_lambda_pion->Fill(lambda_pion->Rho(),lambda_pion->Theta()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);
                h_constrained_lambda_proton->Fill(lambda_proton->Rho(),lambda_proton->Theta()*TMath::RadToDeg(),Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);

                h_kl_momentumacceptance1->Fill(beamrandom,Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);

              } //5
            }   //4
          }     //3
        }       //2
      }         //1
    } // end of large if statement
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// filling histograms
// initial vertex
  h_kaon_plus->Fill(kaon_plus->Rho(),kaon_plus->Theta()*TMath::RadToDeg(), Phasespace_Weight_1);
  h_xi->Fill(xi->Rho(),xi->Theta()*TMath::RadToDeg(),Phasespace_Weight_1);

// first scattering
  h_resc_xi->Fill(resc_xi->Rho(),resc_xi->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);
  h_resc_proton->Fill(resc_proton->Rho(),resc_proton->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2);

// first decay vertex
  h_pi->Fill(pi->Rho(),pi->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3);
  h_lambda->Fill(lambda->Rho(),lambda->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3);

// second decay vertex
  h_lambda_pion->Fill(lambda_pion->Rho(),lambda_pion->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);
  h_lambda_proton->Fill(lambda_proton->Rho(),lambda_proton->Theta()*TMath::RadToDeg(), Phasespace_Weight_1*Phasespace_Weight_2*Phasespace_Weight_3*Phasespace_Weight_4);

// path lengths
  h_l_xy->Fill(phi_xi,l_xy);
  h_l_z->Fill(theta_xi*TMath::RadToDeg(),l_z);
  h_pathlength->Fill(theta_xi*TMath::RadToDeg(),pathlength);
  h_pathlengthdecay->Fill(theta_xi*TMath::RadToDeg(),pathlengthdecay);
  }
// end of loop
//====================================================================================================================================================================================================================
// acceptance histograms
  h_acceptance_kaon_plus->Divide(h_constrained_kaon_plus,h_kaon_plus);
  h_acceptance_resc_proton->Divide(h_constrained_resc_proton,h_resc_proton);
  h_acceptance_pi->Divide(h_constrained_pi,h_pi);
  h_acceptance_lambda_pion->Divide(h_constrained_lambda_pion,h_lambda_pion);
  h_acceptance_lambda_proton->Divide(h_constrained_lambda_proton,h_lambda_proton);

  h_kl_momentumacceptance->Rebin(100);
  h_kl_momentumacceptance->Divide(h_kl_momentumacceptance1,h_kl_momentumacceptance2);

// number of scattered xi particles
  h_scattered_xi->Rebin(100);
  h_scattered_xi = (TH1F*) h_beamflux->Clone("h_scattered_xi");
  h_scattered_xi->Scale((h_filledpathlength->GetMean())*(6.022E-4/2.014)*0.1624*5.25);
  cout<<"Total number of predicted scattered Xi over 100 days: "<< h_scattered_xi->Integral()*100*24*60*60 << endl;

  h_detected_xi->Rebin(100);
  h_detected_xi->Multiply(h_beamflux,h_kl_momentumacceptance);
  h_detected_xi->Scale((h_filledpathlength->GetMean())*(6.022E-4/2.014)*0.1624*5.25);
  cout<<"Total number of detected scattered Xi over 100 days: "<< h_detected_xi->Integral()*100*24*60*60 << endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// drawing plots
//-----------------------------------------------------------------------------------------
TCanvas *c1 = new TCanvas("c1","Angular Dependence of Initial Scattering",800,800);
c1->Divide(2,2);
TCanvas *c2 = new TCanvas("c2","Angular Dependence of Scattered Particle Decay",800,800);
c2->Divide(2,2);
TCanvas *c3 = new TCanvas("c3","Momentum Plots",800,800);
c3->Divide(2,2);
TCanvas *c4 = new TCanvas("c4","Path Length Plots and Decays",1600,800);
c4->Divide(3,2);
TCanvas *c5 = new TCanvas("c5","Acceptance Histograms",1600,800);
c5->Divide(3,2);
TCanvas *c6 = new TCanvas("c6","Number of Scattering Events",800,800);
c6->Divide(1,1);
//------------------------------------------------------------------------------------------
// assigning plot position on print
//------------------------------------------------------------------------------------------
c1->cd(1);
h_xi->Draw("colz");
c1->cd(2);
h_kaon_plus->Draw("colz");
c1->cd(3);
h_resc_xi->Draw("colz");
c1->cd(4);
h_resc_proton->Draw("colz");
//------------------------------------------------------------------------------------------
c2->cd(1);
h_pi->Draw("colz");
c2->cd(2);
h_lambda->Draw("colz");
c2->cd(3);
h_lambda_proton->Draw("colz");
c2->cd(4);
h_lambda_pion->Draw("colz");
//------------------------------------------------------------------------------------------
c3->cd(1);
h_momentum_xi->Draw();
c3->cd(2);
h_momentum_resc_xi->Draw();
c3->cd(3);
h_kl_momentum->Draw();
c3->cd(4);
h_beamflux->Draw();
//------------------------------------------------------------------------------------------
c4->cd(1);
h_l_xy->Draw("colz");
c4->cd(2);
h_l_z->Draw("colz");
c4->cd(3);
h_pathlength->Draw("colz");
c4->cd(4);
h_pathlengthdecay->Draw("colz");
c4->cd(5);
mydecay_plot->Draw();
c4->cd(6);
h_filledpathlength->Draw();
//------------------------------------------------------------------------------------------
c5->cd(1);
h_acceptance_kaon_plus->Draw("colz");
c5->cd(2);
h_acceptance_resc_proton->Draw("colz");
c5->cd(3);
h_acceptance_pi->Draw("colz");
c5->cd(4);
h_acceptance_lambda_pion->Draw("colz");
c5->cd(5);
h_acceptance_lambda_proton->Draw("colz");
c5->cd(6);
h_kl_momentumacceptance->Draw();
//------------------------------------------------------------------------------------------
c6->cd(1);
h_detected_xi->Draw();
} // end of code
