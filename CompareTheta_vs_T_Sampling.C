// CompareTheta_vs_T_Sampling.C
// Compare theta-based vs t-based sampling results

void CompareTheta_vs_T_Sampling()
{
    cout << "\n╔════════════════════════════════════════════════════════╗" << endl;
    cout << "║  Comparing Theta-based vs t-based Sampling            ║" << endl;
    cout << "╚════════════════════════════════════════════════════════╝\n" << endl;
    
    // Open both files
    TFile *f_theta = new TFile("pC_Elas_200MeV_MT_P80_SpinUp_16p2_THETA.root", "READ");
    TFile *f_t = new TFile("pC_Elas_200MeV_MT_P80_SpinUp_16p2_T.root", "READ");
    
    if (!f_theta || f_theta->IsZombie()) {
        cerr << "ERROR: Cannot open THETA file!" << endl;
        return;
    }
    if (!f_t || f_t->IsZombie()) {
        cerr << "ERROR: Cannot open T file!" << endl;
        return;
    }
    
    // Get trees
    TTree *tree_theta = (TTree*)f_theta->Get("data");
    TTree *tree_t = (TTree*)f_t->Get("data");
    
    if (!tree_theta || !tree_t) {
        cerr << "ERROR: Cannot find 'data' tree in files!" << endl;
        return;
    }
    
    cout << "Files opened successfully:" << endl;
    cout << "  THETA file: " << tree_theta->GetEntries() << " events" << endl;
    cout << "  T file:     " << tree_t->GetEntries() << " events" << endl;
    cout << endl;
    
    // ================================================================
    // Setup canvas
    // ================================================================
    TCanvas *c1 = new TCanvas("c1", "Theta vs T Sampling Comparison", 1800, 1200);
    c1->Divide(3, 3);
    
    // ================================================================
    // 1. Theta_Lab distribution
    // ================================================================
    c1->cd(1);
    gPad->SetLogy(0);
    
    TH1D *h_theta_lab_theta = new TH1D("h_theta_lab_theta", "Theta_Lab;#theta_{Lab} [deg];Events", 100, 10, 25);
    TH1D *h_theta_lab_t = new TH1D("h_theta_lab_t", "Theta_Lab;#theta_{Lab} [deg];Events", 100, 10, 25);
    
    tree_theta->Draw("Particles.Theta()*TMath::RadToDeg()>>h_theta_lab_theta", "Particles.ID()==14", "goff");
    tree_t->Draw("Particles.Theta()*TMath::RadToDeg()>>h_theta_lab_t", "Particles.ID()==14", "goff");
    
    h_theta_lab_theta->SetLineColor(kBlue);
    h_theta_lab_theta->SetLineWidth(2);
    h_theta_lab_t->SetLineColor(kRed);
    h_theta_lab_t->SetLineWidth(2);
    h_theta_lab_t->SetLineStyle(2);
    
    h_theta_lab_theta->Draw("HIST");
    h_theta_lab_t->Draw("HIST SAME");
    
    TLegend *leg1 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg1->AddEntry(h_theta_lab_theta, "#theta-based", "l");
    leg1->AddEntry(h_theta_lab_t, "t-based", "l");
    leg1->Draw();
    
    // ================================================================
    // 2. Phi_Lab distribution (CRITICAL TEST!)
    // ================================================================
    c1->cd(2);
    gPad->SetLogy(0);
    
    TH1D *h_phi_lab_theta = new TH1D("h_phi_lab_theta", "Phi_Lab;#phi_{Lab} [deg];Events", 100, -180, 180);
    TH1D *h_phi_lab_t = new TH1D("h_phi_lab_t", "Phi_Lab;#phi_{Lab} [deg];Events", 100, -180, 180);
    
    tree_theta->Draw("Particles.Phi()*TMath::RadToDeg()>>h_phi_lab_theta", "Particles.ID()==14", "goff");
    tree_t->Draw("Particles.Phi()*TMath::RadToDeg()>>h_phi_lab_t", "Particles.ID()==14", "goff");
    
    h_phi_lab_theta->SetLineColor(kBlue);
    h_phi_lab_theta->SetLineWidth(2);
    h_phi_lab_t->SetLineColor(kRed);
    h_phi_lab_t->SetLineWidth(2);
    h_phi_lab_t->SetLineStyle(2);
    
    h_phi_lab_theta->Draw("HIST");
    h_phi_lab_t->Draw("HIST SAME");
    
    TLegend *leg2 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg2->AddEntry(h_phi_lab_theta, "#theta-based", "l");
    leg2->AddEntry(h_phi_lab_t, "t-based", "l");
    leg2->Draw();
    
    // Add text showing counts at 0° and 180°
    Double_t count_0_theta = h_phi_lab_theta->Integral(h_phi_lab_theta->FindBin(-10), h_phi_lab_theta->FindBin(10));
    Double_t count_0_t = h_phi_lab_t->Integral(h_phi_lab_t->FindBin(-10), h_phi_lab_t->FindBin(10));
    Double_t count_180_theta = h_phi_lab_theta->Integral(h_phi_lab_theta->FindBin(170), h_phi_lab_theta->FindBin(180));
    count_180_theta += h_phi_lab_theta->Integral(h_phi_lab_theta->FindBin(-180), h_phi_lab_theta->FindBin(-170));
    Double_t count_180_t = h_phi_lab_t->Integral(h_phi_lab_t->FindBin(170), h_phi_lab_t->FindBin(180));
    count_180_t += h_phi_lab_t->Integral(h_phi_lab_t->FindBin(-180), h_phi_lab_t->FindBin(-170));
    
    TLatex *text_phi = new TLatex();
    text_phi->SetNDC();
    text_phi->SetTextSize(0.03);
    text_phi->DrawLatex(0.15, 0.85, Form("#theta: 0#circ = %.0f, 180#circ = %.0f", count_0_theta, count_180_theta));
    text_phi->DrawLatex(0.15, 0.80, Form("t: 0#circ = %.0f, 180#circ = %.0f", count_0_t, count_180_t));
    
    // ================================================================
    // 3. Proton momentum
    // ================================================================
    c1->cd(3);
    gPad->SetLogy(0);
    
    TH1D *h_p_theta = new TH1D("h_p_theta", "Proton Momentum;p [GeV/c];Events", 100, 0.3, 0.8);
    TH1D *h_p_t = new TH1D("h_p_t", "Proton Momentum;p [GeV/c];Events", 100, 0.3, 0.8);
    
    tree_theta->Draw("Particles.P()>>h_p_theta", "Particles.ID()==14", "goff");
    tree_t->Draw("Particles.P()>>h_p_t", "Particles.ID()==14", "goff");
    
    h_p_theta->SetLineColor(kBlue);
    h_p_theta->SetLineWidth(2);
    h_p_t->SetLineColor(kRed);
    h_p_t->SetLineWidth(2);
    h_p_t->SetLineStyle(2);
    
    h_p_theta->Draw("HIST");
    h_p_t->Draw("HIST SAME");
    
    // ================================================================
    // 4. Proton Pt (transverse momentum)
    // ================================================================
    c1->cd(4);
    gPad->SetLogy(0);
    
    TH1D *h_pt_theta = new TH1D("h_pt_theta", "Proton p_{T};p_{T} [GeV/c];Events", 100, 0, 0.3);
    TH1D *h_pt_t = new TH1D("h_pt_t", "Proton p_{T};p_{T} [GeV/c];Events", 100, 0, 0.3);
    
    tree_theta->Draw("Particles.Perp()>>h_pt_theta", "Particles.ID()==14", "goff");
    tree_t->Draw("Particles.Perp()>>h_pt_t", "Particles.ID()==14", "goff");
    
    h_pt_theta->SetLineColor(kBlue);
    h_pt_theta->SetLineWidth(2);
    h_pt_t->SetLineColor(kRed);
    h_pt_t->SetLineWidth(2);
    h_pt_t->SetLineStyle(2);
    
    h_pt_theta->Draw("HIST");
    h_pt_t->Draw("HIST SAME");
    
    // ================================================================
    // 5. Carbon Theta_Lab
    // ================================================================
    c1->cd(5);
    gPad->SetLogy(0);
    
    TH1D *h_theta_c_theta = new TH1D("h_theta_c_theta", "Carbon Theta_Lab;#theta_{Lab}^{C} [deg];Events", 200, 70, 90);
    TH1D *h_theta_c_t = new TH1D("h_theta_c_t", "Carbon Theta_Lab;#theta_{Lab}^{C} [deg];Events", 200, 70, 90);
    
    tree_theta->Draw("Particles.Theta()*TMath::RadToDeg()>>h_theta_c_theta", "Particles.ID()==614", "goff");
    tree_t->Draw("Particles.Theta()*TMath::RadToDeg()>>h_theta_c_t", "Particles.ID()==614", "goff");
    
    h_theta_c_theta->SetLineColor(kBlue);
    h_theta_c_theta->SetLineWidth(2);
    h_theta_c_t->SetLineColor(kRed);
    h_theta_c_t->SetLineWidth(2);
    h_theta_c_t->SetLineStyle(2);
    
    h_theta_c_theta->Draw("HIST");
    h_theta_c_t->Draw("HIST SAME");
    
    // ================================================================
    // 6. 2D: Theta vs Phi (THETA-based)
    // ================================================================
    c1->cd(6);
    
    TH2D *h2_theta = new TH2D("h2_theta", "#theta-based: #theta vs #phi;#phi_{Lab} [deg];#theta_{Lab} [deg]", 
                              50, -180, 180, 50, 10, 25);
    tree_theta->Draw("Particles.Theta()*TMath::RadToDeg():Particles.Phi()*TMath::RadToDeg()>>h2_theta", 
                     "Particles.ID()==14", "COLZ");
    
    // ================================================================
    // 7. 2D: Theta vs Phi (T-based)
    // ================================================================
    c1->cd(7);
    
    TH2D *h2_t = new TH2D("h2_t", "t-based: #theta vs #phi;#phi_{Lab} [deg];#theta_{Lab} [deg]", 
                          50, -180, 180, 50, 10, 25);
    tree_t->Draw("Particles.Theta()*TMath::RadToDeg():Particles.Phi()*TMath::RadToDeg()>>h2_t", 
                 "Particles.ID()==14", "COLZ");
    
    // ================================================================
    // 8. Asymmetry comparison (RIGHT - LEFT) / (RIGHT + LEFT)
    // ================================================================
    c1->cd(8);
    gPad->SetGridy();
    
    TH1D *h_asym_theta = new TH1D("h_asym_theta", "Azimuthal Asymmetry;Bin;(Right - Left) / (Right + Left)", 10, 0, 10);
    TH1D *h_asym_t = new TH1D("h_asym_t", "Azimuthal Asymmetry;Bin;(Right - Left) / (Right + Left)", 10, 0, 10);
    
    Double_t asym_theta = (count_0_theta - count_180_theta) / (count_0_theta + count_180_theta);
    Double_t asym_t = (count_0_t - count_180_t) / (count_0_t + count_180_t);
    
    h_asym_theta->SetBinContent(1, asym_theta);
    h_asym_t->SetBinContent(1, asym_t);
    
    h_asym_theta->SetLineColor(kBlue);
    h_asym_theta->SetMarkerColor(kBlue);
    h_asym_theta->SetMarkerStyle(20);
    h_asym_theta->SetMarkerSize(2);
    h_asym_t->SetLineColor(kRed);
    h_asym_t->SetMarkerColor(kRed);
    h_asym_t->SetMarkerStyle(21);
    h_asym_t->SetMarkerSize(2);
    
    h_asym_theta->SetMinimum(-1);
    h_asym_theta->SetMaximum(1);
    h_asym_theta->Draw("P");
    h_asym_t->Draw("P SAME");
    
    TLine *line_zero = new TLine(0, 0, 10, 0);
    line_zero->SetLineStyle(2);
    line_zero->Draw();
    
    TLegend *leg8 = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg8->AddEntry(h_asym_theta, Form("#theta: %.3f", asym_theta), "p");
    leg8->AddEntry(h_asym_t, Form("t: %.3f", asym_t), "p");
    leg8->Draw();
    
    // ================================================================
    // 9. Statistical Summary
    // ================================================================
    c1->cd(9);
    gPad->SetFillColor(kWhite);
    
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.04);
    
    pt->AddText("Statistical Comparison");
    pt->AddText(" ");
    pt->AddText(Form("Events: #theta=%d, t=%d", (int)tree_theta->GetEntries(), (int)tree_t->GetEntries()));
    pt->AddText(" ");
    pt->AddText("Phi distribution (±10° bins):");
    pt->AddText(Form("  0° (#theta): %.0f events", count_0_theta));
    pt->AddText(Form("  0° (t):      %.0f events", count_0_t));
    pt->AddText(Form("  180° (#theta): %.0f events", count_180_theta));
    pt->AddText(Form("  180° (t):      %.0f events", count_180_t));
    pt->AddText(" ");
    pt->AddText(Form("Asymmetry A = (R-L)/(R+L):"));
    pt->AddText(Form("  #theta-based: %.3f", asym_theta));
    pt->AddText(Form("  t-based:      %.3f", asym_t));
    pt->AddText(" ");
    
    if (TMath::Abs(asym_theta - asym_t) < 0.1) {
        pt->AddText("Result: ✓ GOOD AGREEMENT");
        pt->GetLine(pt->GetSize()-1)->SetTextColor(kGreen+2);
    } else {
        pt->AddText("Result: ✗ DISAGREEMENT!");
        pt->GetLine(pt->GetSize()-1)->SetTextColor(kRed);
    }
    
    pt->Draw();
    
    // ================================================================
    // Update canvas and print summary
    // ================================================================
    c1->Update();
    
    cout << "\n════════════════════════════════════════════════════════" << endl;
    cout << "SUMMARY:" << endl;
    cout << "════════════════════════════════════════════════════════" << endl;
    cout << "Events:" << endl;
    cout << "  Theta-based: " << tree_theta->GetEntries() << endl;
    cout << "  t-based:     " << tree_t->GetEntries() << endl;
    cout << "\nPhi Distribution (0° = right, 180° = left):" << endl;
    cout << "  Theta-based: 0° = " << count_0_theta << ", 180° = " << count_180_theta << endl;
    cout << "  t-based:     0° = " << count_0_t << ", 180° = " << count_180_t << endl;
    cout << "\nAsymmetry (Right - Left) / (Right + Left):" << endl;
    cout << "  Theta-based: " << asym_theta << endl;
    cout << "  t-based:     " << asym_t << endl;
    cout << "\nDifference: " << TMath::Abs(asym_theta - asym_t) << endl;
    
    if (TMath::Abs(asym_theta - asym_t) < 0.1) {
        cout << "\n✓ DISTRIBUTIONS AGREE!" << endl;
    } else {
        cout << "\n✗ WARNING: Distributions differ significantly!" << endl;
        cout << "  Check azimuthal angle sampling in t-based mode." << endl;
    }
    cout << "════════════════════════════════════════════════════════\n" << endl;
    
    // Save canvas
    c1->SaveAs("Comparison_Theta_vs_T_Sampling.pdf");
    cout << "Plot saved as: Comparison_Theta_vs_T_Sampling.pdf" << endl;
}