{
  TFile* file = new TFile("spectra.root");
  TCanvas *c1 = new TCanvas;
  c1->Divide(3,2);

  c1->cd(1);
  s1_r_matrix->Draw("c");
  s1_r_matrix->GetYaxis()->SetTitle("Differential Cross Section [b/sr]");
  s1_r_matrix->GetYaxis()->CenterTitle();
  s1_r_matrix->GetXaxis()->CenterTitle();
  s1_r_matrix->GetYaxis()->SetTitleOffset(1.4);
  s1_r_matrix->SetLineColor(kRed);
  s1_r_matrix->SetLineWidth(2);
  s1->Draw("same e1");
  s1->SetMarkerStyle(kFullCircle);
  s1->SetMarkerSize(0.6);
  s1->SetLineColor(kBlack);
  
  c1->cd(2);
  s2_r_matrix->Draw("c");
  s2_r_matrix->GetYaxis()->SetTitle("Differential Cross Section [b/sr]");
  s2_r_matrix->GetYaxis()->CenterTitle();
  s2_r_matrix->GetXaxis()->CenterTitle();
  s2_r_matrix->GetYaxis()->SetTitleOffset(1.4);
  s2_r_matrix->SetLineColor(kRed);
  s2_r_matrix->SetLineWidth(2);
  s2->Draw("same e1");
  s2->SetMarkerStyle(kFullCircle);
  s2->SetMarkerSize(0.6);
  s2->SetLineColor(kBlack);
  

  c1->cd(3);
  s3_r_matrix->Draw("c");
  s3_r_matrix->GetYaxis()->SetTitle("Differential Cross Section [b/sr]");
  s3_r_matrix->GetYaxis()->CenterTitle();
  s3_r_matrix->GetXaxis()->CenterTitle();
  s3_r_matrix->GetYaxis()->SetTitleOffset(1.4);
  s3_r_matrix->SetLineColor(kRed);
  s3_r_matrix->SetLineWidth(2);
  s3->Draw("same e1");
  s3->SetMarkerStyle(kFullCircle);
  s3->SetMarkerSize(0.6);
  s3->SetLineColor(kBlack);
  

  c1->cd(4);
  s4_r_matrix->Draw("c");
  s4_r_matrix->GetYaxis()->SetTitle("Differential Cross Section [b/sr]");
  s4_r_matrix->GetYaxis()->CenterTitle();
  s4_r_matrix->GetXaxis()->CenterTitle();
  s4_r_matrix->GetYaxis()->SetTitleOffset(1.4);
  s4_r_matrix->SetLineColor(kRed);
  s4_r_matrix->SetLineWidth(2);
  s4->Draw("same e1");
  s4->SetMarkerStyle(kFullCircle);
  s4->SetMarkerSize(0.6);
  s4->SetLineColor(kBlack);
 

  c1->cd(5);
  s5_r_matrix->Draw("c");
  s5_r_matrix->GetYaxis()->SetTitle("Differential Cross Section [b/sr]");
  s5_r_matrix->GetYaxis()->CenterTitle();
  s5_r_matrix->GetXaxis()->CenterTitle();
  s5_r_matrix->GetYaxis()->SetTitleOffset(1.4);
  s5_r_matrix->SetLineColor(kRed);
  s5_r_matrix->SetLineWidth(2);
  s5->Draw("same e1");
  s5->SetMarkerStyle(kFullCircle);
  s5->SetMarkerSize(0.6);
  s5->SetLineColor(kBlack);
  

  c1->cd(6);
  s6_r_matrix->Draw("c");
  s6_r_matrix->GetYaxis()->SetTitle("Differential Cross Section [b/sr]");
  s6_r_matrix->GetYaxis()->CenterTitle();
  s6_r_matrix->GetXaxis()->CenterTitle();
  s6_r_matrix->GetYaxis()->SetTitleOffset(1.4);
  s6_r_matrix->SetLineColor(kRed);
  s6_r_matrix->SetLineWidth(2);
  s6->Draw("same e1");
  s6->SetMarkerStyle(kFullCircle);
  s6->SetMarkerSize(0.6);
  s6->SetLineColor(kBlack);
  
}
