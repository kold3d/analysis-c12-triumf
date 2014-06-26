void GetMeanWirePos(){
	pos_1->Fit("gaus");
	TF1 *fit_pos_1 = pos_1->GetFunction("gaus");
	Double_t pos_1_mean = fit_pos_1->GetParameter(1);

	pos_2->Fit("gaus");
	TF1 *fit_pos_2 = pos_2->GetFunction("gaus");
	Double_t pos_2_mean = fit_pos_2->GetParameter(1);

	pos_3->Fit("gaus");
	TF1 *fit_pos_3 = pos_3->GetFunction("gaus");
	Double_t pos_3_mean = fit_pos_3->GetParameter(1);

	pos_4->Fit("gaus");
	TF1 *fit_pos_4 = pos_4->GetFunction("gaus");
	Double_t pos_4_mean = fit_pos_4->GetParameter(1);

	pos_5->Fit("gaus");
	TF1 *fit_pos_5 = pos_5->GetFunction("gaus");
	Double_t pos_5_mean = fit_pos_5->GetParameter(1);

	pos_6->Fit("gaus");
	TF1 *fit_pos_6 = pos_6->GetFunction("gaus");
	Double_t pos_6_mean = fit_pos_6->GetParameter(1);

	pos_7->Fit("gaus");
	TF1 *fit_pos_7 = pos_7->GetFunction("gaus");
	Double_t pos_7_mean = fit_pos_7->GetParameter(1);

	pos_8->Fit("gaus");
	TF1 *fit_pos_8 = pos_8->GetFunction("gaus");
	Double_t pos_8_mean = fit_pos_8->GetParameter(1);

	printf("Pos 1 Mean is %f\n",pos_1_mean);
	printf("Pos 2 Mean is %f\n",pos_2_mean);
	printf("Pos 3 Mean is %f\n",pos_3_mean);
	printf("Pos 4 Mean is %f\n",pos_4_mean);
	printf("Pos 5 Mean is %f\n",pos_5_mean);
	printf("Pos 6 Mean is %f\n",pos_6_mean);
	printf("Pos 7 Mean is %f\n",pos_7_mean);
	printf("Pos 8 Mean is %f\n",pos_8_mean);
}