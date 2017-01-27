int readMyTree(){
	TFile* f1 = new TFile("magon1M.root");
	TFile* f2 = new TFile("magon100M.root");

	TTree* t1 = dynamic_cast<TTree*>(f1->Get("mytree"));
	TTree* t2 = dynamic_cast<TTree*>(f2->Get("mytree"));


	TCanvas* c1 = new TCanvas();
	c1->Divide(2,2);

	c1->cd(1);
	t1->Draw("E >> h1"," code==22 && E>0.1 && targetE != -1");

	c1->cd(2);
	t1->Draw("E >> h2"," code==22 && E>0.1 && targetE == -1");

	c1->cd(3);
	t2->Draw("E >> h3"," code==22 && E>0.1 && targetE != -1");

	c1->cd(4);
	t2->Draw("E >> h4"," code==22 && E>0.1 && targetE == -1");

	//h1.SetTitle("hist;Energy;counts");

	return 0;
}
