//make a reco class to hold the coords?

void CMreco(){
  TFile *recoFile=TFile::Open("G4sPHENIX_rcc_g4svtx_eval.root","READ");
  TTree *cTree=(TTree*)recoFile->Get("ntp_cluster");

  TFile *patFile=TFile::Open("~/sphenix/toys/cmpattern/cmDistHitsTree_Event0.root","READ");
  TTree *pTree=(TTree*)patFile->Get("tree");

  
  std::vector<TVector3>pos;
  std::vector<int>layer;
  std::vector<int>i_pair;
  std::vector<float>energy;

  //get our clusters for event 1 (zero doesn't have them for some reason?):
  //eventually, need to do this per-event for comparisons.
  cTree->Draw("x:y:z:e:layer","abs(z)<1 && z>0 && r>25 && event==1","goff");
  double *x=cTree->GetVal(0);
  double *y=cTree->GetVal(1);
  double *z=cTree->GetVal(2);
  double *e=cTree->GetVal(3);
  double *l=cTree->GetVal(4);
  int nClust=cTree->GetSelectedRows();
  int nTpcClust=0;
  
  for (int i=0;i<nClust;i++){
    if (sqrt(x[i]*x[i]+y[i]*y[i])<25) continue; //skip things not in the TPC
    pos.push_back(TVector3(x[i],y[i],z[i]));
    layer.push_back((int)(l[i]));
    if (l[i]-floor(l[i])>0.0001) printf("layer foul!  %d!=%f for i=%d\n",layer.back(),l[i],i);
    i_pair.push_back(-1);
    energy.push_back(e[i]);
    nTpcClust++;
  }

  
  // for each cluster, make some diagnostic plots about nearby clusters
  TH1F *hDist=new TH1F("hDist","3D distance to nearby clusters on same padrow;dist[cm]",100,-1,10);
  TH2F *hDistRow=new TH2F("hDistRow","phi distance to nearby clusters vs (lower)row;dist[rad];padrow",100,-0.001,0.01,60,-0.5,59.5);
  TH1F *hDist2=new TH1F("hDist2","phi distance to nearby clusters on same padrow;dist[rad]",100,-0.001,0.01);
  TH2F *hDistRowAdj=new TH2F("hDistRowAdj","phi distance to nearby clusters vs (lower)row;dist[rad];(lower) padrow",100,-0.001,0.01,60,-0.5,59.5);
  TH1F *hDist2Adj=new TH1F("hDist2Adj","phi distance to nearby clusters on adjacent padrow;dist[rad]",100,-0.001,0.01);
  for (int i=0;i<nTpcClust;i++){
    float maxdist=0.01;
    for (int j=i+1;j<nTpcClust;j++){
      //redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
      if (abs(l[i]-l[j])==0){
	TVector3 delta=pos[i]-pos[j];
	float dphi=abs(pos[i].DeltaPhi(pos[j]));
	hDist->Fill(delta.Mag());
	hDist2->Fill(dphi);
	hDistRow->Fill(dphi,l[i]);
      }
      if (abs(l[i]-l[j])==1){  //match those centers to the known/expected stripe positions

	TVector3 delta=pos[i]-pos[j];
	TVector2 delta2=pos[i].XYvector()-pos[j].XYvector();
	float dphi=abs(pos[i].DeltaPhi(pos[j]));
	float newdist=delta.Mag();
	//hDist->Fill(delta.Mag());
	hDist2Adj->Fill(dphi);
	hDistRowAdj->Fill(dphi,l[i]);
      }
    }
  }

  //now for each cluster, find its nearest partner on an adjacent row:
  const float maxphidist=0.003;//as read off the plots.
  for (int i=0;i<nTpcClust;i++){
    float bestphidist=maxphidist;
    for (int j=0;j<nTpcClust;j++){
      //redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
      if (abs(l[i]-l[j])!=1) continue; //must match to an ADJACENT row.
      //float delta=pos[i].XYvector()-pos[j].XYvector();
      float newphidist=abs(pos[i].DeltaPhi(pos[j]));
      if (newphidist<bestphidist){
	i_pair[i]=j;
	bestphidist=newphidist;
      }
    }  	
  }

  // check to see if the cluster pairs each match each other
  vector<bool>goodPair;
  bool allGood=true;
  int  nGood=0;
  for (int i=0;i<nTpcClust;i++){
    int myPair=i_pair[i];
    int itsPair=myPair<0?-1:i_pair[myPair];
    if (i!=itsPair){
      goodPair.push_back(false);
      //printf("%d doesn't match back (%d matches to %d instead)\n",i,myPair,itsPair);
      allGood=false;
    } else {
      if (i<myPair) nGood++;
      goodPair.push_back(true);
      //good match:
      // if (i<myPair){ // only keep one copy of the pair
    }
  }
  if (allGood){
    printf("All Good!\n");
  } else {
    printf("nGood=%d out of %d\n",nGood,nTpcClust/2);
  }

  TCanvas *c=new TCanvas("c","canvas",1200,1200);


  //build the weighted cluster centers


  //match those centers to the known/expected stripe positions

  //save this to a file sara's code can ingest.
  //do sara's code to measure distortions.


  
  
    TH1F *hClustE[3];
    vector<TVector3> avepos;
    vector<float>avex,avey;
    //eventually the below needs to be turned into building sets of weighted cluster centers, not just computing the cluster total E.
  if (1){//histogram the joined cluster energies
    hClustE[0]= new TH1F("hRawClusterEnergy","Cluster Energy Before Merging;E[?]",100,0,2000);
    hClustE[1] = new TH1F("hMatchedClusterEnergy","Pair Cluster Energy After Merging;E[?]",100,0,2000);
    hClustE[2] = new TH1F("hSoloClusterEnergy","Lone Cluster Energy After Merging;E[?]",100,0,2000);
    for (int i=0;i<nTpcClust;i++){
      hClustE[0]->Fill(energy[i]);
      if (goodPair[i]){
	if (i_pair[i]>i){
	  hClustE[1]->Fill(energy[i]+energy[i_pair[i]]);
	  TVector3 temppos=energy[i]*pos[i];
	  temppos=temppos+(energy[i_pair[i]]*pos[i_pair[i]]);
	  temppos=temppos*(1./(energy[i]+energy[i_pair[i]]));
	  avepos.push_back(temppos);
	  avex.push_back(avepos.back().X());
	  avey.push_back(avepos.back().Y());
	}
      } else {
	hClustE[2]->Fill(energy[i]);
	avepos.push_back(pos[i]);
	avex.push_back(avepos.back().X());
	avey.push_back(avepos.back().Y());
     }
    }
    c->Divide(3,2);
    c->cd(1);
    hClustE[0]->Draw();
    c->cd(2);
    hClustE[1]->SetLineColor(kGreen);
    hClustE[1]->Draw();
    hClustE[2]->SetLineColor(kRed);
    hClustE[2]->Draw("same");
    c->cd(3);
    pTree->Draw("position.Y():position.X()");
    TGraph *recoPos=new TGraph(avex.size(),&(avex[0]),&(avey[0]));
    recoPos->Draw("A*");
    pTree->SetMarkerColor(kRed);
    pTree->Draw("position.Y():position.X()","1","same");


    //for every reco pos, seek the closest real pad
    //this should start at low r and work upward, using the previous offset to guide the next one.
    const float phiseek=0.2;
    const float rseek=2.0;
    float lastphioff=0;
    float lastroff=0;
    TString searchArea;
    TVector3 closest(0,0,0);
    vector<TVector3>matches;

    TH2F *hSample=new TH2F("hS","hS",50,-100,100,50,-100,100);
    TH2F *hR=new TH2F("hR","Radial Matching Distortion vs position",50,-100,100,50,-100,100);
    TH2F *hP=new TH2F("hP","Phi Matching Distortion vs position",50,-100,100,50,-100,100);

    
    for (int i=0;i<avepos.size();i++){
      searchArea=Form("fabs(position.Perp()-%f)<%f && fabs(position.Phi()-%f)<%f",avepos[i].Perp()+lastroff,rseek,avepos[i].Phi()+lastphioff,phiseek);
      pTree->Draw("position.X():position.Y()",searchArea,"goff");
      int nCand=pTree->GetSelectedRows();
      float bestdist=99;
      for (int j=0;j<nCand;j++){
	TVector3 temppos(pTree->GetVal(0)[j],pTree->GetVal(1)[j],0);
	temppos=temppos-avepos[i];
	float dist=temppos.Perp();
	if (dist<bestdist){
	  bestdist=dist;
	  closest.SetXYZ(pTree->GetVal(0)[j],pTree->GetVal(1)[j],0);
	  lastroff=closest.Perp()-avepos[i].Perp();
	  lastphioff=closest.Phi()-avepos[i].Phi();
	  //need to fix the phi wrapping
	}
      }
      matches.push_back(closest);
      hSample->Fill(avepos[i].X(),avepos[i].Y());
      hR->Fill(avepos[i].X(),avepos[i].Y(),matches[i].Perp()-avepos[i].Perp());
      hP->Fill(avepos[i].X(),avepos[i].Y(),matches[i].Phi()-avepos[i].Phi());
    }

    hR->Divide(hSample);
    hP->Divide(hSample);
    c->cd(4);
    hSample->Draw("colz");
   c->cd(5);
    hR->Draw("colz");
   c->cd(6);
    hP->Draw("colz");
    
    return;
  }
  

  
  if (0){
    c->Divide(2,2);
    c->cd(1)->SetLogy();
    hDist2->Draw();
    c->cd(2)->SetLogz();
    hDistRow->Draw("colz");
    c->cd(3)->SetLogy();
    hDist2Adj->Draw();
    c->cd(4)->SetLogz();
    hDistRowAdj->Draw("colz");
    return;
  }
  cTree->SetMarkerStyle(1);
  cTree->SetMarkerColor(kBlack);

  //float xlow=-60, xhigh=-20;ylow=-60,yhigh=-20;
  float xlow=-80, xhigh=80, ylow=-80,yhigh=80;

  cTree->Draw("y:x",Form("abs(z)<1 && z>0  && r>25 && x>%f && y>%f && x<%f && y<%f && event==1",xlow,ylow,xhigh,yhigh));
  TLine *line=new TLine();
  for (int i=0;i<nTpcClust;i++){
    int i2=i_pair[i];
    if (goodPair[i] && i>i2){
      if(pos[i].X()>xlow && pos[i].Y()>ylow && pos[i].X()<xhigh && pos[i].Y()<yhigh){
	line->SetLineColor(kGreen);
	line->DrawLine(pos[i].X(),pos[i].Y(),pos[i2].X(),pos[i2].Y());
      }
    }
    if (!goodPair[i] && i2!=-1){
      if(pos[i].X()>xlow && pos[i].Y()>ylow && pos[i].X()<xhigh && pos[i].Y()<yhigh){
	line->SetLineColor(kRed);
	line->DrawLine(pos[i].X(),pos[i].Y(),pos[i2].X(),pos[i2].Y());
      }
    }
  }
  if (1){//draw module boundaries
    TEllipse *circle=new TEllipse();
    float nominal_radius[]={758.4,583.5,574.9,411.4,402.6,221.0};
    circle->SetLineColorAlpha(kBlue,0.25);
    circle->SetFillColorAlpha(kBlue,0.0);
    for (int i=0;i<6;i++){
      circle->DrawEllipse(0,0,0.1*nominal_radius[i],0,0,360,0);
    }
    circle->SetLineColorAlpha(kMagenta,0.25);
    //for (int i=0;i<6;i++){
    //  float shifted_radius=0.9+0.1*nominal_radius[i]+(0.1*nominal_radius[i]-30)*(1.1/40.);
    //  circle->DrawEllipse(0,0,shifted_radius,0,0,360,0);
    //}
  }

  if (0){//draw row labels
    TLatex *tex=new TLatex();
    tex->SetTextSize(0.015);
    for (int i=0;i<nTpcClust;i++){
      if(pos[i].X()>xlow && pos[i].Y()>ylow && pos[i].X()<xhigh && pos[i].Y()<yhigh){
	if(goodPair[i]) {
	  tex->SetTextColor(kBlack);
	} else {
	  tex->SetTextColor(kRed);
	}
	tex->DrawLatex(pos[i].X(),pos[i].Y(),Form("%d",layer[i]));
      }
    }
  }
  return;
}


