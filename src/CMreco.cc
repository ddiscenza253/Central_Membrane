#include "CMreco.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>
#include <TNtuple.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <string>
#include <tpc/TpcDefs.h>
#include <TVector3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <cmath>

/// Tracking includes
#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterv2.h>
#include <trackbase/TrkrClusterContainer.h>
#include <trackbase/TrkrClusterHitAssocv2.h>
#include <trackbase/TrkrHitTruthAssoc.h>
#include <trackbase_historic/SvtxTrack_v2.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/SvtxVertexMap.h>
#include <trackbase/TrkrHitSetContainerv1.h>
#include <trackbase/TrkrHitSetContainer.h>
#include <trackbase/TrkrHitv2.h>
#include <trackbase/TrkrHitSetv1.h>

#include <trackbase/ActsTrackingGeometry.h>
#include <trackbase/ActsSurfaceMaps.h>
#include <trackbase_historic/ActsTransformations.h>


using namespace std;

//____________________________________________________________________________..
CMreco::CMreco(const std::string &name):
 SubsysReco(name)
{
  /*
  fout = new TFile("readback_ntuple.root","RECREATE");

  ntp = new TNtuple("ntp", "ntp", "pt:x:y:z:dcaxy:dcaz:vtxid:nclus:qual");
  */
  //cout << "CMreco::CMreco(const std::string &name) Calling ctor" << endl;
}

//____________________________________________________________________________..
CMreco::~CMreco()
{

}

//____________________________________________________________________________..
int CMreco::InitRun(PHCompositeNode *topNode)
{
  /*string name;
  string wName;
  for(int layer = 0; layer < 56; ++layer)
    {
       for(int side = 0; side < 2; ++side)
      {
	      name = "h";
	      name += std::to_string(layer);
	
	        name += "_";
	      name += std::to_string(side);
	      wName = name + "_w";
	      
	      h_layer[layer][side] = new TH1F(name.c_str(),"",MAXPAD,0,MAXPAD);
	      h_weight[layer][side] = new TH1F(wName.c_str(),"",MAXPAD,0,MAXPAD);*/
	      

		/*
	      h_layer[layer] = new TH2F(name.c_str(),"",TpcDefs::MAXTBIN,0,TpcDefs::MAXTBIN,TpcDefs::MAXPAD,0,TpcDefs::MAXPAD);
	      h_weight[layer] = new TH2F(wName.c_str(),"",TpcDefs::MAXTBIN,0,TpcDefs::MAXTBIN,TpcDefs::MAXPAD,0,TpcDefs::MAXPAD);
	      */

        //}
//}
  int ret = GetNodes(topNode);
  return ret;
}

//____________________________________________________________________________..
int CMreco::process_event(PHCompositeNode *topNode)
{
  // Loop over all SvtxTracks




  //local coord conversion below
  ActsTrackingGeometry *tgeometry = findNode::getClass<ActsTrackingGeometry>(topNode,"ActsTrackingGeometry");
  std::cout << "got tgeom" << std::endl;
  ActsSurfaceMaps *surfmaps = findNode::getClass<ActsSurfaceMaps>(topNode,"ActsSurfaceMaps");
  std::cout << "got surfmaps" << std::endl;
  if(!tgeometry or !surfmaps)
    {
      std::cout << PHWHERE << "No Acts geometry on node tree. Can't  continue."
		<< std::endl;
    } 
  
  if(Verbosity() > 0) std::cout << std::endl << "original size of cluster map: " << _cluster_map->size() << std::endl;  
  TrkrHitSetContainer::ConstRange hitsetrange = _hitset_map->getHitSets(TrkrDefs::TrkrId::tpcId);
  
  std::vector<TVector3>pos; //position vector in cartesian
  std::vector<int>layer; //cluster layer number
  std::vector<int>i_pair; //vector for pair matching
  std::vector<float>energy;//vector for energy values of clusters
  int nTpcClust = 0;


  for (TrkrHitSetContainer::ConstIterator hitsetitr = hitsetrange.first;
       hitsetitr != hitsetrange.second;
       ++hitsetitr)
    {
      auto clusRange = _cluster_map->getClusters(hitsetitr->first);
      for (auto clusiter = clusRange.first; 
	   clusiter != clusRange.second; ++clusiter)
	{
	  //std::cout << "Getting cluskey" << std::endl;
	  TrkrDefs::cluskey cluskey = clusiter->first;
	  //std::cout << "Cluster Key: " << cluskey << std::endl;
	  TrkrCluster *cluster = clusiter->second;
	  ActsTransformations transformer;
	  auto glob = transformer.getGlobalPosition(cluster,surfmaps,tgeometry);

	  float x = glob(0);
	  float y = glob(1);
	  float z = glob(2);
	  i_pair.push_back(-1);
	  energy.push_back(cluster->getAdc());
	  nTpcClust++;
	  if(Verbosity() > 0)
	    std::cout << "nTpcClust in loop is: " << nTpcClust << std::endl;
	  pos.push_back(TVector3(x,y,z));
	  layer.push_back((int)(TrkrDefs::getLayer(cluskey)));
	  if(Verbosity() > 0)
	    {
	      std::cout << ":\t" << x << "\t" << y << "\t" << z <<std::endl;
	    }
	}
    }
  if(Verbosity() > 0)
    std::cout << "Got out of cluster loop\t nTpcClust =" << nTpcClust << std::endl;
  TVector3 delta;
  float dphi;
  TVector2 delta2;
  //float maxdist=0.01;


  char temp[500];
  sprintf(temp,
  	  "./eval_output/Energy_Histograms_%i.root",process);
  fout = new TFile(temp,"RECREATE");





  TH1F *hDist=new TH1F("hDist","3D distance to nearby clusters on same padrow;dist[cm]",100,-1,10);
  TH2F *hDistRow=new TH2F("hDistRow","phi distance to nearby clusters vs (lower)row;dist[rad];padrow",100,-0.001,0.01,60,-0.5,59.5);
  TH1F *hDist2=new TH1F("hDist2","phi distance to nearby clusters on same padrow;dist[rad]",100,-0.001,0.01);
  TH2F *hDistRowAdj=new TH2F("hDistRowAdj","phi distance to nearby clusters vs (lower)row;dist[rad];(lower) padrow",100,-0.001,0.01,60,-0.5,59.5);
  TH1F *hDist2Adj=new TH1F("hDist2Adj","phi distance to nearby clusters on adjacent padrow;dist[rad]",100,-0.001,0.01);



  
  
  for (int i=0;
       i<nTpcClust;
       i++)
    {
      for (int j=i+1;
	   j<nTpcClust;
	   j++)
	    {      //redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
    	      if (abs(layer[i]-layer[j])==0)
		{
		  delta=pos[i]-pos[j];
		  dphi=abs(pos[i].DeltaPhi(pos[j]));
		  hDist->Fill(delta.Mag());
		  hDist2->Fill(dphi);
		  hDistRow->Fill(dphi,layer[i]);
		}
	      if (abs(layer[i]-layer[j])==1)
		{  //match those centers to the known/expected stripe positions
		  
		  delta=pos[i]-pos[j];
		  dphi=abs(pos[i].DeltaPhi(pos[j]));
		  //float newdist=delta.Mag();
		  //hDist->Fill(delta.Mag());
		  hDist2Adj->Fill(dphi);
		  hDistRowAdj->Fill(dphi,layer[i]);
		}
	    }
	}
     
      //now for each cluster, find its nearest partner on an adjacent row:
      const float maxphidist=0.003;//as read off the plots.
      for (int i=0;i<nTpcClust;i++)
	{
	  float bestphidist=maxphidist;
	  for (int j=0;j<nTpcClust;j++)
	    {
	      //redundant to the 'adjacent row' check:  if (i==j) continue; //can't match yourself.
	      if (abs(layer[i]-layer[j])!=1)
		{
		  continue; //must match to an ADJACENT row.
		}
	      //float delta=pos[i].XYvector()-pos[j].XYvector();
	      float newphidist=abs(pos[i].DeltaPhi(pos[j]));
	      if (newphidist<bestphidist)
		{
		  i_pair[i]=j;
		  bestphidist=newphidist;
		}
	    }  	
	}
     
      // check to see if the cluster pairs each match each other
      vector<bool>goodPair;
      bool allGood=true;
      int  nGood=0;
     
      for (int i=0;i<nTpcClust;i++)
	{
	  //std::cout << "before myPair declaration\t value of i: " << i << std::endl;
	  int myPair=i_pair[i];
	  // std::cout << "Beginning of for loop" << std::endl;
	  int itsPair=myPair<0?-1:i_pair[myPair];
	  if (i!=itsPair){
	    // std::cout << "About to push back goodPair (true)" << std::endl;
	    goodPair.push_back(false);
	    //printf("%d doesn't match back (%d matches to %d instead)\n",i,myPair,itsPair);
	    allGood=false;
	  } else {
	    // std::cout << "got into else" << std::endl;
	    if (i<myPair) nGood++;
	    {
	      // std::cout << "About to push back goodPair (false)" << std::endl;
	      goodPair.push_back(true);
	    }
	      //good match:
	    // if (i<myPair){ // only keep one copy of the pair
	  }
	}
      if (allGood)
	{
	  printf("All Good!\n");
	} 
      else 
	{
	  printf("nGood=%d out of %d\n",nGood,nTpcClust/2);
	}
      //TCanvas *c=new TCanvas("c","canvas",1200,1200); 
     
      //build the weighted cluster centers
      //match those centers to the known/expected stripe positions

      //save this to a file sara's code can ingest.
      //do sara's code to measure distortions.
    


      TH1F *hClustE[3];
      vector<TVector3> avepos;
      vector<float>avex,avey;
      vector<float>aver,avephi;
      //eventually the below needs to be turned into building sets of weighted cluster centers, not just computing the cluster total E.
       if (1)
	{
	  //histogram the joined cluster energies
	  hClustE[0]= new TH1F("hRawClusterEnergy","Cluster Energy Before Merging;E[?]",100,0,2000);
	  hClustE[1] = new TH1F("hMatchedClusterEnergy","Pair Cluster Energy After Merging;E[?]",100,0,2000);
	  hClustE[2] = new TH1F("hSoloClusterEnergy","Lone Cluster Energy After Merging;E[?]",100,0,2000);
	for (int i=0;i<nTpcClust;++i)
	  {
	    // std::cout << "goodPair: " << goodPair[i] << std::endl;
	    //std::cout << "i_pair: " << i_pair[i] << std::endl;
	    hClustE[0]->Fill(energy[i]);
	    //std::cout << "Got in the for loop" << std::endl;
	    if (goodPair[i])
	      {
		//std::cout << "Pair position before matching: \n" << pos[i].X() << std::endl;
		if (i_pair[i]>i)
		  {
		    hClustE[1]->Fill(energy[i]+energy[i_pair[i]]);
		    TVector3 temppos=energy[i]*pos[i];
		    temppos=temppos+(energy[i_pair[i]]*pos[i_pair[i]]);
		    temppos=temppos*(1./(energy[i]+energy[i_pair[i]]));
		    avepos.push_back(temppos);
		    avex.push_back(avepos.back().X());
		    avey.push_back(avepos.back().Y());
		    aver.push_back(pow(avex[i]*avex[i]+avey[i]*avey[i],0.5));
		    avephi.push_back(atan(avey[i]/avex[i]));
		    // std::cout << "Pair position after avg: \n" << aver[i] << std::endl;
		  }
	      } 
	    else 
	      {
		hClustE[2]->Fill(energy[i]);
		//std::cout << "got to the non pair stuff \n" << std::endl;
		avepos.push_back(pos[i]);
		avex.push_back(avepos.back().X());
		avey.push_back(avepos.back().Y());
	      }
	  }
	_corrected_CMcluster_map  = findNode::getClass<TrkrClusterContainer>(topNode, "CORRECTED_TRKR_CLUSTER");
	if(!_corrected_CMcluster_map)
	  {
	    std::cout << "Creating node CORRECTED_CM_CLUSTER" << std::endl;
	    
	    PHNodeIterator iter(topNode);
	    
	    // Looking for the DST node
	    PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
	    if (!dstNode)
	      {
		std::cout << PHWHERE << "DST Node missing, doing nothing." << std::endl;
		return Fun4AllReturnCodes::ABORTRUN;
	      }      
	    PHNodeIterator dstiter(dstNode);
	    PHCompositeNode *DetNode =
	      dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "TRKR"));
	    if (!DetNode)
	      {
		DetNode = new PHCompositeNode("TRKR");
		dstNode->addNode(DetNode);
	      }
	    
	    _corrected_cluster_map = new TrkrClusterContainerv3;
	    PHIODataNode<PHObject> *TrkrClusterContainerNode =
	      new PHIODataNode<PHObject>(_corrected_CMcluster_map, "CORRECTED_CM_CLUSTER", "PHObject");
	    DetNode->addNode(TrkrClusterContainerNode);
	  }    
	
	  /*
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
	
	
	for (int i=0;i<avepos.size();i++)
	  {
	    searchArea=Form("fabs(position.Perp()-%f)<%f && fabs(position.Phi()-%f)<%f",avepos[i].Perp()+lastroff,rseek,avepos[i].Phi()+lastphioff,phiseek);
	    pTree->Draw("position.X():position.Y()",searchArea,"goff");
	    int nCand=pTree->GetSelectedRows();
	    float bestdist=99;
	    for (int j=0;j<nCand;j++)
	      {
		TVector3 temppos(pTree->GetVal(0)[j],pTree->GetVal(1)[j],0);
		temppos=temppos-avepos[i];
		float dist=temppos.Perp();
		if (dist<bestdist)
		  {
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
	hP->Draw("colz");*/
	/*	delete hDist;
	 delete hDistRow;
	 delete hDist2;
	 delete hDistRowAdj;
	 delete hDistRowAdj;
	*/

	



	hClustE[0]->Write();
	hClustE[1]->Write();
	hClustE[2]->Write();
	hDist->Write();
	hDistRow->Write();
	hDist2->Write();
	hDistRowAdj->Write();
	hDist2Adj->Write();
   
	 return Fun4AllReturnCodes::EVENT_OK;
	}
       /*
       delete hDist;
       delete hDistRow;
       delete hDist2;
       delete hDistRowAdj;
       delete hDist2Adj;

       */

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CMreco::EndRun(PHCompositeNode *topNode)
{

  //fout->Write();
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CMreco::End(PHCompositeNode *topNode)
{
  /* for(int layer = 0; layer < 56; ++layer)
    {
       for(int side = 0; side < 2; ++side)
      {
        //for(int sector = 0; sector < 12; ++sector)
	//{
	h_layer[layer][side]->Divide(h_weight[layer][side]);
        //}
       }
       }*/
  //fout->Write();
  //ntp->Write();
  fout->Close();

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..

int  CMreco::GetNodes(PHCompositeNode* topNode)
{
  //---------------------------------
  // Get Objects off of the Node Tree
  //---------------------------------

  _cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!_cluster_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTER" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  _cluster_hit_map = findNode::getClass<TrkrClusterHitAssoc>(topNode, "TRKR_CLUSTERHITASSOC");
  if (!_cluster_hit_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_CLUSTERHITASSOC" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  _hitset_map = findNode::getClass<TrkrHitSetContainer>(topNode, "TRKR_HITSET");
  if (!_hitset_map)
  {
    cerr << PHWHERE << " ERROR: Can't find node TRKR_HITSET" << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }


  /* _vertex_map = findNode::getClass<SvtxVertexMap>(topNode, "SvtxVertexMap");
  if (!_vertex_map)
    {
      cerr << PHWHERE << " ERROR: Can't find SvtxVertexMap." << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  */
    
  /* _track_map = findNode::getClass<SvtxTrackMap>(topNode,  "SvtxTrackMap");
  if (!_track_map)
  {
    cerr << PHWHERE << " ERROR: Can't find SvtxTrackMap: " << endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }
  */
  return Fun4AllReturnCodes::EVENT_OK;
}


