// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef CMRECO_H
#define CMRECO_H

#include <string>

#include <fun4all/SubsysReco.h>

#include <trackbase/TrkrDefs.h>
#include <trackbase/TrkrClusterContainer.h>

class PHCompositeNode;
class SvtxTrackMap;
class SvtxTrack;
class SvtxVertexMap;
class TrkrClusterHitAssoc;
class TrkrHitSetContainer;

class TF1;
class TNtuple;
class TFile;
class TH1F;

class CMreco : public SubsysReco
{
 public:

 CMreco(const std::string &name = "CMreco");

  virtual ~CMreco();

  void set_track_map_name(const std::string &map_name) { _track_map_name = map_name; }
  void set_process(const unsigned int proc)
  {
    process = proc;
  }

 //! run initialization
  int InitRun(PHCompositeNode *topNode);

  //! event processing
  int process_event(PHCompositeNode *topNode);

  int EndRun(PHCompositeNode *topNode);

  //! end of process
int End(PHCompositeNode *topNode);

 protected:
  
 private:

  int GetNodes(PHCompositeNode* topNode);

  std::string _track_map_name;

  SvtxTrackMap *_track_map{nullptr};
  SvtxTrack *_track{nullptr};
  SvtxVertexMap *_vertex_map{nullptr};
  TrkrClusterContainer *_cluster_map;
  TrkrClusterContainer *_corrected_CMcluster_map{nullptr};
  TrkrClusterHitAssoc *_cluster_hit_map;
  TrkrHitSetContainer *_hitset_map;
  const bool clst_track = true;

  TH1F *h_layer[56][2];
  TH1F *h_weight[56][2];
  
  unsigned int MAXPAD = 2400;
  unsigned int process = 0;



  TNtuple *ntp{nullptr};
  TFile *fout;

};

#endif // TPCGEMGAINCALB_H
