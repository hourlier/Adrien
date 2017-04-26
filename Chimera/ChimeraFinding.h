#ifndef CHIMERAFINDING_H
#define CHIMERAFINDING_H

#include <string>
#include <fstream>
#include <vector>

#include "DataFormat/hit.h"
#include "DataFormat/track.h"
#include "DataFormat/storage_manager.h"

#include "ChimeraTrackEvaluator.h"

class ChimeraFinding{
public:
    ChimeraFinding(){
        _TargetFile = "Targets.csv";
        _particleType = "nu";
    }

    virtual ~ChimeraFinding(){}
    virtual bool initialize();
    virtual bool analyze(larlite::storage_manager &storage);
    virtual bool finalize();

    inline void SetTargetFile(std::string targetFile){_TargetFile = targetFile;}
    inline void SetTrackGenerator(std::string trackGenerator){_TrackGenerator = trackGenerator;}
    inline void SetHitGenerator(std::string hitGenerator){_HitGenerator = hitGenerator;}
    inline void SetInputFile(std::string inputFile){_InputFile = inputFile;}
    inline void SetParticleType(std::string particleType){_particleType = particleType;}
    inline void SetSigmaEval(double *sigmaPar){for(int ipar = 0;ipar<3;ipar++){_sigmaPar[ipar] = sigmaPar[ipar];} }

    inline std::vector<double>                     GetScores(){return _max_score;}
    inline std::vector<larlite::track>             GetBestTracks(){return _best_Track;}
    inline std::vector<std::vector<double> >       GetParTarget(){return FullTargetParameters;}
    inline std::vector<std::vector<larlite::hit> > GetBestHitClusters(){return _best_HitCluster;}


protected:
    void ReadTargetFile();

private:
    int _run;
    int _subrun;
    int _event;
    int _track;
    size_t NrequestedEvents;
    double _sigmaPar[3];

    std::vector<std::vector<double> > parTargets;
    std::vector<double> _max_score;
    std::vector<int> _best_run;
    std::vector<int> _best_subrun;
    std::vector<int> _best_event;
    std::vector<int> _best_track_id;
    std::vector<larlite::track> _best_Track;
    std::vector<std::vector<larlite::hit> > _best_HitCluster;

    std::string _TargetFile;
    std::string _InputFile;
    std::string _TrackGenerator;
    std::string _HitGenerator;
    std::string _particleType;

    std::vector<std::vector<double> > FullTargetParameters;
    larlite::ChimeraTrackEvaluator evaluator;
};

#endif
