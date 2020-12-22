#ifndef StudyMTVEfficiencyv2_h
#define StudyMTVEfficiencyv2_h

#include "SDL/Event.h"

#include "Study.h"

#include <vector>
#include <tuple>

#include "TString.h"
#include "trktree.h"
#include "constants.h"

#include "AnalysisConfig.h"

#include "trkCore.h"

class StudyMTVEfficiencyv2 : public Study
{

public:
    const char* studyname;

    std::vector<float> pt_boundaries;

    int pdgid_of_study;

    float dxy_thresh;
    float dz_thresh;
    float pt_thresh;

    float diff_pt;

    std::array<std::vector<float>, 6> tc_matched_track_pt_mtv;
    std::array<std::vector<float>, 6> tc_all_track_pt_mtv;
    std::array<std::vector<float>, 6> tc_matched_track_eta_mtv;
    std::array<std::vector<float>, 6> tc_all_track_eta_mtv;
    std::array<std::vector<float>, 6> tc_matched_track_phi_mtv;
    std::array<std::vector<float>, 6> tc_all_track_phi_mtv;
    std::array<std::vector<float>, 6> tc_matched_track_dxy_mtv;
    std::array<std::vector<float>, 6> tc_all_track_dxy_mtv;

    std::array<std::vector<float>, 6> tl_matched_track_pt_mtv;
    std::array<std::vector<float>, 6> tl_all_track_pt_mtv;
    std::array<std::vector<float>, 6> tl_matched_track_eta_mtv;
    std::array<std::vector<float>, 6> tl_all_track_eta_mtv;
    std::array<std::vector<float>, 6> tl_matched_track_phi_mtv;
    std::array<std::vector<float>, 6> tl_all_track_phi_mtv;
    std::array<std::vector<float>, 6> tl_matched_track_dxy_mtv;
    std::array<std::vector<float>, 6> tl_all_track_dxy_mtv;

    std::array<std::vector<float>, 6> tp_matched_track_pt_mtv;
    std::array<std::vector<float>, 6> tp_all_track_pt_mtv;
    std::array<std::vector<float>, 6> tp_matched_track_eta_mtv;
    std::array<std::vector<float>, 6> tp_all_track_eta_mtv;
    std::array<std::vector<float>, 6> tp_matched_track_phi_mtv;
    std::array<std::vector<float>, 6> tp_all_track_phi_mtv;
    std::array<std::vector<float>, 6> tp_matched_track_dxy_mtv;
    std::array<std::vector<float>, 6> tp_all_track_dxy_mtv;

    StudyMTVEfficiencyv2(
        const char* studyName,
        std::vector<float> ptbounds,
        int pdgid
        );
    virtual void bookStudy();
    virtual void bookStudy_v1();
    virtual void doStudy(SDL::Event& recoevent, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);
    virtual void doStudy_v1(SDL::Event& recoevent, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);
    virtual void doStudy_trackCandidateEfficiency(SDL::Event& recoevent, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);
    virtual void doStudy_trackletEfficiency(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);
    virtual void doStudy_tripletEfficiency(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);

};

#endif
