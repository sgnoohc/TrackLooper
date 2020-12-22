#include "StudyMTVEfficiencyv2.h"

StudyMTVEfficiencyv2::StudyMTVEfficiencyv2(
        const char* studyName,
        std::vector<float> ptbounds,
        int pdgid
        )
{

    studyname = studyName;

    pt_boundaries = ptbounds;

    pdgid_of_study = pdgid;

    dxy_thresh = 2.5;

    dz_thresh = 30;

    pt_thresh = 1.5;

}

void StudyMTVEfficiencyv2::bookStudy()
{
    bookStudy_v1();
}

void StudyMTVEfficiencyv2::bookStudy_v1()
{
    // Book Histograms
    const int nlayers = NLAYERS;

    const float etamax = 2.5;

    const int eta_nbins = 180;

    const float phimax = M_PI;

    const int phi_nbins = 50;

    ana.histograms.addVecHistogram("tc_matched_track_pt_mtv" , pt_boundaries, [&]()                { return tc_matched_track_pt_mtv[0];  } );
    ana.histograms.addVecHistogram("tc_all_track_pt_mtv"     , pt_boundaries, [&]()                { return tc_all_track_pt_mtv[0];      } );
    ana.histograms.addVecHistogram("tc_matched_track_eta_mtv", eta_nbins, -1*etamax, etamax, [&]() { return tc_matched_track_eta_mtv[0]; } );
    ana.histograms.addVecHistogram("tc_all_track_eta_mtv"    , eta_nbins, -1*etamax, etamax, [&]() { return tc_all_track_eta_mtv[0];     } );
    ana.histograms.addVecHistogram("tc_matched_track_dxy_mtv", 50, -10, 10, [&]()                  { return tc_matched_track_dxy_mtv[0]; } );
    ana.histograms.addVecHistogram("tc_all_track_dxy_mtv"    , 50, -10, 10, [&]()                  { return tc_all_track_dxy_mtv[0];     } );

    ana.histograms.addVecHistogram("tc_matched_track_pt_mtv_eta0_0p8" , pt_boundaries, [&]()                { return tc_matched_track_pt_mtv[1];  } );
    ana.histograms.addVecHistogram("tc_all_track_pt_mtv_eta0_0p8"     , pt_boundaries, [&]()                { return tc_all_track_pt_mtv[1];      } );
    ana.histograms.addVecHistogram("tc_matched_track_eta_mtv_eta0_0p8", eta_nbins, -1*etamax, etamax, [&]() { return tc_matched_track_eta_mtv[1]; } );
    ana.histograms.addVecHistogram("tc_all_track_eta_mtv_eta0_0p8"    , eta_nbins, -1*etamax, etamax, [&]() { return tc_all_track_eta_mtv[1];     } );
    ana.histograms.addVecHistogram("tc_matched_track_dxy_mtv_eta0_0p8", 50, -10, 10, [&]()                  { return tc_matched_track_dxy_mtv[1]; } );
    ana.histograms.addVecHistogram("tc_all_track_dxy_mtv_eta0_0p8"    , 50, -10, 10, [&]()                  { return tc_all_track_dxy_mtv[1];     } );

    ana.histograms.addVecHistogram("tc_matched_track_pt_mtv_eta1p0_2p4" , pt_boundaries, [&]()                { return tc_matched_track_pt_mtv[2];  } );
    ana.histograms.addVecHistogram("tc_all_track_pt_mtv_eta1p0_2p4"     , pt_boundaries, [&]()                { return tc_all_track_pt_mtv[2];      } );
    ana.histograms.addVecHistogram("tc_matched_track_eta_mtv_eta1p0_2p4", eta_nbins, -1*etamax, etamax, [&]() { return tc_matched_track_eta_mtv[2]; } );
    ana.histograms.addVecHistogram("tc_all_track_eta_mtv_eta1p0_2p4"    , eta_nbins, -1*etamax, etamax, [&]() { return tc_all_track_eta_mtv[2];     } );
    ana.histograms.addVecHistogram("tc_matched_track_dxy_mtv_eta1p0_2p4", 50, -10, 10, [&]()                  { return tc_matched_track_dxy_mtv[2]; } );
    ana.histograms.addVecHistogram("tc_all_track_dxy_mtv_eta1p0_2p4"    , 50, -10, 10, [&]()                  { return tc_all_track_dxy_mtv[2];     } );

    ana.histograms.addVecHistogram("tc_reco_matched_track_pt_mtv" , pt_boundaries, [&]()                { return tc_matched_track_pt_mtv[3];  } );
    ana.histograms.addVecHistogram("tc_reco_all_track_pt_mtv"     , pt_boundaries, [&]()                { return tc_all_track_pt_mtv[3];      } );
    ana.histograms.addVecHistogram("tc_reco_matched_track_eta_mtv", eta_nbins, -1*etamax, etamax, [&]() { return tc_matched_track_eta_mtv[3]; } );
    ana.histograms.addVecHistogram("tc_reco_all_track_eta_mtv"    , eta_nbins, -1*etamax, etamax, [&]() { return tc_all_track_eta_mtv[3];     } );
    // ana.histograms.addVecHistogram("tc_reco_matched_track_dxy_mtv", 50, -10, 10, [&]()                  { return tc_matched_track_dxy_mtv[3]; } );
    // ana.histograms.addVecHistogram("tc_reco_all_track_dxy_mtv"    , 50, -10, 10, [&]()                  { return tc_all_track_dxy_mtv[3];     } );

    ana.histograms.addVecHistogram("tc_duplicate_matched_track_pt_mtv" , pt_boundaries, [&]()                { return tc_matched_track_pt_mtv[4];  } );
    ana.histograms.addVecHistogram("tc_duplicate_all_track_pt_mtv"     , pt_boundaries, [&]()                { return tc_all_track_pt_mtv[4];      } );
    ana.histograms.addVecHistogram("tc_duplicate_matched_track_eta_mtv", eta_nbins, -1*etamax, etamax, [&]() { return tc_matched_track_eta_mtv[4]; } );
    ana.histograms.addVecHistogram("tc_duplicate_all_track_eta_mtv"    , eta_nbins, -1*etamax, etamax, [&]() { return tc_all_track_eta_mtv[4];     } );
    ana.histograms.addVecHistogram("tc_duplicate_matched_track_dxy_mtv", 50, -10, 10, [&]()                  { return tc_matched_track_dxy_mtv[4]; } );
    ana.histograms.addVecHistogram("tc_duplicate_all_track_dxy_mtv"    , 50, -10, 10, [&]()                  { return tc_all_track_dxy_mtv[4];     } );

    ana.histograms.addVecHistogram("tc_uniqueduplicate_matched_track_pt_mtv" , pt_boundaries, [&]()                { return tc_matched_track_pt_mtv[5];  } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_all_track_pt_mtv"     , pt_boundaries, [&]()                { return tc_all_track_pt_mtv[5];      } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_matched_track_eta_mtv", eta_nbins, -1*etamax, etamax, [&]() { return tc_matched_track_eta_mtv[5]; } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_all_track_eta_mtv"    , eta_nbins, -1*etamax, etamax, [&]() { return tc_all_track_eta_mtv[5];     } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_matched_track_phi_mtv", phi_nbins, -1*phimax, phimax, [&]() { return tc_matched_track_phi_mtv[5]; } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_all_track_phi_mtv"    , phi_nbins, -1*phimax, phimax, [&]() { return tc_all_track_phi_mtv[5];     } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_matched_track_dxy_mtv", 50, -10, 10, [&]()                  { return tc_matched_track_dxy_mtv[5]; } );
    ana.histograms.addVecHistogram("tc_uniqueduplicate_all_track_dxy_mtv"    , 50, -10, 10, [&]()                  { return tc_all_track_dxy_mtv[5];     } );

    ana.histograms.addHistogram("diff_pt" , 180, -0.5, 0.5, [&]() { return diff_pt; } );

    // tracklet

    ana.histograms.addVecHistogram("tl_matched_track_pt_mtv" , pt_boundaries, [&]()                { return tl_matched_track_pt_mtv[0];  } );
    ana.histograms.addVecHistogram("tl_all_track_pt_mtv"     , pt_boundaries, [&]()                { return tl_all_track_pt_mtv[0];      } );
    ana.histograms.addVecHistogram("tl_matched_track_eta_mtv", eta_nbins, -1*etamax, etamax, [&]() { return tl_matched_track_eta_mtv[0]; } );
    ana.histograms.addVecHistogram("tl_all_track_eta_mtv"    , eta_nbins, -1*etamax, etamax, [&]() { return tl_all_track_eta_mtv[0];     } );
    ana.histograms.addVecHistogram("tl_matched_track_dxy_mtv", 50, -10, 10, [&]()                  { return tl_matched_track_dxy_mtv[0]; } );
    ana.histograms.addVecHistogram("tl_all_track_dxy_mtv"    , 50, -10, 10, [&]()                  { return tl_all_track_dxy_mtv[0];     } );

    ana.histograms.addVecHistogram("tl_matched_track_pt_mtv_eta0_0p8" , pt_boundaries, [&]()                { return tl_matched_track_pt_mtv[1];  } );
    ana.histograms.addVecHistogram("tl_all_track_pt_mtv_eta0_0p8"     , pt_boundaries, [&]()                { return tl_all_track_pt_mtv[1];      } );
    ana.histograms.addVecHistogram("tl_matched_track_eta_mtv_eta0_0p8", eta_nbins, -1*etamax, etamax, [&]() { return tl_matched_track_eta_mtv[1]; } );
    ana.histograms.addVecHistogram("tl_all_track_eta_mtv_eta0_0p8"    , eta_nbins, -1*etamax, etamax, [&]() { return tl_all_track_eta_mtv[1];     } );
    ana.histograms.addVecHistogram("tl_matched_track_dxy_mtv_eta0_0p8", 50, -10, 10, [&]()                  { return tl_matched_track_dxy_mtv[1]; } );
    ana.histograms.addVecHistogram("tl_all_track_dxy_mtv_eta0_0p8"    , 50, -10, 10, [&]()                  { return tl_all_track_dxy_mtv[1];     } );

    ana.histograms.addVecHistogram("tl_matched_track_pt_mtv_eta1p0_2p4" , pt_boundaries, [&]()                { return tl_matched_track_pt_mtv[2];  } );
    ana.histograms.addVecHistogram("tl_all_track_pt_mtv_eta1p0_2p4"     , pt_boundaries, [&]()                { return tl_all_track_pt_mtv[2];      } );
    ana.histograms.addVecHistogram("tl_matched_track_eta_mtv_eta1p0_2p4", eta_nbins, -1*etamax, etamax, [&]() { return tl_matched_track_eta_mtv[2]; } );
    ana.histograms.addVecHistogram("tl_all_track_eta_mtv_eta1p0_2p4"    , eta_nbins, -1*etamax, etamax, [&]() { return tl_all_track_eta_mtv[2];     } );
    ana.histograms.addVecHistogram("tl_matched_track_dxy_mtv_eta1p0_2p4", 50, -10, 10, [&]()                  { return tl_matched_track_dxy_mtv[2]; } );
    ana.histograms.addVecHistogram("tl_all_track_dxy_mtv_eta1p0_2p4"    , 50, -10, 10, [&]()                  { return tl_all_track_dxy_mtv[2];     } );

    // triplet

    ana.histograms.addVecHistogram("tp_matched_track_pt_mtv" , pt_boundaries, [&]()                { return tp_matched_track_pt_mtv[0];  } );
    ana.histograms.addVecHistogram("tp_all_track_pt_mtv"     , pt_boundaries, [&]()                { return tp_all_track_pt_mtv[0];      } );
    ana.histograms.addVecHistogram("tp_matched_track_eta_mtv", eta_nbins, -1*etamax, etamax, [&]() { return tp_matched_track_eta_mtv[0]; } );
    ana.histograms.addVecHistogram("tp_all_track_eta_mtv"    , eta_nbins, -1*etamax, etamax, [&]() { return tp_all_track_eta_mtv[0];     } );
    ana.histograms.addVecHistogram("tp_matched_track_dxy_mtv", 50, -10, 10, [&]()                  { return tp_matched_track_dxy_mtv[0]; } );
    ana.histograms.addVecHistogram("tp_all_track_dxy_mtv"    , 50, -10, 10, [&]()                  { return tp_all_track_dxy_mtv[0];     } );

    ana.histograms.addVecHistogram("tp_matched_track_pt_mtv_eta0_0p8" , pt_boundaries, [&]()                { return tp_matched_track_pt_mtv[1];  } );
    ana.histograms.addVecHistogram("tp_all_track_pt_mtv_eta0_0p8"     , pt_boundaries, [&]()                { return tp_all_track_pt_mtv[1];      } );
    ana.histograms.addVecHistogram("tp_matched_track_eta_mtv_eta0_0p8", eta_nbins, -1*etamax, etamax, [&]() { return tp_matched_track_eta_mtv[1]; } );
    ana.histograms.addVecHistogram("tp_all_track_eta_mtv_eta0_0p8"    , eta_nbins, -1*etamax, etamax, [&]() { return tp_all_track_eta_mtv[1];     } );
    ana.histograms.addVecHistogram("tp_matched_track_dxy_mtv_eta0_0p8", 50, -10, 10, [&]()                  { return tp_matched_track_dxy_mtv[1]; } );
    ana.histograms.addVecHistogram("tp_all_track_dxy_mtv_eta0_0p8"    , 50, -10, 10, [&]()                  { return tp_all_track_dxy_mtv[1];     } );

    ana.histograms.addVecHistogram("tp_matched_track_pt_mtv_eta1p0_2p4" , pt_boundaries, [&]()                { return tp_matched_track_pt_mtv[2];  } );
    ana.histograms.addVecHistogram("tp_all_track_pt_mtv_eta1p0_2p4"     , pt_boundaries, [&]()                { return tp_all_track_pt_mtv[2];      } );
    ana.histograms.addVecHistogram("tp_matched_track_eta_mtv_eta1p0_2p4", eta_nbins, -1*etamax, etamax, [&]() { return tp_matched_track_eta_mtv[2]; } );
    ana.histograms.addVecHistogram("tp_all_track_eta_mtv_eta1p0_2p4"    , eta_nbins, -1*etamax, etamax, [&]() { return tp_all_track_eta_mtv[2];     } );
    ana.histograms.addVecHistogram("tp_matched_track_dxy_mtv_eta1p0_2p4", 50, -10, 10, [&]()                  { return tp_matched_track_dxy_mtv[2]; } );
    ana.histograms.addVecHistogram("tp_all_track_dxy_mtv_eta1p0_2p4"    , 50, -10, 10, [&]()                  { return tp_all_track_dxy_mtv[2];     } );

}

void StudyMTVEfficiencyv2::doStudy(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents)
{
    doStudy_trackCandidateEfficiency(event, simtrkevents);
    doStudy_trackletEfficiency(event, simtrkevents);
    doStudy_tripletEfficiency(event, simtrkevents);
}

void StudyMTVEfficiencyv2::doStudy_v1(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents)
{
    // Each do study is performed per event

    // First clear all the output variables that will be used to fill the histograms for this event

    for (unsigned int i = 0; i < 6; ++i)
    {
        tc_matched_track_pt_mtv [i].clear();
        tc_all_track_pt_mtv     [i].clear();
        tc_matched_track_eta_mtv[i].clear();
        tc_all_track_eta_mtv    [i].clear();
        tc_matched_track_phi_mtv[i].clear();
        tc_all_track_phi_mtv    [i].clear();
        tc_matched_track_dxy_mtv[i].clear();
        tc_all_track_dxy_mtv    [i].clear();
    }

    // Barrel layer layer 0
    SDL::Layer barrelLayer0 = event.getLayer(1, SDL::Layer::Barrel);

    //***************************************************************
    // Preselecting Reco Track Candidates that are worth looking into
    //***************************************************************

    // track candidates that match the particle of interest
    std::vector<unsigned int> good_track_candidates;
    std::vector<unsigned int> good_track_candidates_matchable;

    // std::cout <<  " barrelLayer0.getTrackCandidatePtrs().size(): " << barrelLayer0.getTrackCandidatePtrs().size() <<  std::endl;

    // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
    for (unsigned int itc = 0; itc < barrelLayer0.getTrackCandidatePtrs().size(); ++itc)
    {
        SDL::TrackCandidate* tcCandPtr = barrelLayer0.getTrackCandidatePtrs()[itc];
        good_track_candidates.push_back(itc);
        if (matchedSimTrkIdxs(tcCandPtr).size() > 0)
            good_track_candidates_matchable.push_back(itc);
    }

    // std::cout <<  " good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
    // std::cout <<  " good_track_candidates_matchable.size(): " << good_track_candidates_matchable.size() <<  std::endl;

    //***********************
    // Efficiency calculation
    //***********************

    // std::cout <<  " simtrkevents.size(): " << simtrkevents.size() <<  std::endl;

    int nprocessed = 0;

    std::vector<unsigned int> matched_good_track_candidates;

    // Loop over track events
    for (auto& simtrkevent : simtrkevents)
    {

        // std::cout <<  " nprocessed: " << nprocessed <<  std::endl;
        nprocessed++;

        // Unpack the tuple (sim_track_index, SDL::Event containing reco hits only matched to the given sim track)
        unsigned int& isimtrk = std::get<0>(simtrkevent);
        SDL::Event& trackevent = *(std::get<1>(simtrkevent));

        // Create list of booleans
        bool tc_trackmatch_mtv = false;

        // Parse pt and eta of this track
        float pt = std::min((double) trk.sim_pt()[isimtrk], 49.999);
        float eta = trk.sim_eta()[isimtrk];
        float phi = trk.sim_phi()[isimtrk];
        float charge = trk.sim_q()[isimtrk];
        float dxy = trk.sim_pca_dxy()[isimtrk];
        float dz = trk.sim_pca_dz()[isimtrk];

        if (abs(dz) > dz_thresh)
            continue;

        int nmatch = 0;
        int nuniquematch = 0;
        std::vector<int> nbarrels; // unique list of nbarrel in the reco track cand

        std::vector<std::vector<unsigned int>> unique_list_of_module_ids;

        // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        for (auto& igood_tc_idxs : good_track_candidates_matchable)
        {

            SDL::TrackCandidate* tcCandPtr = barrelLayer0.getTrackCandidatePtrs()[igood_tc_idxs];
            std::vector<unsigned int> hitidxs;
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            // std::vector<unsigned int> hitidxs;
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->idx());
            if (isMTVMatch(isimtrk, hitidxs))
            {
                nmatch++;
                tc_trackmatch_mtv = true;
                if (std::find(matched_good_track_candidates.begin(), matched_good_track_candidates.end(), igood_tc_idxs) == matched_good_track_candidates.end())
                    matched_good_track_candidates.push_back(igood_tc_idxs);

                // int nbarrel = getNBarrelTrackCandidate(tcCandPtr);
                // if (std::find(nbarrels.begin(), nbarrels.end(), nbarrel) == nbarrels.end())
                // {
                //     nuniquematch++;
                //     nbarrels.push_back(nbarrel);
                // }
                std::vector<unsigned int> module_ids = getModuleListTrackCandidate(tcCandPtr);
                if (std::find(unique_list_of_module_ids.begin(), unique_list_of_module_ids.end(), module_ids) == unique_list_of_module_ids.end())
                {
                    nuniquematch++;
                    unique_list_of_module_ids.push_back(module_ids);
                }
            }

        }

        // // DEBUG========== CAN ERASE ME LATER
        // if (abs(eta) < 2.4 and pt > 1.0)
        // {
        //     // std::cout << "=============Found a sim track matched to multiple track candidate" << std::endl;
        //     std::cout <<  " nmatch: " << nmatch <<  std::endl;
        //     std::cout <<  " nuniquematch: " << nuniquematch <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " isimtrk: " << isimtrk <<  std::endl;
        //     std::cout <<  " pt: " << pt <<  std::endl;
        //     std::cout <<  " eta: " << eta <<  std::endl;
        //     std::cout <<  " phi: " << phi <<  std::endl;
        //     std::cout <<  " charge: " << charge <<  std::endl;
        // }
        // DEBUG========== CAN ERASE ME LATER

        // // DEBUG========== CAN ERASE ME LATER
        // // if (nmatch > 0)
        // if (pt > 1.5 and abs(eta) < 2.4)
        // {

        //     // std::cout << "=============Found a sim track matched to multiple track candidate" << std::endl;
        //     std::cout <<  " nmatch: " << nmatch <<  std::endl;
        //     std::cout <<  " nuniquematch: " << nuniquematch <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " isimtrk: " << isimtrk <<  std::endl;
        //     std::cout <<  " pt: " << pt <<  std::endl;
        //     std::cout <<  " eta: " << eta <<  std::endl;
        //     std::cout <<  " phi: " << phi <<  std::endl;
        //     std::cout <<  " charge: " << charge <<  std::endl;

        //     for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
        //     {
        //         std::cout <<  " trk.simhit_detId()[trk.sim_simHitIdx()[isimtrk][ith_hit]]: " << trk.simhit_detId()[trk.sim_simHitIdx()[isimtrk][ith_hit]] <<  std::endl;
        //     }

        //     for (auto& module_ids : unique_list_of_module_ids)
        //     {
        //         std::cout << "-----" << std::endl;
        //         for (auto& module_id : module_ids)
        //         {
        //             SDL::Module module(module_id);
        //             std::cout << module;
        //         }
        //     }

        //     // // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        //     // for (auto& igood_tc_idxs : good_track_candidates)
        //     // {

        //     //     SDL::TrackCandidate* tcCandPtr = barrelLayer0.getTrackCandidatePtrs()[igood_tc_idxs];
        //     //     std::vector<unsigned int> hitidxs;
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     if (isMTVMatch(isimtrk, hitidxs))
        //     //     {

        //     //         std::cout << tcCandPtr << std::endl;

        //     //     }

        //     // }

        //     std::cout << "============End printing" << std::endl;

        // }
        // // DEBUG==========




        // if (not tc_trackmatch_mtv and eta > 2.0 and eta < 2.2 and pt > 1.5)
        // {
        //     std::cout << "Found a simulated track with eta > 2.0 that did not match to any track candidate" << std::endl;
        //     printSimTrack(isimtrk);
        //     // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        //     for (auto& igood_tc_idxs : good_track_candidates)
        //     {

        //         SDL::TrackCandidate* tcCandPtr = barrelLayer0.getTrackCandidatePtrs()[igood_tc_idxs];
        //         std::vector<unsigned int> hitidxs;
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //         if (isMTVMatch(isimtrk, hitidxs, true))
        //         {
        //             tc_trackmatch_mtv = true;
        //         }

        //     }
        // }

        // *****************
        // Track Candidate
        // *****************

        float etamaxbound = 2.4;

        if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
        {

            tc_all_track_pt_mtv[0].push_back(pt);
            if (tc_trackmatch_mtv) tc_matched_track_pt_mtv[0].push_back(pt);

            if (abs(eta) < 0.8)
            {
                tc_all_track_pt_mtv[1].push_back(pt);
                if (tc_trackmatch_mtv) tc_matched_track_pt_mtv[1].push_back(pt);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tc_all_track_pt_mtv[2].push_back(pt);
                if (tc_trackmatch_mtv) tc_matched_track_pt_mtv[2].push_back(pt);
            }

            tc_all_track_pt_mtv[4].push_back(pt);
            if (nmatch > 1) tc_matched_track_pt_mtv[4].push_back(pt);

            tc_all_track_pt_mtv[5].push_back(pt);
            if (nuniquematch > 1) tc_matched_track_pt_mtv[5].push_back(pt);
        }

        if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
        {

            tc_all_track_eta_mtv[0].push_back(eta);
            if (tc_trackmatch_mtv) tc_matched_track_eta_mtv[0].push_back(eta);

            if (abs(eta) < 0.8)
            {
                tc_all_track_eta_mtv[1].push_back(eta);
                if (tc_trackmatch_mtv) tc_matched_track_eta_mtv[1].push_back(eta);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tc_all_track_eta_mtv[2].push_back(eta);
                if (tc_trackmatch_mtv) tc_matched_track_eta_mtv[2].push_back(eta);
            }

            tc_all_track_eta_mtv[4].push_back(eta);
            if (nmatch > 1) tc_matched_track_eta_mtv[4].push_back(eta);

            tc_all_track_eta_mtv[5].push_back(eta);
            if (nuniquematch > 1) tc_matched_track_eta_mtv[5].push_back(eta);

            tc_all_track_phi_mtv[5].push_back(phi);
            if (nuniquematch > 1) tc_matched_track_phi_mtv[5].push_back(phi);
        }

        if (pt > pt_thresh and abs(eta) < etamaxbound)
        {

            tc_all_track_dxy_mtv[0].push_back(dxy);
            if (tc_trackmatch_mtv) tc_matched_track_dxy_mtv[0].push_back(dxy);

            if (abs(eta) < 0.8)
            {
                tc_all_track_dxy_mtv[1].push_back(dxy);
                if (tc_trackmatch_mtv) tc_matched_track_dxy_mtv[1].push_back(dxy);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tc_all_track_dxy_mtv[2].push_back(dxy);
                if (tc_trackmatch_mtv) tc_matched_track_dxy_mtv[2].push_back(dxy);
            }

            tc_all_track_dxy_mtv[4].push_back(dxy);
            if (nmatch > 1) tc_matched_track_dxy_mtv[4].push_back(dxy);

            tc_all_track_dxy_mtv[5].push_back(dxy);
            if (nuniquematch > 1) tc_matched_track_dxy_mtv[5].push_back(dxy);
        }

    }

    int n_tc_not_matched = 0;

    for (auto& igood_tc_idxs : good_track_candidates)
    {
        SDL::TrackCandidate* tcCandPtr = barrelLayer0.getTrackCandidatePtrs()[igood_tc_idxs];
        float betaPt_2nd_outer = tcCandPtr->outerTrackletBasePtr()->getRecoVar("betaPt_2nd");
        float betaPt_2nd_inner = tcCandPtr->innerTrackletBasePtr()->getRecoVar("betaPt_2nd");
        if (betaPt_2nd_outer == -999)
            betaPt_2nd_outer = betaPt_2nd_inner;
        float reco_pt = (abs(betaPt_2nd_outer) + abs(betaPt_2nd_inner)) / 2.;
        float reco_eta = tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->eta();

        bool tc_matched = false;
        if (std::find(matched_good_track_candidates.begin(), matched_good_track_candidates.end(), igood_tc_idxs) != matched_good_track_candidates.end())
            tc_matched = true;

        // === DEBUG
        if (not tc_matched)
        {
            n_tc_not_matched++;
        }

        std::vector<unsigned int> module_list = getModuleListTrackCandidate(tcCandPtr);
        std::vector<unsigned int> hitidx_list = getHitIdxListTrackCandidate(tcCandPtr);
        std::vector<unsigned int> simidx_list = allMatchedSimTrkIdxs(tcCandPtr);
        std::vector<         int> gstidx_list = matchedSimTrkIdxs(tcCandPtr);

        // if (not tc_matched)
        // {
        //     std::cout <<  " ==not matched track== " << std::endl;
        //     std::cout <<  " good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
        //     std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
        //     std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
        //     std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
        //     std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
        //     std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
        // }
        // else
        // {
        //     std::cout <<  " ==matched track== " << std::endl;
        //     std::cout <<  " good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
        //     std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
        //     std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
        //     std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
        //     std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
        //     std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
        // }

        // === DEBUG

        float etamaxbound = 2.4;

        // if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
        if (abs(reco_eta) < etamaxbound)
        {
            tc_all_track_pt_mtv[3].push_back(reco_pt);
            if (not tc_matched) tc_matched_track_pt_mtv[3].push_back(reco_pt);
        }

        // if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
        if (reco_pt > pt_thresh and abs(reco_eta) < etamaxbound)
        {
            tc_all_track_eta_mtv[3].push_back(reco_eta);
            if (not tc_matched) tc_matched_track_eta_mtv[3].push_back(reco_eta);
        }

        // if (pt > pt_thresh and abs(eta) < etamaxbound)
        // {
        //     tc_all_track_dxy_mtv[3].push_back(dxy);
        //     if (not tc_matched) tc_matched_track_dxy_mtv[3].push_back(dxy);
        // }
    }

    // std::cout <<  " fakeratestats:n_tc_not_matched: " << n_tc_not_matched <<  std::endl;
    // std::cout <<  " fakeratestats:good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
    // std::cout <<  " fakeratestats:n_tc_not_matched/good_track_candidates.size(): " << float(n_tc_not_matched)/float(good_track_candidates.size()) <<  std::endl;

}

void StudyMTVEfficiencyv2::doStudy_trackCandidateEfficiency(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents)
{
    // Each do study is performed per event

    // First clear all the output variables that will be used to fill the histograms for this event

    for (unsigned int i = 0; i < 6; ++i)
    {
        tc_matched_track_pt_mtv [i].clear();
        tc_all_track_pt_mtv     [i].clear();
        tc_matched_track_eta_mtv[i].clear();
        tc_all_track_eta_mtv    [i].clear();
        tc_matched_track_phi_mtv[i].clear();
        tc_all_track_phi_mtv    [i].clear();
        tc_matched_track_dxy_mtv[i].clear();
        tc_all_track_dxy_mtv    [i].clear();
    }

    // Barrel layer layer 0
    std::vector<SDL::TrackCandidate*> all_trackCandidates;
    const std::vector<SDL::Layer*> layerPtrs = event.getLayerPtrs();
    for (auto& layerPtr : layerPtrs)
    {
        // TrackCandidate ptrs
        const std::vector<SDL::TrackCandidate*>& trackCandidatePtrs = layerPtr->getTrackCandidatePtrs();
        for (auto& trackCandidatePtr : trackCandidatePtrs)
        {
            all_trackCandidates.push_back(trackCandidatePtr);
        }
    }


    //***************************************************************
    // Preselecting Reco Track Candidates that are worth looking into
    //***************************************************************

    // track candidates that match the particle of interest
    std::vector<unsigned int> good_track_candidates;
    std::vector<unsigned int> good_track_candidates_matchable;

    std::cout <<  " all_trackCandidates.size(): " << all_trackCandidates.size() <<  std::endl;

    // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
    for (unsigned int itc = 0; itc < all_trackCandidates.size(); ++itc)
    {
        SDL::TrackCandidate* tcCandPtr = all_trackCandidates[itc];
        good_track_candidates.push_back(itc);
        if (matchedSimTrkIdxs(tcCandPtr).size() > 0)
            good_track_candidates_matchable.push_back(itc);
    }

    // std::cout <<  " good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
    // std::cout <<  " good_track_candidates_matchable.size(): " << good_track_candidates_matchable.size() <<  std::endl;

    //***********************
    // Efficiency calculation
    //***********************

    std::cout <<  " simtrkevents.size(): " << simtrkevents.size() <<  std::endl;

    int nprocessed = 0;

    std::vector<unsigned int> matched_good_track_candidates;

    // DEBUG COUNTING
    int ndenom = 0;
    int ndenom_and = 0;
    int ndenom_anomaly = 0;
    for (auto& simtrkevent : simtrkevents)
    {
        // Unpack the tuple (sim_track_index, SDL::Event containing reco hits only matched to the given sim track)
        unsigned int& isimtrk = std::get<0>(simtrkevent);
        SDL::Event& trackevent = *(std::get<1>(simtrkevent));

        // Parse pt and eta of this track
        float pt = std::min((double) trk.sim_pt()[isimtrk], 49.999);
        float eta = trk.sim_eta()[isimtrk];
        float phi = trk.sim_phi()[isimtrk];
        float charge = trk.sim_q()[isimtrk];
        float dxy = trk.sim_pca_dxy()[isimtrk];
        float dz = trk.sim_pca_dz()[isimtrk];
        int pdgid = trk.sim_pdgId()[isimtrk];
        int bunch = trk.sim_bunchCrossing()[isimtrk];

        bool and_cond = (abs(dz) < 30 and abs(dxy) < 2.5 and abs(eta) < 2.4 and abs(pdgid) == 211 and bunch == 0);
        bool or_cond = (abs(dz) > 30 or abs(dxy) > 2.5 or abs(eta) > 2.4 or abs(pdgid) != 211 or bunch != 0);

        if (and_cond)
            ndenom_and++;

        if ((not or_cond) and (not and_cond))
        {
            ndenom_anomaly++;
            std::cout << "abs(dz): " << abs(dz) << std::endl;
            std::cout << "abs(dxy): " << abs(dxy) << std::endl;
            std::cout << "abs(eta): " << abs(eta) << std::endl;
            std::cout << "abs(pdgid): " << abs(pdgid) << std::endl;
            std::cout << "bunch" << bunch << std::endl;
        }

        if (or_cond)
            continue;

        ndenom++;
    }
    std::cout << "Denom " << ndenom << std::endl;
    std::cout << "Denom_and " << ndenom_and << std::endl;

    // Loop over track events
    for (auto& simtrkevent : simtrkevents)
    {

        // std::cout <<  " nprocessed: " << nprocessed <<  std::endl;
        nprocessed++;

        // Unpack the tuple (sim_track_index, SDL::Event containing reco hits only matched to the given sim track)
        unsigned int& isimtrk = std::get<0>(simtrkevent);
        SDL::Event& trackevent = *(std::get<1>(simtrkevent));

        // Create list of booleans
        bool tc_trackmatch_mtv = false;

        // Parse pt and eta of this track
        float pt = std::min((double) trk.sim_pt()[isimtrk], 49.999);
        float eta = trk.sim_eta()[isimtrk];
        float phi = trk.sim_phi()[isimtrk];
        float charge = trk.sim_q()[isimtrk];
        float dxy = trk.sim_pca_dxy()[isimtrk];
        float dz = trk.sim_pca_dz()[isimtrk];

        if (abs(dz) > dz_thresh)
            continue;

        int nmatch = 0;
        int nuniquematch = 0;
        std::vector<int> nbarrels; // unique list of nbarrel in the reco track cand

        std::vector<std::vector<unsigned int>> unique_list_of_module_ids;

        // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        for (auto& igood_tc_idxs : good_track_candidates_matchable)
        {

            SDL::TrackCandidate* tcCandPtr = all_trackCandidates[igood_tc_idxs];
            std::vector<unsigned int> hitidxs;
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->innerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            // std::vector<unsigned int> hitidxs;
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->idx());
            // hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->idx());
            if (isMTVMatch(isimtrk, hitidxs))
            {
                nmatch++;
                tc_trackmatch_mtv = true;
                if (std::find(matched_good_track_candidates.begin(), matched_good_track_candidates.end(), igood_tc_idxs) == matched_good_track_candidates.end())
                    matched_good_track_candidates.push_back(igood_tc_idxs);

                // int nbarrel = getNBarrelTrackCandidate(tcCandPtr);
                // if (std::find(nbarrels.begin(), nbarrels.end(), nbarrel) == nbarrels.end())
                // {
                //     nuniquematch++;
                //     nbarrels.push_back(nbarrel);
                // }
                std::vector<unsigned int> module_ids = getModuleListTrackCandidate(tcCandPtr);
                if (std::find(unique_list_of_module_ids.begin(), unique_list_of_module_ids.end(), module_ids) == unique_list_of_module_ids.end())
                {
                    nuniquematch++;
                    unique_list_of_module_ids.push_back(module_ids);
                }
            }

        }

        // // DEBUG========== CAN ERASE ME LATER
        // if (abs(eta) < 2.4 and pt > 1.0)
        // {
        //     // std::cout << "=============Found a sim track matched to multiple track candidate" << std::endl;
        //     std::cout <<  " nmatch: " << nmatch <<  std::endl;
        //     std::cout <<  " nuniquematch: " << nuniquematch <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " isimtrk: " << isimtrk <<  std::endl;
        //     std::cout <<  " pt: " << pt <<  std::endl;
        //     std::cout <<  " eta: " << eta <<  std::endl;
        //     std::cout <<  " phi: " << phi <<  std::endl;
        //     std::cout <<  " charge: " << charge <<  std::endl;
        // }
        // DEBUG========== CAN ERASE ME LATER

        // // DEBUG========== CAN ERASE ME LATER
        // // if (nmatch > 0)
        // if (pt > 1.5 and abs(eta) < 2.4)
        // {

        //     // std::cout << "=============Found a sim track matched to multiple track candidate" << std::endl;
        //     std::cout <<  " nmatch: " << nmatch <<  std::endl;
        //     std::cout <<  " nuniquematch: " << nuniquematch <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " isimtrk: " << isimtrk <<  std::endl;
        //     std::cout <<  " pt: " << pt <<  std::endl;
        //     std::cout <<  " eta: " << eta <<  std::endl;
        //     std::cout <<  " phi: " << phi <<  std::endl;
        //     std::cout <<  " charge: " << charge <<  std::endl;

        //     for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
        //     {
        //         std::cout <<  " trk.simhit_detId()[trk.sim_simHitIdx()[isimtrk][ith_hit]]: " << trk.simhit_detId()[trk.sim_simHitIdx()[isimtrk][ith_hit]] <<  std::endl;
        //     }

        //     for (auto& module_ids : unique_list_of_module_ids)
        //     {
        //         std::cout << "-----" << std::endl;
        //         for (auto& module_id : module_ids)
        //         {
        //             SDL::Module module(module_id);
        //             std::cout << module;
        //         }
        //     }

        //     // // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        //     // for (auto& igood_tc_idxs : good_track_candidates)
        //     // {

        //     //     SDL::TrackCandidate* tcCandPtr = all_trackCandidates[igood_tc_idxs];
        //     //     std::vector<unsigned int> hitidxs;
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //     //     hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //     //     if (isMTVMatch(isimtrk, hitidxs))
        //     //     {

        //     //         std::cout << tcCandPtr << std::endl;

        //     //     }

        //     // }

        //     std::cout << "============End printing" << std::endl;

        // }
        // // DEBUG==========




        // if (not tc_trackmatch_mtv and eta > 2.0 and eta < 2.2 and pt > 1.5)
        // {
        //     std::cout << "Found a simulated track with eta > 2.0 that did not match to any track candidate" << std::endl;
        //     printSimTrack(isimtrk);
        //     // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        //     for (auto& igood_tc_idxs : good_track_candidates)
        //     {

        //         SDL::TrackCandidate* tcCandPtr = all_trackCandidates[igood_tc_idxs];
        //         std::vector<unsigned int> hitidxs;
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->innerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
        //         hitidxs.push_back(tcCandPtr->outerTrackletPtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
        //         if (isMTVMatch(isimtrk, hitidxs, true))
        //         {
        //             tc_trackmatch_mtv = true;
        //         }

        //     }
        // }

        // *****************
        // Track Candidate
        // *****************

        float etamaxbound = 2.4;
        // float etamaxbound = 1.5;
        // float etamaxbound = 2.0;

        if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
        {

            tc_all_track_pt_mtv[0].push_back(pt);
            if (tc_trackmatch_mtv) tc_matched_track_pt_mtv[0].push_back(pt);

            if (abs(eta) < 0.8)
            {
                tc_all_track_pt_mtv[1].push_back(pt);
                if (tc_trackmatch_mtv) tc_matched_track_pt_mtv[1].push_back(pt);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tc_all_track_pt_mtv[2].push_back(pt);
                if (tc_trackmatch_mtv) tc_matched_track_pt_mtv[2].push_back(pt);
            }

            tc_all_track_pt_mtv[4].push_back(pt);
            if (nmatch > 1) tc_matched_track_pt_mtv[4].push_back(pt);

            tc_all_track_pt_mtv[5].push_back(pt);
            if (nuniquematch > 1) tc_matched_track_pt_mtv[5].push_back(pt);
        }

        if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
        {

            tc_all_track_eta_mtv[0].push_back(eta);
            if (tc_trackmatch_mtv) tc_matched_track_eta_mtv[0].push_back(eta);

            if (abs(eta) < 0.8)
            {
                tc_all_track_eta_mtv[1].push_back(eta);
                if (tc_trackmatch_mtv) tc_matched_track_eta_mtv[1].push_back(eta);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tc_all_track_eta_mtv[2].push_back(eta);
                if (tc_trackmatch_mtv) tc_matched_track_eta_mtv[2].push_back(eta);
            }

            tc_all_track_eta_mtv[4].push_back(eta);
            if (nmatch > 1) tc_matched_track_eta_mtv[4].push_back(eta);

            tc_all_track_eta_mtv[5].push_back(eta);
            if (nuniquematch > 1) tc_matched_track_eta_mtv[5].push_back(eta);

            tc_all_track_phi_mtv[5].push_back(phi);
            if (nuniquematch > 1) tc_matched_track_phi_mtv[5].push_back(phi);
        }

        if (pt > pt_thresh and abs(eta) < etamaxbound)
        {

            tc_all_track_dxy_mtv[0].push_back(dxy);
            if (tc_trackmatch_mtv) tc_matched_track_dxy_mtv[0].push_back(dxy);

            if (abs(eta) < 0.8)
            {
                tc_all_track_dxy_mtv[1].push_back(dxy);
                if (tc_trackmatch_mtv) tc_matched_track_dxy_mtv[1].push_back(dxy);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tc_all_track_dxy_mtv[2].push_back(dxy);
                if (tc_trackmatch_mtv) tc_matched_track_dxy_mtv[2].push_back(dxy);
            }

            tc_all_track_dxy_mtv[4].push_back(dxy);
            if (nmatch > 1) tc_matched_track_dxy_mtv[4].push_back(dxy);

            tc_all_track_dxy_mtv[5].push_back(dxy);
            if (nuniquematch > 1) tc_matched_track_dxy_mtv[5].push_back(dxy);
        }

    }

    int n_tc_not_matched = 0;

    for (auto& igood_tc_idxs : good_track_candidates)
    {
        SDL::TrackCandidate* tcCandPtr = all_trackCandidates[igood_tc_idxs];
        float betaPt_2nd_outer = tcCandPtr->outerTrackletBasePtr()->getRecoVar("betaPt_2nd");
        float betaPt_2nd_inner = tcCandPtr->innerTrackletBasePtr()->getRecoVar("betaPt_2nd");
        if (betaPt_2nd_outer == -999)
            betaPt_2nd_outer = betaPt_2nd_inner;
        float reco_pt = (abs(betaPt_2nd_outer) + abs(betaPt_2nd_inner)) / 2.;
        float reco_eta = tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->eta();

        bool tc_matched = false;
        if (std::find(matched_good_track_candidates.begin(), matched_good_track_candidates.end(), igood_tc_idxs) != matched_good_track_candidates.end())
            tc_matched = true;

        // === DEBUG
        if (not tc_matched)
        {
            n_tc_not_matched++;
        }

        std::vector<unsigned int> module_list = getModuleListTrackCandidate(tcCandPtr);
        std::vector<unsigned int> hitidx_list = getHitIdxListTrackCandidate(tcCandPtr);
        std::vector<unsigned int> simidx_list = allMatchedSimTrkIdxs(tcCandPtr);
        std::vector<         int> gstidx_list = matchedSimTrkIdxs(tcCandPtr);

        // if (not tc_matched)
        // {
        //     std::cout <<  " ==not matched track== " << std::endl;
        //     std::cout <<  " good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
        //     std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
        //     std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
        //     std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
        //     std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
        //     std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
        // }
        // else
        // {
        //     std::cout <<  " ==matched track== " << std::endl;
        //     std::cout <<  " good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
        //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
        //     std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
        //     std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
        //     std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
        //     std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
        //     std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
        //     std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
        // }

        // === DEBUG

        float etamaxbound = 2.4;

        // if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
        if (abs(reco_eta) < etamaxbound)
        {
            tc_all_track_pt_mtv[3].push_back(reco_pt);
            if (not tc_matched) tc_matched_track_pt_mtv[3].push_back(reco_pt);
        }

        // if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
        if (reco_pt > pt_thresh and abs(reco_eta) < etamaxbound)
        {
            tc_all_track_eta_mtv[3].push_back(reco_eta);
            if (not tc_matched) tc_matched_track_eta_mtv[3].push_back(reco_eta);
        }

        // if (pt > pt_thresh and abs(eta) < etamaxbound)
        // {
        //     tc_all_track_dxy_mtv[3].push_back(dxy);
        //     if (not tc_matched) tc_matched_track_dxy_mtv[3].push_back(dxy);
        // }
    }

    // std::cout <<  " fakeratestats:n_tc_not_matched: " << n_tc_not_matched <<  std::endl;
    // std::cout <<  " fakeratestats:good_track_candidates.size(): " << good_track_candidates.size() <<  std::endl;
    // std::cout <<  " fakeratestats:n_tc_not_matched/good_track_candidates.size(): " << float(n_tc_not_matched)/float(good_track_candidates.size()) <<  std::endl;

}

void StudyMTVEfficiencyv2::doStudy_trackletEfficiency(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents)
{
    // Each do study is performed per event

    // First clear all the output variables that will be used to fill the histograms for this event

    for (unsigned int i = 0; i < 6; ++i)
    {
        tl_matched_track_pt_mtv [i].clear();
        tl_all_track_pt_mtv     [i].clear();
        tl_matched_track_eta_mtv[i].clear();
        tl_all_track_eta_mtv    [i].clear();
        tl_matched_track_phi_mtv[i].clear();
        tl_all_track_phi_mtv    [i].clear();
        tl_matched_track_dxy_mtv[i].clear();
        tl_all_track_dxy_mtv    [i].clear();
    }

    std::vector<SDL::Tracklet*> all_tracklets;
    const std::vector<SDL::Layer*> layerPtrs = event.getLayerPtrs();
    for (auto& layerPtr : layerPtrs)
    {
        // Tracklet ptrs
        const std::vector<SDL::Tracklet*>& trackletPtrs = layerPtr->getTrackletPtrs();
        for (auto& trackletPtr : trackletPtrs)
        {
            all_tracklets.push_back(trackletPtr);
        }
    }

    // // Barrel layer layer 0
    // SDL::Layer barrelLayer0 = event.getLayer(1, SDL::Layer::Barrel);
    // SDL::Layer barrelLayer1 = event.getLayer(2, SDL::Layer::Barrel);
    // SDL::Layer barrelLayer2 = event.getLayer(3, SDL::Layer::Barrel);
    // // SDL::Layer barrelLayer3 = event.getLayer(4, SDL::Layer::Barrel);
    // // SDL::Layer barrelLayer4 = event.getLayer(5, SDL::Layer::Barrel);
    // SDL::Layer endcapLayer1 = event.getLayer(1, SDL::Layer::Barrel);
    // SDL::Layer endcapLayer2 = event.getLayer(2, SDL::Layer::Barrel);

    // std::vector<SDL::Tracklet*> all_tracklets;
    // all_tracklets.insert(all_tracklets.end(), barrelLayer0.getTrackletPtrs().begin(), barrelLayer0.getTrackletPtrs().end());
    // all_tracklets.insert(all_tracklets.end(), barrelLayer1.getTrackletPtrs().begin(), barrelLayer1.getTrackletPtrs().end());
    // all_tracklets.insert(all_tracklets.end(), barrelLayer2.getTrackletPtrs().begin(), barrelLayer2.getTrackletPtrs().end());
    // // all_tracklets.insert(all_tracklets.end(), barrelLayer3.getTrackletPtrs().begin(), barrelLayer3.getTrackletPtrs().end());
    // // all_tracklets.insert(all_tracklets.end(), barrelLayer4.getTrackletPtrs().begin(), barrelLayer4.getTrackletPtrs().end());
    // all_tracklets.insert(all_tracklets.end(), endcapLayer1.getTrackletPtrs().begin(), endcapLayer1.getTrackletPtrs().end());
    // all_tracklets.insert(all_tracklets.end(), endcapLayer2.getTrackletPtrs().begin(), endcapLayer2.getTrackletPtrs().end());

    //***************************************************************
    // Preselecting Reco Track Candidates that are worth looking into
    //***************************************************************

    // track candidates that match the particle of interest
    std::vector<unsigned int> good_tracklets;
    std::vector<unsigned int> good_tracklets_matchable;

    // Loop over the tracklet and ask whether the hits are matched with this track (MTV like)
    // for (unsigned int itl = 0; itl < barrelLayer0.getTrackletPtrs().size(); ++itl)
    for (unsigned int itl = 0; itl < all_tracklets.size(); ++itl)
    {
        SDL::Tracklet* tlPtr = all_tracklets[itl];
        good_tracklets.push_back(itl);
        if (matchedSimTrkIdxs(*tlPtr).size() > 0)
            good_tracklets_matchable.push_back(itl);
    }

    // std::cout <<  " good_tracklets.size(): " << good_tracklets.size() <<  std::endl;
    // std::cout <<  " good_tracklets_matchable.size(): " << good_tracklets_matchable.size() <<  std::endl;

    //***********************
    // Efficiency calculation
    //***********************

    // std::cout <<  " simtrkevents.size(): " << simtrkevents.size() <<  std::endl;

    int nprocessed = 0;

    std::vector<unsigned int> matched_good_tracklets;

    // Loop over track events
    for (auto& simtrkevent : simtrkevents)
    {

        // std::cout <<  " nprocessed: " << nprocessed <<  std::endl;
        nprocessed++;

        // Unpack the tuple (sim_track_index, SDL::Event containing reco hits only matched to the given sim track)
        unsigned int& isimtrk = std::get<0>(simtrkevent);
        SDL::Event& trackevent = *(std::get<1>(simtrkevent));

        // Create list of booleans
        bool tl_trackmatch_mtv = false;

        // Parse pt and eta of this track
        float pt = std::min((double) trk.sim_pt()[isimtrk], 49.999);
        float eta = trk.sim_eta()[isimtrk];
        float phi = trk.sim_phi()[isimtrk];
        float charge = trk.sim_q()[isimtrk];
        float dxy = trk.sim_pca_dxy()[isimtrk];
        float dz = trk.sim_pca_dz()[isimtrk];

        if (abs(dz) > dz_thresh)
            continue;

        int nmatch = 0;
        int nuniquematch = 0;
        std::vector<int> nbarrels; // unique list of nbarrel in the reco track cand

        // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        for (auto& igood_tc_idxs : good_tracklets_matchable)
        {

            SDL::Tracklet* tlPtr = all_tracklets[igood_tc_idxs];
            std::vector<unsigned int> hitidxs;
            hitidxs.push_back(tlPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tlPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            if (isMTVMatch(isimtrk, hitidxs))
            {
                nmatch++;
                tl_trackmatch_mtv = true;
                if (std::find(matched_good_tracklets.begin(), matched_good_tracklets.end(), igood_tc_idxs) == matched_good_tracklets.end())
                    matched_good_tracklets.push_back(igood_tc_idxs);

            }

        }

        // *****************
        // Tracklets
        // *****************

        float etamaxbound = 2.4;

        if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
        {

            tl_all_track_pt_mtv[0].push_back(pt);
            if (tl_trackmatch_mtv) tl_matched_track_pt_mtv[0].push_back(pt);

            if (abs(eta) < 0.8)
            {
                tl_all_track_pt_mtv[1].push_back(pt);
                if (tl_trackmatch_mtv) tl_matched_track_pt_mtv[1].push_back(pt);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tl_all_track_pt_mtv[2].push_back(pt);
                if (tl_trackmatch_mtv) tl_matched_track_pt_mtv[2].push_back(pt);
            }

        }

        if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
        {

            tl_all_track_eta_mtv[0].push_back(eta);
            if (tl_trackmatch_mtv) tl_matched_track_eta_mtv[0].push_back(eta);

            if (abs(eta) < 0.8)
            {
                tl_all_track_eta_mtv[1].push_back(eta);
                if (tl_trackmatch_mtv) tl_matched_track_eta_mtv[1].push_back(eta);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tl_all_track_eta_mtv[2].push_back(eta);
                if (tl_trackmatch_mtv) tl_matched_track_eta_mtv[2].push_back(eta);
            }

        }

        if (pt > pt_thresh and abs(eta) < etamaxbound)
        {

            tl_all_track_dxy_mtv[0].push_back(dxy);
            if (tl_trackmatch_mtv) tl_matched_track_dxy_mtv[0].push_back(dxy);

            if (abs(eta) < 0.8)
            {
                tl_all_track_dxy_mtv[1].push_back(dxy);
                if (tl_trackmatch_mtv) tl_matched_track_dxy_mtv[1].push_back(dxy);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tl_all_track_dxy_mtv[2].push_back(dxy);
                if (tl_trackmatch_mtv) tl_matched_track_dxy_mtv[2].push_back(dxy);
            }

        }

    }

    // int n_tc_not_matched = 0;

    // for (auto& igood_tc_idxs : good_tracklets)
    // {
    //     SDL::Tracklet* tcCandPtr = barrelLayer0.getTrackletPtrs()[igood_tc_idxs];
    //     float betaPt_2nd_outer = tcCandPtr->outerTrackletBasePtr()->getRecoVar("betaPt_2nd");
    //     float betaPt_2nd_inner = tcCandPtr->innerTrackletBasePtr()->getRecoVar("betaPt_2nd");
    //     if (betaPt_2nd_outer == -999)
    //         betaPt_2nd_outer = betaPt_2nd_inner;
    //     float reco_pt = (abs(betaPt_2nd_outer) + abs(betaPt_2nd_inner)) / 2.;
    //     float reco_eta = tcCandPtr->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->eta();

    //     bool tc_matched = false;
    //     if (std::find(matched_good_tracklets.begin(), matched_good_tracklets.end(), igood_tc_idxs) != matched_good_tracklets.end())
    //         tc_matched = true;

    //     // === DEBUG
    //     if (not tc_matched)
    //     {
    //         n_tc_not_matched++;
    //     }

    //     std::vector<unsigned int> module_list = getModuleListTracklet(tcCandPtr);
    //     std::vector<unsigned int> hitidx_list = getHitIdxListTracklet(tcCandPtr);
    //     std::vector<unsigned int> simidx_list = allMatchedSimTrkIdxs(tcCandPtr);
    //     std::vector<         int> gstidx_list = matchedSimTrkIdxs(tcCandPtr);

    //     if (not tc_matched)
    //     {
    //         std::cout <<  " ==not matched track== " << std::endl;
    //         std::cout <<  " good_tracklets.size(): " << good_tracklets.size() <<  std::endl;
    //         std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
    //         std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
    //         std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
    //         std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
    //         std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
    //         std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
    //         std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
    //     }
    //     else
    //     {
    //         std::cout <<  " ==matched track== " << std::endl;
    //         std::cout <<  " good_tracklets.size(): " << good_tracklets.size() <<  std::endl;
    //         std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
    //         std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
    //         std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
    //         std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
    //         std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
    //         std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
    //         std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
    //     }

    //     // === DEBUG

    //     float etamaxbound = 2.4;

    //     // if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
    //     if (abs(reco_eta) < etamaxbound)
    //     {
    //         tc_all_track_pt_mtv[3].push_back(reco_pt);
    //         if (not tc_matched) tc_matched_track_pt_mtv[3].push_back(reco_pt);
    //     }

    //     // if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
    //     if (reco_pt > pt_thresh and abs(reco_eta) < etamaxbound)
    //     {
    //         tc_all_track_eta_mtv[3].push_back(reco_eta);
    //         if (not tc_matched) tc_matched_track_eta_mtv[3].push_back(reco_eta);
    //     }

    //     // if (pt > pt_thresh and abs(eta) < etamaxbound)
    //     // {
    //     //     tc_all_track_dxy_mtv[3].push_back(dxy);
    //     //     if (not tc_matched) tc_matched_track_dxy_mtv[3].push_back(dxy);
    //     // }
    // }

    // std::cout <<  " fakeratestats:n_tc_not_matched: " << n_tc_not_matched <<  std::endl;
    // std::cout <<  " fakeratestats:good_tracklets.size(): " << good_tracklets.size() <<  std::endl;
    // std::cout <<  " fakeratestats:n_tc_not_matched/good_tracklets.size(): " << float(n_tc_not_matched)/float(good_tracklets.size()) <<  std::endl;

}

void StudyMTVEfficiencyv2::doStudy_tripletEfficiency(SDL::Event& event, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents)
{
    // Each do study is performed per event

    // First clear all the output variables that will be used to fill the histograms for this event

    for (unsigned int i = 0; i < 6; ++i)
    {
        tp_matched_track_pt_mtv [i].clear();
        tp_all_track_pt_mtv     [i].clear();
        tp_matched_track_eta_mtv[i].clear();
        tp_all_track_eta_mtv    [i].clear();
        tp_matched_track_phi_mtv[i].clear();
        tp_all_track_phi_mtv    [i].clear();
        tp_matched_track_dxy_mtv[i].clear();
        tp_all_track_dxy_mtv    [i].clear();
    }

    std::vector<SDL::Triplet*> all_triplets;
    const std::vector<SDL::Layer*> layerPtrs = event.getLayerPtrs();
    for (auto& layerPtr : layerPtrs)
    {
        // Triplet ptrs
        const std::vector<SDL::Triplet*>& tripletPtrs = layerPtr->getTripletPtrs();
        for (auto& tripletPtr : tripletPtrs)
        {
            all_triplets.push_back(tripletPtr);
        }
    }

    // // Barrel layer layer 0
    // SDL::Layer barrelLayer0 = event.getLayer(1, SDL::Layer::Barrel);
    // SDL::Layer barrelLayer1 = event.getLayer(2, SDL::Layer::Barrel);
    // SDL::Layer barrelLayer2 = event.getLayer(3, SDL::Layer::Barrel);
    // // SDL::Layer barrelLayer3 = event.getLayer(4, SDL::Layer::Barrel);
    // // SDL::Layer barrelLayer4 = event.getLayer(5, SDL::Layer::Barrel);
    // SDL::Layer endcapLayer1 = event.getLayer(1, SDL::Layer::Barrel);
    // SDL::Layer endcapLayer2 = event.getLayer(2, SDL::Layer::Barrel);

    // std::vector<SDL::Triplet*> all_triplets;
    // all_triplets.insert(all_triplets.end(), barrelLayer0.getTripletPtrs().begin(), barrelLayer0.getTripletPtrs().end());
    // all_triplets.insert(all_triplets.end(), barrelLayer1.getTripletPtrs().begin(), barrelLayer1.getTripletPtrs().end());
    // all_triplets.insert(all_triplets.end(), barrelLayer2.getTripletPtrs().begin(), barrelLayer2.getTripletPtrs().end());
    // // all_triplets.insert(all_triplets.end(), barrelLayer3.getTripletPtrs().begin(), barrelLayer3.getTripletPtrs().end());
    // // all_triplets.insert(all_triplets.end(), barrelLayer4.getTripletPtrs().begin(), barrelLayer4.getTripletPtrs().end());
    // all_triplets.insert(all_triplets.end(), endcapLayer1.getTripletPtrs().begin(), endcapLayer1.getTripletPtrs().end());
    // all_triplets.insert(all_triplets.end(), endcapLayer2.getTripletPtrs().begin(), endcapLayer2.getTripletPtrs().end());

    //***************************************************************
    // Preselecting Reco Track Candidates that are worth looking into
    //***************************************************************

    // track candidates that match the particle of interest
    std::vector<unsigned int> good_triplets;
    std::vector<unsigned int> good_triplets_matchable;

    // Loop over the triplet and ask whether the hits are matched with this track (MTV like)
    // for (unsigned int itl = 0; itl < barrelLayer0.getTripletPtrs().size(); ++itl)
    for (unsigned int itl = 0; itl < all_triplets.size(); ++itl)
    {
        SDL::Triplet* tlPtr = all_triplets[itl];
        good_triplets.push_back(itl);
        if (matchedSimTrkIdxs(tlPtr).size() > 0)
            good_triplets_matchable.push_back(itl);
    }

    // std::cout <<  " good_triplets.size(): " << good_triplets.size() <<  std::endl;
    // std::cout <<  " good_triplets_matchable.size(): " << good_triplets_matchable.size() <<  std::endl;

    //***********************
    // Efficiency calculation
    //***********************

    // std::cout <<  " simtrkevents.size(): " << simtrkevents.size() <<  std::endl;

    int nprocessed = 0;

    std::vector<unsigned int> matched_good_triplets;

    // Loop over track events
    for (auto& simtrkevent : simtrkevents)
    {

        // std::cout <<  " nprocessed: " << nprocessed <<  std::endl;
        nprocessed++;

        // Unpack the tuple (sim_track_index, SDL::Event containing reco hits only matched to the given sim track)
        unsigned int& isimtrk = std::get<0>(simtrkevent);
        SDL::Event& trackevent = *(std::get<1>(simtrkevent));

        // Create list of booleans
        bool tp_trackmatch_mtv = false;

        // Parse pt and eta of this track
        float pt = std::min((double) trk.sim_pt()[isimtrk], 49.999);
        float eta = trk.sim_eta()[isimtrk];
        float phi = trk.sim_phi()[isimtrk];
        float charge = trk.sim_q()[isimtrk];
        float dxy = trk.sim_pca_dxy()[isimtrk];
        float dz = trk.sim_pca_dz()[isimtrk];

        if (abs(dz) > dz_thresh)
            continue;

        int nmatch = 0;
        int nuniquematch = 0;
        std::vector<int> nbarrels; // unique list of nbarrel in the reco track cand

        // Loop over the track candidates and ask whether the hits are matched with this track (MTV like)
        for (auto& igood_tc_idxs : good_triplets_matchable)
        {

            SDL::Triplet* tlPtr = all_triplets[igood_tc_idxs];
            std::vector<unsigned int> hitidxs;
            hitidxs.push_back(tlPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tlPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            hitidxs.push_back(tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
            hitidxs.push_back(tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
            if (isMTVMatch(isimtrk, hitidxs))
            {
                nmatch++;
                tp_trackmatch_mtv = true;
                if (std::find(matched_good_triplets.begin(), matched_good_triplets.end(), igood_tc_idxs) == matched_good_triplets.end())
                    matched_good_triplets.push_back(igood_tc_idxs);

            }

        }

        // *****************
        // Triplets
        // *****************

        float etamaxbound = 2.4;

        if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
        {

            tp_all_track_pt_mtv[0].push_back(pt);
            if (tp_trackmatch_mtv) tp_matched_track_pt_mtv[0].push_back(pt);

            if (abs(eta) < 0.8)
            {
                tp_all_track_pt_mtv[1].push_back(pt);
                if (tp_trackmatch_mtv) tp_matched_track_pt_mtv[1].push_back(pt);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tp_all_track_pt_mtv[2].push_back(pt);
                if (tp_trackmatch_mtv) tp_matched_track_pt_mtv[2].push_back(pt);
            }

        }

        if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
        {

            tp_all_track_eta_mtv[0].push_back(eta);
            if (tp_trackmatch_mtv) tp_matched_track_eta_mtv[0].push_back(eta);

            if (abs(eta) < 0.8)
            {
                tp_all_track_eta_mtv[1].push_back(eta);
                if (tp_trackmatch_mtv) tp_matched_track_eta_mtv[1].push_back(eta);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tp_all_track_eta_mtv[2].push_back(eta);
                if (tp_trackmatch_mtv) tp_matched_track_eta_mtv[2].push_back(eta);
            }

        }

        if (pt > pt_thresh and abs(eta) < etamaxbound)
        {

            tp_all_track_dxy_mtv[0].push_back(dxy);
            if (tp_trackmatch_mtv) tp_matched_track_dxy_mtv[0].push_back(dxy);

            if (abs(eta) < 0.8)
            {
                tp_all_track_dxy_mtv[1].push_back(dxy);
                if (tp_trackmatch_mtv) tp_matched_track_dxy_mtv[1].push_back(dxy);
            }
            if (abs(eta) > 1.0 and abs(eta) < etamaxbound)
            {
                tp_all_track_dxy_mtv[2].push_back(dxy);
                if (tp_trackmatch_mtv) tp_matched_track_dxy_mtv[2].push_back(dxy);
            }

        }

    }

    // int n_tc_not_matched = 0;

    // for (auto& igood_tc_idxs : good_triplets)
    // {
    //     SDL::Triplet* tcCandPtr = barrelLayer0.getTripletPtrs()[igood_tc_idxs];
    //     float betaPt_2nd_outer = tcCandPtr->outerTripletBasePtr()->getRecoVar("betaPt_2nd");
    //     float betaPt_2nd_inner = tcCandPtr->innerTripletBasePtr()->getRecoVar("betaPt_2nd");
    //     if (betaPt_2nd_outer == -999)
    //         betaPt_2nd_outer = betaPt_2nd_inner;
    //     float reco_pt = (abs(betaPt_2nd_outer) + abs(betaPt_2nd_inner)) / 2.;
    //     float reco_eta = tcCandPtr->outerTripletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->eta();

    //     bool tc_matched = false;
    //     if (std::find(matched_good_triplets.begin(), matched_good_triplets.end(), igood_tc_idxs) != matched_good_triplets.end())
    //         tc_matched = true;

    //     // === DEBUG
    //     if (not tc_matched)
    //     {
    //         n_tc_not_matched++;
    //     }

    //     std::vector<unsigned int> module_list = getModuleListTriplet(tcCandPtr);
    //     std::vector<unsigned int> hitidx_list = getHitIdxListTriplet(tcCandPtr);
    //     std::vector<unsigned int> simidx_list = allMatchedSimTrkIdxs(tcCandPtr);
    //     std::vector<         int> gstidx_list = matchedSimTrkIdxs(tcCandPtr);

    //     if (not tc_matched)
    //     {
    //         std::cout <<  " ==not matched track== " << std::endl;
    //         std::cout <<  " good_triplets.size(): " << good_triplets.size() <<  std::endl;
    //         std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
    //         std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
    //         std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
    //         std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
    //         std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
    //         std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
    //         std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
    //     }
    //     else
    //     {
    //         std::cout <<  " ==matched track== " << std::endl;
    //         std::cout <<  " good_triplets.size(): " << good_triplets.size() <<  std::endl;
    //         std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;
    //         std::cout <<  " reco_pt: " << reco_pt <<  std::endl;
    //         std::cout <<  " reco_eta: " << reco_eta <<  std::endl;
    //         std::cout << " module: "; for (auto& module : module_list) std::cout << " " << module; std::cout << std::endl;
    //         std::cout << " hitidx: "; for (auto& hitidx : hitidx_list) std::cout << " " << hitidx; std::cout << std::endl;
    //         std::cout << " simidx: "; for (auto& simidx : simidx_list) std::cout << " " << simidx; std::cout << std::endl;
    //         std::cout << " gstidx: "; for (auto& gstidx : gstidx_list) std::cout << " " << gstidx; std::cout << std::endl;
    //     }

    //     // === DEBUG

    //     float etamaxbound = 2.4;

    //     // if (abs(eta) < etamaxbound and abs(dxy) < dxy_thresh)
    //     if (abs(reco_eta) < etamaxbound)
    //     {
    //         tc_all_track_pt_mtv[3].push_back(reco_pt);
    //         if (not tc_matched) tc_matched_track_pt_mtv[3].push_back(reco_pt);
    //     }

    //     // if (pt > pt_thresh and abs(dxy) < dxy_thresh and abs(eta) < etamaxbound)
    //     if (reco_pt > pt_thresh and abs(reco_eta) < etamaxbound)
    //     {
    //         tc_all_track_eta_mtv[3].push_back(reco_eta);
    //         if (not tc_matched) tc_matched_track_eta_mtv[3].push_back(reco_eta);
    //     }

    //     // if (pt > pt_thresh and abs(eta) < etamaxbound)
    //     // {
    //     //     tc_all_track_dxy_mtv[3].push_back(dxy);
    //     //     if (not tc_matched) tc_matched_track_dxy_mtv[3].push_back(dxy);
    //     // }
    // }

    // std::cout <<  " fakeratestats:n_tc_not_matched: " << n_tc_not_matched <<  std::endl;
    // std::cout <<  " fakeratestats:good_triplets.size(): " << good_triplets.size() <<  std::endl;
    // std::cout <<  " fakeratestats:n_tc_not_matched/good_triplets.size(): " << float(n_tc_not_matched)/float(good_triplets.size()) <<  std::endl;

}

