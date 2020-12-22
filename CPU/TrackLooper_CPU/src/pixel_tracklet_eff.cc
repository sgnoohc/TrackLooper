#include "pixel_tracklet_eff.h"

void pixel_tracklet_eff()
{

    // pt_boundaries
    std::vector<float> pt_boundaries = getPtBounds();

    // histogram variable
    RooUtil::TTreeX var("variable", "Variable");
    var.createBranch<vector<float>>("pt_denom");
    var.createBranch<vector<float>>("pt_denom_P");
    var.createBranch<vector<float>>("pt_numer_P");
    var.createBranch<vector<float>>("pt_numer_PB1B2");
    var.createBranch<vector<float>>("pt_numer_PB2B3");
    var.createBranch<vector<float>>("pt_numer_PB3B4");
    var.createBranch<vector<float>>("pt_numer_PB4B5");
    var.createBranch<vector<float>>("pt_numer_PB5B6");

    ana.histograms.addVecHistogram("pixel_pt_denom" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_denom"); } );
    ana.histograms.addVecHistogram("pixel_pt_denom_P" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_denom_P"); } );
    ana.histograms.addVecHistogram("pixel_pt_numer_P" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_numer_P"); } );
    ana.histograms.addVecHistogram("pixel_pt_numer_PB1B2" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_numer_PB1B2"); } );
    ana.histograms.addVecHistogram("pixel_pt_numer_PB2B3" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_numer_PB2B3"); } );
    ana.histograms.addVecHistogram("pixel_pt_numer_PB3B4" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_numer_PB3B4"); } );
    ana.histograms.addVecHistogram("pixel_pt_numer_PB4B5" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_numer_PB4B5"); } );
    ana.histograms.addVecHistogram("pixel_pt_numer_PB5B6" , pt_boundaries, [&]() { return var.getBranchLazy<vector<float>>("pt_numer_PB5B6"); } );

    var.createBranch<vector<float>>("dBeta_numer_PB1B2");
    var.createBranch<vector<float>>("dBeta_numer_PB2B3");
    var.createBranch<vector<float>>("dBeta_numer_PB3B4");
    var.createBranch<vector<float>>("dBeta_numer_PB4B5");
    var.createBranch<vector<float>>("dBeta_numer_PB5B6");

    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB1B2" , 50, -0.1, 0.1, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB1B2"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB2B3" , 50, -0.1, 0.1, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB2B3"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB3B4" , 50, -0.1, 0.1, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB3B4"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB4B5" , 50, -0.1, 0.1, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB4B5"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB5B6" , 50, -0.1, 0.1, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB5B6"); } );

    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB1B2_zoom" , 50, -0.04, 0.04, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB1B2"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB2B3_zoom" , 50, -0.04, 0.04, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB2B3"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB3B4_zoom" , 50, -0.04, 0.04, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB3B4"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB4B5_zoom" , 50, -0.04, 0.04, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB4B5"); } );
    ana.histograms.addVecHistogram("pixel_dBeta_numer_PB5B6_zoom" , 50, -0.04, 0.04, [&]() { return var.getBranchLazy<vector<float>>("dBeta_numer_PB5B6"); } );

    // Book Histograms
    ana.cutflow.bookHistograms(ana.histograms); // if just want to book everywhere

    // Load various maps used in the SDL reconstruction
    loadMaps();

    // Looping input file
    while (ana.looper.nextEvent())
    {

        if (not goodEvent())
            continue;

        // *********************************************************************
        // SDL::Event from each sim track (using only hits from given sim track)
        // *********************************************************************

        // Each SDL::Event object in simtrkevents will hold single sim-track related hits
        // It will be a vector of tuple of <sim_track_index, SDL::Event*>.
        std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents;

        // Loop over sim-tracks that is from in time (bx = 0) tracks with pdgid matching (against ana.pdg_id) and per sim-track aggregate reco hits
        // and only use those hits, and run SDL on them
        for (auto&& isimtrk : iter::filter([&](int i) { return inTimeTrackWithPdgId(i, ana.pdg_id); }, iter::range(trk.sim_pt().size())))
        {

            // event just for this track
            SDL::Event* trackevent = new SDL::Event();

            // Add reco hits associated with the given sim track to the SDL::Event container
            addOuterTrackerHitsFromSimTrack((*trackevent), isimtrk);

            // Add Pixel Seeds to the event
            addPixelSegments((*trackevent), isimtrk);

            // if (fabs(trk.sim_eta()[isimtrk]) > 0.8 or trk.sim_pt()[isimtrk] < 1.5)
            //     continue;

            if (ana.verbose != 0)
            {
                std::cout << " sim trk pt = " << trk.sim_pt()[isimtrk] << std::endl;
                std::cout << " sim trk eta = " << trk.sim_eta()[isimtrk] << std::endl;
            }

            // Run SDL reconstruction on the event
            runSDLTest_PixelSegments((*trackevent));

            // Push to the vector so we have a data-base of per hit, mini-doublets
            simtrkevents.push_back(std::make_tuple(isimtrk, trackevent));

        }


        // ********************************************************************************************
        // Perform various studies with reco events and sim-track-matched-reco-hits-based mini-doublets
        // ********************************************************************************************

        var.clear();

        // Loop over track events
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

            if (fabs(eta) > 0.8)
                continue;

            if (fabs(dxy) > 2.5)
                continue;

            if (fabs(dz) > 30.)
                continue;

            if (not hasAll12HitsWithNBarrelUsingModuleMap(isimtrk, 6))
                continue;

            // Studying Pixel Segment and Tracklets
            bool pix_sg_trackmatch_by_layer = trackevent.getPixelLayer().getSegmentPtrs().size() > 0;
            if (ana.verbose != 0)
            {
                std::cout <<  " pix_sg_trackmatch_by_layer: " << pix_sg_trackmatch_by_layer <<  std::endl;
                std::cout <<  " trackevent.getPixelLayer().getTrackletPtrs().size(): " << trackevent.getPixelLayer().getTrackletPtrs().size() <<  std::endl;
            }

            var.pushbackToBranch<float>("pt_denom_P", pt);

            if (not pix_sg_trackmatch_by_layer)
                continue;

            var.pushbackToBranch<float>("pt_denom", pt);

            bool has_P = pix_sg_trackmatch_by_layer;
            bool has_PB1B2 = false;
            bool has_PB2B3 = false;
            bool has_PB3B4 = false;
            bool has_PB4B5 = false;
            bool has_PB5B6 = false;

            for (auto& tlPtr : trackevent.getPixelLayer().getTrackletPtrs())
            {

                tlPtr->runTrackletAlgo(SDL::Default_TLAlgo);
                bool pass = tlPtr->passesTrackletAlgo(SDL::Default_TLAlgo);
                int passbit = tlPtr->getPassBitsDefaultAlgo();
                float _pt = trk.sim_pt()[isimtrk];
                float _eta = trk.sim_eta()[isimtrk];
                float _dxy = trk.sim_pca_dxy()[isimtrk];
                float _dz = trk.sim_pca_dz()[isimtrk];
                float _charge = trk.sim_q()[isimtrk];

                // if (not pass)
                //     continue;
                // if (passbit < 63)
                //     continue;

                float dBeta = tlPtr->getRecoVar("dBeta");
                float dBetaCut2 = tlPtr->getRecoVar("dBetaCut2");
                float betaOut = tlPtr->getRecoVar("betaOut");
                float pt_betaIn = tlPtr->getRecoVar("pt_betaIn");

                if (ana.verbose != 0 and pt > 1.5 and fabs(eta) < 0.8)
                {
                    std::cout <<  " _pt: " << _pt <<  " _eta: " << _eta <<  " _dxy: " << _dxy <<  " _dz: " << _dz <<  " _charge: " << _charge <<  std::endl;
                    std::cout <<  " tlPtr->getRecoVar('rawBetaIn'): " << tlPtr->getRecoVar("rawBetaIn") <<  std::endl;
                    std::cout <<  " tlPtr->getRecoVar('betaIn'): " << tlPtr->getRecoVar("betaIn") <<  std::endl;
                    std::cout <<  " tlPtr->getRecoVar('betaOut'): " << tlPtr->getRecoVar("betaOut") <<  std::endl;
                    std::cout <<  " tlPtr->getDeltaBeta(): " << tlPtr->getDeltaBeta() <<  std::endl;
                    std::cout <<  " tlPtr->getDeltaBetaCut(): " << tlPtr->getDeltaBetaCut() <<  std::endl;
                    std::cout <<  " pass: " << pass <<  std::endl;
                    std::bitset<8> x(passbit);
                    std::cout <<  " passbit: " << x <<  std::endl;
                    std::cout <<  " tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule().layer(): " << tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule().layer() <<  std::endl;
                    std::cout <<  " tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule().subdet(): " << tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule().subdet() <<  std::endl;
                    std::cout <<  " tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule().layer(): " << tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule().layer() <<  std::endl;
                    std::cout <<  " tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule().subdet(): " << tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule().subdet() <<  std::endl;
                }

                int layer1 = tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule().layer();
                int subdet1 = tlPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule().subdet();
                int layer2 = tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule().layer();
                int subdet2 = tlPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule().subdet();

                if (layer1 == 1 and layer2 == 2 and subdet1 == 5 and subdet2 == 5) has_PB1B2 = true;
                if (layer1 == 2 and layer2 == 3 and subdet1 == 5 and subdet2 == 5) has_PB2B3 = true;
                if (layer1 == 3 and layer2 == 4 and subdet1 == 5 and subdet2 == 5) has_PB3B4 = true;
                if (layer1 == 4 and layer2 == 5 and subdet1 == 5 and subdet2 == 5) has_PB4B5 = true;
                if (layer1 == 5 and layer2 == 6 and subdet1 == 5 and subdet2 == 5) has_PB5B6 = true;

                if (layer1 == 1 and layer2 == 2 and subdet1 == 5 and subdet2 == 5) var.pushbackToBranch<float>("dBeta_numer_PB1B2", dBeta);
                if (layer1 == 2 and layer2 == 3 and subdet1 == 5 and subdet2 == 5) var.pushbackToBranch<float>("dBeta_numer_PB2B3", dBeta);
                if (layer1 == 3 and layer2 == 4 and subdet1 == 5 and subdet2 == 5) var.pushbackToBranch<float>("dBeta_numer_PB3B4", dBeta);
                if (layer1 == 4 and layer2 == 5 and subdet1 == 5 and subdet2 == 5) var.pushbackToBranch<float>("dBeta_numer_PB4B5", dBeta);
                if (layer1 == 5 and layer2 == 6 and subdet1 == 5 and subdet2 == 5) var.pushbackToBranch<float>("dBeta_numer_PB5B6", dBeta);

            }

            if (has_P) var.pushbackToBranch<float>("pt_numer_P", pt);
            if (has_PB1B2) var.pushbackToBranch<float>("pt_numer_PB1B2", pt);
            if (has_PB2B3) var.pushbackToBranch<float>("pt_numer_PB2B3", pt);
            if (has_PB3B4) var.pushbackToBranch<float>("pt_numer_PB3B4", pt);
            if (has_PB4B5) var.pushbackToBranch<float>("pt_numer_PB4B5", pt);
            if (has_PB5B6) var.pushbackToBranch<float>("pt_numer_PB5B6", pt);

        }

        // ************************************************
        // Now fill all the histograms booked by each study
        // ************************************************

        // Fill all the histograms
        ana.cutflow.fill();

        // <--------------------------
        // <--------------------------
        // <--------------------------
    }

    // Writing output file
    ana.cutflow.saveOutput();

    // Writing ttree output to file
    ana.output_ttree->Write();

    // The below can be sometimes crucial
    delete ana.output_tfile;

}
