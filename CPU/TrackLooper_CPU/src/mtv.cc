#include "mtv.h"

void mtv()
{
    // List of studies to perform
    std::vector<Study*> studies;

    // pt_boundaries
    std::vector<float> pt_boundaries = getPtBounds();

    studies.push_back(new StudyMTVEfficiencyv2("MTVEfficiency", pt_boundaries, ana.pdg_id));

    // book the studies
    for (auto& study : studies)
    {
        study->bookStudy();
    }

    // Book Histograms
    ana.cutflow.bookHistograms(ana.histograms); // if just want to book everywhere

    // Load various maps used in the SDL reconstruction
    loadMaps();

    // Looping input file
    while (ana.looper.nextEvent())
    {

        if (not goodEvent())
            continue;

        // *************************************************
        // Reconstructed hits and formation of mini-doublets
        // *************************************************

        // // === DEBUG where i apply delta R requirement on all pairs of sim tracks with pt > 1.2 and eta < 2.4
        // float mindr = 99999;
        // for (auto&& [itrk, itrkdata] : iter::enumerate(iter::zip(trk.sim_pt(), trk.sim_eta(), trk.sim_phi())))
        // {
        //     auto&& [ipt, ieta, iphi] = itrkdata;
        //     for (auto&& [jtrk, jtrkdata] : iter::enumerate(iter::zip(trk.sim_pt(), trk.sim_eta(), trk.sim_phi())))
        //     {
        //         if (itrk == jtrk)
        //             continue;
        //         auto&& [jpt, jeta, jphi] = jtrkdata;
        //         float deta = (jeta - ieta);
        //         float dphitemp = (jphi - iphi);
        //         float f = dphitemp;
        //         while (f >= M_PI) f -= 2. * M_PI;
        //         while (f < -M_PI) f += 2. * M_PI;
        //         float dphi = f;
        //         float dr = sqrt(deta*deta + dphi*dphi);
        //         // std::cout <<  " deta: " << deta <<  " dphitemp: " << dphitemp <<  " dphi: " << dphi <<  " dr: " << dr <<  std::endl;
        //         if (dr < mindr)
        //         {
        //             mindr = dr;
        //         }
        //     }
        // }
        // if (mindr < 0.055)
        //     continue;
        // // === DEBUG where i apply delta R requirement on all pairs of sim tracks with pt > 1.2 and eta < 2.4

        // Main instance that will hold modules, hits, minidoublets, etc. (i.e. main data structure)
        SDL::Event event;

        // Add hits to the event
        addOuterTrackerHits(event);

        // Run the SDL reconstruction on the event
        // runSDLTest_v1(event);
        runSDLTest_v2(event);

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

            // // Add reco hits associated with the given sim track to the SDL::Event container
            // addOuterTrackerHitsFromSimTrack((*trackevent), isimtrk);

            // // Run SDL reconstruction on the event
            // runSDL((*trackevent));

            // Push to the vector so we have a data-base of per hit, mini-doublets
            simtrkevents.push_back(std::make_tuple(isimtrk, trackevent));

        }


        // ********************************************************************************************
        // Perform various studies with reco events and sim-track-matched-reco-hits-based mini-doublets
        // ********************************************************************************************

        for (auto& study : studies)
        {
            study->doStudy(event, simtrkevents);
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
