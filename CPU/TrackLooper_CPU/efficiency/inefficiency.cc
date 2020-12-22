#include "inefficiency.h"

#include "sdl_types.h"

int main(int argc, char** argv)
{

    // Parse arguments
    parseArguments(argc, argv);

    // Initialize input and output root files
    initializeInputsAndOutputs();

    ana.tx.createBranch<vector<float>>("tc_types");
    ana.tx.createBranch<vector<float>>("tc_uniq_types");
    ana.histograms.addVecHistogram("tc_types", 119, 0, 119, [&]() { return ana.tx.getBranchLazy<vector<float>>("tc_types"); });
    ana.histograms.addVecHistogram("tc_uniq_types", 119, 0, 119, [&]() { return ana.tx.getBranchLazy<vector<float>>("tc_uniq_types"); });

    // Book Histograms
    ana.cutflow.bookHistograms(ana.histograms); // if just want to book everywhere

    // Looping input file
    while (ana.looper.nextEvent())
    {

        // If splitting jobs are requested then determine whether to process the event or not based on remainder
        if (ana.job_index != -1 and ana.nsplit_jobs != -1)
        {
            if (ana.looper.getNEventsProcessed() % ana.nsplit_jobs != (unsigned int) ana.job_index)
                continue;
        }

        // Reset all variables
        ana.tx.clear();

        for (unsigned int isimtrk = 0; isimtrk < sdl.sim_pt().size(); ++isimtrk)
        {

            const float& pt = sdl.sim_pt()[isimtrk];
            const float& eta = sdl.sim_eta()[isimtrk];
            const float& dz = sdl.sim_pca_dz()[isimtrk];
            const float& dxy = sdl.sim_pca_dxy()[isimtrk];
            const float& phi = sdl.sim_phi()[isimtrk];
            const int& bunch = sdl.sim_bunchCrossing()[isimtrk];
            const int& pdgid = sdl.sim_pdgId()[isimtrk];

            if (abs(dz) > 30 or abs(dxy) > 2.5)
                continue;
            if (bunch != 0)
                continue;
            if (ana.pdgid != 0 and abs(pdgid) != abs(ana.pdgid))
                continue;

            if (abs(eta) > 2.5)
                continue;
            if (abs(pt) < 1.)
                continue;

            // if (sdl.sim_tcIdx()[isimtrk].size() == 0)
            // {
            //     std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  " isimtrk: " << isimtrk <<  std::endl;
            // }

            bool found = false;
            for (auto& tcIdx : sdl.sim_tcIdx()[isimtrk])
            {
                const std::vector<int>& layers = sdl.tc_layer()[tcIdx];
                if (TC_set4_types_map.find(layers) != TC_set4_types_map.end())
                {
                    found = true;
                }

            }

            if (not found)
            {

                std::vector<int> tc_uniq_types;
                for (auto& tcIdx : sdl.sim_tcIdx()[isimtrk])
                {
                    const std::vector<int>& layers = sdl.tc_layer()[tcIdx];
                    if ((TC_set4_types_map.find(layers) == TC_set4_types_map.end()) and (TC_types_map.find(layers) != TC_types_map.end()))
                    {

                        int type = TC_types_map[layers];
                        if (std::find(tc_uniq_types.begin(), tc_uniq_types.end(), type) == tc_uniq_types.end())
                        {
                            tc_uniq_types.push_back(type);
                            ana.tx.pushbackToBranch<float>("tc_uniq_types", type);
                        }

                    }

                }

            }

        }

        //Do what you need to do in for each event here
        //To save use the following function
        ana.cutflow.fill();
    }

    // Writing output file
    ana.cutflow.saveOutput();

    // The below can be sometimes crucial
    delete ana.output_tfile;
}
