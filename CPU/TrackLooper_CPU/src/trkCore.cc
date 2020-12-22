#include "trkCore.h"

//__________________________________________________________________________________________
float simhit_p(unsigned int simhitidx)
{

    // |momentum| calculation

    float px = trk.simhit_px()[simhitidx];
    float py = trk.simhit_py()[simhitidx];
    float pz = trk.simhit_pz()[simhitidx];
    return sqrt(px*px + py*py + pz*pz);
}

//__________________________________________________________________________________________
float hitAngle(unsigned int simhitidx)
{

    // This is the angle calculation between position vector and the momentum vector

    float x = trk.simhit_x()[simhitidx];
    float y = trk.simhit_y()[simhitidx];
    float z = trk.simhit_z()[simhitidx];
    float r3 = sqrt(x*x + y*y + z*z);
    float px = trk.simhit_px()[simhitidx];
    float py = trk.simhit_py()[simhitidx];
    float pz = trk.simhit_pz()[simhitidx];
    float p = sqrt(px*px + py*py + pz*pz);
    float rdotp = x*px + y*py + z*pz;
    rdotp = rdotp / r3;
    rdotp = rdotp / p;
    float angle = acos(rdotp);
    return angle;
}

//__________________________________________________________________________________________
bool isMuonCurlingHit(unsigned int isimtrk, unsigned int ith_hit)
{

    // To assess whether the ith_hit for isimtrk is a "curling" hit

    // Retrieve the sim hit idx
    unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

    // We're only concerned about hits in the outer tracker
    // This is more of a sanity check
    if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
        return false;

    // Retrieve the sim hit pdgId
    int simhit_particle = trk.simhit_particle()[simhitidx];

    // If the hit is not muon then we can't tell anything
    if (abs(simhit_particle) != 13)
        return false;

    // Get the angle
    float angle = hitAngle(simhitidx);

    // If the angle is quite different then it's the last hit
    if (abs(angle) > 1.6)
        return true;

    // Afterwards, we check the energy loss
    //
    // If it is the first hit without any previous hit present,
    // we can't tell whether it is last hit or not
    // so we just say false to be conservative
    if (ith_hit == 0)
        return false;

    // Retrieve the module where the hit lies
    int detid = trk.simhit_detId()[simhitidx];
    SDL::Module module = SDL::Module(detid);

    // Calculate the momentum
    float p = simhit_p(simhitidx);

    // Find the previous simhit that is on the lower side of the module
    int simhitidx_previous = -999;
    for (unsigned int ith_back = 1; ith_back <= ith_hit; ith_back++)
    {
        // Retrieve the hit idx of ith_hit - ith_back;
        unsigned int simhitidx_previous_candidate = trk.sim_simHitIdx()[isimtrk][ith_hit-ith_back];

        if (not (trk.simhit_subdet()[simhitidx_previous_candidate] == 4 or trk.simhit_subdet()[simhitidx_previous_candidate] == 5))
            continue;

        if (not (trk.simhit_particle()[simhitidx_previous_candidate] == 13))
            continue;

        // Retrieve the module where the previous candidate hit lies
        int detid = trk.simhit_detId()[simhitidx_previous_candidate];
        SDL::Module module = SDL::Module(detid);

        // If the module is lower, then we get the index
        if (module.isLower())
        {
            simhitidx_previous = simhitidx_previous_candidate;
            break;
        }

    }

    // If no previous layer is found then can't do much
    if (simhitidx_previous == -999)
        return false;

    // Get the previous layer momentum
    float p_previous = simhit_p(simhitidx_previous);

    // Calculate the momentum loss
    float loss = fabs(p_previous - p) / p_previous;

    // If the momentum loss is large it is the last hit
    if (loss > 0.35)
        return true;

    // If it reaches this point again, we're not sure what this hit is
    // So we return false
    return false;

}


//__________________________________________________________________________________________
// Currently the denominator tracks for computing algorithmic efficiency are defined as sim track
// with a sequence of 12 reco hits that is expected to be able to create track candidate with
// current algorithm This definition will evolve as the project develops
std::tuple<bool, bool, int, int> getDenomSimTrackTypeForAlgoEff(float isimtrk)
{

    // Read track parameters
    float vx = trk.simvtx_x()[0];
    float vy = trk.simvtx_y()[0];
    float vz = trk.simvtx_z()[0];
    float pt = trk.sim_pt()[isimtrk];
    float eta = trk.sim_eta()[isimtrk];
    float phi = trk.sim_phi()[isimtrk];
    float charge = trk.sim_q()[isimtrk];

    // Construct helix object
    SDLMath::Helix helix(pt, eta, phi, vx, vy, vz, charge);

    // std::cout <<  " pt: " << pt <<  " eta: " << eta <<  " phi: " << phi <<  " vx: " << vx <<  " vy: " << vy <<  " vz: " << vz <<  " charge: " << charge <<  std::endl;

    // List of detids
    std::vector<unsigned int> lower_module_detids;
    std::vector<unsigned int> upper_module_detids;

    // List of pdgids
    std::vector<int> lower_module_pdgids;
    std::vector<int> upper_module_pdgids;

    // List of layers
    std::vector<unsigned int> lower_module_layers;
    std::vector<unsigned int> upper_module_layers;

    // List of n reco hits
    std::vector<int> lower_module_nrecos;
    std::vector<int> upper_module_nrecos;

    // List of n reco hits
    std::vector<float> lower_module_hit_xs;
    std::vector<float> upper_module_hit_xs;
    std::vector<float> lower_module_hit_ys;
    std::vector<float> upper_module_hit_ys;
    std::vector<float> lower_module_hit_zs;
    std::vector<float> upper_module_hit_zs;
    std::vector<float> lower_module_hit_rs;
    std::vector<float> upper_module_hit_rs;

    // Array to store the modules
    std::array<std::vector<SDL::Module>, 6> layers_modules_barrel; // Watch out for duplicates in this vector, do not count with this for unique count.
    std::array<std::vector<SDL::Module>, 6> layers_modules_endcap; // Watch out for duplicates in this vector, do not count with this for unique count.

    // sim hit selection
    // - Must be same pdgid as the parent particle
    // - expected radius based on helix is < 2% away
    for (unsigned int isimhit = 0; isimhit < trk.sim_simHitIdx()[isimtrk].size(); ++isimhit)
    {

        unsigned int isimhitidx = trk.sim_simHitIdx()[isimtrk][isimhit];

        // if (not (trk.simhit_particle()[isimhitidx] == trk.sim_pdgId()[isimtrk]))
        //     continue;

        if (not (trk.simhit_subdet()[isimhitidx] == 4 or trk.simhit_subdet()[isimhitidx] == 5))
            continue;

        // Sim hit vector
        std::vector<float> point = {trk.simhit_x()[isimhitidx], trk.simhit_y()[isimhitidx], trk.simhit_z()[isimhitidx]};

        // Inferring parameter t of helix parametric function via z position
        float t = helix.infer_t(point);

        // If the best fit is more than pi parameter away then it's a returning hit and should be ignored
        if (not (t <= M_PI))
            continue;

        // Expected hit position with given z
        auto [x, y, z, r] = helix.get_helix_point(t);

        // ( expected_r - simhit_r ) / expected_r
        float drfrac = abs(helix.compare_radius(point)) / r;

        // Require that the simhit is within 2% of the expected radius
        if (not (drfrac < 0.02))
            continue;

        // Require that the simhit contains a matched reco
        if (not (trk.simhit_hitIdx()[isimhitidx].size() > 0))
            continue;

        // Sim hit detid
        unsigned int simhit_detid = trk.simhit_detId()[isimhitidx];
        int simhit_pdgid = trk.simhit_particle()[isimhitidx];
        int simhit_layer = 6 * (trk.simhit_subdet()[isimhitidx] == 4) + trk.simhit_layer()[isimhitidx];
        int simhit_nreco = trk.simhit_hitIdx()[isimhitidx].size();
        float simhit_hit_x = trk.simhit_x()[isimhitidx];
        float simhit_hit_y = trk.simhit_y()[isimhitidx];
        float simhit_hit_z = trk.simhit_z()[isimhitidx];
        float simhit_hit_r = sqrt(simhit_hit_x * simhit_hit_x + simhit_hit_y * simhit_hit_y);

        SDL::Module module(simhit_detid);

        if (module.isLower())
        {
            lower_module_detids.push_back(simhit_detid);
            lower_module_pdgids.push_back(simhit_pdgid);
            lower_module_layers.push_back(simhit_layer);
            lower_module_nrecos.push_back(simhit_nreco);
            lower_module_hit_xs.push_back(simhit_hit_x);
            lower_module_hit_ys.push_back(simhit_hit_y);
            lower_module_hit_zs.push_back(simhit_hit_z);
            lower_module_hit_rs.push_back(simhit_hit_r);
        }
        else
        {
            upper_module_detids.push_back(simhit_detid);
            upper_module_pdgids.push_back(simhit_pdgid);
            upper_module_layers.push_back(simhit_layer);
            upper_module_nrecos.push_back(simhit_nreco);
            upper_module_hit_xs.push_back(simhit_hit_x);
            upper_module_hit_ys.push_back(simhit_hit_y);
            upper_module_hit_zs.push_back(simhit_hit_z);
            upper_module_hit_rs.push_back(simhit_hit_r);
        }


        // // list of reco hit matched to this sim hit
        // for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[isimhitidx].size(); ++irecohit)
        // {
        //     // Get the recohit type
        //     int recohittype = trk.simhit_hitType()[isimhitidx][irecohit];

        //     // Consider only ph2 hits (i.e. outer tracker hits)
        //     if (recohittype == 4)
        //     {

        //         int ihit = trk.simhit_hitIdx()[isimhitidx][irecohit];

        //         unsigned int detid = trk.ph2_detId()[ihit];

        //         SDL::Module module(detid);

        //         if (module.isLower())
        //         {
        //             lower_module_detids.push_back(detid);
        //             lower_module_pdgids.push_back(simhit_pdgid);
        //             lower_module_layers.push_back(module.layer());
        //         }
        //         else
        //         {
        //             upper_module_detids.push_back(detid);
        //             upper_module_pdgids.push_back(simhit_pdgid);
        //             upper_module_layers.push_back(module.layer());
        //         }

        //         // if (trk.ph2_subdet()[ihit] == 5)
        //         // {
        //         //     layers_modules_barrel[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
        //         // }
        //         // if (trk.ph2_subdet()[ihit] == 4)
        //         // {
        //         //     layers_modules_endcap[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
        //         // }

        //     }

        // }

    }

    std::vector<unsigned int> matched_lower_module_detids;

    for (auto& lower_module_detid : lower_module_detids)
    {
        SDL::Module lower_module(lower_module_detid);

        if (std::find(upper_module_detids.begin(), upper_module_detids.end(), lower_module.partnerDetId()) != upper_module_detids.end())
        {
            matched_lower_module_detids.push_back(lower_module_detid);
        }
    }

    std::vector<unsigned int> matched_layers;

    for (auto& matched_lower_module_detid : matched_lower_module_detids)
    {
        SDL::Module matched_lower_module(matched_lower_module_detid);
        unsigned int layer = matched_lower_module.layer();
        unsigned int subdet = matched_lower_module.subdet();
        matched_layers.push_back(layer + 6 * (subdet == 4));
    }

    std::vector<int> matched_unique_layers;
    std::vector<int> matched_unique_barrel_layers;
    std::vector<int> matched_unique_endcap_layers;
    for (auto&& i : iter::unique_justseen(matched_layers))
    {
        matched_unique_layers.push_back(i);
        if (i < 7)
            matched_unique_barrel_layers.push_back(i);
        else
            matched_unique_endcap_layers.push_back(i);
    }

    bool has_missing_hits = false;

    if (matched_unique_barrel_layers.size() > 0)
    {
        // std::cout <<  " matched_unique_barrel_layers.back(): " << matched_unique_barrel_layers.back() <<  std::endl;
        // std::cout <<  " matched_unique_barrel_layers.size(): " << matched_unique_barrel_layers.size() <<  std::endl;

        if (matched_unique_barrel_layers.back() != matched_unique_barrel_layers.size())
        {
            has_missing_hits = true;
        }
    }

    if (matched_unique_endcap_layers.size() > 0)
    {

        // std::cout <<  " (matched_unique_endcap_layers.back()-6): " << (matched_unique_endcap_layers.back()-6) <<  std::endl;
        // std::cout <<  " matched_unique_endcap_layers.size(): " << matched_unique_endcap_layers.size() <<  std::endl;

        if ((matched_unique_endcap_layers.back() - 6) != matched_unique_endcap_layers.size())
        {
            has_missing_hits = true;
        }
    }

    int nbarrel = matched_unique_barrel_layers.size();
    int nendcap = matched_unique_endcap_layers.size();

    bool is_clean_punch_through = has_missing_hits ? false : nbarrel + nendcap >= 6;

    bool verbose = has_missing_hits;

    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 8) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;
    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 9) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;
    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(),10) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;
    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(),11) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;

    if (verbose)
    {

        std::cout << "lower_module_detids: ";
        for (auto& lower_module_detid : lower_module_detids) std::cout << std::setw(10) << lower_module_detid << " ";
        std::cout << std::endl;

        std::cout << "lower_module_pdgids: ";
        for (auto& lower_module_pdgid : lower_module_pdgids) std::cout << std::setw(10) << lower_module_pdgid << " ";
        std::cout << std::endl;

        std::cout << "lower_module_layers: ";
        for (auto& lower_module_layer : lower_module_layers) std::cout << std::setw(10) << lower_module_layer << " ";
        std::cout << std::endl;

        std::cout << "lower_module_nrecos: ";
        for (auto& lower_module_nreco : lower_module_nrecos) std::cout << std::setw(10) << lower_module_nreco << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_xs: ";
        for (auto& lower_module_hit_x : lower_module_hit_xs) std::cout << std::setw(10) << lower_module_hit_x << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_ys: ";
        for (auto& lower_module_hit_y : lower_module_hit_ys) std::cout << std::setw(10) << lower_module_hit_y << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_zs: ";
        for (auto& lower_module_hit_z : lower_module_hit_zs) std::cout << std::setw(10) << lower_module_hit_z << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_rs: ";
        for (auto& lower_module_hit_r : lower_module_hit_rs) std::cout << std::setw(10) << lower_module_hit_r << " ";
        std::cout << std::endl;

        std::cout << "upper_module_detids: ";
        for (auto& upper_module_detid : upper_module_detids) std::cout << std::setw(10) << upper_module_detid << " ";
        std::cout << std::endl;

        std::cout << "upper_module_pdgids: ";
        for (auto& upper_module_pdgid : upper_module_pdgids) std::cout << std::setw(10) << upper_module_pdgid << " ";
        std::cout << std::endl;

        std::cout << "upper_module_layers: ";
        for (auto& upper_module_layer : upper_module_layers) std::cout << std::setw(10) << upper_module_layer << " ";
        std::cout << std::endl;

        std::cout << "upper_module_nrecos: ";
        for (auto& upper_module_nreco : upper_module_nrecos) std::cout << std::setw(10) << upper_module_nreco << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_xs: ";
        for (auto& upper_module_hit_x : upper_module_hit_xs) std::cout << std::setw(10) << upper_module_hit_x << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_ys: ";
        for (auto& upper_module_hit_y : upper_module_hit_ys) std::cout << std::setw(10) << upper_module_hit_y << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_zs: ";
        for (auto& upper_module_hit_z : upper_module_hit_zs) std::cout << std::setw(10) << upper_module_hit_z << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_rs: ";
        for (auto& upper_module_hit_r : upper_module_hit_rs) std::cout << std::setw(10) << upper_module_hit_r << " ";
        std::cout << std::endl;

        std::cout << "matched_lower_module_detids: ";
        for (auto& matched_lower_module_detid : matched_lower_module_detids) std::cout << matched_lower_module_detid << " ";
        std::cout << std::endl;

        std::cout << "matched_layers: ";
        for (auto&& i : matched_layers) std::cout << i << " ";
        std::cout << std::endl;

        std::cout << "matched_unique_layers: ";
        for (auto&& i : matched_unique_layers) std::cout << i << " ";
        std::cout << std::endl;

    }

    return std::make_tuple(has_missing_hits, is_clean_punch_through, nbarrel, nendcap);

}

//__________________________________________________________________________________________
bool hasAll12HitsWithNBarrelUsingModuleMap(unsigned int isimtrk, int nbarrel, bool usesimhits)
{

    // Select only tracks that left all 12 hits in the barrel
    std::array<std::vector<SDL::Module>, 6> layers_modules_barrel; // Watch out for duplicates in this vector, do not count with this for unique count.
    std::array<std::vector<SDL::Module>, 6> layers_modules_endcap; // Watch out for duplicates in this vector, do not count with this for unique count.

    std::vector<float> ps;

    for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
    {

        // Retrieve the sim hit idx
        unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

        // Select only the hits in the outer tracker
        if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
            continue;

        if (not (trk.simhit_particle()[simhitidx] == trk.sim_pdgId()[isimtrk]))
            continue;

        if (isMuonCurlingHit(isimtrk, ith_hit))
            break;

        if (not usesimhits)
        {

            // list of reco hit matched to this sim hit
            for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
            {
                // Get the recohit type
                int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

                // Consider only ph2 hits (i.e. outer tracker hits)
                if (recohittype == 4)
                {

                    int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];

                    if (trk.ph2_subdet()[ihit] == 5)
                    {
                        layers_modules_barrel[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
                    }
                    if (trk.ph2_subdet()[ihit] == 4)
                    {
                        layers_modules_endcap[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
                    }

                }

            }

        }
        else
        {

            if (trk.simhit_hitIdx()[simhitidx].size() > 0)
            {

                if (trk.simhit_subdet()[simhitidx] == 5)
                {
                    layers_modules_barrel[trk.simhit_layer()[simhitidx] - 1].push_back(SDL::Module(trk.simhit_detId()[simhitidx]));
                }
                if (trk.simhit_subdet()[simhitidx] == 4)
                {
                    layers_modules_endcap[trk.simhit_layer()[simhitidx] - 1].push_back(SDL::Module(trk.simhit_detId()[simhitidx]));
                }
            }

            // // list of reco hit matched to this sim hit
            // for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
            // {
            //     // Get the recohit type
            //     int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

            //     // Consider only ph2 hits (i.e. outer tracker hits)
            //     if (recohittype == 4)
            //     {

            //         int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];

            //         if (trk.ph2_subdet()[ihit] == 5)
            //         {
            //             layers_modules_barrel[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
            //         }
            //         if (trk.ph2_subdet()[ihit] == 4)
            //         {
            //             layers_modules_endcap[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
            //         }

            //     }

            // }

        }

    }

    // Aggregating good pair modules (i.e. a double module with each side having non-zero reco hits.)
    std::array<std::vector<unsigned int>, 6> layers_good_paired_modules; // Watch out for duplicates in this vector, do not count with this for unique count.

    for (int logical_ilayer = 0; logical_ilayer < 6; ++logical_ilayer)
    {
        // Raw layer number
        bool is_ilayer_barrel = logical_ilayer < nbarrel;
        int raw_ilayer = is_ilayer_barrel ? logical_ilayer: logical_ilayer - nbarrel;
        const std::array<std::vector<SDL::Module>, 6>& layers_modules = is_ilayer_barrel ? layers_modules_barrel : layers_modules_endcap;

        // Then get the module in the ilayer and check that it has a good module pair
        // Loop over modules in the given raw_ilayer and match the pairs and if a good pair of modules have hits in each module
        // then save the lower modules to layers_good_paired_modules.
        // NOTE there may be duplicates being pushed to layers_good_paired_modules
        // Do not count with these vectors
        for (unsigned imod = 0; imod < layers_modules[raw_ilayer].size(); ++imod)
        {
            for (unsigned jmod = imod + 1; jmod < layers_modules[raw_ilayer].size(); ++jmod)
            {
                // if two modules are a good pair
                if (layers_modules[raw_ilayer][imod].partnerDetId() == layers_modules[raw_ilayer][jmod].detId())
                {
                    // add the lower module one to the good_paired_modules list
                    if (layers_modules[raw_ilayer][imod].isLower())
                    {
                        if (std::find(
                                    layers_good_paired_modules[logical_ilayer].begin(),
                                    layers_good_paired_modules[logical_ilayer].end(),
                                    layers_modules[raw_ilayer][imod].detId()) == layers_good_paired_modules[logical_ilayer].end())
                            layers_good_paired_modules[logical_ilayer].push_back(layers_modules[raw_ilayer][imod].detId());
                    }
                    else
                    {
                        if (std::find(
                                    layers_good_paired_modules[logical_ilayer].begin(),
                                    layers_good_paired_modules[logical_ilayer].end(),
                                    layers_modules[raw_ilayer][imod].partnerDetId()) == layers_good_paired_modules[logical_ilayer].end())
                        layers_good_paired_modules[logical_ilayer].push_back(layers_modules[raw_ilayer][imod].partnerDetId());
                    }
                }
            }
        }
    }

    return checkModuleConnectionsAreGood(layers_good_paired_modules);

}

//__________________________________________________________________________________________
bool checkModuleConnectionsAreGood(std::array<std::vector<unsigned int>, 6>& layers_good_paired_modules)
{
    // Dumbest possible solution
    for (auto& module0 : layers_good_paired_modules[0])
    {
        const std::vector<unsigned int>& connectedModule1s = ana.moduleConnectiongMapLoose.getConnectedModuleDetIds(module0);
        for (auto& module1 : layers_good_paired_modules[1])
        {
            if (std::find(connectedModule1s.begin(), connectedModule1s.end(), module1) == connectedModule1s.end())
                break;
            const std::vector<unsigned int>& connectedModule2s = ana.moduleConnectiongMapLoose.getConnectedModuleDetIds(module1);
            for (auto& module2 : layers_good_paired_modules[2])
            {
                if (std::find(connectedModule2s.begin(), connectedModule2s.end(), module2) == connectedModule2s.end())
                    break;
                const std::vector<unsigned int>& connectedModule3s = ana.moduleConnectiongMapLoose.getConnectedModuleDetIds(module2);
                for (auto& module3 : layers_good_paired_modules[3])
                {
                    if (std::find(connectedModule3s.begin(), connectedModule3s.end(), module3) == connectedModule3s.end())
                        break;
                    const std::vector<unsigned int>& connectedModule4s = ana.moduleConnectiongMapLoose.getConnectedModuleDetIds(module3);
                    for (auto& module4 : layers_good_paired_modules[4])
                    {
                        if (std::find(connectedModule4s.begin(), connectedModule4s.end(), module4) == connectedModule4s.end())
                            break;
                        const std::vector<unsigned int>& connectedModule5s = ana.moduleConnectiongMapLoose.getConnectedModuleDetIds(module4);
                        for (auto& module5 : layers_good_paired_modules[5])
                        {
                            if (std::find(connectedModule5s.begin(), connectedModule5s.end(), module5) == connectedModule5s.end())
                                break;
                            return true;
                        }
                    }
                }
            }
        }
    }
    return false;
}

//__________________________________________________________________________________________
bool hasAll12HitsWithNBarrel(unsigned int isimtrk, int nbarrel)
{

    // Select only tracks that left all 12 hits in the barrel
    std::array<std::vector<SDL::Module>, 6> layers_modules_barrel;
    std::array<std::vector<SDL::Module>, 6> layers_modules_endcap;

    std::vector<float> ps;

    for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
    {

        // Retrieve the sim hit idx
        unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

        // Select only the hits in the outer tracker
        if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
            continue;

        // if (not (trk.simhit_particle()[simhitidx] == trk.sim_pdgId()[isimtrk]))
        //     continue;

        if (isMuonCurlingHit(isimtrk, ith_hit))
            break;

        // list of reco hit matched to this sim hit
        for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
        {
            // Get the recohit type
            int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

            // Consider only ph2 hits (i.e. outer tracker hits)
            if (recohittype == 4)
            {

                int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];

                if (trk.ph2_subdet()[ihit] == 5)
                {
                    layers_modules_barrel[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
                }
                if (trk.ph2_subdet()[ihit] == 4)
                {
                    layers_modules_endcap[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
                }

            }

        }

    }

    std::array<bool, 6> has_good_pair_by_layer_barrel;
    has_good_pair_by_layer_barrel[0] = false;
    has_good_pair_by_layer_barrel[1] = false;
    has_good_pair_by_layer_barrel[2] = false;
    has_good_pair_by_layer_barrel[3] = false;
    has_good_pair_by_layer_barrel[4] = false;
    has_good_pair_by_layer_barrel[5] = false;

    std::array<bool, 6> has_good_pair_by_layer_endcap;
    has_good_pair_by_layer_endcap[0] = false;
    has_good_pair_by_layer_endcap[1] = false;
    has_good_pair_by_layer_endcap[2] = false;
    has_good_pair_by_layer_endcap[3] = false;
    has_good_pair_by_layer_endcap[4] = false;
    has_good_pair_by_layer_endcap[5] = false;

    bool has_good_pair_all_layer = true;

    for (int ilayer = 0; ilayer < 6; ++ilayer)
    {

        bool has_good_pair = false;

        if (ilayer < nbarrel)
        {

            for (unsigned imod = 0; imod < layers_modules_barrel[ilayer].size(); ++imod)
            {
                for (unsigned jmod = imod + 1; jmod < layers_modules_barrel[ilayer].size(); ++jmod)
                {
                    if (layers_modules_barrel[ilayer][imod].partnerDetId() == layers_modules_barrel[ilayer][jmod].detId())
                        has_good_pair = true;
                }
            }

        }
        else
        {

            int endcap_ilayer = ilayer - nbarrel;

            for (unsigned imod = 0; imod < layers_modules_endcap[endcap_ilayer].size(); ++imod)
            {
                for (unsigned jmod = imod + 1; jmod < layers_modules_endcap[endcap_ilayer].size(); ++jmod)
                {
                    if (layers_modules_endcap[endcap_ilayer][imod].partnerDetId() == layers_modules_endcap[endcap_ilayer][jmod].detId())
                        has_good_pair = true;
                }
            }

        }

        if (not has_good_pair)
            has_good_pair_all_layer = false;

        if (ilayer < nbarrel)
            has_good_pair_by_layer_barrel[ilayer] = has_good_pair;
        else
            has_good_pair_by_layer_endcap[ilayer] = has_good_pair;

    }


    float pt = trk.sim_pt()[isimtrk];
    float eta = trk.sim_eta()[isimtrk];

    // if (abs((trk.sim_pt()[isimtrk] - 0.71710)) < 0.00001)
    // {
    //     std::cout << std::endl;
    //     std::cout <<  " has_good_pair_by_layer[0]: " << has_good_pair_by_layer[0] <<  " has_good_pair_by_layer[1]: " << has_good_pair_by_layer[1] <<  " has_good_pair_by_layer[2]: " << has_good_pair_by_layer[2] <<  " has_good_pair_by_layer[3]: " << has_good_pair_by_layer[3] <<  " has_good_pair_by_layer[4]: " << has_good_pair_by_layer[4] <<  " has_good_pair_by_layer[5]: " << has_good_pair_by_layer[5] <<  " pt: " << pt <<  " eta: " << eta <<  std::endl;
    // }

    return has_good_pair_all_layer;

}

//__________________________________________________________________________________________
// Check for at least one sim hit in each layer and nothing more
bool goodBarrelTrack(unsigned int isimtrk, int pdgid)
{

    std::vector<int> layers;

    for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
    {

        // Retrieve the sim hit idx
        unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

        // Select only the hits in barrel
        if (not (trk.simhit_subdet()[simhitidx] == 5))
            continue;

        // Select only sim hits matching the particle pdgid
        if (not (trk.simhit_particle()[simhitidx] == trk.sim_pdgId()[isimtrk]))
            continue;

        // if pdgid is provided then check that the pdgid
        if (pdgid != 0)
            if (not (trk.sim_pdgId()[isimtrk] == abs(pdgid)))
                continue;

        // add to layers
        layers.push_back(trk.simhit_layer()[simhitidx]);

    }

    if (not (std::find(layers.begin(), layers.end(), 1) != layers.end())) return false;
    if (not (std::find(layers.begin(), layers.end(), 2) != layers.end())) return false;
    if (not (std::find(layers.begin(), layers.end(), 3) != layers.end())) return false;
    if (not (std::find(layers.begin(), layers.end(), 4) != layers.end())) return false;
    if (not (std::find(layers.begin(), layers.end(), 5) != layers.end())) return false;
    if (not (std::find(layers.begin(), layers.end(), 6) != layers.end())) return false;
    return true;

}

//__________________________________________________________________________________________
bool hasAll12HitsInBarrel(unsigned int isimtrk)
{

    // Select only tracks that left all 12 hits in the barrel
    std::array<std::vector<SDL::Module>, 6> layers_modules;

    std::vector<float> ps;

    for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
    {

        // Retrieve the sim hit idx
        unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

        // Select only the hits in the outer tracker
        if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
            continue;

        // if (not (trk.simhit_particle()[simhitidx] == trk.sim_pdgId()[isimtrk]))
        //     continue;

        if (isMuonCurlingHit(isimtrk, ith_hit))
            break;

        // list of reco hit matched to this sim hit
        for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
        {
            // Get the recohit type
            int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

            // Consider only ph2 hits (i.e. outer tracker hits)
            if (recohittype == 4)
            {

                int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];

                if (trk.ph2_subdet()[ihit] == 5)
                {
                    layers_modules[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
                }

            }

        }

    }

    std::array<bool, 6> has_good_pair_by_layer;
    has_good_pair_by_layer[0] = false;
    has_good_pair_by_layer[1] = false;
    has_good_pair_by_layer[2] = false;
    has_good_pair_by_layer[3] = false;
    has_good_pair_by_layer[4] = false;
    has_good_pair_by_layer[5] = false;

    bool has_good_pair_all_layer = true;

    for (int ilayer = 0; ilayer < 6; ++ilayer)
    {

        bool has_good_pair = false;

        for (unsigned imod = 0; imod < layers_modules[ilayer].size(); ++imod)
        {
            for (unsigned jmod = imod + 1; jmod < layers_modules[ilayer].size(); ++jmod)
            {
                if (layers_modules[ilayer][imod].partnerDetId() == layers_modules[ilayer][jmod].detId())
                    has_good_pair = true;
            }
        }

        if (not has_good_pair)
            has_good_pair_all_layer = false;

        has_good_pair_by_layer[ilayer] = has_good_pair;

    }

    float pt = trk.sim_pt()[isimtrk];
    float eta = trk.sim_eta()[isimtrk];

    // if (abs((trk.sim_pt()[isimtrk] - 0.71710)) < 0.00001)
    // {
    //     std::cout << std::endl;
    //     std::cout <<  " has_good_pair_by_layer[0]: " << has_good_pair_by_layer[0] <<  " has_good_pair_by_layer[1]: " << has_good_pair_by_layer[1] <<  " has_good_pair_by_layer[2]: " << has_good_pair_by_layer[2] <<  " has_good_pair_by_layer[3]: " << has_good_pair_by_layer[3] <<  " has_good_pair_by_layer[4]: " << has_good_pair_by_layer[4] <<  " has_good_pair_by_layer[5]: " << has_good_pair_by_layer[5] <<  " pt: " << pt <<  " eta: " << eta <<  std::endl;
    // }

    return has_good_pair_all_layer;

}



bool hasAtLeastOneHitPairinEndcapLikeTiltedModule(unsigned short layer, unsigned int isimtrk)
{
    //Checking done only for the layer specified by "layer", otherwise function returns true always
    std::vector<SDL::Module> layer_modules;


    for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
    {

        // Retrieve the sim hit idx
        unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

        // Select only the hits in the outer tracker
        if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
            continue;

        if (isMuonCurlingHit(isimtrk, ith_hit))
            break;

        // list of reco hit matched to this sim hit
        for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
        {
            // Get the recohit type
            int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

            // Consider only ph2 hits (i.e. outer tracker hits)
            if (recohittype == 4)
            {

                int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];


                if(trk.ph2_layer()[ihit] == layer and trk.ph2_subdet()[ihit] == 5 and (not SDL::MiniDoublet::useBarrelLogic(trk.ph2_detId()[ihit])))
                {
                    layer_modules.push_back(SDL::Module(trk.ph2_detId()[ihit]));

                }

            }

        }
    }
    //Check if there is at least one pair in an endcap like tilted module
    //in the layer denoted by the parameter "layer"

    bool hasTiltedPairAtLayer = false;

    for (unsigned imod = 0; imod < layer_modules.size(); ++imod)
    {

        for (unsigned jmod = imod + 1; jmod < layer_modules.size(); ++jmod)
        {
            if (layer_modules[imod].partnerDetId() == layer_modules[jmod].detId())
            {
                hasTiltedPairAtLayer = true;
            }

        }
    }

    return hasTiltedPairAtLayer;
}


//__________________________________________________________________________________________
bool isMTVMatch(unsigned int isimtrk, std::vector<unsigned int> hit_idxs, bool verbose)
{
    std::vector<unsigned int> sim_trk_ihits;
    for (auto& i_simhit_idx : trk.sim_simHitIdx()[isimtrk])
    {
        for (auto& ihit : trk.simhit_hitIdx()[i_simhit_idx])
        {
            sim_trk_ihits.push_back(ihit);
        }
    }

    std::sort(sim_trk_ihits.begin(), sim_trk_ihits.end());
    std::sort(hit_idxs.begin(), hit_idxs.end());

    std::vector<unsigned int> v_intersection;

    std::set_intersection(sim_trk_ihits.begin(), sim_trk_ihits.end(),
                          hit_idxs.begin(), hit_idxs.end(),
                          std::back_inserter(v_intersection));

    if (verbose)
    {
        if (v_intersection.size() > ana.nmatch_threshold)
        {
            std::cout << "Matched" << std::endl;
        }
        else
        {
            std::cout << "Not matched" << std::endl;
        }
        std::cout << "sim_trk_ihits: ";
        for (auto& i_simhit_idx : sim_trk_ihits)
            std::cout << i_simhit_idx << " ";
        std::cout << std::endl;

        std::cout << "     hit_idxs: ";
        for (auto& i_hit_idx : hit_idxs)
            std::cout << i_hit_idx << " ";
        std::cout << std::endl;
    }

    int nhits = hit_idxs.size();

    float factor = nhits / 12.;

    // If 75% of 12 hits have been found than it is matched
    return (v_intersection.size() > ana.nmatch_threshold * factor);
}

//__________________________________________________________________________________________
bool is75percentFromSimMatchedHits(std::vector<unsigned int> hitidxs, int pdgid)
{
    std::vector<unsigned int> hitidxs_w_matched_pdgid;
    for (auto& i_reco_hit : hitidxs)
    {
        bool matched = false;
        for (auto& i_simhit_idx : trk.ph2_simHitIdx()[i_reco_hit])
        {
            matched = true;
            break;
            // int particle_id = trk.simhit_particle()[i_simhit_idx];
            // if (abs(pdgid) == abs(particle_id))
            // {
            //     matched = true;
            //     break;
            // }
        }
        if (matched)
            hitidxs_w_matched_pdgid.push_back(i_reco_hit);
    }

    std::vector<unsigned int> v_intersection;

    std::set_intersection(hitidxs_w_matched_pdgid.begin(), hitidxs_w_matched_pdgid.end(),
                          hitidxs.begin(), hitidxs.end(),
                          std::back_inserter(v_intersection));

    // If 75% of 12 hits have been found than it is matched
    return (v_intersection.size() >= ana.nmatch_threshold);

}

//__________________________________________________________________________________________
TVector3 linePropagateR(const TVector3& r3, const TVector3& p3, double rDest, int& status, bool useClosest, bool verbose)
{
  double rt = r3.Pt();
  double d = rDest - rt;

  double dotPR2D = r3.x()*p3.x() + r3.y()*p3.y();

  double pt = p3.Pt();
  double p =  p3.Mag();

  // r3 + p3/p*x*|d| = dest : dest.t = rt + d <=> rt^2 + 2*dotPR2D/p*x*|d| + pt^2/p^2*x^2*d^2 = rt^2 + 2*rt*d + d^2
  // 2*dotPR2D/p*|d|* x + pt^2/p^2*d^2* x^2 - ( 2*rt*d + d^2) = 0
  // 2*dotPR2D/p/|d|* x + pt^2/p^2* x^2 - ( 2*rt/d + 1) = 0
  // x^2 + 2*dotPR2D/p/|d|*(p/pt)^2* x  - ( 2*rt/d + 1)*(p/pt)^2 = 0
  // - dotPR2D/p/|d|*(p/pt)^2  +/- sqrt( (dotPR2D/p/|d|*(p/pt)^2)^2 + ( 2*rt/d + 1)*(p/pt)^2 )
  // (p/pt)*( - dotPR2D/pt/|d|  +/- sqrt( (dotPR2D/pt/|d| )^2 + ( 2*rt/d + 1) ) )
  double bb = dotPR2D/pt/std::abs(d);
  double disc = bb*bb + (2.*rt/d + 1.);
  status = 0;
  if (disc < 0){
    status = 1;
    return r3;
  }
  double dSign = useClosest ? 1. : -1.;
  double xxP = (p/pt)*( sqrt(bb*bb + (2.*rt/d + 1.)) - bb);
  double xxM = (p/pt)*( - sqrt(bb*bb + (2.*rt/d + 1.)) - bb);
  double xx;
  if (useClosest){
    xx = std::abs(xxP) < std::abs(xxM) ? xxP : xxM;
  } else {
    xx = std::abs(xxP) < std::abs(xxM) ? xxM : xxP;
  }
  TVector3 dest = r3 + p3*(std::abs(d)/p)*xx;
  if (verbose || std::abs(dest.Pt() - rDest)>0.001){
    std::cout<<"linePropagateR "<<r3.Pt()<<" "<<r3.Phi()<<" "<<r3.z()<<" "<<pt<<" "<<p
	     <<" "<<d<<" "<<r3.x()*p3.x()<<" "<<r3.y()*p3.y()<<" "<<dotPR2D<<" "<<bb<<" "<<(2.*rt/d + 1.)<<" "<<bb*bb + (2.*rt/d + 1.)
	     <<" => "<<rDest
	     <<" => "<<dest.Pt()<<" "<<dest.Phi()<<" "<<dest.z()
	     <<std::endl;
  }
  return dest;

}

std::pair<TVector3,TVector3> helixPropagateApproxR(const TVector3& r3, const TVector3& p3, double rDest, int q, int& status, bool useClosest, bool verbose)
{
  double epsilon = 0.001;
  double p = p3.Mag();
  double kap = (2.99792458e-3*3.8*q/p);

  auto lastR3 = r3;
  auto lastT3 = p3.Unit();
  int nIts = 7;

  while (std::abs(lastR3.Perp() - rDest) > epsilon && nIts >= 0){
    auto lineEst = linePropagateR(lastR3, lastT3*p, rDest, status, useClosest, verbose);
    if (status){
      if (verbose) std::cout<<" failed with status "<<status<<std::endl;
      return { lineEst, lastT3*p};
    }
    if (q==0) return {lineEst, lastT3*p};

    double dir = (lineEst.x() - lastR3.x())*lastT3.x() + (lineEst.y() - lastR3.y())*lastT3.y() > 0 ? 1. : -1;
    double dS = (lineEst - lastR3).Mag()*dir;
    double phi = kap*dS;
    if (std::abs(phi) > 1) {
      if (verbose) std::cout<<" return line for very large angle "<<status<<std::endl;
      return { lineEst, lastT3*p};
    }
    double alpha = 1 - sin(phi)/phi;
    double beta = (1 - cos(phi))/phi;

    TVector3 tau = lastT3;

    TVector3 hEstR3(tau.x()*(1.-alpha) + tau.y()*beta, tau.y()*(1.-alpha) - tau.x()*beta, tau.z());
    hEstR3 *= dS;
    hEstR3 += lastR3;
    lastR3 = hEstR3;

    TVector3 hEstT3(tau.x()*cos(phi) + tau.y()*sin(phi), tau.y()*cos(phi) - tau.x()*sin(phi), tau.z());
    lastT3 = hEstT3;
    --nIts;
    if (verbose){
      std::cout<<"nIts "<<nIts<<" rDest "<<rDest<<" dS "<<dS<<" phi "<<phi
	       <<" r3In ("<<r3.Pt()<<", "<<r3.Eta()<<", "<<r3.Phi()<<")"
	       <<" p3In ("<<p3.Pt()<<", "<<p3.Eta()<<", "<<p3.Phi()<<")"
	       <<" r3out ("<<lastR3.Pt()<<", "<<lastR3.Eta()<<", "<<lastR3.Phi()<<")"
	       <<" p3Out ("<<lastT3.Pt()*p<<", "<<lastT3.Eta()<<", "<<lastT3.Phi()<<")"
	       <<std::endl;
    }
  }
  status = (std::abs(lastR3.Perp() - rDest) > epsilon);
  return {lastR3, lastT3*p};

}


void fitCircle(std::vector<float> xs, std::vector<float> ys)
{

    TCanvas *c1 = new TCanvas("c1","c1",600,600);
    c1->SetGrid();

    // Generate graph that contains the data
    TGraph* gr = new TGraph(xs.size());
    for (unsigned int i = 0; i < xs.size(); ++i)
    {
        gr->SetPoint(i,xs[i],ys[i]);
    }

    c1->DrawFrame(-120,-120,120,120);
    gr->Draw("p");

    auto chi2Function = [&](const Double_t *par) {
        // minimisation function computing the sum of squares of residuals
        // looping at the graph points
        Int_t np = gr->GetN();
        Double_t f = 0;
        Double_t *x = gr->GetX();
        Double_t *y = gr->GetY();
        for (Int_t i=0;i<np;i++) {
            Double_t u = x[i] - par[0];
            Double_t v = y[i] - par[1];
            Double_t dr = par[2] - std::sqrt(u*u+v*v);
            f += dr*dr;
        }
        return f;
    };

    // wrap chi2 funciton in a function object for the fit
    // 3 is the number of fit parameters (size of array par)
    ROOT::Math::Functor fcn(chi2Function,3);
    ROOT::Fit::Fitter  fitter;
    double pStart[3] = {0,0,1};
    fitter.SetFCN(fcn, pStart);
    fitter.Config().ParSettings(0).SetName("x0");
    fitter.Config().ParSettings(1).SetName("y0");
    fitter.Config().ParSettings(2).SetName("R");

    // do the fit
    bool ok = fitter.FitFCN();
    if (!ok) {
        Error("line3Dfit","Line3D Fit failed");
    }

    const ROOT::Fit::FitResult & result = fitter.Result();

    result.Print(std::cout);

    //Draw the circle on top of the points
    TArc *arc = new TArc(result.Parameter(0),result.Parameter(1),result.Parameter(2));
    arc->SetLineColor(kRed);
    arc->SetLineWidth(4);
    arc->SetFillColorAlpha(0, 0.35);
    arc->Draw();
    c1->SaveAs("result.pdf");

}

//__________________________________________________________________________________________
void printMiniDoubletConnectionMultiplicitiesBarrel(SDL::Event& event, int layer, int depth, bool goinside)
{

    std::vector<int> multiplicities;
    int total_nmult = 0;
    float avg_mult = 0;

    if (not goinside)
    {
        // ----------------
        multiplicities.clear();
        total_nmult = 0;
        for (auto& miniDoubletPtr : event.getLayer(layer, SDL::Layer::Barrel).getMiniDoubletPtrs())
        {
            int nmult = 0;
            for (auto& seg1 : miniDoubletPtr->getListOfOutwardSegmentPtrs())
            {
                if (depth == 1)
                    nmult++;
                for (auto& seg2 : seg1->outerMiniDoubletPtr()->getListOfOutwardSegmentPtrs())
                {
                    if (depth == 2)
                        nmult++;
                    for (auto& seg3 : seg2->outerMiniDoubletPtr()->getListOfOutwardSegmentPtrs())
                    {
                        if (depth == 3)
                            nmult++;
                        for (auto& seg4 : seg3->outerMiniDoubletPtr()->getListOfOutwardSegmentPtrs())
                        {
                            if (depth == 4)
                                nmult++;
                            for (auto& seg5 : seg4->outerMiniDoubletPtr()->getListOfOutwardSegmentPtrs())
                            {
                                if (depth == 5)
                                    nmult++;
                            }
                        }
                    }
                }
            }
            multiplicities.push_back(nmult);
            total_nmult += nmult;
        }
        avg_mult = ((float) total_nmult) / ((float) multiplicities.size());
        std::cout <<  " layer: " << layer <<  " depth: " << depth <<  " total_nmult: " << total_nmult <<  " avg_mult: " << avg_mult <<  " goinside: " << goinside <<  std::endl;
    }
    else
    {

        // ----------------
        multiplicities.clear();
        total_nmult = 0;
        for (auto& miniDoubletPtr : event.getLayer(layer, SDL::Layer::Barrel).getMiniDoubletPtrs())
        {
            int nmult = 0;
            for (auto& seg1 : miniDoubletPtr->getListOfInwardSegmentPtrs())
            {
                if (depth == 1)
                    nmult++;
                for (auto& seg2 : seg1->innerMiniDoubletPtr()->getListOfInwardSegmentPtrs())
                {
                    if (depth == 2)
                        nmult++;
                    for (auto& seg3 : seg2->innerMiniDoubletPtr()->getListOfInwardSegmentPtrs())
                    {
                        if (depth == 3)
                            nmult++;
                        for (auto& seg4 : seg3->innerMiniDoubletPtr()->getListOfInwardSegmentPtrs())
                        {
                            if (depth == 4)
                                nmult++;
                            for (auto& seg5 : seg4->innerMiniDoubletPtr()->getListOfInwardSegmentPtrs())
                            {
                                if (depth == 5)
                                    nmult++;
                            }
                        }
                    }
                }
            }
            multiplicities.push_back(nmult);
            total_nmult += nmult;
        }
        avg_mult = ((float) total_nmult) / ((float) multiplicities.size());
        std::cout <<  " layer: " << layer <<  " depth: " << depth <<  " total_nmult: " << total_nmult <<  " avg_mult: " << avg_mult <<  " goinside: " << goinside <<  std::endl;

    }

}

//__________________________________________________________________________________________
std::vector<unsigned int> getHitIdxListTrackCandidate(SDL::TrackCandidate* tc)
{
    std::vector<unsigned int> hit_idxs;
    hit_idxs.push_back(tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
    hit_idxs.push_back(tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
    hit_idxs.push_back(tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
    hit_idxs.push_back(tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx());
    hit_idxs.push_back(tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx());
    return hit_idxs;
}

//__________________________________________________________________________________________
std::vector<unsigned int> getModuleListTrackCandidate(SDL::TrackCandidate* tc)
{
    std::vector<unsigned int> module_ids;
    module_ids.push_back(tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
    module_ids.push_back(tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
    module_ids.push_back(tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
    module_ids.push_back(tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
    module_ids.push_back(tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
    module_ids.push_back(tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
    return module_ids;
}

//__________________________________________________________________________________________
int getNBarrelTrackCandidate(SDL::TrackCandidate* tc)
{

    int nbarrel = 0;
    if ( (tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).subdet() == SDL::Module::Barrel) nbarrel++;
    if ( (tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).subdet() == SDL::Module::Barrel) nbarrel++;
    if ( (tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).subdet() == SDL::Module::Barrel) nbarrel++;
    if ( (tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).subdet() == SDL::Module::Barrel) nbarrel++;
    if ( (tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).subdet() == SDL::Module::Barrel) nbarrel++;
    if ( (tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).subdet() == SDL::Module::Barrel) nbarrel++;
    return nbarrel;

}

//__________________________________________________________________________________________
vector<unsigned int> allMatchedSimTrkIdxs(SDL::TrackCandidate* tc)
{

    std::vector<int> hitidxs = {
        tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx()
        };

    std::vector<unsigned int> unique_idxs; // to aggregate which ones to count and test

    for (auto& hitidx : hitidxs)
    {
        for (auto& simhit_idx : trk.ph2_simHitIdx()[hitidx])
        {
            int simtrk_idx = trk.simhit_simTrkIdx()[simhit_idx];
            // simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }
    }

    return unique_idxs;

}

//__________________________________________________________________________________________
vector<int> matchedSimTrkIdxs(SDL::MiniDoublet* md)
{
    std::vector<int> hitidxs = {md->lowerHitPtr()->idx(), md->upperHitPtr()->idx()};

    std::vector<vector<int>> simtrk_idxs;
    std::vector<int> unique_idxs; // to aggregate which ones to count and test

    for (auto& hitidx : hitidxs)
    {
        std::vector<int> simtrk_idxs_per_hit;
        for (auto& simhit_idx : trk.ph2_simHitIdx()[hitidx])
        {
            int simtrk_idx = trk.simhit_simTrkIdx()[simhit_idx];
            simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }
        if (simtrk_idxs_per_hit.size() == 0)
        {
            simtrk_idxs_per_hit.push_back(-1);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), -1) == unique_idxs.end())
                unique_idxs.push_back(-1);
        }
        simtrk_idxs.push_back(simtrk_idxs_per_hit);
    }

    // Compute all permutations
    std::function<void(vector<vector<int>>&, vector<int>, size_t, vector<vector<int>>&)> perm = [&](vector<vector<int>>& result, vector<int> intermediate, size_t n, vector<vector<int>>& va)
    {
        if (va.size() > n)
        {
            for (auto x : va[n])
            {
                intermediate.push_back(x);
                perm(result, intermediate, n+1, va);
            }
        }
        else
        {
            result.push_back(intermediate);
        }
    };

    vector<vector<int>> allperms;
    perm(allperms, vector<int>(), 0, simtrk_idxs);

    // for (auto& iperm : allperms)
    // {
    //     for (auto& idx : iperm)
    //         std::cout << idx << " ";
    //     std::cout << std::endl;
    // }

    std::vector<int> matched_sim_trk_idxs;
    for (auto& trkidx_perm : allperms)
    {
        std::vector<int> counts;
        for (auto& unique_idx : unique_idxs)
        {
            int cnt = std::count(trkidx_perm.begin(), trkidx_perm.end(), unique_idx);
            counts.push_back(cnt);
        }
        auto result = std::max_element(counts.begin(), counts.end());
        int rawidx = std::distance(counts.begin(), result);
        int trkidx = unique_idxs[rawidx];
        if (trkidx < 0)
            continue;
        if (counts[rawidx] > 1.5)
            matched_sim_trk_idxs.push_back(trkidx);
    }

    return matched_sim_trk_idxs;

}

//__________________________________________________________________________________________
vector<int> matchedSimTrkIdxs(SDL::Segment* sg, bool matchOnlyAnchor)
{
    std::vector<int> hitidxs_all = {
        sg->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        sg->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        sg->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        sg->outerMiniDoubletPtr()->upperHitPtr()->idx()
    };

    std::vector<int> hitidxs_onlyAnchor = {
        sg->innerMiniDoubletPtr()->anchorHitPtr()->idx(),
        sg->outerMiniDoubletPtr()->anchorHitPtr()->idx()
    };

    std::vector<int> hitidxs;
    if (matchOnlyAnchor)
    {
        hitidxs = hitidxs_onlyAnchor;
    }
    else
    {
        hitidxs = hitidxs_all;
    }

    int threshold = matchOnlyAnchor ? 1 : 3;

    std::vector<vector<int>> simtrk_idxs;
    std::vector<int> unique_idxs; // to aggregate which ones to count and test

    for (auto& hitidx : hitidxs)
    {
        std::vector<int> simtrk_idxs_per_hit;
        for (auto& simhit_idx : trk.ph2_simHitIdx()[hitidx])
        {
            int simtrk_idx = trk.simhit_simTrkIdx()[simhit_idx];
            simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }
        if (simtrk_idxs_per_hit.size() == 0)
        {
            simtrk_idxs_per_hit.push_back(-1);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), -1) == unique_idxs.end())
                unique_idxs.push_back(-1);
        }
        simtrk_idxs.push_back(simtrk_idxs_per_hit);
    }

    // Compute all permutations
    std::function<void(vector<vector<int>>&, vector<int>, size_t, vector<vector<int>>&)> perm = [&](vector<vector<int>>& result, vector<int> intermediate, size_t n, vector<vector<int>>& va)
    {
        if (va.size() > n)
        {
            for (auto x : va[n])
            {
                intermediate.push_back(x);
                perm(result, intermediate, n+1, va);
            }
        }
        else
        {
            result.push_back(intermediate);
        }
    };

    vector<vector<int>> allperms;
    perm(allperms, vector<int>(), 0, simtrk_idxs);

    // for (auto& iperm : allperms)
    // {
    //     for (auto& idx : iperm)
    //         std::cout << idx << " ";
    //     std::cout << std::endl;
    // }

    std::vector<int> matched_sim_trk_idxs;
    for (auto& trkidx_perm : allperms)
    {
        std::vector<int> counts;
        for (auto& unique_idx : unique_idxs)
        {
            int cnt = std::count(trkidx_perm.begin(), trkidx_perm.end(), unique_idx);
            counts.push_back(cnt);
        }
        auto result = std::max_element(counts.begin(), counts.end());
        int rawidx = std::distance(counts.begin(), result);
        int trkidx = unique_idxs[rawidx];
        if (trkidx < 0)
            continue;
        if (counts[rawidx] > threshold)
            matched_sim_trk_idxs.push_back(trkidx);
    }

    return matched_sim_trk_idxs;

}

//__________________________________________________________________________________________
vector<int> matchedSimTrkIdxs(SDL::TrackCandidate* tc)
{

    std::vector<int> hitidxs = {
        tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->innerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->innerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tc->outerTrackletBasePtr()->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx()
        };

    // Resize the hitidxs so that it is different for different track candidate type
    // For some cases it has 6 minidoublets, some cases it has 5 minidoublets
    std::vector<int>::iterator ip;

    // Using std::unique
    ip = std::unique(hitidxs.begin(), hitidxs.end());
    // Now v becomes {1 3 10 1 3 7 8 * * * * *}
    // * means undefined

    // Resizing the vector so as to remove the undefined terms
    hitidxs.resize(std::distance(hitidxs.begin(), ip));

    std::vector<vector<int>> simtrk_idxs;
    std::vector<int> unique_idxs; // to aggregate which ones to count and test

    for (auto& hitidx : hitidxs)
    {
        std::vector<int> simtrk_idxs_per_hit;
        for (auto& simhit_idx : trk.ph2_simHitIdx()[hitidx])
        {
            int simtrk_idx = trk.simhit_simTrkIdx()[simhit_idx];
            simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }
        if (simtrk_idxs_per_hit.size() == 0)
        {
            simtrk_idxs_per_hit.push_back(-1);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), -1) == unique_idxs.end())
                unique_idxs.push_back(-1);
        }
        simtrk_idxs.push_back(simtrk_idxs_per_hit);
    }

    // // print
    // std::cout << "va print" << std::endl;
    // for (auto& vec : simtrk_idxs)
    // {
    //     for (auto& idx : vec)
    //     {
    //         std::cout << idx << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "va print end" << std::endl;

    // Compute all permutations
    std::function<void(vector<vector<int>>&, vector<int>, size_t, vector<vector<int>>&)> perm = [&](vector<vector<int>>& result, vector<int> intermediate, size_t n, vector<vector<int>>& va)
    {
        if (va.size() > n)
        {
            for (auto x : va[n])
            {
                intermediate.push_back(x);
                perm(result, intermediate, n+1, va);
            }
        }
        else
        {
            result.push_back(intermediate);
        }
    };

    vector<vector<int>> allperms;
    perm(allperms, vector<int>(), 0, simtrk_idxs);

    // for (auto& iperm : allperms)
    // {
    //     for (auto& idx : iperm)
    //         std::cout << idx << " ";
    //     std::cout << std::endl;
    // }

    int nhits = hitidxs.size();

    float factor = nhits / 12.;

    std::vector<int> matched_sim_trk_idxs;
    for (auto& trkidx_perm : allperms)
    {
        std::vector<int> counts;
        for (auto& unique_idx : unique_idxs)
        {
            int cnt = std::count(trkidx_perm.begin(), trkidx_perm.end(), unique_idx);
            counts.push_back(cnt);
        }
        auto result = std::max_element(counts.begin(), counts.end());
        int rawidx = std::distance(counts.begin(), result);
        int trkidx = unique_idxs[rawidx];
        if (trkidx < 0)
            continue;
        // if (counts[rawidx] > ana.nmatch_threshold)
        if (counts[rawidx] > ana.nmatch_threshold)
            matched_sim_trk_idxs.push_back(trkidx);
    }

    return matched_sim_trk_idxs;

}

//__________________________________________________________________________________________
std::vector<float> getPtBounds()
{
    std::vector<float> pt_boundaries;
    if (ana.ptbound_mode == 0)
        pt_boundaries = {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};
    else if (ana.ptbound_mode == 1)
        pt_boundaries = {0.988, 0.99, 0.992, 0.994, 0.996, 0.998, 1.0, 1.002, 1.004, 1.006, 1.008, 1.01, 1.012}; // lowpt
    else if (ana.ptbound_mode == 2)
        pt_boundaries = {0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.00, 1.005, 1.01, 1.015, 1.02, 1.025, 1.03, 1.035, 1.04, 1.045, 1.05}; // pt 0p95 1p05
    else if (ana.ptbound_mode == 3)
        pt_boundaries = {0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.2, 1.5}; // lowpt
    else if (ana.ptbound_mode == 4)
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}; // lowpt
    else if (ana.ptbound_mode == 5)
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.24, 1.28, 1.32, 1.36, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0}; // lowpt
    else if (ana.ptbound_mode == 6)
        pt_boundaries = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 3.0, 4.0, 5.0}; // lowpt
    else if (ana.ptbound_mode == 7)
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50}; // lowpt
    else if (ana.ptbound_mode == 8)
        pt_boundaries = {0, 0.5, 1.0, 3.0, 5.0, 10, 15., 25, 50};
    return pt_boundaries;
}

//__________________________________________________________________________________________
bool goodEvent()
{
    if (ana.specific_event_index >= 0)
    {
        if ((int)ana.looper.getCurrentEventIndex() != ana.specific_event_index)
            return false;
    }

    // If splitting jobs are requested then determine whether to process the event or not based on remainder
    if (ana.nsplit_jobs >= 0 and ana.job_index >= 0)
    {
        if (ana.looper.getNEventsProcessed() % ana.nsplit_jobs != (unsigned int) ana.job_index)
            return false;
    }

    if (ana.verbose) std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;

    return true;
}

//__________________________________________________________________________________________
bool inTimeTrackWithPdgId(int isimtrk, int pdgid)
{
    // Then select all charged particle
    if (pdgid == 0)
    {
        // Select all charged particle tracks
        if (abs(trk.sim_q()[isimtrk]) == 0)
            return false;
    }
    else
    {
        // Select tracks with given pdgid
        if (abs(trk.sim_pdgId()[isimtrk]) != pdgid)
            return false;
    }

    // Select in time only
    if (abs(trk.sim_bunchCrossing()[isimtrk]) != 0)
        return false;

    return true;
}

//__________________________________________________________________________________________
TrackletType getTrackletCategory(SDL::Tracklet& tl)
{
    const SDL::Module& innerSgInnerMDAnchorHitModule = tl.innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& outerSgInnerMDAnchorHitModule = tl.outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& innerSgOuterMDAnchorHitModule = tl.innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& outerSgOuterMDAnchorHitModule = tl.outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule();

    const int innerLayerIdx = innerSgInnerMDAnchorHitModule.layer();
    const int outerLayerIdx = outerSgInnerMDAnchorHitModule.layer();

    const bool l1_Barrel = (innerSgInnerMDAnchorHitModule.subdet() == SDL::Module::Barrel);
    const bool l2_Barrel = (innerSgOuterMDAnchorHitModule.subdet() == SDL::Module::Barrel);
    const bool l3_Barrel = (outerSgInnerMDAnchorHitModule.subdet() == SDL::Module::Barrel);
    const bool l4_Barrel = (outerSgOuterMDAnchorHitModule.subdet() == SDL::Module::Barrel);
    const bool l1_Endcap = (innerSgInnerMDAnchorHitModule.subdet() == SDL::Module::Endcap);
    const bool l2_Endcap = (innerSgOuterMDAnchorHitModule.subdet() == SDL::Module::Endcap);
    const bool l3_Endcap = (outerSgInnerMDAnchorHitModule.subdet() == SDL::Module::Endcap);
    const bool l4_Endcap = (outerSgOuterMDAnchorHitModule.subdet() == SDL::Module::Endcap);

    if (innerLayerIdx == 1 and outerLayerIdx == 3 and l1_Barrel and l2_Barrel and l3_Barrel and l4_Barrel) return BB1BB3;
    if (innerLayerIdx == 2 and outerLayerIdx == 4 and l1_Barrel and l2_Barrel and l3_Barrel and l4_Barrel) return BB2BB4;
    if (innerLayerIdx == 3 and outerLayerIdx == 5 and l1_Barrel and l2_Barrel and l3_Barrel and l4_Barrel) return BB3BB5;
    if (innerLayerIdx == 1 and outerLayerIdx == 3 and l1_Barrel and l2_Barrel and l3_Barrel and l4_Endcap) return BB1BE3;
    if (innerLayerIdx == 2 and outerLayerIdx == 4 and l1_Barrel and l2_Barrel and l3_Barrel and l4_Endcap) return BB2BE4;
    if (innerLayerIdx == 3 and outerLayerIdx == 5 and l1_Barrel and l2_Barrel and l3_Barrel and l4_Endcap) return BB3BE5;
    if (innerLayerIdx == 1 and outerLayerIdx == 1 and l1_Barrel and l2_Barrel and l3_Endcap and l4_Endcap) return BB1EE3;
    if (innerLayerIdx == 2 and outerLayerIdx == 1 and l1_Barrel and l2_Barrel and l3_Endcap and l4_Endcap) return BB2EE4;
    if (innerLayerIdx == 3 and outerLayerIdx == 1 and l1_Barrel and l2_Barrel and l3_Endcap and l4_Endcap) return BB3EE5;
    if (innerLayerIdx == 1 and outerLayerIdx == 2 and l1_Barrel and l2_Endcap and l3_Endcap and l4_Endcap) return BE1EE3;
    if (innerLayerIdx == 2 and outerLayerIdx == 2 and l1_Barrel and l2_Endcap and l3_Endcap and l4_Endcap) return BE2EE4;
    if (innerLayerIdx == 3 and outerLayerIdx == 2 and l1_Barrel and l2_Endcap and l3_Endcap and l4_Endcap) return BE3EE5;
    if (innerLayerIdx == 1 and outerLayerIdx == 3 and l1_Endcap and l2_Endcap and l3_Endcap and l4_Endcap) return EE1EE3;
    if (innerLayerIdx == 2 and outerLayerIdx == 4 and l1_Endcap and l2_Endcap and l3_Endcap and l4_Endcap) return EE2EE4;
}

//__________________________________________________________________________________________
int getNPSModules(SDL::Tracklet& tl)
{
    const SDL::Module& innerSgInnerMDAnchorHitModule = tl.innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& outerSgInnerMDAnchorHitModule = tl.outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& innerSgOuterMDAnchorHitModule = tl.innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& outerSgOuterMDAnchorHitModule = tl.outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule();

    int nPS = 0;

    if (innerSgInnerMDAnchorHitModule.moduleType() == SDL::Module::PS)
        nPS++;
    if (innerSgOuterMDAnchorHitModule.moduleType() == SDL::Module::PS)
        nPS++;
    if (outerSgInnerMDAnchorHitModule.moduleType() == SDL::Module::PS)
        nPS++;
    if (outerSgOuterMDAnchorHitModule.moduleType() == SDL::Module::PS)
        nPS++;

    return nPS;
}

//__________________________________________________________________________________________
std::vector<int> matchedSimTrkIdxs(SDL::Tracklet& tl)
{
    std::vector<int> hitidxs = {
        tl.innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tl.innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tl.innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tl.innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx(),
        tl.outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tl.outerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tl.outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tl.outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx()
        };

    std::vector<vector<int>> simtrk_idxs;
    std::vector<int> unique_idxs; // to aggregate which ones to count and test

    for (auto& hitidx : hitidxs)
    {
        std::vector<int> simtrk_idxs_per_hit;
        for (auto& simhit_idx : trk.ph2_simHitIdx()[hitidx])
        {
            int simtrk_idx = trk.simhit_simTrkIdx()[simhit_idx];
            simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }
        if (simtrk_idxs_per_hit.size() == 0)
        {
            simtrk_idxs_per_hit.push_back(-1);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), -1) == unique_idxs.end())
                unique_idxs.push_back(-1);
        }
        simtrk_idxs.push_back(simtrk_idxs_per_hit);
    }

    // // print
    // std::cout << "va print" << std::endl;
    // for (auto& vec : simtrk_idxs)
    // {
    //     for (auto& idx : vec)
    //     {
    //         std::cout << idx << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "va print end" << std::endl;

    // Compute all permutations
    std::function<void(vector<vector<int>>&, vector<int>, size_t, vector<vector<int>>&)> perm = [&](vector<vector<int>>& result, vector<int> intermediate, size_t n, vector<vector<int>>& va)
    {
        if (va.size() > n)
        {
            for (auto x : va[n])
            {
                intermediate.push_back(x);
                perm(result, intermediate, n+1, va);
            }
        }
        else
        {
            result.push_back(intermediate);
        }
    };

    vector<vector<int>> allperms;
    perm(allperms, vector<int>(), 0, simtrk_idxs);

    // for (auto& iperm : allperms)
    // {
    //     for (auto& idx : iperm)
    //         std::cout << idx << " ";
    //     std::cout << std::endl;
    // }

    std::vector<int> matched_sim_trk_idxs;
    for (auto& trkidx_perm : allperms)
    {
        std::vector<int> counts;
        for (auto& unique_idx : unique_idxs)
        {
            int cnt = std::count(trkidx_perm.begin(), trkidx_perm.end(), unique_idx);
            counts.push_back(cnt);
        }
        auto result = std::max_element(counts.begin(), counts.end());
        int rawidx = std::distance(counts.begin(), result);
        int trkidx = unique_idxs[rawidx];
        if (trkidx < 0)
            continue;
        if (counts[rawidx] > 6)
            matched_sim_trk_idxs.push_back(trkidx);
    }

    return matched_sim_trk_idxs;

}

//__________________________________________________________________________________________
std::vector<int> matchedSimTrkIdxs(SDL::Triplet* tp)
{
    std::vector<int> hitidxs = {
        tp->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tp->innerSegmentPtr()->innerMiniDoubletPtr()->upperHitPtr()->idx(),
        tp->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tp->innerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx(),
        tp->outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->idx(),
        tp->outerSegmentPtr()->outerMiniDoubletPtr()->upperHitPtr()->idx()
        };

    std::vector<vector<int>> simtrk_idxs;
    std::vector<int> unique_idxs; // to aggregate which ones to count and test

    for (auto& hitidx : hitidxs)
    {
        std::vector<int> simtrk_idxs_per_hit;
        for (auto& simhit_idx : trk.ph2_simHitIdx()[hitidx])
        {
            int simtrk_idx = trk.simhit_simTrkIdx()[simhit_idx];
            simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }
        if (simtrk_idxs_per_hit.size() == 0)
        {
            simtrk_idxs_per_hit.push_back(-1);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), -1) == unique_idxs.end())
                unique_idxs.push_back(-1);
        }
        simtrk_idxs.push_back(simtrk_idxs_per_hit);
    }

    // // print
    // std::cout << "va print" << std::endl;
    // for (auto& vec : simtrk_idxs)
    // {
    //     for (auto& idx : vec)
    //     {
    //         std::cout << idx << " ";
    //     }
    //     std::cout << std::endl;
    // }
    // std::cout << "va print end" << std::endl;

    // Compute all permutations
    std::function<void(vector<vector<int>>&, vector<int>, size_t, vector<vector<int>>&)> perm = [&](vector<vector<int>>& result, vector<int> intermediate, size_t n, vector<vector<int>>& va)
    {
        if (va.size() > n)
        {
            for (auto x : va[n])
            {
                intermediate.push_back(x);
                perm(result, intermediate, n+1, va);
            }
        }
        else
        {
            result.push_back(intermediate);
        }
    };

    vector<vector<int>> allperms;
    perm(allperms, vector<int>(), 0, simtrk_idxs);

    // for (auto& iperm : allperms)
    // {
    //     for (auto& idx : iperm)
    //         std::cout << idx << " ";
    //     std::cout << std::endl;
    // }

    std::vector<int> matched_sim_trk_idxs;
    for (auto& trkidx_perm : allperms)
    {
        std::vector<int> counts;
        for (auto& unique_idx : unique_idxs)
        {
            int cnt = std::count(trkidx_perm.begin(), trkidx_perm.end(), unique_idx);
            counts.push_back(cnt);
        }
        auto result = std::max_element(counts.begin(), counts.end());
        int rawidx = std::distance(counts.begin(), result);
        int trkidx = unique_idxs[rawidx];
        if (trkidx < 0)
            continue;
        if (counts[rawidx] > 4.5)
            matched_sim_trk_idxs.push_back(trkidx);
    }

    return matched_sim_trk_idxs;

}

//__________________________________________________________________________________________
std::vector<int> matchedSimTrkIdxs(std::vector<int> hitidxs, std::vector<int> hittypes)
{
    if (hitidxs.size() != hittypes.size())
    {
        std::cout << "Error: matched_sim_trk_idxs()   hitidxs and hittypes have different lengths" << std::endl;
        std::cout << "hitidxs.size(): " << hitidxs.size() << std::endl;
        std::cout << "hittypes.size(): " << hittypes.size() << std::endl;
    }

    std::vector<std::pair<int, int>> to_check_duplicate;
    for (auto&& [ihit, ihitdata] : iter::enumerate(iter::zip(hitidxs, hittypes)))
    {
        auto&& [hitidx, hittype] = ihitdata;
        auto item = std::make_pair(hitidx, hittype);
        if (std::find(to_check_duplicate.begin(), to_check_duplicate.end(), item) == to_check_duplicate.end())
        {
            to_check_duplicate.push_back(item);
        }
    }

    int nhits_input = to_check_duplicate.size();

    std::vector<vector<int>> simtrk_idxs;
    std::vector<int> unique_idxs; // to aggregate which ones to count and test

    for (auto&& [ihit, ihitdata] : iter::enumerate(to_check_duplicate))
    {
        auto&& [hitidx, hittype] = ihitdata;

        std::vector<int> simtrk_idxs_per_hit;

        const std::vector<vector<int>>* simHitIdxs;

        if (hittype == 4)
            simHitIdxs = &trk.ph2_simHitIdx();
        else
            simHitIdxs = &trk.pix_simHitIdx();

        if ( (*simHitIdxs).size() <= hitidx)
        {
                std::cout << (*simHitIdxs).size() << " " << hittype << std::endl;
                std::cout << hitidx << " " << hittype << std::endl;
        }

        for (auto& simhit_idx : (*simHitIdxs).at(hitidx))
        {
            // std::cout << "  " << trk.simhit_simTrkIdx().size() << std::endl;
            // std::cout << " " << simhit_idx << std::endl;
            if (trk.simhit_simTrkIdx().size() <= simhit_idx)
            {
                std::cout << (*simHitIdxs).size() << " " << hittype << std::endl;
                std::cout << hitidx << " " << hittype << std::endl;
                std::cout << trk.simhit_simTrkIdx().size() << " " << simhit_idx << std::endl;
            }
            int simtrk_idx = trk.simhit_simTrkIdx().at(simhit_idx);
            simtrk_idxs_per_hit.push_back(simtrk_idx);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), simtrk_idx) == unique_idxs.end())
                unique_idxs.push_back(simtrk_idx);
        }

        if (simtrk_idxs_per_hit.size() == 0)
        {
            simtrk_idxs_per_hit.push_back(-1);
            if (std::find(unique_idxs.begin(), unique_idxs.end(), -1) == unique_idxs.end())
                unique_idxs.push_back(-1);
        }

        simtrk_idxs.push_back(simtrk_idxs_per_hit);
    }

    // // print
    // if (ana.verbose != 0)
    // {
    //     std::cout << "va print" << std::endl;
    //     for (auto& vec : simtrk_idxs)
    //     {
    //         for (auto& idx : vec)
    //         {
    //             std::cout << idx << " ";
    //         }
    //         std::cout << std::endl;
    //     }
    //     std::cout << "va print end" << std::endl;
    // }

    // Compute all permutations
    std::function<void(vector<vector<int>>&, vector<int>, size_t, vector<vector<int>>&)> perm =
        [&](vector<vector<int>>& result, vector<int> intermediate, size_t n, vector<vector<int>>& va)
    {
        if (va.size() > n)
        {
            for (auto x : va[n])
            {
                intermediate.push_back(x);
                perm(result, intermediate, n+1, va);
            }
        }
        else
        {
            result.push_back(intermediate);
        }
    };

    vector<vector<int>> allperms;
    perm(allperms, vector<int>(), 0, simtrk_idxs);

    std::vector<int> matched_sim_trk_idxs;
    for (auto& trkidx_perm : allperms)
    {
        std::vector<int> counts;
        for (auto& unique_idx : unique_idxs)
        {
            int cnt = std::count(trkidx_perm.begin(), trkidx_perm.end(), unique_idx);
            counts.push_back(cnt);
        }
        auto result = std::max_element(counts.begin(), counts.end());
        int rawidx = std::distance(counts.begin(), result);
        int trkidx = unique_idxs[rawidx];
        if (trkidx < 0)
            continue;
        if (counts[rawidx] > (((float)nhits_input) * 0.75))
            matched_sim_trk_idxs.push_back(trkidx);
    }

    return matched_sim_trk_idxs;
}

//__________________________________________________________________________________________
void loadMaps()
{
    SDL::endcapGeometry.load("scripts/endcap_orientation_data_v2.txt"); // centroid values added to the map
    SDL::tiltedGeometry.load("scripts/tilted_orientation_data.txt");

    SDL::moduleConnectionMap.load("data/module_connection_combined_2020_0520_helixray.txt");
    ana.moduleConnectiongMapLoose.load("data/module_connection_combined_2020_0520_helixray.txt");

    // SDL::moduleConnectionMap.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/module_connection_map_data_10_e0_200_100_pt0p8_2p0_400_pt0p8_2p0_nolossers_dxy35cm_endcaplayer2.txt");
    // ana.moduleConnectiongMapLoose.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/module_connection_map_data_10_e0_200_100_pt0p8_2p0_400_pt0p8_2p0_nolossers_dxy35cm_endcaplayer2.txt");

    // SDL::moduleConnectionMap.load("data/module_connection_2020_0429.txt");
    // ana.moduleConnectiongMapLoose.load("data/module_connection_2020_0429.txt");

    // SDL::moduleConnectionMap.load("data/module_connection_tracing_2020_0514.txt");
    // SDL::moduleConnectionMap.load("data/module_connection_combined_2020_0518_helixray.txt");

}

//__________________________________________________________________________________________
void addOuterTrackerHits(SDL::Event& event)
{
    // Adding hits to modules
    for (auto&& [ihit, data] : iter::enumerate(iter::zip(trk.ph2_x(), trk.ph2_y(), trk.ph2_z(), trk.ph2_subdet(), trk.ph2_detId())))
    {

        auto&& [x, y, z, subdet, detid] = data;

        if (not (subdet == 5 or subdet == 4))
            continue;

        if (abs(ana.pdg_id) != 0)
        {
            bool found = false;
            for (auto& isimhit : trk.ph2_simHitIdx()[ihit])
            {
                int isimtrk = trk.simhit_simTrkIdx()[isimhit];

                if (abs(trk.sim_pdgId()[isimtrk]) == ana.pdg_id)
                {
                    found = true;
                    break;
                }

            }

            if (not found)
                continue;

        }

        // Takes two arguments, SDL::Hit, and detId
        // SDL::Event internally will structure whether we already have the module instance or we need to create a new one.
        event.addHitToModule(
                // a hit
                SDL::Hit(x, y, z, ihit),
                // add to module with "detId"
                detid
                );

    }
}

//__________________________________________________________________________________________
TVector3 r3FromPCA(const TVector3& p3, const float dxy, const float dz){
  const float pt = p3.Pt();
  const float p = p3.Mag();
  const float vz = dz*pt*pt/p/p;

  const float vx = -dxy*p3.y()/pt - p3.x()/p*p3.z()/p*dz;
  const float vy =  dxy*p3.x()/pt - p3.y()/p*p3.z()/p*dz;
  return TVector3(vx, vy, vz);
}


//__________________________________________________________________________________________
void addPixelSegments(SDL::Event& event, int isimtrk)
{
    for (auto&& [iSeed, _] : iter::enumerate(trk.see_stateTrajGlbPx()))
    {

        if (isimtrk >= 0)
        {
            bool match = false;
            for (auto& seed_simtrkidx : trk.see_simTrkIdx()[iSeed])
            {
                if (seed_simtrkidx == isimtrk)
                {
                    match = true;
                }
            }
            if (not match)
                continue;
        }

        if (abs(ana.pdg_id) != 0)
        {
            bool match = false;
            for (auto& seed_simtrkidx : trk.see_simTrkIdx()[iSeed])
            {
                if (abs(trk.sim_pdgId()[seed_simtrkidx]) == ana.pdg_id)
                {
                    match = true;
                }
            }
            if (not match)
                continue;
        }

        TVector3 p3LH(trk.see_stateTrajGlbPx()[iSeed], trk.see_stateTrajGlbPy()[iSeed], trk.see_stateTrajGlbPz()[iSeed]);
        TVector3 r3LH(trk.see_stateTrajGlbX()[iSeed], trk.see_stateTrajGlbY()[iSeed], trk.see_stateTrajGlbZ()[iSeed]);
        TVector3 p3PCA(trk.see_px()[iSeed], trk.see_py()[iSeed], trk.see_pz()[iSeed]);
        TVector3 r3PCA(r3FromPCA(p3PCA, trk.see_dxy()[iSeed], trk.see_dz()[iSeed]));
        auto const& seedHitsV = trk.see_hitIdx()[iSeed];
        auto const& seedHitTypesV = trk.see_hitType()[iSeed];

        // /// track algorithm; partial copy from TrackBase.h
        // enum class TrackAlgorithm {
        //     undefAlgorithm = 0,
        //     ctf = 1,
        //     duplicateMerge = 2,
        //     cosmics = 3,
        //     initialStep = 4,
        //     lowPtTripletStep = 5,
        //     pixelPairStep = 6,
        //     detachedTripletStep = 7,
        //     mixedTripletStep = 8,
        //     pixelLessStep = 9,
        //     tobTecStep = 10,
        //     jetCoreRegionalStep = 11,
        //     conversionStep = 12,
        //     muonSeededStepInOut = 13,
        //     muonSeededStepOutIn = 14,
        //     outInEcalSeededConv = 15, inOutEcalSeededConv = 16,
        //     nuclInter = 17,
        //     standAloneMuon = 18, globalMuon = 19, cosmicStandAloneMuon = 20, cosmicGlobalMuon = 21,
        //     // Phase1
        //     highPtTripletStep = 22, lowPtQuadStep = 23, detachedQuadStep = 24,
        //     reservedForUpgrades1 = 25, reservedForUpgrades2 = 26,
        //     bTagGhostTracks = 27,
        //     beamhalo = 28,
        //     gsf = 29
        // };
        bool good_seed_type = false;
        if (trk.see_algo()[iSeed] == 4) good_seed_type = true;
        if (trk.see_algo()[iSeed] == 5) good_seed_type = true;
        if (trk.see_algo()[iSeed] == 7) good_seed_type = true;
        if (trk.see_algo()[iSeed] == 22) good_seed_type = true;
        if (trk.see_algo()[iSeed] == 23) good_seed_type = true;
        if (trk.see_algo()[iSeed] == 24) good_seed_type = true;
        if (not good_seed_type)
            continue;
        // if (trk.see_algo()[iSeed] != 4 and trk.see_algo()[iSeed] != 5)
        //     continue;
        // if (trk.see_algo()[iSeed] < 4 and trk.see_algo()[iSeed] > 8)
        //     continue;

        int nHits = seedHitsV.size();

        //assert(nHits == 4);
        //for (int iH = 0; iH < nHits; ++iH){
        //    //FIXME: make this configurable
        //    assert(seedHitTypesV[iH] == 0);
        //}

        // float seedSD_mdRef_pixL = HitIndexWithType(trk.see_hitIdx()[iSeed][0], HitType(trk.see_hitType()[iSeed][0])).indexWithType;
        // float seedSD_mdRef_pixU = HitIndexWithType(trk.see_hitIdx()[iSeed][1], HitType(trk.see_hitType()[iSeed][1])).indexWithType;
        int seedSD_mdRef_pixL = trk.see_hitIdx()[iSeed][0];
        int seedSD_mdRef_pixU = trk.see_hitIdx()[iSeed][1];
        TVector3 seedSD_mdRef_r3 = r3PCA;
        float seedSD_mdRef_rt = r3PCA.Pt();
        float seedSD_mdRef_z = r3PCA.Z();
        float seedSD_mdRef_r = r3PCA.Mag();
        float seedSD_mdRef_phi = r3PCA.Phi();
        float seedSD_mdRef_alpha = r3PCA.DeltaPhi(p3PCA);
        // const int itpRL = simsPerHitAll(seedSD_mdRef_pixL); // TODO: best sim trk idx SO PERHAPS NOT NEEDED
        // const int itpRU = simsPerHitAll(seedSD_mdRef_pixU); // TODO: best sim trk idx
        // float seedSD_mdRef_itp = itpRL;
        // float seedSD_mdRef_ntp = 1;
        // float seedSD_mdRef_itpLL = itpRL;

        // if (itpRL >= 0 && itpRU >= 0)
        // {
        //     if (itpRL == itpRU)
        //     {
        //         seedSD_mdRef_ntp = 2;
        //     }
        // }
        // else if (itpRL == -1 && itpRU == -1)
        // {
        //     seedSD_mdRef_ntp = 0;
        // }
        // else if (itpRU >= 0)
        // {
        //     seedSD_mdRef_itp = itpRU;
        // }

        int seedSD_mdOut_pixL = trk.see_hitIdx()[iSeed][2];
        int seedSD_mdOut_pixU = trk.see_hitIdx()[iSeed][3];
        // if (nPix >= 4)
        //     seedSD_mdOut_pixU = trk.see_hitIdx()[iSeed][3];
        TVector3 seedSD_mdOut_r3 = r3LH;
        float seedSD_mdOut_rt = r3LH.Pt();
        float seedSD_mdOut_z = r3LH.Z();
        float seedSD_mdOut_r = r3LH.Mag();
        float seedSD_mdOut_phi = r3LH.Phi();
        float seedSD_mdOut_alpha = r3LH.DeltaPhi(p3LH);
        // const int itpOL = simsPerHitAll(seedSD_mdOut_pixL);
        // const int itpOU = simsPerHitAll(seedSD_mdOut_pixU);
        // float seedSD_mdOut_itp = itpOL;
        // float seedSD_mdOut_ntp = 1;
        // float seedSD_mdOut_itpLL = itpOL;
        // if (itpOL >= 0 && itpOU >= 0)
        // {
        //     if (itpOL == itpOU)
        //     {
        //         seedSD_mdOut_ntp = 2;
        //     }
        // }
        // else if (itpOL == -1 && itpOU == -1)
        // {
        //     seedSD_mdOut_ntp = 0;
        // }
        // else if (itpOU >= 0)
        // {
        //     seedSD_mdOut_itp = itpOU;
        // }

        float seedSD_iRef = iSeed;
        float seedSD_iOut = iSeed;
        TVector3 seedSD_r3 = r3LH;
        float seedSD_rt = r3LH.Pt();
        float seedSD_rtInv = 1.f / seedSD_rt;
        float seedSD_z = seedSD_r3.Z();
        TVector3 seedSD_p3 = p3LH;
        float seedSD_alpha = r3LH.DeltaPhi(p3LH);
        float seedSD_dr = (r3LH - r3PCA).Pt();
        float seedSD_d = seedSD_rt - r3PCA.Pt();
        float seedSD_zeta = seedSD_p3.Pt() / seedSD_p3.Z();

        // std::map<int, int> tps;
        // int seedSD_itp = -1;
        // int seedSD_ntp = 0;
        // tps[itpRL]++;
        // tps[itpRU]++;
        // tps[itpOL]++;
        // tps[itpOU]++;

        // for (auto m : tps)
        // {
        //     if (m.first >= 0 && m.second > seedSD_ntp)
        //     {
        //         seedSD_itp = m.first;
        //         seedSD_ntp = m.second;
        //     }
        // }

        // // LL indexing is not very useful for seeds; to it anyways for "full" coverage
        // int seedSD_itpLL = -1;
        // int seedSD_ntpLL = 0;
        // tps.clear();
        // tps[itpRL]++;
        // tps[itpOL]++;
        // for (auto m : tps)
        // {
        //     if (m.first >= 0 && m.second > seedSD_ntpLL)
        //     {
        //         seedSD_itpLL = m.first;
        //         seedSD_ntpLL = m.second;
        //     }
        // }

        // For this segment in this layer, I need to add this many
        // I suppose for now... I can have a SDL::Module of ginormous size
        // with a special detId of ... 1? (does this require Module.h modification?)
        // SDL::Hit
        // SDL::Hit
        // SDL::MiniDoublet
        // SDL::Hit
        // SDL::Hit
        // SDL::MiniDoublet
        // SDL::Segment
        // std::cout <<  " seedSD_mdRef_rt: " << seedSD_mdRef_rt <<  std::endl;
        // std::cout <<  " seedSD_mdRef_z: " << seedSD_mdRef_z <<  std::endl;
        // std::cout <<  " seedSD_mdRef_alpha: " << seedSD_mdRef_alpha <<  std::endl;
        // std::cout <<  " seedSD_mdOut_rt: " << seedSD_mdOut_rt <<  std::endl;
        // std::cout <<  " seedSD_mdOut_z: " << seedSD_mdOut_z <<  std::endl;
        // std::cout <<  " seedSD_mdOut_alpha: " << seedSD_mdOut_alpha <<  std::endl;
        // std::cout <<  " trk.pix_x()[trk.see_hitIdx()[iSeed][0]]: " << trk.pix_x()[trk.see_hitIdx()[iSeed][0]] <<  std::endl;
        // std::cout <<  " trk.pix_y()[trk.see_hitIdx()[iSeed][0]]: " << trk.pix_y()[trk.see_hitIdx()[iSeed][0]] <<  std::endl;
        // std::cout <<  " sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][0]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][0]],2)): " << sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][0]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][0]],2)) <<  std::endl;
        // std::cout <<  " trk.pix_z()[trk.see_hitIdx()[iSeed][0]]: " << trk.pix_z()[trk.see_hitIdx()[iSeed][0]] <<  std::endl;
        // std::cout <<  " trk.pix_x()[trk.see_hitIdx()[iSeed][1]]: " << trk.pix_x()[trk.see_hitIdx()[iSeed][1]] <<  std::endl;
        // std::cout <<  " trk.pix_y()[trk.see_hitIdx()[iSeed][1]]: " << trk.pix_y()[trk.see_hitIdx()[iSeed][1]] <<  std::endl;
        // std::cout <<  " sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][1]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][1]],2)): " << sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][1]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][1]],2)) <<  std::endl;
        // std::cout <<  " trk.pix_z()[trk.see_hitIdx()[iSeed][1]]: " << trk.pix_z()[trk.see_hitIdx()[iSeed][1]] <<  std::endl;
        // std::cout <<  " trk.pix_x()[trk.see_hitIdx()[iSeed][2]]: " << trk.pix_x()[trk.see_hitIdx()[iSeed][2]] <<  std::endl;
        // std::cout <<  " trk.pix_y()[trk.see_hitIdx()[iSeed][2]]: " << trk.pix_y()[trk.see_hitIdx()[iSeed][2]] <<  std::endl;
        // std::cout <<  " sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][2]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][2]],2)): " << sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][2]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][2]],2)) <<  std::endl;
        // std::cout <<  " trk.pix_z()[trk.see_hitIdx()[iSeed][2]]: " << trk.pix_z()[trk.see_hitIdx()[iSeed][2]] <<  std::endl;
        // std::cout <<  " trk.pix_x()[trk.see_hitIdx()[iSeed][3]]: " << trk.pix_x()[trk.see_hitIdx()[iSeed][3]] <<  std::endl;
        // std::cout <<  " trk.pix_y()[trk.see_hitIdx()[iSeed][3]]: " << trk.pix_y()[trk.see_hitIdx()[iSeed][3]] <<  std::endl;
        // std::cout <<  " sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][3]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][3]],2)): " << sqrt(pow(trk.pix_x()[trk.see_hitIdx()[iSeed][3]],2)+pow(trk.pix_y()[trk.see_hitIdx()[iSeed][3]],2)) <<  std::endl;
        // std::cout <<  " trk.pix_z()[trk.see_hitIdx()[iSeed][3]]: " << trk.pix_z()[trk.see_hitIdx()[iSeed][3]] <<  std::endl;

        // Inner most hit
        std::vector<SDL::Hit> hits;
        int hitidx0 = trk.see_hitIdx()[iSeed][0];
        int hittype0 = trk.see_hitType()[iSeed][0];
        // hits.push_back(SDL::Hit(trk.pix_x()[hitidx0], trk.pix_y()[hitidx0], trk.pix_z()[hitidx0], hitidx0));
        hits.push_back(SDL::Hit(r3PCA.X(), r3PCA.Y(), r3PCA.Z(), hitidx0));
        int hitidx1 = trk.see_hitIdx()[iSeed][1];
        int hittype1 = trk.see_hitType()[iSeed][1];
        // hits.push_back(SDL::Hit(trk.pix_x()[hitidx1], trk.pix_y()[hitidx1], trk.pix_z()[hitidx1], hitidx1));
        hits.push_back(SDL::Hit(r3PCA.X(), r3PCA.Y(), r3PCA.Z(), hitidx1));
        int hitidx2 = trk.see_hitIdx()[iSeed][2];
        int hittype2 = trk.see_hitType()[iSeed][2];
        // hits.push_back(SDL::Hit(trk.pix_x()[hitidx2], trk.pix_y()[hitidx2], trk.pix_z()[hitidx2], hitidx2));
        hits.push_back(SDL::Hit(r3LH.X(), r3LH.Y(), r3LH.Z(), hitidx2));
        int hitidx3 = trk.see_hitIdx()[iSeed].size() > 3 ? trk.see_hitIdx()[iSeed][3] : trk.see_hitIdx()[iSeed][2]; // repeat last one if triplet
        int hittype3 = trk.see_hitIdx()[iSeed].size() > 3 ? trk.see_hitType()[iSeed][3] : trk.see_hitIdx()[iSeed][2]; // repeat last one if triplet
        // hits.push_back(SDL::Hit(trk.pix_x()[hitidx3], trk.pix_y()[hitidx3], trk.pix_z()[hitidx3], hitidx3));
        hits.push_back(SDL::Hit(r3LH.X(), r3LH.Y(), r3LH.Z(), hitidx3));

        float pixelSegmentDeltaPhiChange = r3LH.DeltaPhi(p3LH);
        float ptIn = p3LH.Pt();
        float ptErr = trk.see_ptErr()[iSeed];
        float etaErr = trk.see_etaErr()[iSeed];
        float px = p3LH.X();
        float py = p3LH.Y();
        float pz = p3LH.Z();

        if ((ptIn > 0.7) and (fabs(p3LH.Eta()) < 3))
            event.addPixelSegmentsToEvent(hits, pixelSegmentDeltaPhiChange, ptIn, ptErr, px, py, pz, etaErr, iSeed);

    }
}

//__________________________________________________________________________________________
void addOuterTrackerHitsFromSimTrack(SDL::Event& event, int isimtrk)
{

    // loop over the simulated hits
    for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
    {

        // Retrieve the sim hit idx
        unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

        // To clean up low pt junk in mini-doublet efficiency plots
        if (isMuonCurlingHit(isimtrk, ith_hit))
            break;

        // list of reco hit matched to this sim hit
        for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
        {

            // Get the recohit type
            int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

            // Consider only ph2 hits (i.e. outer tracker hits)
            if (recohittype == 4)
            {

                int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];

                event.addHitToModule(
                        // a hit
                        SDL::Hit(trk.ph2_x()[ihit], trk.ph2_y()[ihit], trk.ph2_z()[ihit], ihit),
                        // add to module with "detId"
                        trk.ph2_detId()[ihit]
                        );

            }

        }

    }
}

//__________________________________________________________________________________________
void printHitSummary(SDL::Event& event)
{
    if (ana.verbose != 0) std::cout << "Summary of hits" << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits: " << event.getNumberOfHits() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in barrel layer 1: " << event.getNumberOfHitsByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in barrel layer 2: " << event.getNumberOfHitsByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in barrel layer 3: " << event.getNumberOfHitsByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in barrel layer 4: " << event.getNumberOfHitsByLayerBarrel(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in barrel layer 5: " << event.getNumberOfHitsByLayerBarrel(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in barrel layer 6: " << event.getNumberOfHitsByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in endcap layer 1: " << event.getNumberOfHitsByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in endcap layer 2: " << event.getNumberOfHitsByLayerEndcap(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in endcap layer 3: " << event.getNumberOfHitsByLayerEndcap(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in endcap layer 4: " << event.getNumberOfHitsByLayerEndcap(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits in endcap layer 5: " << event.getNumberOfHitsByLayerEndcap(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 1: " << event.getNumberOfHitsByLayerBarrelUpperModule(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 2: " << event.getNumberOfHitsByLayerBarrelUpperModule(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 3: " << event.getNumberOfHitsByLayerBarrelUpperModule(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 4: " << event.getNumberOfHitsByLayerBarrelUpperModule(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 5: " << event.getNumberOfHitsByLayerBarrelUpperModule(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 6: " << event.getNumberOfHitsByLayerBarrelUpperModule(5) << std::endl;
}

//__________________________________________________________________________________________
void printMiniDoubletSummary(SDL::Event& event)
{
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced: " << event.getNumberOfMiniDoublets() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced barrel layer 1: " << event.getNumberOfMiniDoubletsByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced barrel layer 2: " << event.getNumberOfMiniDoubletsByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced barrel layer 3: " << event.getNumberOfMiniDoubletsByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced barrel layer 4: " << event.getNumberOfMiniDoubletsByLayerBarrel(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced barrel layer 5: " << event.getNumberOfMiniDoubletsByLayerBarrel(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced barrel layer 6: " << event.getNumberOfMiniDoubletsByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced endcap layer 1: " << event.getNumberOfMiniDoubletsByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced endcap layer 2: " << event.getNumberOfMiniDoubletsByLayerEndcap(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced endcap layer 3: " << event.getNumberOfMiniDoubletsByLayerEndcap(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced endcap layer 4: " << event.getNumberOfMiniDoubletsByLayerEndcap(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets produced endcap layer 5: " << event.getNumberOfMiniDoubletsByLayerEndcap(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered: " << event.getNumberOfMiniDoubletCandidates() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered barrel layer 1: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered barrel layer 2: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered barrel layer 3: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered barrel layer 4: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered barrel layer 5: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered barrel layer 6: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered endcap layer 1: " << event.getNumberOfMiniDoubletCandidatesByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered endcap layer 2: " << event.getNumberOfMiniDoubletCandidatesByLayerEndcap(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered endcap layer 3: " << event.getNumberOfMiniDoubletCandidatesByLayerEndcap(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered endcap layer 4: " << event.getNumberOfMiniDoubletCandidatesByLayerEndcap(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Mini-doublets considered endcap layer 5: " << event.getNumberOfMiniDoubletCandidatesByLayerEndcap(4) << std::endl;
}

//__________________________________________________________________________________________
void printSegmentSummary(SDL::Event& event)
{
    if (ana.verbose != 0) std::cout << "# of Pixel Segments: " << event.getPixelLayer().getSegmentPtrs().size() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced: " << event.getNumberOfSegments() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced barrel layer 1-2: " << event.getNumberOfSegmentsByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced barrel layer 2-3: " << event.getNumberOfSegmentsByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced barrel layer 3-4: " << event.getNumberOfSegmentsByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced barrel layer 4-5: " << event.getNumberOfSegmentsByLayerBarrel(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced barrel layer 5-6: " << event.getNumberOfSegmentsByLayerBarrel(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced endcap layer 1-2: " << event.getNumberOfSegmentsByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced endcap layer 2-3: " << event.getNumberOfSegmentsByLayerEndcap(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced endcap layer 3-4: " << event.getNumberOfSegmentsByLayerEndcap(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments produced endcap layer 4-5: " << event.getNumberOfSegmentsByLayerEndcap(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Segments produced layer 6: " << event.getNumberOfSegmentsByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered: " << event.getNumberOfSegmentCandidates() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 1-2: " << event.getNumberOfSegmentCandidatesByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 2-3: " << event.getNumberOfSegmentCandidatesByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 3-4: " << event.getNumberOfSegmentCandidatesByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 4-5: " << event.getNumberOfSegmentCandidatesByLayerBarrel(3) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 5-6: " << event.getNumberOfSegmentCandidatesByLayerBarrel(4) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 1-2: " << event.getNumberOfSegmentCandidatesByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 2-3: " << event.getNumberOfSegmentCandidatesByLayerEndcap(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 3-4: " << event.getNumberOfSegmentCandidatesByLayerEndcap(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Segments considered layer 4-5: " << event.getNumberOfSegmentCandidatesByLayerEndcap(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Segments considered layer 6: " << event.getNumberOfSegmentCandidatesByLayerBarrel(5) << std::endl;
}

//__________________________________________________________________________________________
void printMiniDoubletConnectionMultiplicitiesSummary(SDL::Event& event)
{
    if (ana.verbose != 0) std::cout << "Printing connection information" << std::endl;
    if (ana.verbose != 0)
    {
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 1, 1);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 1, 2);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 1, 3);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 1, 4);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 1, 5);
        std::cout << "--------" << std::endl;

        printMiniDoubletConnectionMultiplicitiesBarrel(event, 2, 1, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 2, 1);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 2, 2);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 2, 3);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 2, 4);
        std::cout << "--------" << std::endl;

        printMiniDoubletConnectionMultiplicitiesBarrel(event, 3, 2, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 3, 1, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 3, 1);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 3, 2);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 3, 3);
        std::cout << "--------" << std::endl;

        printMiniDoubletConnectionMultiplicitiesBarrel(event, 4, 3, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 4, 2, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 4, 1, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 4, 1);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 4, 2);
        std::cout << "--------" << std::endl;

        printMiniDoubletConnectionMultiplicitiesBarrel(event, 5, 4, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 5, 3, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 5, 2, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 5, 1, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 5, 1);
        std::cout << "--------" << std::endl;

        printMiniDoubletConnectionMultiplicitiesBarrel(event, 6, 5, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 6, 4, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 6, 3, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 6, 2, true);
        printMiniDoubletConnectionMultiplicitiesBarrel(event, 6, 1, true);
        std::cout << "--------" << std::endl;
    }
}

//__________________________________________________________________________________________
void printTripletSummary(SDL::Event& event)
{
    // ----------------
    if (ana.verbose != 0) std::cout << "# of Triplets produced: " << event.getNumberOfTriplets() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets produced layer 1-2-3: " << event.getNumberOfTripletsByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets produced layer 2-3-4: " << event.getNumberOfTripletsByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets produced layer 3-4-5: " << event.getNumberOfTripletsByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets produced layer 4-5-6: " << event.getNumberOfTripletsByLayerBarrel(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Triplets produced layer 5: " << event.getNumberOfTripletsByLayerBarrel(4) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Triplets produced layer 6: " << event.getNumberOfTripletsByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets considered: " << event.getNumberOfTripletCandidates() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets considered layer 1-2-3: " << event.getNumberOfTripletCandidatesByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets considered layer 2-3-4: " << event.getNumberOfTripletCandidatesByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets considered layer 3-4-5: " << event.getNumberOfTripletCandidatesByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Triplets considered layer 4-5-6: " << event.getNumberOfTripletCandidatesByLayerBarrel(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Triplets considered layer 5: " << event.getNumberOfTripletCandidatesByLayerBarrel(4) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Triplets considered layer 6: " << event.getNumberOfTripletCandidatesByLayerBarrel(5) << std::endl;
    // ----------------
}

//__________________________________________________________________________________________
void printTrackletSummary(SDL::Event& event)
{
    // ----------------
    if (ana.verbose != 0) std::cout << "# of Tracklets produced with pixel segment: " << event.getPixelLayer().getTrackletPtrs().size() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets produced: " << event.getNumberOfTracklets() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 1-2-3-4: " << event.getNumberOfTrackletsByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 2-3-4-5: " << event.getNumberOfTrackletsByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 3-4-5-6: " << event.getNumberOfTrackletsByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 1-2-3-4: " << event.getNumberOfTrackletsByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 2-3-4-5: " << event.getNumberOfTrackletsByLayerEndcap(1) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 4: " << event.getNumberOfTrackletsByLayerBarrel(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 5: " << event.getNumberOfTrackletsByLayerBarrel(4) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 6: " << event.getNumberOfTrackletsByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets considered: " << event.getNumberOfTrackletCandidates() << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 1-2-3-4: " << event.getNumberOfTrackletCandidatesByLayerBarrel(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 2-3-4-5: " << event.getNumberOfTrackletCandidatesByLayerBarrel(1) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 3-4-5-6: " << event.getNumberOfTrackletCandidatesByLayerBarrel(2) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 1-2-3-4: " << event.getNumberOfTrackletCandidatesByLayerEndcap(0) << std::endl;
    if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 2-3-4-5: " << event.getNumberOfTrackletCandidatesByLayerEndcap(1) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 4: " << event.getNumberOfTrackletCandidatesByLayerBarrel(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 5: " << event.getNumberOfTrackletCandidatesByLayerBarrel(4) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 6: " << event.getNumberOfTrackletCandidatesByLayerBarrel(5) << std::endl;
    // ----------------
}

//__________________________________________________________________________________________
void printTrackCandidateSummary(SDL::Event& event)
{
    // ----------------
    if (ana.verbose != 0) std::cout << "# of TrackCandidates produced: " << event.getNumberOfTrackCandidates() << std::endl;
    if (ana.verbose != 0) std::cout << "# of TrackCandidates produced layer 1-2-3-4-5-6: " << event.getNumberOfTrackCandidatesByLayerBarrel(0) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates produced layer 2: " << event.getNumberOfTrackCandidatesByLayerBarrel(1) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates produced layer 3: " << event.getNumberOfTrackCandidatesByLayerBarrel(2) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates produced layer 4: " << event.getNumberOfTrackCandidatesByLayerBarrel(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates produced layer 5: " << event.getNumberOfTrackCandidatesByLayerBarrel(4) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates produced layer 6: " << event.getNumberOfTrackCandidatesByLayerBarrel(5) << std::endl;
    if (ana.verbose != 0) std::cout << "# of TrackCandidates considered: " << event.getNumberOfTrackCandidateCandidates() << std::endl;
    if (ana.verbose != 0) std::cout << "# of TrackCandidates considered layer 1-2-3-4-5-6: " << event.getNumberOfTrackCandidateCandidatesByLayerBarrel(0) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates considered layer 2: " << event.getNumberOfTrackCandidateCandidatesByLayerBarrel(1) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates considered layer 3: " << event.getNumberOfTrackCandidateCandidatesByLayerBarrel(2) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates considered layer 4: " << event.getNumberOfTrackCandidateCandidatesByLayerBarrel(3) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates considered layer 5: " << event.getNumberOfTrackCandidateCandidatesByLayerBarrel(4) << std::endl;
    // if (ana.verbose != 0) std::cout << "# of TrackCandidates considered layer 6: " << event.getNumberOfTrackCandidateCandidatesByLayerBarrel(5) << std::endl;
    // ----------------

}

//__________________________________________________________________________________________
void runMiniDoublet(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Mini-Doublet start" << std::endl;
    my_timer.Start();
    event.createMiniDoublets();
    float md_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Mini-doublet processing time: " << md_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runSegment(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Segment start" << std::endl;
    my_timer.Start(kFALSE);
    event.createSegmentsWithModuleMap();
    float sg_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Segment processing time: " << sg_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTriplet(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Triplet start" << std::endl;
    my_timer.Start(kFALSE);
    event.createTriplets();
    float tp_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Triplet processing time: " << tp_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTracklet(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Tracklet Without A Gap start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTracklets();
    event.createTrackletsWithModuleMap();
    // event.createTrackletsWithAGapWithModuleMap();
    // event.createTrackletsWithTwoGapsWithModuleMap();
    // event.createTrackletsViaNavigation();
    float tl_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Tracklet Without A Gap processing time: " << tl_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackletWithAGap(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Tracklet With A Gap start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTracklets();
    event.createTrackletsWithAGapWithModuleMap();
    // event.createTrackletsWithTwoGapsWithModuleMap();
    // event.createTrackletsViaNavigation();
    float tl_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Tracklet With A Gap processing time: " << tl_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackletTest_PixelSegment_v1(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Tracklet with pixel segment start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTracklets();
    // event.createTrackletsWithModuleMap();
    // event.createTrackletsWithAGapWithModuleMap();
    if (ana.verbose != 0) event.setLogLevel(SDL::Log_Debug);
    event.createTrackletsWithPixelAndBarrel();
    if (ana.verbose != 0) event.setLogLevel(SDL::Log_Nothing);
    // event.createTrackletsWithPixelAndBarrel(SDL::AllComb_TLAlgo);
    // event.createTrackletsWithTwoGapsWithModuleMap();
    // event.createTrackletsViaNavigation();
    float tl_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Tracklet with pixel segment processing time: " << tl_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackletTest_v1(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco Tracklet start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTracklets();
    event.createTrackletsWithModuleMap();
    event.createTrackletsWithAGapWithModuleMap();
    // event.createTrackletsWithTwoGapsWithModuleMap();
    // event.createTrackletsViaNavigation();
    float tl_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco Tracklet processing time: " << tl_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackCandidate(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTrackCandidatesFromTriplets();
    // event.createTrackCandidates();
    event.createTrackCandidatesFromTracklets();
    float tc_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate processing time: " << tc_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackCandidateTest_v1(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTrackCandidatesFromTriplets();
    // event.createTrackCandidates();
    event.createTrackCandidatesTest_v1();
    float tc_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate processing time: " << tc_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackCandidateTest_v2(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTrackCandidatesFromTriplets();
    // event.createTrackCandidates();
    event.createTrackCandidatesTest_v2();
    float tc_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate processing time: " << tc_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runTrackCandidateTest_v3(SDL::Event& event)
{
    TStopwatch my_timer;
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate start" << std::endl;
    my_timer.Start(kFALSE);
    // event.createTrackCandidatesFromTriplets();
    // event.createTrackCandidates();
    event.createTrackCandidatesTest_v3();
    float tc_elapsed = my_timer.RealTime();
    if (ana.verbose != 0) std::cout << "Reco TrackCandidate processing time: " << tc_elapsed << " secs" << std::endl;
}

//__________________________________________________________________________________________
void runSDL(SDL::Event& event)
{

    printHitSummary(event);
    runMiniDoublet(event);
    printMiniDoubletSummary(event);
    runSegment(event);
    printSegmentSummary(event);
    runTriplet(event);
    printTripletSummary(event);
    runTracklet(event);
    printTrackletSummary(event);
    runTrackCandidate(event);
    printTrackCandidateSummary(event);

}

//__________________________________________________________________________________________
void runSDLTest_v1(SDL::Event& event)
{

    printHitSummary(event);
    runMiniDoublet(event);
    printMiniDoubletSummary(event);
    runSegment(event);
    printSegmentSummary(event);
    runTriplet(event);
    printTripletSummary(event);
    runTrackletTest_v1(event);
    printTrackletSummary(event);
    runTrackCandidateTest_v1(event);
    printTrackCandidateSummary(event);

}

//__________________________________________________________________________________________
void runSDLTest_v2(SDL::Event& event)
{

    printHitSummary(event);
    runMiniDoublet(event);
    printMiniDoubletSummary(event);
    runSegment(event);
    printSegmentSummary(event);
    runTriplet(event);
    printTripletSummary(event);
    runTracklet(event);
    printTrackletSummary(event);
    runTrackCandidateTest_v2(event);
    printTrackCandidateSummary(event);

}

//__________________________________________________________________________________________
void runSDLTest_PixelSegments(SDL::Event& event)
{

    printHitSummary(event);
    runMiniDoublet(event);
    printMiniDoubletSummary(event);
    runSegment(event);
    printSegmentSummary(event);
    runTriplet(event);
    printTripletSummary(event);
    runTrackletTest_PixelSegment_v1(event);
    printTrackletSummary(event);
    // runTrackCandidateTest_v2(event);
    // printTrackCandidateSummary(event);

}

//__________________________________________________________________________________________
void printSimTrack(int isimtrk)
{

    float vx = trk.simvtx_x()[0];
    float vy = trk.simvtx_y()[0];
    float vz = trk.simvtx_z()[0];
    float pt = trk.sim_pt()[isimtrk];
    float eta = trk.sim_eta()[isimtrk];
    float phi = trk.sim_phi()[isimtrk];
    float charge = trk.sim_q()[isimtrk];

    // Construct helix object
    SDLMath::Helix helix(pt, eta, phi, vx, vy, vz, charge);

    std::cout <<  " pt: " << pt <<  " eta: " << eta <<  " phi: " << phi <<  " vx: " << vx <<  " vy: " << vy <<  " vz: " << vz <<  " charge: " << charge <<  std::endl;

    // List of detids
    std::vector<unsigned int> lower_module_detids;
    std::vector<unsigned int> upper_module_detids;

    // List of pdgids
    std::vector<int> lower_module_pdgids;
    std::vector<int> upper_module_pdgids;

    // List of layers
    std::vector<unsigned int> lower_module_layers;
    std::vector<unsigned int> upper_module_layers;

    // List of n reco hits
    std::vector<int> lower_module_nrecos;
    std::vector<int> upper_module_nrecos;

    // List of n reco hits
    std::vector<float> lower_module_hit_xs;
    std::vector<float> upper_module_hit_xs;
    std::vector<float> lower_module_hit_ys;
    std::vector<float> upper_module_hit_ys;
    std::vector<float> lower_module_hit_zs;
    std::vector<float> upper_module_hit_zs;
    std::vector<float> lower_module_hit_rs;
    std::vector<float> upper_module_hit_rs;

    // sim hit selection
    // - Must be same pdgid as the parent particle
    // - expected radius based on helix is < 2% away
    for (unsigned int isimhit = 0; isimhit < trk.sim_simHitIdx()[isimtrk].size(); ++isimhit)
    {

        unsigned int isimhitidx = trk.sim_simHitIdx()[isimtrk][isimhit];

        // if (not (trk.simhit_particle()[isimhitidx] == trk.sim_pdgId()[isimtrk]))
        //     continue;

        if (not (trk.simhit_subdet()[isimhitidx] == 4 or trk.simhit_subdet()[isimhitidx] == 5))
            continue;

        // Sim hit vector
        std::vector<float> point = {trk.simhit_x()[isimhitidx], trk.simhit_y()[isimhitidx], trk.simhit_z()[isimhitidx]};

        // Inferring parameter t of helix parametric function via z position
        float t = helix.infer_t(point);

        // If the best fit is more than pi parameter away then it's a returning hit and should be ignored
        if (not (t <= M_PI))
            continue;

        // Expected hit position with given z
        auto [x, y, z, r] = helix.get_helix_point(t);

        // ( expected_r - simhit_r ) / expected_r
        float drfrac = abs(helix.compare_radius(point)) / r;

        // Require that the simhit is within 2% of the expected radius
        if (not (drfrac < 0.02))
            continue;

        // Require that the simhit contains a matched reco
        if (not (trk.simhit_hitIdx()[isimhitidx].size() > 0))
            continue;

        // Sim hit detid
        unsigned int simhit_detid = trk.simhit_detId()[isimhitidx];
        int simhit_pdgid = trk.simhit_particle()[isimhitidx];
        int simhit_layer = 6 * (trk.simhit_subdet()[isimhitidx] == 4) + trk.simhit_layer()[isimhitidx];
        int simhit_nreco = trk.simhit_hitIdx()[isimhitidx].size();
        float simhit_hit_x = trk.simhit_x()[isimhitidx];
        float simhit_hit_y = trk.simhit_y()[isimhitidx];
        float simhit_hit_z = trk.simhit_z()[isimhitidx];
        float simhit_hit_r = sqrt(simhit_hit_x * simhit_hit_x + simhit_hit_y * simhit_hit_y);

        SDL::Module module(simhit_detid);

        if (module.isLower())
        {
            lower_module_detids.push_back(simhit_detid);
            lower_module_pdgids.push_back(simhit_pdgid);
            lower_module_layers.push_back(simhit_layer);
            lower_module_nrecos.push_back(simhit_nreco);
            lower_module_hit_xs.push_back(simhit_hit_x);
            lower_module_hit_ys.push_back(simhit_hit_y);
            lower_module_hit_zs.push_back(simhit_hit_z);
            lower_module_hit_rs.push_back(simhit_hit_r);
        }
        else
        {
            upper_module_detids.push_back(simhit_detid);
            upper_module_pdgids.push_back(simhit_pdgid);
            upper_module_layers.push_back(simhit_layer);
            upper_module_nrecos.push_back(simhit_nreco);
            upper_module_hit_xs.push_back(simhit_hit_x);
            upper_module_hit_ys.push_back(simhit_hit_y);
            upper_module_hit_zs.push_back(simhit_hit_z);
            upper_module_hit_rs.push_back(simhit_hit_r);
        }

    }

    std::vector<unsigned int> matched_lower_module_detids;

    for (auto& lower_module_detid : lower_module_detids)
    {
        SDL::Module lower_module(lower_module_detid);

        if (std::find(upper_module_detids.begin(), upper_module_detids.end(), lower_module.partnerDetId()) != upper_module_detids.end())
        {
            matched_lower_module_detids.push_back(lower_module_detid);
        }
    }

    std::vector<unsigned int> matched_layers;

    for (auto& matched_lower_module_detid : matched_lower_module_detids)
    {
        SDL::Module matched_lower_module(matched_lower_module_detid);
        unsigned int layer = matched_lower_module.layer();
        unsigned int subdet = matched_lower_module.subdet();
        matched_layers.push_back(layer + 6 * (subdet == 4));
    }

    std::vector<int> matched_unique_layers;
    std::vector<int> matched_unique_barrel_layers;
    std::vector<int> matched_unique_endcap_layers;
    for (auto&& i : iter::unique_justseen(matched_layers))
    {
        matched_unique_layers.push_back(i);
        if (i < 7)
            matched_unique_barrel_layers.push_back(i);
        else
            matched_unique_endcap_layers.push_back(i);
    }

    bool has_missing_hits = false;

    if (matched_unique_barrel_layers.size() > 0)
    {
        // std::cout <<  " matched_unique_barrel_layers.back(): " << matched_unique_barrel_layers.back() <<  std::endl;
        // std::cout <<  " matched_unique_barrel_layers.size(): " << matched_unique_barrel_layers.size() <<  std::endl;

        if (matched_unique_barrel_layers.back() != matched_unique_barrel_layers.size())
        {
            has_missing_hits = true;
        }
    }

    if (matched_unique_endcap_layers.size() > 0)
    {

        // std::cout <<  " (matched_unique_endcap_layers.back()-6): " << (matched_unique_endcap_layers.back()-6) <<  std::endl;
        // std::cout <<  " matched_unique_endcap_layers.size(): " << matched_unique_endcap_layers.size() <<  std::endl;

        if ((matched_unique_endcap_layers.back() - 6) != matched_unique_endcap_layers.size())
        {
            has_missing_hits = true;
        }
    }

    int nbarrel = matched_unique_barrel_layers.size();
    int nendcap = matched_unique_endcap_layers.size();

    bool is_clean_punch_through = has_missing_hits ? false : nbarrel + nendcap >= 6;

    bool verbose = true;

    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 8) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;
    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 9) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;
    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(),10) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;
    // if (std::find(matched_unique_layers.begin(), matched_unique_layers.end(),11) != matched_unique_layers.end() and std::find(matched_unique_layers.begin(), matched_unique_layers.end(), 7) == matched_unique_layers.end()) verbose = true;

    if (verbose)
    {

        std::cout << "lower_module_detids: ";
        for (auto& lower_module_detid : lower_module_detids) std::cout << std::setw(10) << lower_module_detid << " ";
        std::cout << std::endl;

        std::cout << "lower_module_pdgids: ";
        for (auto& lower_module_pdgid : lower_module_pdgids) std::cout << std::setw(10) << lower_module_pdgid << " ";
        std::cout << std::endl;

        std::cout << "lower_module_layers: ";
        for (auto& lower_module_layer : lower_module_layers) std::cout << std::setw(10) << lower_module_layer << " ";
        std::cout << std::endl;

        std::cout << "lower_module_nrecos: ";
        for (auto& lower_module_nreco : lower_module_nrecos) std::cout << std::setw(10) << lower_module_nreco << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_xs: ";
        for (auto& lower_module_hit_x : lower_module_hit_xs) std::cout << std::setw(10) << lower_module_hit_x << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_ys: ";
        for (auto& lower_module_hit_y : lower_module_hit_ys) std::cout << std::setw(10) << lower_module_hit_y << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_zs: ";
        for (auto& lower_module_hit_z : lower_module_hit_zs) std::cout << std::setw(10) << lower_module_hit_z << " ";
        std::cout << std::endl;

        std::cout << "lower_module_hit_rs: ";
        for (auto& lower_module_hit_r : lower_module_hit_rs) std::cout << std::setw(10) << lower_module_hit_r << " ";
        std::cout << std::endl;

        std::cout << "upper_module_detids: ";
        for (auto& upper_module_detid : upper_module_detids) std::cout << std::setw(10) << upper_module_detid << " ";
        std::cout << std::endl;

        std::cout << "upper_module_pdgids: ";
        for (auto& upper_module_pdgid : upper_module_pdgids) std::cout << std::setw(10) << upper_module_pdgid << " ";
        std::cout << std::endl;

        std::cout << "upper_module_layers: ";
        for (auto& upper_module_layer : upper_module_layers) std::cout << std::setw(10) << upper_module_layer << " ";
        std::cout << std::endl;

        std::cout << "upper_module_nrecos: ";
        for (auto& upper_module_nreco : upper_module_nrecos) std::cout << std::setw(10) << upper_module_nreco << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_xs: ";
        for (auto& upper_module_hit_x : upper_module_hit_xs) std::cout << std::setw(10) << upper_module_hit_x << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_ys: ";
        for (auto& upper_module_hit_y : upper_module_hit_ys) std::cout << std::setw(10) << upper_module_hit_y << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_zs: ";
        for (auto& upper_module_hit_z : upper_module_hit_zs) std::cout << std::setw(10) << upper_module_hit_z << " ";
        std::cout << std::endl;

        std::cout << "upper_module_hit_rs: ";
        for (auto& upper_module_hit_r : upper_module_hit_rs) std::cout << std::setw(10) << upper_module_hit_r << " ";
        std::cout << std::endl;

        std::cout << "matched_lower_module_detids: ";
        for (auto& matched_lower_module_detid : matched_lower_module_detids) std::cout << matched_lower_module_detid << " ";
        std::cout << std::endl;

        std::cout << "matched_layers: ";
        for (auto&& i : matched_layers) std::cout << i << " ";
        std::cout << std::endl;

        std::cout << "matched_unique_layers: ";
        for (auto&& i : matched_unique_layers) std::cout << i << " ";
        std::cout << std::endl;

    }

}

//__________________________________________________________________________________________
float drfracSimHitConsistentWithHelix(int isimtrk, int isimhitidx)
{
    // Read track parameters
    float vx = trk.simvtx_x()[0];
    float vy = trk.simvtx_y()[0];
    float vz = trk.simvtx_z()[0];
    float pt = trk.sim_pt()[isimtrk];
    float eta = trk.sim_eta()[isimtrk];
    float phi = trk.sim_phi()[isimtrk];
    float charge = trk.sim_q()[isimtrk];

    // Construct helix object
    SDLMath::Helix helix(pt, eta, phi, vx, vy, vz, charge);

    return drfracSimHitConsistentWithHelix(helix, isimhitidx);

}

//__________________________________________________________________________________________
float drfracSimHitConsistentWithHelix(SDLMath::Helix& helix, int isimhitidx)
{

    // Sim hit vector
    std::vector<float> point = {trk.simhit_x()[isimhitidx], trk.simhit_y()[isimhitidx], trk.simhit_z()[isimhitidx]};

    // Inferring parameter t of helix parametric function via z position
    float t = helix.infer_t(point);

    // If the best fit is more than pi parameter away then it's a returning hit and should be ignored
    if (not (t <= M_PI))
        return 999;

    // Expected hit position with given z
    auto [x, y, z, r] = helix.get_helix_point(t);

    // ( expected_r - simhit_r ) / expected_r
    float drfrac = abs(helix.compare_radius(point)) / r;

    return drfrac;

}

//__________________________________________________________________________________________
float distxySimHitConsistentWithHelix(int isimtrk, int isimhitidx)
{
    // Read track parameters
    float vx = trk.simvtx_x()[0];
    float vy = trk.simvtx_y()[0];
    float vz = trk.simvtx_z()[0];
    float pt = trk.sim_pt()[isimtrk];
    float eta = trk.sim_eta()[isimtrk];
    float phi = trk.sim_phi()[isimtrk];
    float charge = trk.sim_q()[isimtrk];

    // Construct helix object
    SDLMath::Helix helix(pt, eta, phi, vx, vy, vz, charge);

    return distxySimHitConsistentWithHelix(helix, isimhitidx);

}

//__________________________________________________________________________________________
float distxySimHitConsistentWithHelix(SDLMath::Helix& helix, int isimhitidx)
{

    // Sim hit vector
    std::vector<float> point = {trk.simhit_x()[isimhitidx], trk.simhit_y()[isimhitidx], trk.simhit_z()[isimhitidx]};

    // Inferring parameter t of helix parametric function via z position
    float t = helix.infer_t(point);

    // If the best fit is more than pi parameter away then it's a returning hit and should be ignored
    if (not (t <= M_PI))
        return 999;

    // Expected hit position with given z
    auto [x, y, z, r] = helix.get_helix_point(t);

    // ( expected_r - simhit_r ) / expected_r
    float distxy = helix.compare_xy(point);

    return distxy;

}
