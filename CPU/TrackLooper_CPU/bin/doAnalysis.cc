#include "doAnalysis.h"

// ./process INPUTFILEPATH OUTPUTFILE [NEVENTS]
int main(int argc, char** argv)
{

//********************************************************************************
//
// 1. Parsing options
//
//********************************************************************************

    // cxxopts is just a tool to parse argc, and argv easily

    // Grand option setting
    cxxopts::Options options("\n  $ doAnalysis",  "\n         **********************\n         *                    *\n         *       Looper       *\n         *                    *\n         **********************\n");

    // Read the options
    options.add_options()
        ("i,input"          , "Comma separated input file list OR if just a directory is provided it will glob all in the directory BUT must end with '/' for the path", cxxopts::value<std::string>())
        ("t,tree"           , "Name of the tree in the root file to open and loop over"                                             , cxxopts::value<std::string>())
        ("o,output"         , "Output file name"                                                                                    , cxxopts::value<std::string>())
        ("n,nevents"        , "N events to loop over"                                                                               , cxxopts::value<int>()->default_value("-1"))
        ("x,event_index"    , "specific event index to process"                                                                     , cxxopts::value<int>()->default_value("-1"))
        ("p,ptbound_mode"   , "Pt bound mode (i.e. 0 = default, 1 = pt~1, 2 = pt~0.95-1.5, 3 = pt~0.5-1.5, 4 = pt~0.5-2.0"          , cxxopts::value<int>()->default_value("0"))
        ("g,pdg_id"         , "The simhit pdgId match option (default = 13)"                                                        , cxxopts::value<int>()->default_value("13"))
        ("v,verbose"        , "Verbose mode"                                                                                        , cxxopts::value<int>()->default_value("0"))
        ("l,run_ineff_study", "Write debug ntuples 0 = MDs, 1 = SGs, 2 = TLs, 3 = TCs"                                              , cxxopts::value<int>()->default_value("-1"))
        ("d,debug"          , "Run debug job. i.e. overrides output option to 'debug.root' and 'recreate's the file.")
        ("c,print_conn"     , "Print module connections")
        ("b,print_boundary" , "Print module boundaries")
        ("r,centroid"       , "Print centroid information")
        ("e,run_eff_study"  , "Run efficiency study")
        ("m,run_mtv_study"  , "Run MTV study")
        ("j,nsplit_jobs"    , "Enable splitting jobs by N blocks (--job_index must be set)", cxxopts::value<int>())
        ("I,job_index"      , "job_index of split jobs (--nsplit_jobs must be set. index starts from 0. i.e. 0, 1, 2, 3, etc...)", cxxopts::value<int>())
        ("h,help"           , "Print help")
        ;

    auto result = options.parse(argc, argv);

    // NOTE: When an option was provided (e.g. -i or --input), then the result.count("<option name>") is more than 0
    // Therefore, the option can be parsed easily by asking the condition if (result.count("<option name>");
    // That's how the several options are parsed below

    //_______________________________________________________________________________
    // --help
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(1);
    }

    //_______________________________________________________________________________
    // --input
    if (result.count("input"))
    {
        ana.input_file_list_tstring = result["input"].as<std::string>();
    }
    else
    {
        std::cout << options.help() << std::endl;
        std::cout << "ERROR: Input list is not provided! Check your arguments" << std::endl;
        exit(1);
    }

    //_______________________________________________________________________________
    // --tree
    if (result.count("tree"))
    {
        ana.input_tree_name = result["tree"].as<std::string>();
    }
    else
    {
        std::cout << options.help() << std::endl;
        std::cout << "ERROR: Input tree name is not provided! Check your arguments" << std::endl;
        exit(1);
    }

    //_______________________________________________________________________________
    // --debug
    if (result.count("debug"))
    {
        ana.output_tfile = new TFile("debug.root", "recreate");
    }
    else
    {
        //_______________________________________________________________________________
        // --output
        if (result.count("output"))
        {
            ana.output_tfile = new TFile(result["output"].as<std::string>().c_str(), "create");
            if (not ana.output_tfile->IsOpen())
            {
                std::cout << options.help() << std::endl;
                std::cout << "ERROR: output already exists! provide new output name or delete old file. OUTPUTFILE=" << result["output"].as<std::string>() << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cout << options.help() << std::endl;
            std::cout << "ERROR: Output file name is not provided! Check your arguments" << std::endl;
            exit(1);
        }
    }

    //_______________________________________________________________________________
    // --print_conn
    if (result.count("print_conn"))
    {
        ana.print_conn = true;
    }
    else
    {
        ana.print_conn = false;
    }

    //_______________________________________________________________________________
    // --print_boundary
    if (result.count("print_boundary"))
    {
        ana.print_boundary = true;
    }
    else
    {
        ana.print_boundary = false;
    }

    //_______________________________________________________________________________
    // --run_eff_study
    if (result.count("run_eff_study"))
    {
        ana.run_eff_study = true;
    }
    else
    {
        ana.run_eff_study = false;
    }

    //_______________________________________________________________________________
    // --run_ineff_study
    if (result.count("run_ineff_study"))
    {
        ana.mode_write_ineff_study_debug_ntuple = result["run_ineff_study"].as<int>();
        if (result["run_ineff_study"].as<int>() >= 0)
        {
            ana.run_ineff_study = true;
        }
        else
        {
            ana.run_ineff_study = false;
        }
    }

    //_______________________________________________________________________________
    // --run_mtv_study
    if (result.count("run_mtv_study"))
    {
        ana.run_mtv_study = true;
    }
    else
    {
        ana.run_mtv_study = false;
    }

    //_______________________________________________________________________________
    // --centroid
    if (result.count("centroid"))
    {
        ana.print_centroid = true;
    }
    else
    {
        ana.print_centroid = false;
    }

    //_______________________________________________________________________________
    // --nevents
    ana.n_events = result["nevents"].as<int>();
    ana.specific_event_index = result["event_index"].as<int>();

    // -1 upto mini-doublet is all-comb
    // -2 upto segment is all-comb
    // -3 upto tracklet is all-comb NOTE: MEMORY WILL BLOW UP FOR HIGH PU
    // -4 upto trackcandidate is all-comb NOTE: MEMORY WILL BLOW UP FOR HIGH PU
    //  0 nothing
    //  1 upto mini-doublet is all-comb
    //  2 upto mini-doublet is default segment is all-comb
    //  3 upto segment is default tracklet is all-comb
    //  4 upto tracklet is default trackcandidate is all-comb
    ana.ptbound_mode = result["ptbound_mode"].as<int>();

    //_______________________________________________________________________________
    // --pdg_id
    ana.pdg_id = result["pdg_id"].as<int>();

    //_______________________________________________________________________________
    // --nsplit_jobs
    if (result.count("nsplit_jobs"))
    {
        ana.nsplit_jobs = result["nsplit_jobs"].as<int>();
        if (ana.nsplit_jobs <= 0)
        {
            std::cout << options.help() << std::endl;
            std::cout << "ERROR: option string --nsplit_jobs" << ana.nsplit_jobs << " has zero or negative value!" << std::endl;
            std::cout << "I am not sure what this means..." << std::endl;
            exit(1);
        }
    }
    else
    {
        ana.nsplit_jobs = -1;
    }

    //_______________________________________________________________________________
    // --job_index
    if (result.count("job_index"))
    {
        ana.job_index = result["job_index"].as<int>();
        if (ana.job_index < 0)
        {
            std::cout << options.help() << std::endl;
            std::cout << "ERROR: option string --job_index" << ana.job_index << " has negative value!" << std::endl;
            std::cout << "I am not sure what this means..." << std::endl;
            exit(1);
        }
    }
    else
    {
        ana.job_index = -1;
    }

    // Sanity check for split jobs (if one is set the other must be set too)
    if (result.count("job_index") or result.count("nsplit_jobs"))
    {
        // If one is not provided then throw error
        if ( not (result.count("job_index") and result.count("nsplit_jobs")))
        {
            std::cout << options.help() << std::endl;
            std::cout << "ERROR: option string --job_index and --nsplit_jobs must be set at the same time!" << std::endl;
            exit(1);
        }
        // If it is set then check for sanity
        else
        {
            if (ana.job_index >= ana.nsplit_jobs)
            {
                std::cout << options.help() << std::endl;
                std::cout << "ERROR: --job_index >= --nsplit_jobs ! This does not make sense..." << std::endl;
                exit(1);
            }
        }
    }

    //_______________________________________________________________________________
    // --verbose
    ana.verbose = result["verbose"].as<int>();

    //
    // Printing out the option settings overview
    //
    std::cout <<  "=========================================================" << std::endl;
    std::cout <<  " Setting of the analysis job based on provided arguments " << std::endl;
    std::cout <<  "---------------------------------------------------------" << std::endl;
    std::cout <<  " ana.input_file_list_tstring: " << ana.input_file_list_tstring <<  std::endl;
    std::cout <<  " ana.output_tfile: " << ana.output_tfile->GetName() <<  std::endl;
    std::cout <<  " ana.n_events: " << ana.n_events <<  std::endl;
    std::cout <<  " ana.run_eff_study: " << ana.run_eff_study <<  std::endl;
    std::cout <<  " ana.run_ineff_study: " << ana.run_ineff_study <<  std::endl;
    std::cout <<  " ana.run_mtv_study: " << ana.run_mtv_study <<  std::endl;
    std::cout <<  " ana.print_centroid: " << ana.print_centroid <<  std::endl;
    std::cout <<  " ana.print_conn: " << ana.print_conn <<  std::endl;
    std::cout <<  " ana.print_boundary: " << ana.print_boundary <<  std::endl;
    std::cout <<  " ana.ptbound_mode: " << ana.ptbound_mode <<  std::endl;
    std::cout <<  " ana.nsplit_jobs: " << ana.nsplit_jobs <<  std::endl;
    std::cout <<  " ana.job_index: " << ana.job_index <<  std::endl;
    std::cout <<  "=========================================================" << std::endl;

    // Consistency check
    if ((ana.run_ineff_study and ana.run_eff_study) or
        (ana.run_ineff_study and ana.run_mtv_study) or
        (ana.run_eff_study and ana.run_mtv_study)
       )
    {
        RooUtil::error("More than one of -e, -l, or -m option is set! Please only set one of them");
    }


//********************************************************************************
//
// 2. Opening input baby files
//
//********************************************************************************

    // Create the TChain that holds the TTree's of the baby ntuples
    ana.events_tchain = RooUtil::FileUtil::createTChain(ana.input_tree_name, ana.input_file_list_tstring);

    // Create a Looper object to loop over input files
    // the "www" object is defined in the wwwtree.h/cc
    // This is an instance which helps read variables in the WWW baby TTree
    // It is a giant wrapper that facilitates reading TBranch values.
    // e.g. if there is a TBranch named "lep_pt" which is a std::vector<float> then, one can access the branch via
    //
    //    std::vector<float> lep_pt = www.lep_pt();
    //
    // and no need for "SetBranchAddress" and declaring variable shenanigans necessary
    // This is a standard thing SNT does pretty much every looper we use
    ana.looper.init(ana.events_tchain, &trk, ana.n_events);

//********************************************************************************
//
// Interlude... notes on RooUtil framework
//
//********************************************************************************

    // Set the cutflow object output file
    ana.cutflow.setTFile(ana.output_tfile);
    // ana.cutflow.addCut("CutWeight", UNITY, UNITY);

    // Create ttree to output to the ana.output_tfile
    ana.output_ttree = new TTree("tree", "tree");

    // Create TTreeX instance that will take care of the interface part of TTree
    ana.tx = new RooUtil::TTreeX(ana.output_ttree);

    // Print cut structure
    ana.cutflow.printCuts();

    // pt_boundaries
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
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0}; // lowpt
    else if (ana.ptbound_mode == 6)
        pt_boundaries = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 3.0, 4.0, 5.0}; // lowpt
        // pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.5, 5.0}; // lowpt

    // List of studies to perform
    std::vector<Study*> studies;
    if (ana.run_eff_study)
    {
        // studies.push_back(new StudySDLEfficiency("efficiency",
        //             StudySDLEfficiency::kStudySDLMDEffBarrel,
        //             StudySDLEfficiency::kStudySDLSGEffBB,
        //             StudySDLEfficiency::kStudySDLTLEffBBBB,
        //             StudySDLEfficiency::kStudySDLTCEffBBBBBB,
        //             pt_boundaries));
        studies.push_back(new StudySDLEfficiencyv2("eff_bbbbbb", StudySDLEfficiencyv2::kStudySDLBBBBBB, pt_boundaries));
        studies.push_back(new StudySDLEfficiencyv2("eff_bbbbbe", StudySDLEfficiencyv2::kStudySDLBBBBBE, pt_boundaries));
        studies.push_back(new StudySDLEfficiencyv2("eff_bbbbee", StudySDLEfficiencyv2::kStudySDLBBBBEE, pt_boundaries));
        studies.push_back(new StudySDLEfficiencyv2("eff_bbbeee", StudySDLEfficiencyv2::kStudySDLBBBEEE, pt_boundaries));
        studies.push_back(new StudySDLEfficiencyv2("eff_bbeeee", StudySDLEfficiencyv2::kStudySDLBBEEEE, pt_boundaries));
        studies.push_back(new StudySDLEfficiencyv2("eff_beeeee", StudySDLEfficiencyv2::kStudySDLBEEEEE, pt_boundaries));
    }
    else if (ana.run_ineff_study)
    {
        switch (ana.mode_write_ineff_study_debug_ntuple)
        {
            case 0: studies.push_back(new StudySDLMiniDoubletDebugNtupleWriter()); break;
            case 1: studies.push_back(new StudySDLSegmentDebugNtupleWriter()); break;
            case 2: studies.push_back(new StudySDLTrackletDebugNtupleWriter()); break;
            case 3: studies.push_back(new StudySDLTrackCandidateDebugNtupleWriter()); break;
            default:
                std::cout << options.help() << std::endl;
                std::cout << "ERROR: ana.mode_write_ineff_study_debug_ntuple not recognized! value = " << ana.mode_write_ineff_study_debug_ntuple << std::endl;
                exit(1);
                break;
        }
    }
    else if (ana.run_mtv_study)
    {
        studies.push_back(new StudyMTVEfficiency("MTVEfficiency",
                    pt_boundaries,
                    ana.pdg_id
                    ));
        studies.push_back(new StudyConditionalHitEfficiency("ConditionalHitEfficiency", pt_boundaries, ana.pdg_id));
    }
    else
    {
        studies.push_back(new StudyOccupancy("studyOccupancy"));
        studies.push_back(new StudyMDOccupancy("studyMDOccupancy"));
        studies.push_back(new StudyLinkedModule("studyLinkedModule"));
        studies.push_back(new StudySegmentOccupancy("studySegmentOccupancy"));
        studies.push_back(new StudyTrackletSelection("studySelTlBB1BB3", StudyTrackletSelection::kStudySelBB1BB3));
        studies.push_back(new StudyTrackletSelection("studySelTlBB2BB4", StudyTrackletSelection::kStudySelBB2BB4));
        studies.push_back(new StudyTrackletSelection("studySelTlBB3BB5", StudyTrackletSelection::kStudySelBB3BB5));
        studies.push_back(new StudyTrackletSelection("studySelTlBB1BE3", StudyTrackletSelection::kStudySelBB1BE3));
        studies.push_back(new StudyTrackletSelection("studySelTlBB2BE4", StudyTrackletSelection::kStudySelBB2BE4));
        studies.push_back(new StudyTrackletSelection("studySelTlBB3BE5", StudyTrackletSelection::kStudySelBB3BE5));
        studies.push_back(new StudyTrackletSelection("studySelTlBB1EE3", StudyTrackletSelection::kStudySelBB1EE3));
        studies.push_back(new StudyTrackletSelection("studySelTlEE1EE3", StudyTrackletSelection::kStudySelEE1EE3));
        // studies.push_back(new StudyTripletSelection("studySelTPBB1BB2", StudyTripletSelection::kStudySelBB1BB2));
        // studies.push_back(new StudyTripletSelection("studySelTPBB2BB3", StudyTripletSelection::kStudySelBB2BB3));
        // studies.push_back(new StudyTripletSelection("studySelTPBB3BB4", StudyTripletSelection::kStudySelBB3BB4));
        // studies.push_back(new StudyTripletSelection("studySelTPBB4BB5", StudyTripletSelection::kStudySelBB4BB5));
        // studies.push_back(new StudyTrackCandidateSelection("studySelTCAll", StudyTrackCandidateSelection::kStudySelAll, pt_boundaries));
    }

    // book the studies
    for (auto& study : studies)
    {
        study->bookStudy();
    }

    // Book Histograms
    ana.cutflow.bookHistograms(ana.histograms); // if just want to book everywhere

    // SDL::endcapGeometry.load("scripts/endcap_orientation_data.txt");
    SDL::endcapGeometry.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/endcap_orientation_data_v2.txt"); // centroid values added to the map
    SDL::tiltedGeometry.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/tilted_orientation_data.txt");
    SDL::moduleConnectionMap.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/module_connection_map_data_10_e0_200_100_pt0p8_2p0_400_pt0p8_2p0_nolossers_dxy35cm_endcaplayer2.txt");

    // // Following maps to compute centroid of each modules
    std::map<unsigned int, std::vector<float>> module_xs;
    std::map<unsigned int, std::vector<float>> module_ys;
    std::map<unsigned int, std::vector<float>> module_zs;

    // connection information
    std::ofstream module_connection_log_output;
    if (ana.print_conn)
        module_connection_log_output.open("conn.txt");

    // module boundary information to be written out in case module boundary info is asked to be printed
    std::ofstream module_boundary_output_info;
    if (ana.print_boundary)
        module_boundary_output_info.open("module_boundary.txt");

    // Write the simhits in a given module to the output TTree
    if (ana.print_boundary)
    {
        ana.tx->createBranch<int>("detId");
        ana.tx->createBranch<int>("subdet");
        ana.tx->createBranch<int>("side");
        ana.tx->createBranch<int>("layer");
        ana.tx->createBranch<int>("rod");
        ana.tx->createBranch<int>("module");
        ana.tx->createBranch<int>("ring");
        ana.tx->createBranch<int>("isPS");
        ana.tx->createBranch<int>("isStrip");
        ana.tx->createBranch<vector<float>>("x");
        ana.tx->createBranch<vector<float>>("y");
        ana.tx->createBranch<vector<float>>("z");
    }

    // Looping input file
    while (ana.looper.nextEvent())
    {

        if (ana.specific_event_index >= 0)
        {
            if ((int)ana.looper.getCurrentEventIndex() != ana.specific_event_index)
                continue;
        }

        // If splitting jobs are requested then determine whether to process the event or not based on remainder
        if (result.count("job_index") and result.count("nsplit_jobs"))
        {
            if (ana.looper.getNEventsProcessed() % ana.nsplit_jobs != (unsigned int) ana.job_index)
                continue;
        }

        if (ana.verbose) std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;

        // *****************************************************
        // Print module -> module connection info from sim track
        // *****************************************************
        if (ana.print_conn)
        {
            // Print the module connection info and do nothing else on the event
            printModuleConnectionInfo(module_connection_log_output);
            continue;
        }

        // *****************************************************
        // Print module boundaries
        // *****************************************************
        if (ana.print_boundary)
        {
            // Print the module connection info and do nothing else on the event
            processModuleBoundaryInfo();
            continue;
        }

        if (ana.print_centroid)
        {
            // Adding hits to modules
            for (unsigned int ihit = 0; ihit < trk.ph2_x().size(); ++ihit)
            {

                // To print the reco hits per module to create a table of centroids of each module
                SDL::Module module = SDL::Module(trk.ph2_detId()[ihit]);
                if (module.isLower())
                {
                    module_xs[trk.ph2_detId()[ihit]].push_back(trk.ph2_x()[ihit]);
                    module_ys[trk.ph2_detId()[ihit]].push_back(trk.ph2_y()[ihit]);
                    module_zs[trk.ph2_detId()[ihit]].push_back(trk.ph2_z()[ihit]);
                }
                else
                {
                    module_xs[module.partnerDetId()].push_back(trk.ph2_x()[ihit]);
                    module_ys[module.partnerDetId()].push_back(trk.ph2_y()[ihit]);
                    module_zs[module.partnerDetId()].push_back(trk.ph2_z()[ihit]);
                }

            }

            continue;

        }

        // <--------------------------
        // <--------------------------
        // <--------------------------
        //
        // ***************************
        // Testing SDL Implementations
        // ***************************
        //

        // *************************************************
        // Reconstructed hits and formation of mini-doublets
        // *************************************************

        // Main instance that will hold modules, hits, minidoublets, etc. (i.e. main data structure)
        SDL::Event event;

        // Each SDL::Event object in simtrkevents will hold single sim-track related hits
        // It will be a vector of tuple of <sim_track_index, SDL::Event*>.
        std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents;

        TStopwatch my_timer;

        // run_eff_study == 0 then run all the reconstruction
        if ((ana.run_eff_study == 0 and ana.run_ineff_study == 0) or ana.run_mtv_study)
        {

            // Adding hits to modules
            for (unsigned int ihit = 0; ihit < trk.ph2_x().size(); ++ihit)
            {


//                if (trk.ph2_subdet()[ihit] != 5)
//                    continue;

                // if (trk.ph2_subdet()[ihit] != 5)
                //     continue;

                // if (trk.sim_pt()[trk.simhit_simTrkIdx()[trk.ph2_simHitIdx()[ihit][0]]] < 2.5)
                //     continue;

                for (auto& isimhit : trk.ph2_simHitIdx()[ihit])
                {
                    int isimtrk = trk.simhit_simTrkIdx()[isimhit];
                    if (trk.sim_pt()[isimtrk] < 2.5)
                        continue;
                }


                // Takes two arguments, SDL::Hit, and detId
                // SDL::Event internally will structure whether we already have the module instance or we need to create a new one.
                event.addHitToModule(
                        // a hit
                        SDL::Hit(trk.ph2_x()[ihit], trk.ph2_y()[ihit], trk.ph2_z()[ihit], ihit),
                        // add to module with "detId"
                        trk.ph2_detId()[ihit]
                        );

            }

            float elapsed = 0;

            // ----------------
            if (ana.verbose != 0) std::cout << "Summary of hits" << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits: " << event.getNumberOfHits() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits in layer 1: " << event.getNumberOfHitsByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits in layer 2: " << event.getNumberOfHitsByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits in layer 3: " << event.getNumberOfHitsByLayerBarrel(2) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits in layer 4: " << event.getNumberOfHitsByLayerBarrel(3) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits in layer 5: " << event.getNumberOfHitsByLayerBarrel(4) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits in layer 6: " << event.getNumberOfHitsByLayerBarrel(5) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 1: " << event.getNumberOfHitsByLayerBarrelUpperModule(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 2: " << event.getNumberOfHitsByLayerBarrelUpperModule(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 3: " << event.getNumberOfHitsByLayerBarrelUpperModule(2) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 4: " << event.getNumberOfHitsByLayerBarrelUpperModule(3) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 5: " << event.getNumberOfHitsByLayerBarrelUpperModule(4) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Hits Upper Module in layer 6: " << event.getNumberOfHitsByLayerBarrelUpperModule(5) << std::endl;
            // ----------------


            // ----------------
            if (ana.verbose != 0) std::cout << "Reco Mini-Doublet start" << std::endl;
            my_timer.Start();
            event.createMiniDoublets();
            // event.createPseudoMiniDoubletsFromAnchorModule(); // Useless.....
            float md_elapsed = my_timer.RealTime();
            if (ana.verbose != 0) std::cout << "Reco Mini-doublet processing time: " << md_elapsed << " secs" << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced: " << event.getNumberOfMiniDoublets() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced layer 1: " << event.getNumberOfMiniDoubletsByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced layer 2: " << event.getNumberOfMiniDoubletsByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced layer 3: " << event.getNumberOfMiniDoubletsByLayerBarrel(2) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced layer 4: " << event.getNumberOfMiniDoubletsByLayerBarrel(3) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced layer 5: " << event.getNumberOfMiniDoubletsByLayerBarrel(4) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets produced layer 6: " << event.getNumberOfMiniDoubletsByLayerBarrel(5) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered: " << event.getNumberOfMiniDoubletCandidates() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered layer 1: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered layer 2: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered layer 3: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(2) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered layer 4: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(3) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered layer 5: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(4) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Mini-doublets considered layer 6: " << event.getNumberOfMiniDoubletCandidatesByLayerBarrel(5) << std::endl;
            // ----------------

            // ----------------
            if (ana.verbose != 0) std::cout << "Reco Segment start" << std::endl;
            my_timer.Start(kFALSE);
            event.createSegmentsWithModuleMap();
            float sg_elapsed = my_timer.RealTime();
            if (ana.verbose != 0) std::cout << "Reco Segment processing time: " << sg_elapsed - md_elapsed << " secs" << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments produced: " << event.getNumberOfSegments() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments produced layer 1-2: " << event.getNumberOfSegmentsByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments produced layer 2-3: " << event.getNumberOfSegmentsByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments produced layer 3-4: " << event.getNumberOfSegmentsByLayerBarrel(2) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments produced layer 4-5: " << event.getNumberOfSegmentsByLayerBarrel(3) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments produced layer 5-6: " << event.getNumberOfSegmentsByLayerBarrel(4) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Segments produced layer 6: " << event.getNumberOfSegmentsByLayerBarrel(5) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments considered: " << event.getNumberOfSegmentCandidates() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments considered layer 1-2: " << event.getNumberOfSegmentCandidatesByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments considered layer 2-3: " << event.getNumberOfSegmentCandidatesByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments considered layer 3-4: " << event.getNumberOfSegmentCandidatesByLayerBarrel(2) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments considered layer 4-5: " << event.getNumberOfSegmentCandidatesByLayerBarrel(3) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Segments considered layer 5-6: " << event.getNumberOfSegmentCandidatesByLayerBarrel(4) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Segments considered layer 6: " << event.getNumberOfSegmentCandidatesByLayerBarrel(5) << std::endl;
            // ----------------

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

            // ----------------
            if (ana.verbose != 0) std::cout << "Reco Triplet start" << std::endl;
            my_timer.Start(kFALSE);
            event.createTriplets();
            float tp_elapsed = my_timer.RealTime();
            if (ana.verbose != 0) std::cout << "Reco Triplet processing time: " << tp_elapsed - sg_elapsed << " secs" << std::endl;
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

            // ----------------
            if (ana.verbose != 0) std::cout << "Reco Tracklet start" << std::endl;
            my_timer.Start(kFALSE);
            // event.createTracklets();
            // event.createTrackletsWithModuleMap();
            event.createTrackletsViaNavigation();
            float tl_elapsed = my_timer.RealTime();
            if (ana.verbose != 0) std::cout << "Reco Tracklet processing time: " << tl_elapsed - tp_elapsed << " secs" << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets produced: " << event.getNumberOfTracklets() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 1-2-3-4: " << event.getNumberOfTrackletsByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 2-3-4-5: " << event.getNumberOfTrackletsByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 3-4-5-6: " << event.getNumberOfTrackletsByLayerBarrel(2) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 4: " << event.getNumberOfTrackletsByLayerBarrel(3) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 5: " << event.getNumberOfTrackletsByLayerBarrel(4) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Tracklets produced layer 6: " << event.getNumberOfTrackletsByLayerBarrel(5) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets considered: " << event.getNumberOfTrackletCandidates() << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 1-2-3-4: " << event.getNumberOfTrackletCandidatesByLayerBarrel(0) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 2-3-4-5: " << event.getNumberOfTrackletCandidatesByLayerBarrel(1) << std::endl;
            if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 3-4-5-6: " << event.getNumberOfTrackletCandidatesByLayerBarrel(2) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 4: " << event.getNumberOfTrackletCandidatesByLayerBarrel(3) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 5: " << event.getNumberOfTrackletCandidatesByLayerBarrel(4) << std::endl;
            // if (ana.verbose != 0) std::cout << "# of Tracklets considered layer 6: " << event.getNumberOfTrackletCandidatesByLayerBarrel(5) << std::endl;
            // ----------------

            // // ----------------
            // multiplicities.clear();
            // total_nmult = 0;
            // for (auto& segmentPtr : event.getLayer(1, SDL::Layer::Barrel).getSegmentPtrs())
            // {
            //     int nmult = 0;
            //     for (auto& tkl1 : segmentPtr->getListOfOutwardTrackletPtrs())
            //     {
            //         for (auto& tkl2 : tkl1->outerSegmentPtr()->getListOfOutwardTrackletPtrs())
            //         {
            //             nmult++;
            //         }
            //     }
            //     multiplicities.push_back(nmult);
            //     total_nmult += nmult;
            // }
            // std::cout <<  " total_nmult: " << total_nmult <<  std::endl;
            // avg_mult = ((float) total_nmult) / ((float) multiplicities.size());
            // std::cout <<  " avg_mult: " << avg_mult <<  std::endl;
            // // ----------------

            // // ----------------
            // multiplicities.clear();
            // total_nmult = 0;
            // for (auto& segmentPtr : event.getLayer(5, SDL::Layer::Barrel).getSegmentPtrs())
            // {
            //     int nmult = 0;
            //     for (auto& tkl1 : segmentPtr->getListOfInwardTrackletPtrs())
            //     {
            //         for (auto& tkl2 : tkl1->innerSegmentPtr()->getListOfInwardTrackletPtrs())
            //         {
            //             nmult++;
            //         }
            //     }
            //     multiplicities.push_back(nmult);
            //     total_nmult += nmult;
            // }
            // std::cout <<  " total_nmult: " << total_nmult <<  std::endl;
            // avg_mult = ((float) total_nmult) / ((float) multiplicities.size());
            // std::cout <<  " avg_mult: " << avg_mult <<  std::endl;
            // // ----------------

            // ----------------
            if (ana.verbose != 0) std::cout << "Reco TrackCandidate start" << std::endl;
            my_timer.Start(kFALSE);
            // event.createTrackCandidatesFromTriplets();
            // event.createTrackCandidates();
            event.createTrackCandidatesFromTracklets();
            float tc_elapsed = my_timer.RealTime();
            if (ana.verbose != 0) std::cout << "Reco TrackCandidate processing time: " << tc_elapsed - tl_elapsed << " secs" << std::endl;
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

            if (ana.verbose != 0) std::cout << "Reco SDL end" << std::endl;

        }
        // If efficiency is to be calculated
        if (ana.run_eff_study or ana.run_ineff_study or ana.run_mtv_study)
        {

            // *******************************************************
            // Formation of mini-doublets "candidates" from sim-tracks
            // *******************************************************

            // Loop over sim-tracks and per sim-track aggregate good hits (i.e. matched with particle ID)
            // and only use those hits, and run mini-doublet reco algorithm on the sim-track-matched-reco-hits
            for (unsigned int isimtrk = 0; isimtrk < trk.sim_q().size(); ++isimtrk)
            {

                // Then select all charged particle
                if (ana.pdg_id == 1)
                {
                    // Select only muon tracks
                    if (abs(trk.sim_q()[isimtrk]) == 0)
                        continue;
                }
                else
                {
                    // Select only muon tracks
                    if (abs(trk.sim_pdgId()[isimtrk]) != ana.pdg_id)
                        continue;
                }

                // // Select hard scatter only
                // if (abs(trk.sim_event()[isimtrk]) != 0)
                //     continue;

                // Select in time only
                if (abs(trk.sim_bunchCrossing()[isimtrk]) != 0)
                    continue;

                // // Select only muon with pt > 1 GeV
                // if (trk.sim_pt()[isimtrk] < 1)
                //     continue;

                // if (not (hasAll12HitsWithNBarrel(isimtrk, 2)))
                //     continue;
                // if (not (hasAll12HitsWithNBarrelUsingModuleMap(isimtrk, 6) or hasAll12HitsWithNBarrelUsingModuleMap(isimtrk, 5)))
                //     continue;
                // if (not hasAll12HitsInBarrel(isimtrk))
                //     continue;

                // event just for this track
                SDL::Event* trackevent = new SDL::Event();

                // loop over the simulated hits
                for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
                {

                    // Retrieve the sim hit idx
                    unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

                    // Select only the hits in the outer tracker
                    // if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
                    //     continue;
                    // if (not (trk.simhit_subdet()[simhitidx] == 5))
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

                            trackevent->addHitToModule(
                                    // a hit
                                    SDL::Hit(trk.ph2_x()[ihit], trk.ph2_y()[ihit], trk.ph2_z()[ihit], ihit),
                                    // add to module with "detId"
                                    trk.ph2_detId()[ihit]
                                    );

                        }

                    }

                }

                if (ana.run_ineff_study)
                {
                    switch (ana.mode_write_ineff_study_debug_ntuple)
                    {
                        case 0: // MD
                            if (ana.verbose != 0) std::cout << "Sim Mini-Doublet start" << std::endl;
                            trackevent->createMiniDoublets(SDL::AllComb_MDAlgo);
                            break;
                        case 1:
                            if (ana.verbose != 0) std::cout << "Sim Mini-Doublet start" << std::endl;
                            trackevent->createMiniDoublets();
                            if (ana.verbose != 0) std::cout << "Sim Segment start" << std::endl;
                            trackevent->createSegmentsWithModuleMap(SDL::AllComb_SGAlgo);
                            break;
                        case 2:
                            if (ana.verbose != 0) std::cout << "Sim Mini-Doublet start" << std::endl;
                            trackevent->createMiniDoublets();
                            if (ana.verbose != 0) std::cout << "Sim Segment start" << std::endl;
                            trackevent->createSegmentsWithModuleMap();
                            if (ana.verbose != 0) std::cout << "Sim Tracklet start" << std::endl;
                            trackevent->createTrackletsWithModuleMap(SDL::AllComb_TLAlgo);
                            break;
                        case 3:
                            if (ana.verbose != 0) std::cout << "Sim Mini-Doublet start" << std::endl;
                            trackevent->createMiniDoublets();
                            if (ana.verbose != 0) std::cout << "Sim Segment start" << std::endl;
                            trackevent->createSegmentsWithModuleMap();
                            if (ana.verbose != 0) std::cout << "Sim Tracklet start" << std::endl;
                            trackevent->createTrackletsWithModuleMap();
                            if (ana.verbose != 0) std::cout << "Sim TrackCandidate start" << std::endl;
                            trackevent->createTrackCandidatesFromTracklets(SDL::AllComb_TCAlgo);
                            break;
                        default:
                            std::cout << options.help() << std::endl;
                            std::cout << "ERROR: ana.mode_write_ineff_study_debug_ntuple not recognized! value = " << ana.mode_write_ineff_study_debug_ntuple << std::endl;
                            exit(1);
                            break;
                    }
                }
                else
                {
                    if (ana.verbose != 0) std::cout << "Sim Mini-Doublet start" << std::endl;
                    trackevent->createMiniDoublets();
                    if (ana.verbose != 0) std::cout << "Sim Segment start" << std::endl;
                    trackevent->createSegmentsWithModuleMap();
                    if (ana.verbose != 0) std::cout << "Sim Tracklet start" << std::endl;
                    // trackevent->createTrackletsWithModuleMap();
                    trackevent->createTrackletsWithModuleMap();
                    if (ana.verbose != 0) std::cout << "Sim Triplet start" << std::endl;
                    trackevent->createTriplets();
                    if (ana.verbose != 0) std::cout << "Sim TrackCandidate start" << std::endl;
                    trackevent->createTrackCandidatesFromTracklets();
                    if (ana.verbose != 0) std::cout << "Sim SDL end" << std::endl;
                }


                // Push to the vector so we have a data-base of per hit, mini-doublets
                simtrkevents.push_back(std::make_tuple(isimtrk, trackevent));

            }


        }

        // Print content in the event
        // (SDL::cout is a modified version of std::cout where each line is prefixed by SDL::)
        // if (ana.looper.getCurrentEventIndex() < 2) // Print for the first 10 events only
        // {
        //     SDL::cout << event;
        // }


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

    if (ana.print_boundary)
    {
        for (auto mb : ana.moduleBoundaries)
        {
            int detid = mb.first;
            std::array<float, 6> bounds = mb.second;
            float zmin = bounds[0];
            float zmax = bounds[1];
            float phimin = bounds[2];
            float phimax = bounds[3];
            float rtmin = bounds[4];
            float rtmax = bounds[5];
            float zlen = abs(zmax - zmin);
            float philen = abs(phimax - phimin);
            float rtlen = abs(rtmax - rtmin);
            module_boundary_output_info <<  " detid: " << detid <<  " zlen: " << zlen <<  " philen: " << philen <<  " rtlen: " << rtlen <<  " zmax: " << zmax <<  " zmin: " << zmin <<  " phimax: " << phimax <<  " phimin: " << phimin <<  " rtmax: " << rtmax <<  " rtmin: " << rtmin <<  std::endl;
            // std::cout <<  " detid: " << detid <<  " zlen: " << zlen <<  " philen: " << philen <<  " rtlen: " << rtlen <<  " zmax: " << zmax <<  " zmin: " << zmin <<  " phimax: " << phimax <<  " phimin: " << phimin <<  " rtmax: " << rtmax <<  " rtmin: " << rtmin <<  std::endl;
        }
        for (auto ms : ana.moduleSimHits)
        {
            int detid = ms.first;
            std::vector<std::vector<float>>& hit_coords = ms.second;
            SDL::Module module(detid);
            ana.tx->setBranch<int>("detId", detid);
            ana.tx->setBranch<int>("subdet", module.subdet());
            ana.tx->setBranch<int>("side", module.side());
            ana.tx->setBranch<int>("layer", module.layer());
            ana.tx->setBranch<int>("rod", module.rod());
            ana.tx->setBranch<int>("module", module.module());
            ana.tx->setBranch<int>("ring", module.ring());
            ana.tx->setBranch<int>("isPS", module.moduleType() == SDL::Module::PS);
            ana.tx->setBranch<int>("isStrip", module.moduleLayerType() == SDL::Module::Strip);
            for (auto& hit_coord : hit_coords)
            {
                ana.tx->pushbackToBranch<float>("x", hit_coord[0]);
                ana.tx->pushbackToBranch<float>("y", hit_coord[1]);
                ana.tx->pushbackToBranch<float>("z", hit_coord[2]);
            }
            ana.tx->fill();
            ana.tx->clear();
        }
    }

    if (ana.print_centroid)
    {

        std::map<unsigned int, float> module_center_x;
        std::map<unsigned int, float> module_center_y;
        std::map<unsigned int, float> module_center_z;

        for (auto& key : module_xs)
        {
            unsigned int detid = key.first;
            std::vector<float> fs = key.second;
            float sumf = 0;
            for (auto& f : fs)
                sumf += f;
            if (fs.size() > 0)
                module_center_x[detid] = sumf / fs.size();
        }

        for (auto& key : module_ys)
        {
            unsigned int detid = key.first;
            std::vector<float> fs = key.second;
            float sumf = 0;
            for (auto& f : fs)
                sumf += f;
            if (fs.size() > 0)
                module_center_y[detid] = sumf / fs.size();
        }

        for (auto& key : module_zs)
        {
            unsigned int detid = key.first;
            std::vector<float> fs = key.second;
            float sumf = 0;
            for (auto& f : fs)
                sumf += f;
            if (fs.size() > 0)
                module_center_z[detid] = sumf / fs.size();
        }

        for (auto& key : module_center_x)
        {
            unsigned int detid = key.first;
            std::cout <<  " detid: " << detid <<  " module_center_x[detid]: " << module_center_x[detid] <<  " module_center_y[detid]: " << module_center_y[detid] <<  " module_center_z[detid]: " << module_center_z[detid] <<  std::endl;
        }

    }

    // Writing output file
    ana.cutflow.saveOutput();

    // Writing ttree output to file
    ana.output_ttree->Write();

    // The below can be sometimes crucial
    delete ana.output_tfile;
}

int layer(int lay, int det)
{
    //try to restore the SLHC indexing:
    // barrel: 5-10 OT
    // endcap: 11-15 OT disks
    // IT is not handled: "4" layers are filled from seed layer counting (no geometry ref)
    if (det == 4) return 10 + lay; //OT endcap
    if (det == 5) return 4 + lay; //OT barrel
    return lay;
}

void processModuleBoundaryInfo()
{
    int iidOld = -1;

    constexpr int nBoundPoints = 6;
    constexpr int nLayersB = 10; //barrel layers
    std::array<float, nBoundPoints>* cbound;

    int dpop = 0;
    int* cpop = nullptr;
    std::array<float, nBoundPoints> dbound {999, -999, 999, -999, 999, -999}; //zmin, zmax, phimin, phimax, rtmin, rtmax

    std::array<int, nLayersB + 1> hitsBarrelLayer {};
    auto nPh2 = trk.ph2_isBarrel().size();
    for (auto iph2 = 0U; iph2 < nPh2; ++iph2)
    {
        bool isBarrel = trk.ph2_isBarrel()[iph2];
        int lay = layer(trk.ph2_layer()[iph2], trk.ph2_detId()[iph2]);
        hitsBarrelLayer[lay]++;
        int iid = trk.ph2_detId()[iph2];

        SDL::Module module = SDL::Module(iid);

        // if (not (module.layer() == 1 and module.subdet() == SDL::Module::Barrel and module.side() == SDL::Module::Center and module.rod() == 1 and module.module() == 1))
        //     continue;

        if (iidOld != iid)
        {
            iidOld = iid;
            if (ana.modulePopulation.find(iid) == ana.modulePopulation.end())
            {
                ana.modulePopulation[iid] = dpop;
                ana.moduleBoundaries[iid] = dbound;
            }
            cpop = &ana.modulePopulation[iid];
            cbound = &ana.moduleBoundaries[iid];
        }

        //look at all associated simhits, ignoring double-counting
        auto const& ph2shV = trk.ph2_simHitIdx()[iph2];
        auto nPh2sh = ph2shV.size();
        for (auto iph2sh = 0U; iph2sh < nPh2sh; ++iph2sh)
        {
            const float x = trk.simhit_x()[ph2shV[iph2sh]];
            const float y = trk.simhit_y()[ph2shV[iph2sh]];
            const float rt = sqrt(x * x + y * y);
            // const float z = isBarrel ? trk.simhit_z()[ph2shV[iph2sh]] : rt;
            const float z = trk.simhit_z()[ph2shV[iph2sh]];

            if (z == 0)
                continue;

            float phi = atan2(y, x);

            std::vector<float> simhit_coords = {x, y, z};
            ana.moduleSimHits[iid].push_back(simhit_coords);

            if ((*cbound)[0] > z)
                (*cbound)[0] = z;
            if ((*cbound)[1] < z)
                (*cbound)[1] = z;

            if (*cpop == 0)
            {
                (*cbound)[2] = phi;
                (*cbound)[3] = phi;
                (*cbound)[4] = rt;
                (*cbound)[5] = rt;
            }
            else
            {
                if (sin((*cbound)[2] - phi) > 0)
                    (*cbound)[2] = phi;
                if (sin((*cbound)[3] - phi) < 0)
                    (*cbound)[3] = phi;

                if ((*cbound)[4] > rt)
                    (*cbound)[4] = rt;
                if ((*cbound)[5] < rt)
                    (*cbound)[5] = rt;
            }

            (*cpop)++;
        }
    }

}

void printModuleConnectionInfo(std::ofstream& ostrm)
{
    // *****************************************************
    // Print module -> module connection info from sim track
    // *****************************************************

    // Loop over sim-tracks and per sim-track aggregate good hits (i.e. matched with particle ID)
    // and only use those hits, and run mini-doublet reco algorithm on the sim-track-matched-reco-hits
    for (unsigned int isimtrk = 0; isimtrk < trk.sim_q().size(); ++isimtrk)
    {

        // Select only muon tracks
        if (abs(trk.sim_pdgId()[isimtrk]) != 13)
            continue;

        if (trk.sim_pt()[isimtrk] < 1)
            continue;

        // if (trk.sim_pt()[isimtrk] > 1.05)
        //     continue;

        // // Select only 1 cm from center tracks
        // if (fabs(trk.sim_pca_dz()[isimtrk]) > 0.1)
        //     continue;

        std::vector<int> layers;
        std::vector<int> subdets;
        std::vector<int> detids;
        std::vector<float> ps;

        // std::cout <<  " trk.event(): " << trk.event() <<  std::endl;
        // std::cout <<  " isimtrk: " << isimtrk <<  std::endl;

        // loop over the simulated hits
        for (auto& simhitidx : trk.sim_simHitIdx()[isimtrk])
        {

            // Select only the hits in the outer tracker
            if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
                continue;

            int simhit_particle = trk.simhit_particle()[simhitidx];
            int detid = trk.simhit_detId()[simhitidx];
            int layer = trk.simhit_layer()[simhitidx];
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
            int subdet = trk.simhit_subdet()[simhitidx];
            int trkidx = trk.simhit_simTrkIdx()[simhitidx];
            SDL::Module module = SDL::Module(detid);

            // Select only the sim hits that is matched to the muon
            if (abs(simhit_particle) != 13)
                continue;

            // std::cout <<  " layer: " << layer <<  " subdet: " << subdet <<  " simhit_particle: " << simhit_particle <<  " x: " << x <<  " y: " << y <<  " z: " << z <<  " r3: " << r3 <<  " px: " << px <<  " py: " << py <<  " pz: " << pz <<  " p: " << p <<  " angle: " << angle <<  " trkidx: " << trkidx <<  " module.isLower(): " << module.isLower() <<  " detid: " << detid <<  std::endl;

            // if (abs(angle) > 1.6)
            //     break;

            if (ps.size() > 0)
            {
                float loss = fabs(ps.back() - p) / ps.back();
                if (loss > 0.35)
                    break;
            }

            // Access hits on the S side of the PS modules in the endcaps and get three numbers, (detId, x, y)
            if (module.isLower())
            {
                layers.push_back(layer);
                subdets.push_back(subdet);
                detids.push_back(detid);
                ps.push_back(p);
            }

        }

        // std::cout << "momentums: ";
        // for (auto& p : ps)
        // {
        //     std::cout << p << " ";
        // }
        // std::cout << std::endl;

        if (layers.size() > 0)
        {

            // std::cout <<  " trk.sim_pt()[isimtrk]: " << trk.sim_pt()[isimtrk] <<  " trk.sim_phi()[isimtrk]: " << trk.sim_phi()[isimtrk] <<  " trk.sim_eta()[isimtrk]: " << trk.sim_eta()[isimtrk] <<  " trk.sim_pca_dz()[isimtrk]: " << trk.sim_pca_dz()[isimtrk] <<  " trk.sim_pca_dxy()[isimtrk]: " << trk.sim_pca_dxy()[isimtrk] <<  " trk.sim_pca_lambda()[isimtrk]: " << trk.sim_pca_lambda()[isimtrk] <<  " trk.sim_pca_cotTheta()[isimtrk]: " << trk.sim_pca_cotTheta()[isimtrk] <<  " trk.sim_pca_phi()[isimtrk]: " << trk.sim_pca_phi()[isimtrk] <<  std::endl;
            ostrm << "moduleconnection: ";
            for (unsigned int i = 0; i < layers.size(); ++i)
            {
                ostrm << "(" << layers[i] << "," << subdets[i] << "," << detids[i] << "," << SDL::Module(detids[i]).partnerDetId() << ");";
            }
            ostrm << std::endl;
            // ostrm << trk.event() << " " << isimtrk << " ";
            // for (unsigned int i = 0; i < layers.size(); ++i)
            // {
            //     ostrm << detids[i] << "," << SDL::Module(detids[i]).partnerDetId();
            //     if (i != layers.size() - 1)
            //         ostrm << ",";
            // }
            // ostrm << std::endl;

            // // std::cout <<  " trk.sim_pt()[isimtrk]: " << trk.sim_pt()[isimtrk] <<  " trk.sim_phi()[isimtrk]: " << trk.sim_phi()[isimtrk] <<  " trk.sim_eta()[isimtrk]: " << trk.sim_eta()[isimtrk] <<  " trk.sim_pca_dz()[isimtrk]: " << trk.sim_pca_dz()[isimtrk] <<  " trk.sim_pca_dxy()[isimtrk]: " << trk.sim_pca_dxy()[isimtrk] <<  " trk.sim_pca_lambda()[isimtrk]: " << trk.sim_pca_lambda()[isimtrk] <<  " trk.sim_pca_cotTheta()[isimtrk]: " << trk.sim_pca_cotTheta()[isimtrk] <<  " trk.sim_pca_phi()[isimtrk]: " << trk.sim_pca_phi()[isimtrk] <<  std::endl;
            // std::cout << "moduleconnection: ";
            // for (unsigned int i = 0; i < layers.size(); ++i)
            // {
            //     std::cout << "(" << layers[i] << "," << subdets[i] << "," << detids[i] << "," << SDL::Module(detids[i]).partnerDetId() << ");";
            // }
            // std::cout << std::endl;
            // std::cout << trk.event() << " " << isimtrk << " ";
            // for (unsigned int i = 0; i < layers.size(); ++i)
            // {
            //     std::cout << detids[i] << "," << SDL::Module(detids[i]).partnerDetId();
            //     if (i != layers.size() - 1)
            //         std::cout << ",";
            // }
            // std::cout << std::endl;

        }

    }
}

// bool isSDLDenominator(unsigned int isimtrk)
// {

//     int pdgid = trk.sim_pdgId()[isimtrk];
//     int abspdgid = trk.sim_pdgId()[isimtrk];

//     // Muon SDL denominator object definition
//     if (pdgid == 13)
//     {
//     }

// }

// bool hasAll12HitsInBarrel(unsigned int isimtrk)
// {

//     // Select only tracks that left all 12 hits in the barrel
//     std::array<std::vector<SDL::Module>, 6> layers_modules;

//     std::vector<float> ps;

//     for (unsigned int ith_hit = 0; ith_hit < trk.sim_simHitIdx()[isimtrk].size(); ++ith_hit)
//     {

//         // Retrieve the sim hit idx
//         unsigned int simhitidx = trk.sim_simHitIdx()[isimtrk][ith_hit];

//         // Select only the hits in the outer tracker
//         if (not (trk.simhit_subdet()[simhitidx] == 4 or trk.simhit_subdet()[simhitidx] == 5))
//             continue;

//         if (isMuonCurlingHit(isimtrk, ith_hit))
//             break;

//         // list of reco hit matched to this sim hit
//         for (unsigned int irecohit = 0; irecohit < trk.simhit_hitIdx()[simhitidx].size(); ++irecohit)
//         {
//             // Get the recohit type
//             int recohittype = trk.simhit_hitType()[simhitidx][irecohit];

//             // Consider only ph2 hits (i.e. outer tracker hits)
//             if (recohittype == 4)
//             {

//                 int ihit = trk.simhit_hitIdx()[simhitidx][irecohit];

//                 if (trk.ph2_subdet()[ihit] == 5)
//                 {
//                     layers_modules[trk.ph2_layer()[ihit] - 1].push_back(SDL::Module(trk.ph2_detId()[ihit]));
//                 }

//             }

//         }

//     }

//     std::array<bool, 6> has_good_pair_by_layer;
//     has_good_pair_by_layer[0] = false;
//     has_good_pair_by_layer[1] = false;
//     has_good_pair_by_layer[2] = false;
//     has_good_pair_by_layer[3] = false;
//     has_good_pair_by_layer[4] = false;
//     has_good_pair_by_layer[5] = false;

//     bool has_good_pair_all_layer = true;

//     for (int ilayer = 0; ilayer < 6; ++ilayer)
//     {

//         bool has_good_pair = false;

//         for (unsigned imod = 0; imod < layers_modules[ilayer].size(); ++imod)
//         {
//             for (unsigned jmod = imod + 1; jmod < layers_modules[ilayer].size(); ++jmod)
//             {
//                 if (layers_modules[ilayer][imod].partnerDetId() == layers_modules[ilayer][jmod].detId())
//                     has_good_pair = true;
//             }
//         }

//         if (not has_good_pair)
//             has_good_pair_all_layer = false;

//         has_good_pair_by_layer[ilayer] = has_good_pair;

//     }

//     float pt = trk.sim_pt()[isimtrk];
//     float eta = trk.sim_eta()[isimtrk];

//     // if (abs((trk.sim_pt()[isimtrk] - 0.71710)) < 0.00001)
//     // {
//     //     std::cout << std::endl;
//     //     std::cout <<  " has_good_pair_by_layer[0]: " << has_good_pair_by_layer[0] <<  " has_good_pair_by_layer[1]: " << has_good_pair_by_layer[1] <<  " has_good_pair_by_layer[2]: " << has_good_pair_by_layer[2] <<  " has_good_pair_by_layer[3]: " << has_good_pair_by_layer[3] <<  " has_good_pair_by_layer[4]: " << has_good_pair_by_layer[4] <<  " has_good_pair_by_layer[5]: " << has_good_pair_by_layer[5] <<  " pt: " << pt <<  " eta: " << eta <<  std::endl;
//     // }

//     return has_good_pair_all_layer;

// }

