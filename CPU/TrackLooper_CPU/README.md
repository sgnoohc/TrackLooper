# TrackLooper

## Install

    git clone --recurse-submodules git@github.com:sgnoohc/TrackLooper.git
    cd TrackLooper
    source rooutil/thisrooutil.sh # set rooutil
    source rooutil/root.sh # set ROOT
    make clean;
    make -j;

## Obtaining detector geometry via the centroids

    ./bin/sdl -i "/nfs-7/userdata/phchang/trackingNtuple/trackingNtuple_*.root" -t trackingNtuple/tree -d -m 1 # creates data/all_sim_hits.root
    python scripts/module_bound_fit.py # creates data/phase2.txt

## Building map

Folloing command creates ```data/conn.txt```.  Afterwards, copy it to a ```conn.py``` format.

    ./bin/sdl -m 0 -i /nfs-7/userdata/phchang/trackingNtuple/trackingNtuple_400_pt_0p8_2p0.root,/nfs-7/userdata/phchang/trackingNtuple/trackingNtuple_10MuGun.root,/nfs-7/userdata/phchang/trackingNtuple/trackingNtuple_hundred_pt_0p8_2p0.root

Then run the following command. (Watch out for the import statements)

    python SDLpython/build_module_connection.py

## Running SDL on CPU and creating "SDL tracking ntuple" and creating standard efficiency plots in one go

Following are list of available sample choices

    SAMPLEID = 1  : pt 0.5 to 2.0 hundred-mu-gun sample
               2  : e 0 to 200 ten-mu-gun sample
               3  : pu200 ttbar (252 events without sim info on pileup tracks) with sdl option --pdg_id == 13
               4  : pu200 ttbar (252 events without sim info on pileup tracks) with sdl option --pdg_id == 211
               5  : pu200 ttbar (252 events without sim info on pileup tracks) with sdl option --pdg_id == 11
               6  : DO NOT USE: displaced SUSY stuff (not sure what)
               7  : pu200 ctau100 of some SUSY? with sdl option --pdg_id == 13
               8  : pu200 ctau100 of some SUSY? with sdl option --pdg_id == 211
               9  : pu200 ctau100 of some SUSY? with sdl option --pdg_id == 11
               10 : pt 0.5 to 3.0 hundred-mu-gun sample
               11 : pt 0.5 to 5.0 ten-mu-gun sample
               12 : pu200 ttbar (500 evt) with sdl option --pdg_id == 13
               13 : pu200 ttbar (500 evt) with sdl option --pdg_id == 211
               14 : pu200 ttbar (500 evt) with sdl option --pdg_id == 11
               15 : pu200 ttbar (500 evt) with sdl option --pdg_id == 0 (0 means all charged particle)
               17 : 5 cm 'cube' with pt 0.5 to 50 ten-mu-gun
               18 : 50 cm 'cube' with pt 0.5 to 50 ten-mu-gun
               19 : ttbar (pu 200?)
               20 : 100 pion gun pt 0.5 to 50 GeV

To run the 100 mu gun 0.5 to 2.0 GeV sample setup

    sh shell/sdlv2.sh -s 1 -t MyTag -d "JobTag"  -w

## Running standard efficiency plots only

The setup assumes one has obtained a "SDL tracking ntuple"

    cd efficiency/
    make -j
    sh run.sh your_SDL_tracking_ntuple.root ${PTBOUNDARYOPTION} # This creates efficiency.root file with histograms
    sh plot.sh ${SAMPLETAG} $(git rev-parse --short HEAD) > plot.log 2>&1 # This actually creates the plots
    cp index.html plots/ # helper index.html for easy presentation

where ${PTBOUNDARYOPTION} is an integer option with one of the following to set the pt efficiency plot bin boundaries

     0:
        pt_boundaries = {0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 3.0, 5.0, 10, 15., 25, 50};
     1:
        pt_boundaries = {0.988, 0.99, 0.992, 0.994, 0.996, 0.998, 1.0, 1.002, 1.004, 1.006, 1.008, 1.01, 1.012}; // lowpt
     2:
        pt_boundaries = {0.955, 0.96, 0.965, 0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 1.00, 1.005, 1.01, 1.015, 1.02, 1.025, 1.03, 1.035, 1.04, 1.045, 1.05}; // pt 0p95 1p05
     3:
        pt_boundaries = {0.5, 0.6, 0.7, 0.8, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.2, 1.5}; // lowpt
     4:
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0}; // lowpt
     5:
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.24, 1.28, 1.32, 1.36, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0}; // lowpt
     6:
        pt_boundaries = {0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 2.0, 3.0, 4.0, 5.0}; // lowpt
     7:
        pt_boundaries = {0.5, 0.52, 0.54, 0.56, 0.58, 0.6, 0.62, 0.64, 0.66, 0.68, 0.7, 0.72, 0.74, 0.76, 0.78, 0.8, 0.82, 0.84, 0.86, 0.88, 0.9, 0.92, 0.94, 0.96, 0.98, 1.0, 1.02, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.22, 1.24, 1.26, 1.28, 1.3, 1.32, 1.34, 1.36, 1.38, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20, 25, 30, 35, 40, 45, 50}; // lowpt
     8:
        pt_boundaries = {0, 0.5, 1.0, 3.0, 5.0, 10, 15., 25, 50};

and ${SAMPLETAG} can be something like ```pt0p5_2p0```, or ```pt0p5_2p0_5cm_cube``` to indicate some information about your sample.  
Keep it short! Otherwise your plot will contain a very long string that will look ugly.

