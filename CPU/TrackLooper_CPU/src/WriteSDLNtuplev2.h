#ifndef WriteSDLNtuplev2_h
#define WriteSDLNtuplev2_h

#include "SDL/Event.h"

#include "Study.h"

#include <vector>
#include <tuple>

#include "TString.h"
#include "trktree.h"
#include "constants.h"

#include "AnalysisConfig.h"

#include "trkCore.h"

#include "SDLMath.h"

#include "SDL/MathUtil.h"

#include <cppitertools/itertools.hpp>

class WriteSDLNtuplev2 : public Study
{

public:
    const char* studyname;

    WriteSDLNtuplev2(const char* studyName);

    virtual void bookStudy();
    virtual void doStudy(SDL::Event& recoevent, std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents);

    void createHitsSimHitsSimTracksBranches();
    void createPixelSeedBranches();
    void createMiniDoubletBranches();
    void createSegmentBranches();
    void createPixelSegmentBranches();
    void createTripletBranches();
    void createQuadrupletBranches();
    void createPixelQuadrupletBranches();
    void createTrackCandidateBranches();

    void setHitsSimHitsSimTracksBranches();
    void setPixelSeedBranches();
    void setMiniDoubletBranches(SDL::Event& recoevent);
    void setSegmentBranches(SDL::Event& recoevent);
    void setPixelSegmentBranches(SDL::Event& recoevent);
    void setTripletBranches(SDL::Event& recoevent);
    void setQuadrupletBranches(SDL::Event& recoevent);
    void setPixelQuadrupletBranches(SDL::Event& recoevent);
    void setTrackCandidateBranches(SDL::Event& recoevent);

};

#endif
