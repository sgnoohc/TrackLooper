#include "Event.h"

SDL::Event::Event() : logLevel_(SDL::Log_Nothing), pixelLayer_(0, 0)
{
    createLayers();
    n_hits_by_layer_barrel_.fill(0);
    n_hits_by_layer_endcap_.fill(0);
    n_hits_by_layer_barrel_upper_.fill(0);
    n_hits_by_layer_endcap_upper_.fill(0);
    n_miniDoublet_candidates_by_layer_barrel_.fill(0);
    n_segment_candidates_by_layer_barrel_.fill(0);
    n_tracklet_candidates_by_layer_barrel_.fill(0);
    n_triplet_candidates_by_layer_barrel_.fill(0);
    n_trackcandidate_candidates_by_layer_barrel_.fill(0);
    n_miniDoublet_by_layer_barrel_.fill(0);
    n_segment_by_layer_barrel_.fill(0);
    n_tracklet_by_layer_barrel_.fill(0);
    n_triplet_by_layer_barrel_.fill(0);
    n_trackcandidate_by_layer_barrel_.fill(0);
    n_miniDoublet_candidates_by_layer_endcap_.fill(0);
    n_segment_candidates_by_layer_endcap_.fill(0);
    n_tracklet_candidates_by_layer_endcap_.fill(0);
    n_triplet_candidates_by_layer_endcap_.fill(0);
    n_trackcandidate_candidates_by_layer_endcap_.fill(0);
    n_miniDoublet_by_layer_endcap_.fill(0);
    n_segment_by_layer_endcap_.fill(0);
    n_tracklet_by_layer_endcap_.fill(0);
    n_triplet_by_layer_endcap_.fill(0);
    n_trackcandidate_by_layer_endcap_.fill(0);

}

SDL::Event::~Event()
{
}

bool SDL::Event::hasModule(unsigned int detId)
{
    if (modulesMapByDetId_.find(detId) == modulesMapByDetId_.end())
    {
        return false;
    }
    else
    {
        return true;
    }
}

void SDL::Event::setLogLevel(SDL::LogLevel logLevel)
{
    logLevel_ = logLevel;
}

SDL::Module& SDL::Event::getModule(unsigned int detId)
{
    // using std::map::emplace
    std::pair<std::map<unsigned int, Module>::iterator, bool> emplace_result = modulesMapByDetId_.emplace(detId, detId);

    // Retreive the module
    auto& inserted_or_existing = (*(emplace_result.first)).second;

    // If new was inserted, then insert to modulePtrs_ pointer list
    if (emplace_result.second) // if true, new was inserted
    {

        // The pointer to be added
        Module* module_ptr = &((*(emplace_result.first)).second);

        // Add the module pointer to the list of modules
        modulePtrs_.push_back(module_ptr);

        // If the module is lower module then add to list of lower modules
        if (module_ptr->isLower())
            lowerModulePtrs_.push_back(module_ptr);
    }

    return inserted_or_existing;
}

const std::vector<SDL::Module*> SDL::Event::getModulePtrs() const
{
    return modulePtrs_;
}

const std::vector<SDL::Module*> SDL::Event::getLowerModulePtrs() const
{
    return lowerModulePtrs_;
}

void SDL::Event::createLayers()
{
    // Create barrel layers
    for (int ilayer = SDL::Layer::BarrelLayer0; ilayer < SDL::Layer::nBarrelLayer; ++ilayer)
    {
        barrelLayers_[ilayer] = SDL::Layer(ilayer, SDL::Layer::Barrel);
        layerPtrs_.push_back(&(barrelLayers_[ilayer]));
    }

    // Create endcap layers
    for (int ilayer = SDL::Layer::EndcapLayer0; ilayer < SDL::Layer::nEndcapLayer; ++ilayer)
    {
        endcapLayers_[ilayer] = SDL::Layer(ilayer, SDL::Layer::Endcap);
        layerPtrs_.push_back(&(endcapLayers_[ilayer]));
    }
}

SDL::Layer& SDL::Event::getLayer(int ilayer, SDL::Layer::SubDet subdet)
{
    if (subdet == SDL::Layer::Barrel)
        return barrelLayers_[ilayer];
    else // if (subdet == SDL::Layer::Endcap)
        return endcapLayers_[ilayer];
}

SDL::Layer& SDL::Event::getPixelLayer()
{
    return pixelLayer_;
}

const std::vector<SDL::Layer*> SDL::Event::getLayerPtrs() const
{
    return layerPtrs_;
}

void SDL::Event::addHitToModule(SDL::Hit hit, unsigned int detId)
{
    // Add to global list of hits, where it will hold the object's instance
    hits_.push_back(hit);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addHit(&(hits_.back()));

    // Count number of hits in the event
    incrementNumberOfHits(getModule(detId));

    // If the hit is 2S in the endcap then the hit boundary needs to be set
    if (getModule(detId).subdet() == SDL::Module::Endcap and getModule(detId).moduleType() == SDL::Module::TwoS)
    {
        hits_2s_edges_.push_back(GeometryUtil::stripHighEdgeHit(hits_.back()));
        hits_.back().setHitHighEdgePtr(&(hits_2s_edges_.back()));
        hits_2s_edges_.push_back(GeometryUtil::stripLowEdgeHit(hits_.back()));
        hits_.back().setHitLowEdgePtr(&(hits_2s_edges_.back()));
    }
}

void SDL::Event::addMiniDoubletToEvent(SDL::MiniDoublet md, unsigned int detId, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of mini doublets, where it will hold the object's instance
    miniDoublets_.push_back(md);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addMiniDoublet(&(miniDoublets_.back()));

    // And get the layer
    getLayer(layerIdx, subdet).addMiniDoublet(&(miniDoublets_.back()));
}

[[deprecated("SDL:: addMiniDoubletToLowerModule() is deprecated. Use addMiniDoubletToEvent")]]
void SDL::Event::addMiniDoubletToLowerModule(SDL::MiniDoublet md, unsigned int detId)
{
    // Add to global list of mini doublets, where it will hold the object's instance
    miniDoublets_.push_back(md);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addMiniDoublet(&(miniDoublets_.back()));
}

void SDL::Event::addSegmentToEvent(SDL::Segment sg, unsigned int detId, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of segments, where it will hold the object's instance
    segments_.push_back(sg);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addSegment(&(segments_.back()));

    // And get the layer andd the segment to it
    getLayer(layerIdx, subdet).addSegment(&(segments_.back()));

    // Link segments to mini-doublets
    segments_.back().addSelfPtrToMiniDoublets();

}

void SDL::Event::addTrackletToEvent(SDL::Tracklet tl, unsigned int detId, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of segments, where it will hold the object's instance
    tracklets_.push_back(tl);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addTracklet(&(tracklets_.back()));

    // And get the layer andd the segment to it
    if (layerIdx == 0)
        getPixelLayer().addTracklet(&(tracklets_.back()));
    else
        getLayer(layerIdx, subdet).addTracklet(&(tracklets_.back()));

    // Link segments to mini-doublets
    tracklets_.back().addSelfPtrToSegments();

}

void SDL::Event::addTripletToEvent(SDL::Triplet tp, unsigned int detId, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of segments, where it will hold the object's instance
    triplets_.push_back(tp);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addTriplet(&(triplets_.back()));

    // And get the layer andd the triplet to it
    getLayer(layerIdx, subdet).addTriplet(&(triplets_.back()));
}

[[deprecated("SDL:: addSegmentToLowerModule() is deprecated. Use addSegmentToEvent")]]
void SDL::Event::addSegmentToLowerModule(SDL::Segment sg, unsigned int detId)
{
    // Add to global list of segments, where it will hold the object's instance
    segments_.push_back(sg);

    // And get the module (if not exists, then create), and add the address to Module.hits_
    getModule(detId).addSegment(&(segments_.back()));
}

[[deprecated("SDL:: addSegmentToLowerLayer() is deprecated. Use addSegmentToEvent")]]
void SDL::Event::addSegmentToLowerLayer(SDL::Segment sg, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of segments, where it will hold the object's instance
    segments_.push_back(sg);

    // And get the layer
    getLayer(layerIdx, subdet).addSegment(&(segments_.back()));
}

void SDL::Event::addTrackletToLowerLayer(SDL::Tracklet tl, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of tracklets, where it will hold the object's instance
    tracklets_.push_back(tl);

    // And get the layer
    getLayer(layerIdx, subdet).addTracklet(&(tracklets_.back()));
}

void SDL::Event::addTrackCandidateToLowerLayer(SDL::TrackCandidate tc, int layerIdx, SDL::Layer::SubDet subdet)
{
    // Add to global list of trackcandidates, where it will hold the object's instance
    trackcandidates_.push_back(tc);

    // And get the layer
    if (layerIdx == 0)
        getPixelLayer().addTrackCandidate(&(trackcandidates_.back()));
    else
        getLayer(layerIdx, subdet).addTrackCandidate(&(trackcandidates_.back()));
}

void SDL::Event::addPixelSegmentsToEvent(std::vector<SDL::Hit> hits, float dPhiChange, float ptIn, float ptErr, float px, float py, float pz, float etaErr, int iSeed)
{
    // detId = 1 module means grand "pixel module" where everything related to pixel hits/md/segments will be stored to
    Module& pixelModule = getModule(1);

    // Assert that we provided a quadruplet pixel segment
    assert(hits.size() == 4);
    pixel_hits_.push_back(hits[0]);
    SDL::Hit* hit0_ptr = &pixel_hits_.back();
    pixel_hits_.push_back(hits[1]);
    SDL::Hit* hit1_ptr = &pixel_hits_.back();
    pixel_hits_.push_back(hits[2]);
    SDL::Hit* hit2_ptr = &pixel_hits_.back();
    pixel_hits_.push_back(hits[3]);
    SDL::Hit* hit3_ptr = &pixel_hits_.back();

    // Add hits to the pixel module
    pixelModule.addHit(hit0_ptr);
    pixelModule.addHit(hit1_ptr);
    pixelModule.addHit(hit2_ptr);
    pixelModule.addHit(hit3_ptr);

    // Create MiniDoublets
    SDL::MiniDoublet innerMiniDoublet(hit0_ptr, hit1_ptr);
    SDL::MiniDoublet outerMiniDoublet(hit2_ptr, hit3_ptr);
    // innerMiniDoublet.runMiniDoubletAllCombAlgo(); // Setting "all combination" pass in order to flag the pass bool flag
    // outerMiniDoublet.runMiniDoubletAllCombAlgo(); // Setting "all combination" pass in order to flag the pass bool flag

    pixel_miniDoublets_.push_back(innerMiniDoublet);
    SDL::MiniDoublet* innerMiniDoubletPtr = &pixel_miniDoublets_.back();
    pixel_miniDoublets_.push_back(outerMiniDoublet);
    SDL::MiniDoublet* outerMiniDoubletPtr = &pixel_miniDoublets_.back();

    // Create Segments
    segments_.push_back(SDL::Segment(innerMiniDoubletPtr, outerMiniDoubletPtr));
    SDL::Segment* pixelSegmentPtr = &segments_.back();

    // Set the deltaPhiChange
    pixelSegmentPtr->setDeltaPhiChange(dPhiChange);
    pixelSegmentPtr->setRecoVars("ptIn", ptIn);
    pixelSegmentPtr->setRecoVars("ptErr", ptErr);
    pixelSegmentPtr->setRecoVars("px", px);
    pixelSegmentPtr->setRecoVars("py", py);
    pixelSegmentPtr->setRecoVars("pz", pz);
    pixelSegmentPtr->setRecoVars("etaErr", etaErr);
    pixelSegmentPtr->setRecoVars("iSeed", iSeed);

    getPixelLayer().addSegment(pixelSegmentPtr);

}

void SDL::Event::createMiniDoublets(MDAlgo algo)
{
    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        // Create mini doublets
        createMiniDoubletsFromLowerModule(lowerModulePtr->detId(), algo);

    }
}

void SDL::Event::createMiniDoubletsFromLowerModule(unsigned int detId, SDL::MDAlgo algo)
{
    // Get reference to the lower Module
    Module& lowerModule = getModule(detId);

    // Get reference to the upper Module
    Module& upperModule = getModule(lowerModule.partnerDetId());

    // Double nested loops
    // Loop over lower module hits
    for (auto& lowerHitPtr : lowerModule.getHitPtrs())
    {

        // Get reference to lower Hit
        SDL::Hit& lowerHit = *lowerHitPtr;

        // Loop over upper module hits
        for (auto& upperHitPtr : upperModule.getHitPtrs())
        {

            // Get reference to upper Hit
            SDL::Hit& upperHit = *upperHitPtr;

            // Create a mini-doublet candidate
            SDL::MiniDoublet mdCand(lowerHitPtr, upperHitPtr);

            // Count the number of mdCand considered
            incrementNumberOfMiniDoubletCandidates(lowerModule);

            // Run mini-doublet algorithm on mdCand (mini-doublet candidate)
            mdCand.runMiniDoubletAlgo(algo, logLevel_);

            if (mdCand.passesMiniDoubletAlgo(algo))
            {

                // Count the number of md formed
                incrementNumberOfMiniDoublets(lowerModule);

                if (lowerModule.subdet() == SDL::Module::Barrel)
                    addMiniDoubletToEvent(mdCand, lowerModule.detId(), lowerModule.layer(), SDL::Layer::Barrel);
                else
                    addMiniDoubletToEvent(mdCand, lowerModule.detId(), lowerModule.layer(), SDL::Layer::Endcap);
            }

        }

    }

}

void SDL::Event::createPseudoMiniDoubletsFromAnchorModule(SDL::MDAlgo algo)
{

    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        unsigned int detId = lowerModulePtr->detId();

        // Get reference to the lower Module
        Module& lowerModule = getModule(detId);

        // Assign anchor hit pointers based on their hit type
        bool loopLower = true;
        if (lowerModule.moduleType() == SDL::Module::PS)
        {
            if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
            {
                loopLower = true;
            }
            else
            {
                loopLower = false;
            }
        }
        else
        {
            loopLower = true;
        }

        // Get reference to the upper Module
        Module& upperModule = getModule(lowerModule.partnerDetId());

        if (loopLower)
        {
            // Loop over lower module hits
            for (auto& lowerHitPtr : lowerModule.getHitPtrs())
            {
                // Get reference to lower Hit
                SDL::Hit& lowerHit = *lowerHitPtr;

                // Loop over upper module hits
                for (auto& upperHitPtr : upperModule.getHitPtrs())
                {

                    // Get reference to upper Hit
                    SDL::Hit& upperHit = *upperHitPtr;

                    // Create a mini-doublet candidate
                    SDL::MiniDoublet mdCand(lowerHitPtr, upperHitPtr);

                    // Count the number of mdCand considered
                    incrementNumberOfMiniDoubletCandidates(lowerModule);

                    // Run mini-doublet algorithm on mdCand (mini-doublet candidate)
                    mdCand.runMiniDoubletAlgo(SDL::AllComb_MDAlgo, logLevel_);

                    if (mdCand.passesMiniDoubletAlgo(SDL::AllComb_MDAlgo))
                    {

                        // Count the number of md formed
                        incrementNumberOfMiniDoublets(lowerModule);

                        if (lowerModule.subdet() == SDL::Module::Barrel)
                            addMiniDoubletToEvent(mdCand, lowerModule.detId(), lowerModule.layer(), SDL::Layer::Barrel);
                        else
                            addMiniDoubletToEvent(mdCand, lowerModule.detId(), lowerModule.layer(), SDL::Layer::Endcap);

                        // Break to exit on first pseudo mini-doublet
                        break;
                    }

                }

            }

        }
        else
        {
            // Loop over lower module hits
            for (auto& upperHitPtr : upperModule.getHitPtrs())
            {
                // Get reference to upper Hit
                SDL::Hit& upperHit = *upperHitPtr;

                // Loop over upper module hits
                for (auto& lowerHitPtr : lowerModule.getHitPtrs())
                {

                    // Get reference to lower Hit
                    SDL::Hit& lowerHit = *lowerHitPtr;

                    // Create a mini-doublet candidate
                    SDL::MiniDoublet mdCand(lowerHitPtr, upperHitPtr);

                    // Count the number of mdCand considered
                    incrementNumberOfMiniDoubletCandidates(lowerModule);

                    // Run mini-doublet algorithm on mdCand (mini-doublet candidate)
                    mdCand.runMiniDoubletAlgo(SDL::AllComb_MDAlgo, logLevel_);

                    if (mdCand.passesMiniDoubletAlgo(SDL::AllComb_MDAlgo))
                    {

                        // Count the number of md formed
                        incrementNumberOfMiniDoublets(lowerModule);

                        if (lowerModule.subdet() == SDL::Module::Barrel)
                            addMiniDoubletToEvent(mdCand, lowerModule.detId(), lowerModule.layer(), SDL::Layer::Barrel);
                        else
                            addMiniDoubletToEvent(mdCand, lowerModule.detId(), lowerModule.layer(), SDL::Layer::Endcap);

                        // Break to exit on first pseudo mini-doublet
                        break;
                    }

                }

            }

        }
    }

}

void SDL::Event::createSegments(SGAlgo algo)
{

    for (auto& segment_compatible_layer_pair : SDL::Layer::getListOfSegmentCompatibleLayerPairs())
    {
        int innerLayerIdx = segment_compatible_layer_pair.first.first;
        SDL::Layer::SubDet innerLayerSubDet = segment_compatible_layer_pair.first.second;
        int outerLayerIdx = segment_compatible_layer_pair.second.first;
        SDL::Layer::SubDet outerLayerSubDet = segment_compatible_layer_pair.second.second;
        createSegmentsFromTwoLayers(innerLayerIdx, innerLayerSubDet, outerLayerIdx, outerLayerSubDet, algo);
    }
}

void SDL::Event::createSegmentsWithModuleMap(SGAlgo algo)
{
    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        // Create mini doublets
        createSegmentsFromInnerLowerModule(lowerModulePtr->detId(), algo);

    }
}

void SDL::Event::createSegmentsFromInnerLowerModule(unsigned int detId, SDL::SGAlgo algo)
{

    // x's and y's are mini doublets
    // -------x--------
    // --------x------- <--- outer lower module
    //
    // --------y-------
    // -------y-------- <--- inner lower module

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module mini-doublets
    for (auto& innerMiniDoubletPtr : innerLowerModule.getMiniDoubletPtrs())
    {

        // Get reference to mini-doublet in inner lower module
        SDL::MiniDoublet& innerMiniDoublet = *innerMiniDoubletPtr;

        // Get connected outer lower module detids
        const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(detId);

        // Loop over connected outer lower modules
        for (auto& outerLowerModuleDetId : connectedModuleDetIds)
        {

            if (not hasModule(outerLowerModuleDetId))
                continue;

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(outerLowerModuleDetId);

            // Loop over outer lower module mini-doublets
            for (auto& outerMiniDoubletPtr : outerLowerModule.getMiniDoubletPtrs())
            {

                // Get reference to mini-doublet in outer lower module
                SDL::MiniDoublet& outerMiniDoublet = *outerMiniDoubletPtr;

                // Create a segment candidate
                SDL::Segment sgCand(innerMiniDoubletPtr, outerMiniDoubletPtr);

                // Run segment algorithm on sgCand (segment candidate)
                sgCand.runSegmentAlgo(algo, logLevel_);

                // Count the # of sgCands considered by layer
                incrementNumberOfSegmentCandidates(innerLowerModule);

                if (sgCand.passesSegmentAlgo(algo))
                {

                    // Count the # of sg formed by layer
                    incrementNumberOfSegments(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addSegmentToEvent(sgCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addSegmentToEvent(sgCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }
        }
    }
}

void SDL::Event::createTriplets(TPAlgo algo)
{
    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        // if (lowerModulePtr->layer() != 4 and lowerModulePtr->layer() != 3)
        //     continue;

        // Create mini doublets
        createTripletsFromInnerLowerModule(lowerModulePtr->detId(), algo);

    }
}

void SDL::Event::createTripletsFromInnerLowerModule(unsigned int detId, SDL::TPAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerSegmentPtr : innerLowerModule.getSegmentPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Segment& innerSegment = *innerSegmentPtr;

        // Get connected outer lower module detids
        const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(detId);

        // Loop over connected outer lower modules
        for (auto& outerLowerModuleDetId : connectedModuleDetIds)
        {

            if (not hasModule(outerLowerModuleDetId))
                continue;

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(outerLowerModuleDetId);

            // Loop over outer lower module mini-doublets
            for (auto& outerSegmentPtr : outerLowerModule.getSegmentPtrs())
            {

                // Get reference to mini-doublet in outer lower module
                SDL::Segment& outerSegment = *outerSegmentPtr;

                // Create a segment candidate
                SDL::Triplet tpCand(innerSegmentPtr, outerSegmentPtr);

                // Run segment algorithm on tpCand (segment candidate)
                tpCand.runTripletAlgo(algo, logLevel_);

                // Count the # of tpCands considered by layer
                incrementNumberOfTripletCandidates(innerLowerModule);

                if (tpCand.passesTripletAlgo(algo))
                {

                    // Count the # of sg formed by layer
                    incrementNumberOfTriplets(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTripletToEvent(tpCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTripletToEvent(tpCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }
        }
    }
}

void SDL::Event::createTracklets(TLAlgo algo)
{
    for (auto& tracklet_compatible_layer_pair : SDL::Layer::getListOfTrackletCompatibleLayerPairs())
    {
        int innerLayerIdx = tracklet_compatible_layer_pair.first.first;
        SDL::Layer::SubDet innerLayerSubDet = tracklet_compatible_layer_pair.first.second;
        int outerLayerIdx = tracklet_compatible_layer_pair.second.first;
        SDL::Layer::SubDet outerLayerSubDet = tracklet_compatible_layer_pair.second.second;
        createTrackletsFromTwoLayers(innerLayerIdx, innerLayerSubDet, outerLayerIdx, outerLayerSubDet, algo);
    }
}

void SDL::Event::createTrackletsWithModuleMap(TLAlgo algo)
{
    if (logLevel_ == SDL::Log_Debug)
        SDL::cout << "SDL::Event::createTrackletsWithModuleMap()" << std::endl;

    // Loop over lower modules
    int nModuleProcessed = 0;
    int nTotalLowerModule = getLowerModulePtrs().size();

    if (logLevel_ == SDL::Log_Debug)
        SDL::cout <<  " nTotalLowerModule: " << nTotalLowerModule <<  std::endl;

    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        if (logLevel_ == SDL::Log_Debug)
            if (nModuleProcessed % 1000 == 0)
                SDL::cout <<  "    nModuleProcessed: " << nModuleProcessed <<  std::endl;

        if (logLevel_ == SDL::Log_Debug)
        {
            std::cout <<  " lowerModulePtr->subdet(): " << lowerModulePtr->subdet() <<  std::endl;
            std::cout <<  " lowerModulePtr->layer(): " << lowerModulePtr->layer() <<  std::endl;
            std::cout <<  " lowerModulePtr->getSegmentPtrs().size(): " << lowerModulePtr->getSegmentPtrs().size() <<  std::endl;
        }

        // if (lowerModulePtr->layer() != 1)
        //     continue;

        // Create mini doublets
        createTrackletsFromInnerLowerModule(lowerModulePtr->detId(), algo);

        nModuleProcessed++;

    }
}

// Create tracklets from inner modules
void SDL::Event::createTrackletsFromInnerLowerModule(unsigned int detId, SDL::TLAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerSegmentPtr : innerLowerModule.getSegmentPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Segment& innerSegment = *innerSegmentPtr;

        // Get the outer mini-doublet module detId
        const SDL::Module& innerSegmentOuterModule = innerSegment.outerMiniDoubletPtr()->lowerHitPtr()->getModule();

        unsigned int innerSegmentOuterModuleDetId = innerSegmentOuterModule.detId();

        // Get connected outer lower module detids
        const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(innerSegmentOuterModuleDetId);

        // Loop over connected outer lower modules
        for (auto& outerLowerModuleDetId : connectedModuleDetIds)
        {

            if (not hasModule(outerLowerModuleDetId))
                continue;

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(outerLowerModuleDetId);

            // Loop over outer lower module mini-doublets
            for (auto& outerSegmentPtr : outerLowerModule.getSegmentPtrs())
            {

                // Count the # of tlCands considered by layer
                incrementNumberOfTrackletCandidates(innerLowerModule);

                // Get reference to mini-doublet in outer lower module
                SDL::Segment& outerSegment = *outerSegmentPtr;

                // Create a tracklet candidate
                SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

                // Run segment algorithm on tlCand (tracklet candidate)
                tlCand.runTrackletAlgo(algo, logLevel_);

                if (tlCand.passesTrackletAlgo(algo))
                {

                    // Count the # of sg formed by layer
                    incrementNumberOfTracklets(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

}

// Create tracklets
void SDL::Event::createTrackletsWithAGapWithModuleMap(TLAlgo algo)
{
    if (logLevel_ == SDL::Log_Debug)
        SDL::cout << "SDL::Event::createTrackletsWithAGapWithModuleMap()" << std::endl;

    // Loop over lower modules
    int nModuleProcessed = 0;
    int nTotalLowerModule = getLowerModulePtrs().size();

    if (logLevel_ == SDL::Log_Debug)
        SDL::cout <<  " nTotalLowerModule: " << nTotalLowerModule <<  std::endl;

    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        if (logLevel_ == SDL::Log_Debug)
            if (nModuleProcessed % 1000 == 0)
                SDL::cout <<  "    nModuleProcessed: " << nModuleProcessed <<  std::endl;

        if (logLevel_ == SDL::Log_Debug)
        {
            std::cout <<  " lowerModulePtr->subdet(): " << lowerModulePtr->subdet() <<  std::endl;
            std::cout <<  " lowerModulePtr->layer(): " << lowerModulePtr->layer() <<  std::endl;
            std::cout <<  " lowerModulePtr->getSegmentPtrs().size(): " << lowerModulePtr->getSegmentPtrs().size() <<  std::endl;
        }

        // if (lowerModulePtr->layer() != 1)
        //     continue;

        // Create mini doublets
        createTrackletsWithAGapFromInnerLowerModule(lowerModulePtr->detId(), algo);

        nModuleProcessed++;

    }
}

// Create tracklets from inner modules
void SDL::Event::createTrackletsWithAGapFromInnerLowerModule(unsigned int detId, SDL::TLAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerSegmentPtr : innerLowerModule.getSegmentPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Segment& innerSegment = *innerSegmentPtr;

        // Get the outer mini-doublet module detId
        const SDL::Module& innerSegmentOuterModule = innerSegment.outerMiniDoubletPtr()->lowerHitPtr()->getModule();

        unsigned int innerSegmentOuterModuleDetId = innerSegmentOuterModule.detId();

        // Get connected middle module detids
        const std::vector<unsigned int>& connectedMiddleModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(innerSegmentOuterModuleDetId);

        for (auto& middleLowerModuleDetId : connectedMiddleModuleDetIds)
        {

            // Get connected outer lower module detids
            const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(middleLowerModuleDetId);

            // Loop over connected outer lower modules
            for (auto& outerLowerModuleDetId : connectedModuleDetIds)
            {

                if (not hasModule(outerLowerModuleDetId))
                    continue;

                // Get reference to the outer lower module
                Module& outerLowerModule = getModule(outerLowerModuleDetId);

                // Loop over outer lower module mini-doublets
                for (auto& outerSegmentPtr : outerLowerModule.getSegmentPtrs())
                {

                    // Count the # of tlCands considered by layer
                    incrementNumberOfTrackletCandidates(innerLowerModule);

                    // Get reference to mini-doublet in outer lower module
                    SDL::Segment& outerSegment = *outerSegmentPtr;

                    // Create a tracklet candidate
                    SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

                    // Run segment algorithm on tlCand (tracklet candidate)
                    tlCand.runTrackletAlgo(algo, logLevel_);

                    if (tlCand.passesTrackletAlgo(algo))
                    {

                        // Count the # of sg formed by layer
                        incrementNumberOfTracklets(innerLowerModule);

                        if (innerLowerModule.subdet() == SDL::Module::Barrel)
                            addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                        else
                            addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
                    }

                }

            }

        }

    }

}

// Create tracklets
void SDL::Event::createTrackletsWithTwoGapsWithModuleMap(TLAlgo algo)
{
    if (logLevel_ == SDL::Log_Debug)
        SDL::cout << "SDL::Event::createTrackletsWithTwoGapsWithModuleMap()" << std::endl;

    // Loop over lower modules
    int nModuleProcessed = 0;
    int nTotalLowerModule = getLowerModulePtrs().size();

    if (logLevel_ == SDL::Log_Debug)
        SDL::cout <<  " nTotalLowerModule: " << nTotalLowerModule <<  std::endl;

    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        if (logLevel_ == SDL::Log_Debug)
            if (nModuleProcessed % 1000 == 0)
                SDL::cout <<  "    nModuleProcessed: " << nModuleProcessed <<  std::endl;

        if (logLevel_ == SDL::Log_Debug)
        {
            std::cout <<  " lowerModulePtr->subdet(): " << lowerModulePtr->subdet() <<  std::endl;
            std::cout <<  " lowerModulePtr->layer(): " << lowerModulePtr->layer() <<  std::endl;
            std::cout <<  " lowerModulePtr->getSegmentPtrs().size(): " << lowerModulePtr->getSegmentPtrs().size() <<  std::endl;
        }

        // if (lowerModulePtr->layer() != 1)
        //     continue;

        // Create mini doublets
        createTrackletsWithTwoGapsFromInnerLowerModule(lowerModulePtr->detId(), algo);

        nModuleProcessed++;

    }
}

// Create tracklets from inner modules
void SDL::Event::createTrackletsWithTwoGapsFromInnerLowerModule(unsigned int detId, SDL::TLAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerSegmentPtr : innerLowerModule.getSegmentPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Segment& innerSegment = *innerSegmentPtr;

        // Get the outer mini-doublet module detId
        const SDL::Module& innerSegmentOuterModule = innerSegment.outerMiniDoubletPtr()->lowerHitPtr()->getModule();

        unsigned int innerSegmentOuterModuleDetId = innerSegmentOuterModule.detId();

        // Get connected middle module detids
        const std::vector<unsigned int>& connectedMiddleModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(innerSegmentOuterModuleDetId);

        for (auto& middleLowerModuleDetId : connectedMiddleModuleDetIds)
        {

            // Get connected second middle module detids
            const std::vector<unsigned int>& connectedMiddle2ModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(middleLowerModuleDetId);

            for (auto& middle2LowerModuleDetId : connectedMiddle2ModuleDetIds)
            {

                // Get connected outer lower module detids
                const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(middle2LowerModuleDetId);

                // Loop over connected outer lower modules
                for (auto& outerLowerModuleDetId : connectedModuleDetIds)
                {

                    if (not hasModule(outerLowerModuleDetId))
                        continue;

                    // Get reference to the outer lower module
                    Module& outerLowerModule = getModule(outerLowerModuleDetId);

                    // Loop over outer lower module mini-doublets
                    for (auto& outerSegmentPtr : outerLowerModule.getSegmentPtrs())
                    {

                        // Count the # of tlCands considered by layer
                        incrementNumberOfTrackletCandidates(innerLowerModule);

                        // Get reference to mini-doublet in outer lower module
                        SDL::Segment& outerSegment = *outerSegmentPtr;

                        // Create a tracklet candidate
                        SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

                        // Run segment algorithm on tlCand (tracklet candidate)
                        tlCand.runTrackletAlgo(algo, logLevel_);

                        if (tlCand.passesTrackletAlgo(algo))
                        {

                            // Count the # of sg formed by layer
                            incrementNumberOfTracklets(innerLowerModule);

                            if (innerLowerModule.subdet() == SDL::Module::Barrel)
                                addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                            else
                                addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
                        }

                    }

                }

            }

        }

    }

}

// Create tracklets via navigation
void SDL::Event::createTrackletsViaNavigation(SDL::TLAlgo algo)
{
    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        // Get reference to the inner lower Module
        Module& innerLowerModule = getModule(lowerModulePtr->detId());

        // Triple nested loops
        // Loop over inner lower module for segments
        for (auto& innerSegmentPtr : innerLowerModule.getSegmentPtrs())
        {

            // Get reference to segment in inner lower module
            SDL::Segment& innerSegment = *innerSegmentPtr;

            // Get the connecting segment ptrs
            for (auto& connectingSegmentPtr : innerSegmentPtr->outerMiniDoubletPtr()->getListOfOutwardSegmentPtrs())
            {

                for (auto& outerSegmentPtr : connectingSegmentPtr->outerMiniDoubletPtr()->getListOfOutwardSegmentPtrs())
                {

                    // Count the # of tlCands considered by layer
                    incrementNumberOfTrackletCandidates(innerLowerModule);

                    // Get reference to mini-doublet in outer lower module
                    SDL::Segment& outerSegment = *outerSegmentPtr;

                    // Create a tracklet candidate
                    SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

                    // Run segment algorithm on tlCand (tracklet candidate)
                    tlCand.runTrackletAlgo(algo, logLevel_);

                    if (tlCand.passesTrackletAlgo(algo))
                    {

                        // Count the # of sg formed by layer
                        incrementNumberOfTracklets(innerLowerModule);

                        if (innerLowerModule.subdet() == SDL::Module::Barrel)
                            addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                        else
                            addTrackletToEvent(tlCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
                    }

                }

            }

        }
    }

}


// Create tracklets with pixel to barrel
void SDL::Event::createTrackletsWithPixelAndBarrel(TLAlgo algo)
{
    if (logLevel_ == SDL::Log_Debug)
        SDL::cout << "SDL::Event::createTrackletsWithPixelAndBarrel()" << std::endl;

    // Loop over lower modules
    int nModuleProcessed = 0;
    int nTotalLowerModule = getLowerModulePtrs().size();

    if (logLevel_ == SDL::Log_Debug)
        SDL::cout <<  " nTotalLowerModule: " << nTotalLowerModule <<  std::endl;

    int nCombinations = 0;


    // Loop over inner lower module for segments
    for (auto& innerSegmentPtr : getPixelLayer().getSegmentPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Segment& innerSegment = *innerSegmentPtr;

        // Get reference to the inner lower Module
        Module& pixelModule = getModule(1);
        Module& innerLowerModule = pixelModule;

        for (auto& lowerModulePtr : getLowerModulePtrs())
        {

            if (lowerModulePtr->getSegmentPtrs().size() == 0)
                continue;

            // if (lowerModulePtr->moduleType() != SDL::Module::PS)
            //     continue;

            // if (logLevel_ == SDL::Log_Debug)
            //     if (nModuleProcessed % 1000 == 0)
            //         SDL::cout <<  "    nModuleProcessed: " << nModuleProcessed <<  std::endl;

            // if (logLevel_ == SDL::Log_Debug)
            // {
            //     std::cout <<  " lowerModulePtr->subdet(): " << lowerModulePtr->subdet() <<  std::endl;
            //     std::cout <<  " lowerModulePtr->layer(): " << lowerModulePtr->layer() <<  std::endl;
            //     std::cout <<  " lowerModulePtr->getSegmentPtrs().size(): " << lowerModulePtr->getSegmentPtrs().size() <<  std::endl;
            // }

            // Get reference to the outer lower module
            Module& outerLowerModule = *lowerModulePtr;

            // Loop over outer lower module mini-doublets
            for (auto& outerSegmentPtr : outerLowerModule.getSegmentPtrs())
            {

                // // Count the # of tlCands considered by layer
                // incrementNumberOfTrackletCandidates(innerLowerModule);

                // Get reference to mini-doublet in outer lower module
                SDL::Segment& outerSegment = *outerSegmentPtr;

                // Create a tracklet candidate
                SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

                // Run segment algorithm on tlCand (tracklet candidate)
                tlCand.runTrackletAlgo(algo, logLevel_);

                if (logLevel_ == SDL::Log_Debug)
                {
                    // int passbit = tlCand.getPassBitsDefaultAlgo();
                    // std::bitset<8> x(passbit);
                    // std::cout <<  " passbit: " << x <<  std::endl;
                }

                nCombinations++;

                if (tlCand.passesTrackletAlgo(algo))
                {

                    // Count the # of sg formed by layer
                    incrementNumberOfTracklets(innerLowerModule);

                    addTrackletToEvent(tlCand, 1/*pixel module=1*/, 0/*pixel is layer=0*/, SDL::Layer::Barrel);
                }

            }

            nModuleProcessed++;

        }

    }

    if (logLevel_ == SDL::Log_Debug)
        std::cout << "SDL::Event::createTrackletsWithPixelAndBarrel(): nCombinations = " << nCombinations << std::endl;


    // for (auto& lowerModulePtr : getLowerModulePtrs())
    // {

    //     if (lowerModulePtr->getSegmentPtrs().size() == 0)
    //         continue;

    //     if (logLevel_ == SDL::Log_Debug)
    //         if (nModuleProcessed % 1000 == 0)
    //             SDL::cout <<  "    nModuleProcessed: " << nModuleProcessed <<  std::endl;

    //     if (logLevel_ == SDL::Log_Debug)
    //     {
    //         std::cout <<  " lowerModulePtr->subdet(): " << lowerModulePtr->subdet() <<  std::endl;
    //         std::cout <<  " lowerModulePtr->layer(): " << lowerModulePtr->layer() <<  std::endl;
    //         std::cout <<  " lowerModulePtr->getSegmentPtrs().size(): " << lowerModulePtr->getSegmentPtrs().size() <<  std::endl;
    //     }

    //     // Get reference to the inner lower Module
    //     Module& pixelModule = getModule(1);
    //     Module& innerLowerModule = pixelModule;

    //     // Triple nested loops
    //     // Loop over inner lower module for segments
    //     for (auto& innerSegmentPtr : getPixelLayer().getSegmentPtrs())
    //     {

    //         // Get reference to segment in inner lower module
    //         SDL::Segment& innerSegment = *innerSegmentPtr;

    //         // Get reference to the outer lower module
    //         Module& outerLowerModule = *lowerModulePtr;

    //         // Loop over outer lower module mini-doublets
    //         for (auto& outerSegmentPtr : outerLowerModule.getSegmentPtrs())
    //         {

    //             // // Count the # of tlCands considered by layer
    //             // incrementNumberOfTrackletCandidates(innerLowerModule);

    //             // Get reference to mini-doublet in outer lower module
    //             SDL::Segment& outerSegment = *outerSegmentPtr;

    //             // Create a tracklet candidate
    //             SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

    //             // Run segment algorithm on tlCand (tracklet candidate)
    //             tlCand.runTrackletAlgo(algo, logLevel_);

    //             nCombinations++;

    //             if (tlCand.passesTrackletAlgo(algo))
    //             {

    //                 // Count the # of sg formed by layer
    //                 incrementNumberOfTracklets(innerLowerModule);

    //                 addTrackletToEvent(tlCand, 1/*pixel module=1*/, 0/*pixel is layer=0*/, SDL::Layer::Barrel);
    //             }

    //         }

    //     }

    //     nModuleProcessed++;

    // }

    // if (logLevel_ == SDL::Log_Debug)
    //     std::cout << "SDL::Event::createTrackletsWithPixelAndBarrel(): nCombinations = " << nCombinations << std::endl;
}


// Create tracklets from two layers (inefficient way)
void SDL::Event::createTrackletsFromTwoLayers(int innerLayerIdx, SDL::Layer::SubDet innerLayerSubDet, int outerLayerIdx, SDL::Layer::SubDet outerLayerSubDet, TLAlgo algo)
{
    Layer& innerLayer = getLayer(innerLayerIdx, innerLayerSubDet);
    Layer& outerLayer = getLayer(outerLayerIdx, outerLayerSubDet);

    for (auto& innerSegmentPtr : innerLayer.getSegmentPtrs())
    {
        SDL::Segment& innerSegment = *innerSegmentPtr;
        for (auto& outerSegmentPtr : outerLayer.getSegmentPtrs())
        {
            // SDL::Segment& outerSegment = *outerSegmentPtr;

            // if (SDL::Tracklet::isSegmentPairATracklet(innerSegment, outerSegment, algo, logLevel_))
            //     addTrackletToLowerLayer(SDL::Tracklet(innerSegmentPtr, outerSegmentPtr), innerLayerIdx, innerLayerSubDet);

            SDL::Segment& outerSegment = *outerSegmentPtr;

            SDL::Tracklet tlCand(innerSegmentPtr, outerSegmentPtr);

            tlCand.runTrackletAlgo(algo, logLevel_);

            // Count the # of tracklet candidate considered by layer
            incrementNumberOfTrackletCandidates(innerLayer);

            if (tlCand.passesTrackletAlgo(algo))
            {

                // Count the # of tracklet formed by layer
                incrementNumberOfTracklets(innerLayer);

                addTrackletToLowerLayer(tlCand, innerLayerIdx, innerLayerSubDet);
            }

        }
    }
}

// Create segments from two layers (inefficient way)
void SDL::Event::createSegmentsFromTwoLayers(int innerLayerIdx, SDL::Layer::SubDet innerLayerSubDet, int outerLayerIdx, SDL::Layer::SubDet outerLayerSubDet, SGAlgo algo)
{
    Layer& innerLayer = getLayer(innerLayerIdx, innerLayerSubDet);
    Layer& outerLayer = getLayer(outerLayerIdx, outerLayerSubDet);

    for (auto& innerMiniDoubletPtr : innerLayer.getMiniDoubletPtrs())
    {
        SDL::MiniDoublet& innerMiniDoublet = *innerMiniDoubletPtr;

        for (auto& outerMiniDoubletPtr : outerLayer.getMiniDoubletPtrs())
        {
            SDL::MiniDoublet& outerMiniDoublet = *outerMiniDoubletPtr;

            SDL::Segment sgCand(innerMiniDoubletPtr, outerMiniDoubletPtr);

            sgCand.runSegmentAlgo(algo, logLevel_);

            if (sgCand.passesSegmentAlgo(algo))
            {
                const SDL::Module& innerLowerModule = innerMiniDoubletPtr->lowerHitPtr()->getModule();
                if (innerLowerModule.subdet() == SDL::Module::Barrel)
                    addSegmentToEvent(sgCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Barrel);
                else
                    addSegmentToEvent(sgCand, innerLowerModule.detId(), innerLowerModule.layer(), SDL::Layer::Endcap);
            }

        }
    }
}

void SDL::Event::createTrackCandidates(TCAlgo algo)
{
    // TODO Implement some structure for Track Candidates
    // for (auto& trackCandidate_compatible_layer_pair : SDL::Layer::getListOfTrackCandidateCompatibleLayerPairs())
    // {
    //     int innerLayerIdx = trackCandidate_compatible_layer_pair.first.first;
    //     SDL::Layer::SubDet innerLayerSubDet = trackCandidate_compatible_layer_pair.first.second;
    //     int outerLayerIdx = trackCandidate_compatible_layer_pair.second.first;
    //     SDL::Layer::SubDet outerLayerSubDet = trackCandidate_compatible_layer_pair.second.second;
    //     createTrackCandidatesFromTwoLayers(innerLayerIdx, innerLayerSubDet, outerLayerIdx, outerLayerSubDet, algo);
    // }

    createTrackCandidatesFromTwoLayers(1, SDL::Layer::Barrel, 3, SDL::Layer::Barrel, algo);

}

// Create trackCandidates from two layers (inefficient way)
void SDL::Event::createTrackCandidatesFromTwoLayers(int innerLayerIdx, SDL::Layer::SubDet innerLayerSubDet, int outerLayerIdx, SDL::Layer::SubDet outerLayerSubDet, TCAlgo algo)
{
    Layer& innerLayer = getLayer(innerLayerIdx, innerLayerSubDet);
    Layer& outerLayer = getLayer(outerLayerIdx, outerLayerSubDet);

    for (auto& innerTrackletPtr : innerLayer.getTrackletPtrs())
    {
        SDL::Tracklet& innerTracklet = *innerTrackletPtr;

        for (auto& outerTrackletPtr : outerLayer.getTrackletPtrs())
        {

            SDL::Tracklet& outerTracklet = *outerTrackletPtr;

            SDL::TrackCandidate tcCand(innerTrackletPtr, outerTrackletPtr);

            tcCand.runTrackCandidateAlgo(algo, logLevel_);

            // Count the # of track candidates considered
            incrementNumberOfTrackCandidateCandidates(innerLayer);

            if (tcCand.passesTrackCandidateAlgo(algo))
            {

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidates(innerLayer);

                addTrackCandidateToLowerLayer(tcCand, innerLayerIdx, innerLayerSubDet);
            }

        }
    }
}

void SDL::Event::createTrackCandidatesFromTriplets(TCAlgo algo)
{
    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        // if (not (lowerModulePtr->layer() == 1))
        //     continue;

        // Create mini doublets
        createTrackCandidatesFromInnerModulesFromTriplets(lowerModulePtr->detId(), algo);

    }
}

void SDL::Event::createTrackCandidatesFromInnerModulesFromTriplets(unsigned int detId, SDL::TCAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerTripletPtr : innerLowerModule.getTripletPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Triplet& innerTriplet = *innerTripletPtr;

        // Get the outer mini-doublet module detId
        const SDL::Module& innerTripletOutermostModule = innerTriplet.outerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule();

        unsigned int innerTripletOutermostModuleDetId = innerTripletOutermostModule.detId();

        // Get connected outer lower module detids
        const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(innerTripletOutermostModuleDetId);

        // Loop over connected outer lower modules
        for (auto& outerLowerModuleDetId : connectedModuleDetIds)
        {

            if (not hasModule(outerLowerModuleDetId))
                continue;

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(outerLowerModuleDetId);

            // Loop over outer lower module mini-doublets
            for (auto& outerTripletPtr : outerLowerModule.getTripletPtrs())
            {

                // Count the # of tlCands considered by layer
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                // Segment between innerSgOuterMD - outerSgInnerMD
                SDL::Segment sgCand(innerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr(),outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr());

                // Run the segment algo (supposedly is fast)
                sgCand.runSegmentAlgo(SDL::Default_SGAlgo, logLevel_);

                if (not (sgCand.passesSegmentAlgo(SDL::Default_SGAlgo)))
                {
                    continue;
                }

                // SDL::Tracklet tlCand(innerTripletPtr->innerSegmentPtr(), &sgCand);

                // // Run the segment algo (supposedly is fast)
                // tlCand.runTrackletAlgo(SDL::Default_TLAlgo, logLevel_);

                // if (not (tlCand.passesTrackletAlgo(SDL::Default_TLAlgo)))
                // {
                //     continue;
                // }

                SDL::Tracklet tlCandOuter(&sgCand, outerTripletPtr->outerSegmentPtr());

                // Run the segment algo (supposedly is fast)
                tlCandOuter.runTrackletAlgo(SDL::Default_TLAlgo, logLevel_);

                if (not (tlCandOuter.passesTrackletAlgo(SDL::Default_TLAlgo)))
                {
                    continue;
                }

                SDL::TrackCandidate tcCand(innerTripletPtr, outerTripletPtr);

                // if (tcCand.passesTrackCandidateAlgo(algo))
                // {

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidates(innerLowerModule);

                addTrackCandidateToLowerLayer(tcCand, 1, SDL::Layer::Barrel);
                // if (innerLowerModule.subdet() == SDL::Module::Barrel)
                //     addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                // else
                //     addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);

                // }

            }

        }

    }


}

void SDL::Event::createTrackCandidatesFromTracklets(TCAlgo algo)
{
    // Loop over lower modules
    for (auto& lowerModulePtr : getLowerModulePtrs())
    {

        // if (not (lowerModulePtr->layer() == 1))
        //     continue;

        // Create mini doublets
        createTrackCandidatesFromInnerModulesFromTracklets(lowerModulePtr->detId(), algo);
        // createTrackCandidatesFromInnerModulesFromTrackletsToTriplets(lowerModulePtr->detId(), algo);

    }
}

void SDL::Event::createTrackCandidatesFromInnerModulesFromTracklets(unsigned int detId, SDL::TCAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerTrackletPtr : innerLowerModule.getTrackletPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Tracklet& innerTracklet = *innerTrackletPtr;

        // Get the outer mini-doublet module detId
        const SDL::Module& innerTrackletSecondModule = innerTracklet.innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule();

        unsigned int innerTrackletSecondModuleDetId = innerTrackletSecondModule.detId();

        // Get connected outer lower module detids
        const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(innerTrackletSecondModuleDetId);

        // Loop over connected outer lower modules
        for (auto& outerLowerModuleDetId : connectedModuleDetIds)
        {

            if (not hasModule(outerLowerModuleDetId))
                continue;

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(outerLowerModuleDetId);

            // Loop over outer lower module mini-doublets
            for (auto& outerTrackletPtr : outerLowerModule.getTrackletPtrs())
            {

                SDL::Tracklet& outerTracklet = *outerTrackletPtr;

                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateAlgo(algo, logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

}

void SDL::Event::createTrackCandidatesTest_v1(SDL::TCAlgo algo)
{

    // August 31, 2020 Considering only the following kinds of track candidates in barrel
    // - 1 2 3 4 + 3 4 5 6  : 0 missing
    // - 1 2 3 4 + 3 4 5    : 6 missing
    // - 1 2 4 5 + 4 5 6    : 3 missing
    // - 2 3 4 5 + 4 5 6    : 1 missing
    // - 1 2 3 + 2 3 5 6    : 4 missing
    // 2 and 5 missing not done at this point

    const SDL::Layer& barrelLayer1 = getLayer(1, SDL::Layer::Barrel);
    const SDL::Layer& barrelLayer2 = getLayer(2, SDL::Layer::Barrel);
    const SDL::Layer& barrelLayer3 = getLayer(3, SDL::Layer::Barrel);
    const SDL::Layer& barrelLayer4 = getLayer(4, SDL::Layer::Barrel);
    const SDL::Layer& barrelLayer5 = getLayer(5, SDL::Layer::Barrel);
    const SDL::Layer& barrelLayer6 = getLayer(6, SDL::Layer::Barrel);

    for (auto& innerTrackletPtr : barrelLayer1.getTrackletPtrs())
    {

        // Check if it is a barrel only tracklet
        bool isBBBB =
        ((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

        bool is1234 = 
        ((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 1) and
        ((innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 2) and
        ((innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
        ((innerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4);

        bool is1245 = 
        ((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 1) and
        ((innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 2) and
        ((innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
        ((innerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5);

        if (isBBBB and is1234)
        {

            // Get reference to the inner lower Module
            Module& innerLowerModule = getModule((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId());

            unsigned int innerTrackletModule2 = (innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();
            unsigned int innerTrackletModule3 = (innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(innerTrackletModule3);

            // Loop over outer lower module Tracklets
            for (auto& outerTrackletPtr : outerLowerModule.getTrackletPtrs())
            {

                // Check if it is a barrel only tracklet
                bool isOuterBBBB =
                    ((outerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

                // Check if it is a barrel only tracklet
                bool isOuter3456 =
                    ((outerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
                    ((outerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
                    ((outerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5) and
                    ((outerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 6);

                if (not (isOuterBBBB and isOuter3456))
                {
                    continue;
                }

                SDL::Tracklet& outerTracklet = *outerTrackletPtr;

                //============================================================
                //
                //
                // Type 1: 1 2 3 4 + 3 4 5 6 : no layer missing
                //
                //
                //============================================================
                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateAlgo(algo, logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

            // Loop over outer lower triplets
            for (auto& outerTripletPtr : outerLowerModule.getTripletPtrs())
            {

                // Check if it is a barrel only triplet
                bool isOuterBBB =
                    ((outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

                // Check if it is a barrel only triplet
                bool isOuter345 =
                    ((outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
                    ((outerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
                    ((outerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5);

                if (not (isOuterBBB and isOuter345))
                {
                    continue;
                }

                SDL::Triplet& outerTriplet = *outerTripletPtr;

                //============================================================
                //
                //
                // Type 1: 1 2 3 4 + 3 4 5   : layer 6 is missing
                //
                //
                //============================================================
                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

                tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }


        }
        else if (isBBBB and is1245)
        {

            // Get reference to the inner lower Module
            Module& innerLowerModule = getModule((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId());

            unsigned int innerTrackletModule2 = (innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();
            unsigned int innerTrackletModule3 = (innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(innerTrackletModule3);

            // Loop over outer lower triplets
            for (auto& outerTripletPtr : outerLowerModule.getTripletPtrs())
            {

                // Check if it is a barrel only triplet
                bool isOuterBBB =
                    ((outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

                // Check if it is a barrel only triplet
                bool isOuter456 =
                    ((outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
                    ((outerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5) and
                    ((outerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 6);

                if (not (isOuterBBB and isOuter456))
                {
                    continue;
                }

                SDL::Triplet& outerTriplet = *outerTripletPtr;

                //============================================================
                //
                //
                // Type 1: 1 2 4 5 + 4 5 6   : layer 3 is missing
                //
                //
                //============================================================
                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

                tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

    for (auto& innerTrackletPtr : barrelLayer2.getTrackletPtrs())
    {

        // Check if it is a barrel only tracklet
        bool isBBBB =
        ((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

        bool is2345 = 
        ((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 2) and
        ((innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
        ((innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
        ((innerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5);

        if (isBBBB and is2345)
        {

            // Get reference to the inner lower Module
            Module& innerLowerModule = getModule((innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId());

            unsigned int innerTrackletModule2 = (innerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();
            unsigned int innerTrackletModule3 = (innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(innerTrackletModule3);

            // Loop over outer lower triplets
            for (auto& outerTripletPtr : outerLowerModule.getTripletPtrs())
            {

                // Check if it is a barrel only triplet
                bool isOuterBBB =
                    ((outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

                // Check if it is a barrel only triplet
                bool isOuter456 =
                    ((outerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
                    ((outerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5) and
                    ((outerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 6);

                if (not (isOuterBBB and isOuter456))
                {
                    continue;
                }

                SDL::Triplet& outerTriplet = *outerTripletPtr;

                //============================================================
                //
                //
                // Type 1: 2 3 4 5 + 4 5 6   : layer 1 is missing
                //
                //
                //============================================================
                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

                tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

    for (auto& innerTripletPtr : barrelLayer1.getTripletPtrs())
    {

        // Check if it is a barrel only tracklet
        bool isBBB =
        ((innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

        bool is123 = 
        ((innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 1) and
        ((innerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 2) and
        ((innerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3);

        if (isBBB and is123)
        {

            // Get reference to the inner lower Module
            Module& innerLowerModule = getModule((innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId());

            // unsigned int innerTripletModule2 = (innerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();
            unsigned int innerTripletModule2 = (innerTripletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(innerTripletModule2);

            // Loop over outer lower triplets
            for (auto& outerTrackletPtr : outerLowerModule.getTrackletPtrs())
            {

                // Check if it is a barrel only triplet
                bool isOuterBBBB =
                    ((outerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

                // Check if it is a barrel only triplet
                bool isOuter2356 =
                    ((outerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 2) and
                    ((outerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
                    ((outerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5) and
                    ((outerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 6);

                if (not (isOuterBBBB and isOuter2356))
                {
                    continue;
                }

                SDL::Tracklet& outerTracklet = *outerTrackletPtr;

                //============================================================
                //
                //
                // Type 1: 2 3 4 5 + 4 5 6   : layer 1 is missing
                //
                //
                //============================================================
                SDL::TrackCandidate tcCand(innerTripletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateInnerTripletToOuterTracklet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

    for (auto& innerTripletPtr : barrelLayer2.getTripletPtrs())
    {

        // Check if it is a barrel only tracklet
        bool isBBB =
        ((innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
        ((innerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

        bool is234 = 
        ((innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 2) and
        ((innerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
        ((innerTripletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4);

        if (isBBB and is234)
        {

            // Get reference to the inner lower Module
            Module& innerLowerModule = getModule((innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId());

            // unsigned int innerTripletModule2 = (innerTripletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();
            unsigned int innerTripletModule2 = (innerTripletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule()).detId();

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(innerTripletModule2);

            // Loop over outer lower triplets
            for (auto& outerTrackletPtr : outerLowerModule.getTrackletPtrs())
            {

                // Check if it is a barrel only triplet
                bool isOuterBBBB =
                    ((outerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel) and
                    ((outerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Barrel);

                // Check if it is a barrel only triplet
                bool isOuter3456 =
                    ((outerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 3) and
                    ((outerTrackletPtr->innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 4) and
                    ((outerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 5) and
                    ((outerTrackletPtr->outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).layer() == 6);

                if (not (isOuterBBBB and isOuter3456))
                {
                    continue;
                }

                SDL::Tracklet& outerTracklet = *outerTrackletPtr;

                //============================================================
                //
                //
                // Type 1: 2 3 4 5 + 4 5 6   : layer 1 is missing
                //
                //
                //============================================================
                SDL::TrackCandidate tcCand(innerTripletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateInnerTripletToOuterTracklet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

}

void SDL::Event::createTrackCandidatesTest_v2(SDL::TCAlgo algo)
{

    // September 10, 2020 Consider ALL cases
    // Loop over all tracklets, a-b-c-d go to c then get tracklets or triplets and ask whether segment is shared
    // Ditto for Triplet -> Tracklet
    for (auto& layerPtr : getLayerPtrs())
    {

        for (auto& innerTrackletPtr : layerPtr->getTrackletPtrs())
        {
            SDL::Module& innerLowerModule = getModule(innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
            const SDL::Module& commonModule = innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule();

            for (auto& outerTrackletPtr : commonModule.getTrackletPtrs())
            {

                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateAlgo(algo, logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

            for (auto& outerTripletPtr : commonModule.getTripletPtrs())
            {

                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

                tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

        for (auto& innerTripletPtr : layerPtr->getTripletPtrs())
        {
            SDL::Module& innerLowerModule = getModule(innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
            const SDL::Module& commonModule = innerTripletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule();

            for (auto& outerTrackletPtr : commonModule.getTrackletPtrs())
            {

                SDL::TrackCandidate tcCand(innerTripletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateInnerTripletToOuterTracklet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

}

void SDL::Event::createTrackCandidatesTest_v3(SDL::TCAlgo algo)
{

    // September 10, 2020 Consider ALL cases
    // Loop over all tracklets, a-b-c-d go to c then get tracklets or triplets and ask whether segment is shared
    // Ditto for Triplet -> Tracklet
    for (auto& layerPtr : getLayerPtrs())
    {

        for (auto& innerTrackletPtr : layerPtr->getTrackletPtrs())
        {
            SDL::Module& innerLowerModule = getModule(innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
            const SDL::Module& commonModule = innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule();

            for (auto& outerTrackletPtr : commonModule.getTrackletPtrs())
            {

                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateAlgo(algo, logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

            for (auto& outerTripletPtr : commonModule.getTripletPtrs())
            {

                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

                tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

        for (auto& innerTripletPtr : layerPtr->getTripletPtrs())
        {
            SDL::Module& innerLowerModule = getModule(innerTripletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
            const SDL::Module& commonModule = innerTripletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule();

            for (auto& outerTrackletPtr : commonModule.getTrackletPtrs())
            {

                SDL::TrackCandidate tcCand(innerTripletPtr, outerTrackletPtr);

                tcCand.runTrackCandidateInnerTripletToOuterTracklet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

    for (auto& innerTrackletPtr : getPixelLayer().getTrackletPtrs())
    {
        SDL::Module& innerLowerModule = getModule(innerTrackletPtr->innerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule().detId());
        const SDL::Module& commonModule = innerTrackletPtr->outerSegmentPtr()->innerMiniDoubletPtr()->lowerHitPtr()->getModule();

        for (auto& outerTrackletPtr : commonModule.getTrackletPtrs())
        {

            SDL::TrackCandidate tcCand(innerTrackletPtr, outerTrackletPtr);

            tcCand.runTrackCandidateAlgo(algo, logLevel_);

            if (tcCand.passesTrackCandidateAlgo(algo))
            {
                addTrackCandidateToLowerLayer(tcCand, 0, SDL::Layer::Barrel);
            }

        }

        for (auto& outerTripletPtr : commonModule.getTripletPtrs())
        {

            SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

            tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

            if (tcCand.passesTrackCandidateAlgo(algo))
            {
                addTrackCandidateToLowerLayer(tcCand, 0, SDL::Layer::Barrel);
            }

        }

    }

}


void SDL::Event::createTrackCandidatesFromInnerModulesFromTrackletsToTriplets(unsigned int detId, SDL::TCAlgo algo)
{

    // Get reference to the inner lower Module
    Module& innerLowerModule = getModule(detId);

    // Triple nested loops
    // Loop over inner lower module for segments
    for (auto& innerTrackletPtr : innerLowerModule.getTrackletPtrs())
    {

        // Get reference to segment in inner lower module
        SDL::Tracklet& innerTracklet = *innerTrackletPtr;

        // Get the outer mini-doublet module detId
        const SDL::Module& innerTrackletSecondModule = innerTracklet.innerSegmentPtr()->outerMiniDoubletPtr()->lowerHitPtr()->getModule();

        unsigned int innerTrackletSecondModuleDetId = innerTrackletSecondModule.detId();

        // Get connected outer lower module detids
        const std::vector<unsigned int>& connectedModuleDetIds = moduleConnectionMap.getConnectedModuleDetIds(innerTrackletSecondModuleDetId);

        // Loop over connected outer lower modules
        for (auto& outerLowerModuleDetId : connectedModuleDetIds)
        {

            if (not hasModule(outerLowerModuleDetId))
                continue;

            // Get reference to the outer lower module
            Module& outerLowerModule = getModule(outerLowerModuleDetId);

            // Loop over outer lower module mini-doublets
            for (auto& outerTripletPtr : outerLowerModule.getTripletPtrs())
            {

                SDL::Triplet& outerTriplet = *outerTripletPtr;

                SDL::TrackCandidate tcCand(innerTrackletPtr, outerTripletPtr);

                tcCand.runTrackCandidateInnerTrackletToOuterTriplet(logLevel_);

                // Count the # of track candidates considered
                incrementNumberOfTrackCandidateCandidates(innerLowerModule);

                if (tcCand.passesTrackCandidateAlgo(algo))
                {

                    // Count the # of track candidates considered
                    incrementNumberOfTrackCandidates(innerLowerModule);

                    if (innerLowerModule.subdet() == SDL::Module::Barrel)
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Barrel);
                    else
                        addTrackCandidateToLowerLayer(tcCand, innerLowerModule.layer(), SDL::Layer::Endcap);
                }

            }

        }

    }

}

// Multiplicity of mini-doublets
unsigned int SDL::Event::getNumberOfHits() { return hits_.size(); }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfHitsByLayerBarrel(unsigned int ilayer) { return n_hits_by_layer_barrel_[ilayer]; }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfHitsByLayerEndcap(unsigned int ilayer) { return n_hits_by_layer_endcap_[ilayer]; }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfHitsByLayerBarrelUpperModule(unsigned int ilayer) { return n_hits_by_layer_barrel_upper_[ilayer]; }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfHitsByLayerEndcapUpperModule(unsigned int ilayer) { return n_hits_by_layer_endcap_upper_[ilayer]; }

// Multiplicity of mini-doublets
unsigned int SDL::Event::getNumberOfMiniDoublets() { return miniDoublets_.size(); }

// Multiplicity of segments
unsigned int SDL::Event::getNumberOfSegments() { return segments_.size(); }

// Multiplicity of tracklets
unsigned int SDL::Event::getNumberOfTracklets() { return tracklets_.size(); }

// Multiplicity of triplets
unsigned int SDL::Event::getNumberOfTriplets() { return triplets_.size(); }

// Multiplicity of track candidates
unsigned int SDL::Event::getNumberOfTrackCandidates() { return trackcandidates_.size(); }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfMiniDoubletCandidates() { unsigned int n = 0; for (unsigned int i = 0; i < 6; ++i) {n += n_miniDoublet_candidates_by_layer_barrel_[i];} for (unsigned int i = 0; i < 5; ++i) {n += n_miniDoublet_candidates_by_layer_endcap_[i];} return n; }

// Multiplicity of segment candidates considered in this event
unsigned int SDL::Event::getNumberOfSegmentCandidates() { unsigned int n = 0; for (unsigned int i = 0; i < 6; ++i) {n += n_segment_candidates_by_layer_barrel_[i];} for (unsigned int i = 0; i < 5; ++i) {n += n_segment_candidates_by_layer_endcap_[i];} return n; }

// Multiplicity of tracklet candidates considered in this event
unsigned int SDL::Event::getNumberOfTrackletCandidates() { unsigned int n = 0; for (unsigned int i = 0; i < 6; ++i) {n += n_tracklet_candidates_by_layer_barrel_[i];} for (unsigned int i = 0; i < 5; ++i) {n += n_tracklet_candidates_by_layer_endcap_[i];} return n; }

// Multiplicity of triplet candidates considered in this event
unsigned int SDL::Event::getNumberOfTripletCandidates() { unsigned int n = 0; for (unsigned int i = 0; i < 6; ++i) {n += n_triplet_candidates_by_layer_barrel_[i];} for (unsigned int i = 0; i < 5; ++i) {n += n_triplet_candidates_by_layer_endcap_[i];} return n; }

// Multiplicity of track candidate candidates considered in this event
unsigned int SDL::Event::getNumberOfTrackCandidateCandidates() { unsigned int n = 0; for (unsigned int i = 0; i < 6; ++i) {n += n_trackcandidate_candidates_by_layer_barrel_[i];} for (unsigned int i = 0; i < 5; ++i) {n += n_trackcandidate_candidates_by_layer_endcap_[i];} return n; }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfMiniDoubletCandidatesByLayerBarrel(unsigned int ilayer) { return n_miniDoublet_candidates_by_layer_barrel_[ilayer]; }

// Multiplicity of segment candidates considered in this event
unsigned int SDL::Event::getNumberOfSegmentCandidatesByLayerBarrel(unsigned int ilayer) { return n_segment_candidates_by_layer_barrel_[ilayer]; }

// Multiplicity of tracklet candidates considered in this event
unsigned int SDL::Event::getNumberOfTrackletCandidatesByLayerBarrel(unsigned int ilayer) { return n_tracklet_candidates_by_layer_barrel_[ilayer]; }

// Multiplicity of triplet candidates considered in this event
unsigned int SDL::Event::getNumberOfTripletCandidatesByLayerBarrel(unsigned int ilayer) { return n_triplet_candidates_by_layer_barrel_[ilayer]; }

// Multiplicity of track candidate candidates considered in this event
unsigned int SDL::Event::getNumberOfTrackCandidateCandidatesByLayerBarrel(unsigned int ilayer) { return n_trackcandidate_candidates_by_layer_barrel_[ilayer]; }

// Multiplicity of mini-doublet candidates considered in this event
unsigned int SDL::Event::getNumberOfMiniDoubletCandidatesByLayerEndcap(unsigned int ilayer) { return n_miniDoublet_candidates_by_layer_endcap_[ilayer]; }

// Multiplicity of segment candidates considered in this event
unsigned int SDL::Event::getNumberOfSegmentCandidatesByLayerEndcap(unsigned int ilayer) { return n_segment_candidates_by_layer_endcap_[ilayer]; }

// Multiplicity of tracklet candidates considered in this event
unsigned int SDL::Event::getNumberOfTrackletCandidatesByLayerEndcap(unsigned int ilayer) { return n_tracklet_candidates_by_layer_endcap_[ilayer]; }

// Multiplicity of triplet candidates considered in this event
unsigned int SDL::Event::getNumberOfTripletCandidatesByLayerEndcap(unsigned int ilayer) { return n_triplet_candidates_by_layer_endcap_[ilayer]; }

// Multiplicity of track candidate candidates considered in this event
unsigned int SDL::Event::getNumberOfTrackCandidateCandidatesByLayerEndcap(unsigned int ilayer) { return n_trackcandidate_candidates_by_layer_endcap_[ilayer]; }

// Multiplicity of mini-doublet formed in this event
unsigned int SDL::Event::getNumberOfMiniDoubletsByLayerBarrel(unsigned int ilayer) { return n_miniDoublet_by_layer_barrel_[ilayer]; }

// Multiplicity of segment formed in this event
unsigned int SDL::Event::getNumberOfSegmentsByLayerBarrel(unsigned int ilayer) { return n_segment_by_layer_barrel_[ilayer]; }

// Multiplicity of tracklet formed in this event
unsigned int SDL::Event::getNumberOfTrackletsByLayerBarrel(unsigned int ilayer) { return n_tracklet_by_layer_barrel_[ilayer]; }

// Multiplicity of triplet formed in this event
unsigned int SDL::Event::getNumberOfTripletsByLayerBarrel(unsigned int ilayer) { return n_triplet_by_layer_barrel_[ilayer]; }

// Multiplicity of track candidate formed in this event
unsigned int SDL::Event::getNumberOfTrackCandidatesByLayerBarrel(unsigned int ilayer) { return n_trackcandidate_by_layer_barrel_[ilayer]; }

// Multiplicity of mini-doublet formed in this event
unsigned int SDL::Event::getNumberOfMiniDoubletsByLayerEndcap(unsigned int ilayer) { return n_miniDoublet_by_layer_endcap_[ilayer]; }

// Multiplicity of segment formed in this event
unsigned int SDL::Event::getNumberOfSegmentsByLayerEndcap(unsigned int ilayer) { return n_segment_by_layer_endcap_[ilayer]; }

// Multiplicity of tracklet formed in this event
unsigned int SDL::Event::getNumberOfTrackletsByLayerEndcap(unsigned int ilayer) { return n_tracklet_by_layer_endcap_[ilayer]; }

// Multiplicity of triplet formed in this event
unsigned int SDL::Event::getNumberOfTripletsByLayerEndcap(unsigned int ilayer) { return n_triplet_by_layer_endcap_[ilayer]; }

// Multiplicity of track candidate formed in this event
unsigned int SDL::Event::getNumberOfTrackCandidatesByLayerEndcap(unsigned int ilayer) { return n_trackcandidate_by_layer_endcap_[ilayer]; }

// Multiplicity of hits in this event
void SDL::Event::incrementNumberOfHits(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);

    // Only count hits in lower module
    if (not module.isLower())
    {
        if (isbarrel)
            n_hits_by_layer_barrel_upper_[layer-1]++;
        else
            n_hits_by_layer_endcap_upper_[layer-1]++;
    }
    else
    {
        if (isbarrel)
            n_hits_by_layer_barrel_[layer-1]++;
        else
            n_hits_by_layer_endcap_[layer-1]++;
    }
}

// Multiplicity of mini-doublet candidates considered in this event
void SDL::Event::incrementNumberOfMiniDoubletCandidates(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_miniDoublet_candidates_by_layer_barrel_[layer-1]++;
    else
        n_miniDoublet_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of segment candidates considered in this event
void SDL::Event::incrementNumberOfSegmentCandidates(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_segment_candidates_by_layer_barrel_[layer-1]++;
    else
        n_segment_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of triplet candidates considered in this event
void SDL::Event::incrementNumberOfTripletCandidates(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_triplet_candidates_by_layer_barrel_[layer-1]++;
    else
        n_triplet_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of tracklet candidates considered in this event
void SDL::Event::incrementNumberOfTrackletCandidates(SDL::Layer& _layer)
{
    int layer = _layer.layerIdx();
    int isbarrel = (_layer.subdet() == SDL::Layer::Barrel);
    if (isbarrel)
        n_tracklet_candidates_by_layer_barrel_[layer-1]++;
    else
        n_tracklet_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of tracklet candidates considered in this event
void SDL::Event::incrementNumberOfTrackletCandidates(SDL::Module& _module)
{
    int layer = _module.layer();
    int isbarrel = (_module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_tracklet_candidates_by_layer_barrel_[layer-1]++;
    else
        n_tracklet_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of track candidate candidates considered in this event
void SDL::Event::incrementNumberOfTrackCandidateCandidates(SDL::Layer& _layer)
{
    int layer = _layer.layerIdx();
    int isbarrel = (_layer.subdet() == SDL::Layer::Barrel);
    if (isbarrel)
        n_trackcandidate_candidates_by_layer_barrel_[layer-1]++;
    else
        n_trackcandidate_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of track candidate candidates considered in this event
void SDL::Event::incrementNumberOfTrackCandidateCandidates(SDL::Module& _module)
{
    int layer = _module.layer();
    int isbarrel = (_module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_trackcandidate_candidates_by_layer_barrel_[layer-1]++;
    else
        n_trackcandidate_candidates_by_layer_endcap_[layer-1]++;
}

// Multiplicity of mini-doublet formed in this event
void SDL::Event::incrementNumberOfMiniDoublets(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_miniDoublet_by_layer_barrel_[layer-1]++;
    else
        n_miniDoublet_by_layer_endcap_[layer-1]++;
}

// Multiplicity of segment formed in this event
void SDL::Event::incrementNumberOfSegments(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_segment_by_layer_barrel_[layer-1]++;
    else
        n_segment_by_layer_endcap_[layer-1]++;
}

// Multiplicity of triplet formed in this event
void SDL::Event::incrementNumberOfTriplets(SDL::Module& module)
{
    int layer = module.layer();
    int isbarrel = (module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_triplet_by_layer_barrel_[layer-1]++;
    else
        n_triplet_by_layer_endcap_[layer-1]++;
}

// Multiplicity of tracklet formed in this event
void SDL::Event::incrementNumberOfTracklets(SDL::Layer& _layer)
{
    int layer = _layer.layerIdx();
    int isbarrel = (_layer.subdet() == SDL::Layer::Barrel);
    if (isbarrel)
        n_tracklet_by_layer_barrel_[layer-1]++;
    else
        n_tracklet_by_layer_endcap_[layer-1]++;
}

// Multiplicity of tracklet formed in this event
void SDL::Event::incrementNumberOfTracklets(SDL::Module& _module)
{
    int layer = _module.layer();
    int isbarrel = (_module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_tracklet_by_layer_barrel_[layer-1]++;
    else
        n_tracklet_by_layer_endcap_[layer-1]++;
}

// Multiplicity of track candidate formed in this event
void SDL::Event::incrementNumberOfTrackCandidates(SDL::Layer& _layer)
{
    int layer = _layer.layerIdx();
    int isbarrel = (_layer.subdet() == SDL::Layer::Barrel);
    if (isbarrel)
        n_trackcandidate_by_layer_barrel_[layer-1]++;
    else
        n_trackcandidate_by_layer_endcap_[layer-1]++;
}

// Multiplicity of track candidate formed in this event
void SDL::Event::incrementNumberOfTrackCandidates(SDL::Module& _module)
{
    int layer = _module.layer();
    int isbarrel = (_module.subdet() == SDL::Module::Barrel);
    if (isbarrel)
        n_trackcandidate_by_layer_barrel_[layer-1]++;
    else
        n_trackcandidate_by_layer_endcap_[layer-1]++;
}

namespace SDL
{
    std::ostream& operator<<(std::ostream& out, const Event& event)
    {

        out << "" << std::endl;
        out << "==============" << std::endl;
        out << "Printing Event" << std::endl;
        out << "==============" << std::endl;
        out << "" << std::endl;

        for (auto& modulePtr : event.modulePtrs_)
        {
            out << modulePtr;
        }

        for (auto& layerPtr : event.layerPtrs_)
        {
            out << layerPtr;
        }

        return out;
    }

    std::ostream& operator<<(std::ostream& out, const Event* event)
    {
        out << *event;
        return out;
    }

}
