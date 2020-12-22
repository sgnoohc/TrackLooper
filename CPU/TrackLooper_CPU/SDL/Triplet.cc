#include "Triplet.h"

SDL::Triplet::Triplet()
{
}

SDL::Triplet::~Triplet()
{
}

SDL::Triplet::Triplet(const Triplet& tl) :
    TrackletBase(tl),
    tlCand(tl.tlCand)
{
}

SDL::Triplet::Triplet(SDL::Segment* innerSegmentPtr, SDL::Segment* outerSegmentPtr) :
    TrackletBase(innerSegmentPtr, outerSegmentPtr)
{
}

bool SDL::Triplet::passesTripletAlgo(SDL::TPAlgo algo) const
{
    // Each algorithm is an enum shift it by its value and check against the flag
    return passAlgo_ & (1 << algo);
}

void SDL::Triplet::runTripletAlgo(SDL::TPAlgo algo, SDL::LogLevel logLevel)
{
    if (algo == SDL::AllComb_TPAlgo)
    {
        runTripletAllCombAlgo();
    }
    else if (algo == SDL::Default_TPAlgo)
    {

        setRecoVars("betaPt_2nd", -999);

        runTripletDefaultAlgo(logLevel);
    }
    else
    {
        SDL::cout << "Warning: Unrecognized segment algorithm!" << algo << std::endl;
        return;
    }
}

void SDL::Triplet::runTripletAllCombAlgo()
{
    passAlgo_ |= (1 << SDL::AllComb_TPAlgo);
}

void SDL::Triplet::runTripletDefaultAlgo(SDL::LogLevel logLevel)
{

    passAlgo_ &= (0 << SDL::Default_TPAlgo);

    if (not (innerSegmentPtr()->hasCommonMiniDoublet(*(outerSegmentPtr()))))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #1 in " << __FUNCTION__ << std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TPAlgo);
        return;
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::commonSegment);

    // If it does not pass pointing constraint
    if (not (passPointingConstraint(logLevel)))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #2 in " << __FUNCTION__ << std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TPAlgo);
        return;
    }

    //====================================================
    //
    // Running Tracklet algo within triplet
    //
    // if (false)
    {
        // Check tracklet algo on triplet
        tlCand = SDL::Tracklet(innerSegmentPtr(), outerSegmentPtr());

        tlCand.runTrackletDefaultAlgo(logLevel);

        if (not (tlCand.passesTrackletAlgo(SDL::Default_TLAlgo)))
        {
            if (logLevel >= SDL::Log_Debug3)
            {
                SDL::cout << "Failed Cut #3 in " << __FUNCTION__ << std::endl;
            }
            passAlgo_ &= (0 << SDL::Default_TPAlgo);
            return;
        }
        // Flag the pass bit
        passBitsDefaultAlgo_ |= (1 << TripletSelection::tracklet);
    }
    //
    //====================================================

    //====================================================
    //
    // Compute momentum of Triplet
    //
    if (false)
    {
        SDL::Hit& HitA = (*innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr());
        SDL::Hit& HitB = (*innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr());
        SDL::Hit& HitC = (*outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr());
        SDL::Hit center = SDL::MathUtil::getCenterFromThreePoints(HitA, HitB, HitC);
        float ptEst = SDL::MathUtil::ptEstimateFromRadius((HitA - center).rt());

        if (not (ptEst > 1.))
        {
            if (logLevel >= SDL::Log_Debug3)
            {
                SDL::cout << "Failed Cut #4 in " << __FUNCTION__ << std::endl;
            }
            passAlgo_ &= (0 << SDL::Default_TPAlgo);
            return;
        }
        // Flag the pass bit
        passBitsDefaultAlgo_ |= (1 << TripletSelection::tracklet);
    }
    //
    //====================================================

    passAlgo_ |= (1 << SDL::Default_TPAlgo);
}

bool SDL::Triplet::passPointingConstraint(SDL::LogLevel logLevel)
{
    // SDL::cout << innerSegmentPtr();
    // SDL::cout << outerSegmentPtr();
    // return false;
    
    const SDL::Module& ModuleA = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& ModuleB = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule();
    const SDL::Module& ModuleC = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule();

    if (ModuleA.subdet() == SDL::Module::Barrel
    and ModuleB.subdet() == SDL::Module::Barrel
    and ModuleC.subdet() == SDL::Module::Barrel)
    {
        return passPointingConstraintBarrelBarrelBarrel(logLevel);
    }
    else if (ModuleA.subdet() == SDL::Module::Barrel
         and ModuleB.subdet() == SDL::Module::Barrel
         and ModuleC.subdet() == SDL::Module::Endcap)
    {
        return passPointingConstraintBarrelBarrelEndcap(logLevel);
    }
    else if (ModuleA.subdet() == SDL::Module::Barrel
         and ModuleB.subdet() == SDL::Module::Endcap
         and ModuleC.subdet() == SDL::Module::Endcap)
    {
        return passPointingConstraintBarrelEndcapEndcap(logLevel);
    }
    else if (ModuleA.subdet() == SDL::Module::Endcap
         and ModuleB.subdet() == SDL::Module::Endcap
         and ModuleC.subdet() == SDL::Module::Endcap)
    {
        return passPointingConstraintEndcapEndcapEndcap(logLevel);
        // return false;
    }
    else
    {
        SDL::cout << ModuleA.subdet() << std::endl;
        SDL::cout << ModuleB.subdet() << std::endl;
        SDL::cout << ModuleC.subdet() << std::endl;
        SDL::cout << "WHY AM I HERE?" << std::endl;
        return false;
    }
}

bool SDL::Triplet::passPointingConstraintBarrelBarrelBarrel(SDL::LogLevel logLevel)
{


    // Z pointing between inner most MD to outer most MD

    const float deltaZLum = 15.f;
    const float kRinv1GeVf = (2.99792458e-3 * 3.8);
    const float k2Rinv1GeVf = kRinv1GeVf / 2.;
    const float ptCut = PTCUT;
    const float sinAlphaMax = 0.95; // 54.43 degrees
    const float pixelPSZpitch = 0.15;
    const float strip2SZpitch = 5.0;

    // Nomenclature
    // outer segment inner mini-doublet will be referred to as "OutLo"
    // inner segment inner mini-doublet will be referred to as "InLo"
    // outer segment outer mini-doublet will be referred to as "OutUp"
    // inner segment outer mini-doublet will be referred to as "InUp"
    // Pair of MDs in oSG will be referred to as "OutSeg"
    // Pair of MDs in iSG will be referred to as "InSeg"

    // NOTE in triplet InUp == OutLo

    const bool isPS_InLo = ((innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS);
    const bool isPS_OutLo = ((outerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS);
    const bool isPS_OutUp = ((outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS);
    const bool isEC_OutUp = ((outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).subdet() == SDL::Module::Endcap);

    const float rt_InLo = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float rt_OutUp = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float z_InLo = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->z();
    const float z_OutUp = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->z();

    const float alpha1GeV_OutUp = std::asin(std::min(rt_OutUp * k2Rinv1GeVf / ptCut, sinAlphaMax));

    const float rtRatio_OutUpInLo = rt_OutUp / rt_InLo; // Outer segment beginning rt divided by inner segment beginning rt;
    const float dzDrtScale = tan(alpha1GeV_OutUp) / alpha1GeV_OutUp; // The track can bend in r-z plane slightly
    const float zpitch_InLo = (isPS_InLo ? pixelPSZpitch : strip2SZpitch);
    const float zpitch_OutUp = (isPS_OutUp ? pixelPSZpitch : strip2SZpitch);

    // Tracklet is two segments
    // So physically it will look like the following:
    //
    // Below, the pair of x's are one segment and the pair of y's are another segment
    //
    // The x's are outer segment
    // The y's are inner segment
    //
    // rt
    //  ^
    //  |    --------------x-- ----------- --x--------------
    //  |    ---------------x- ---------|- -x--|------------ <- z_OutUp, rt_OutUp
    //  |                               <------>
    //  |                              zLo     zHi
    //  |    ----------------- -y-------y- -----------------
    //  |    ----------------- --y-----y-- -----------------
    //  |
    //  |    ----------------- ---y---y--- -----------------
    //  |    ----------------- ----y-y---- ----------------- <- z_InLo, rt_InLo
    //  |
    //  |                           *
    //  |                     <----> <----> <-- deltaZLum
    //  |
    //  |-----------------------------------------------------------------> z
    //

    // From the picture above, let's say we want to guess the zHi for example.
    //
    // Then we write down the equation
    //
    //    (zHi - z_InLo)           (z_InLo + deltaZLum)
    // ------------------------   = ---------------------------
    // (rt_OutUp - rt_InLo)            rt_InLo
    //
    // Then, if you solve for zHi you get most of the the below equation
    // Then there are two corrections.
    // dzDrtScale to account for the forward-bending of the helix tracks
    // z_pitch terms to additonally add error

    const float zHi = z_InLo + (z_InLo + deltaZLum) * (rtRatio_OutUpInLo - 1.f) * (z_InLo < 0.f ? 1.f : dzDrtScale) + (zpitch_InLo + zpitch_OutUp);
    const float zLo = z_InLo + (z_InLo - deltaZLum) * (rtRatio_OutUpInLo - 1.f) * (z_InLo > 0.f ? 1.f : dzDrtScale) - (zpitch_InLo + zpitch_OutUp); //slope-correction only on outer end

    // Reset passBitsDefaultAlgo_;
    passBitsDefaultAlgo_ = 0;

    //==========================================================================
    //
    // Cut #1: Z compatibility
    //
    //==========================================================================
    setRecoVars("z_OutUp", z_OutUp);
    setRecoVars("zLo", zLo);
    setRecoVars("zHi", zHi);
    if (not (z_OutUp >= zLo and z_OutUp <= zHi))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #2 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " zLo: " << zLo <<  " z_OutUp: " << z_OutUp <<  " zHi: " << zHi <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::deltaZ);
    //--------------------------------------------------------------------------

    const float drt_OutUp_InLo = (rt_OutUp - rt_InLo);
    const float invRt_InLo = 1. / rt_InLo;
    const float r3_InLo = sqrt(z_InLo * z_InLo + rt_InLo * rt_InLo);
    const float drt_InSeg = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->rt() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float dz_InSeg  = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->z() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->z();
    const float dr3_InSeg = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->r3() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->r3();

    // The x's are outer segment
    // The y's are inner segment
    //
    // rt
    //  ^
    //  |    --------------x-- ----------- --x--------------
    //  |    ---------------x- ---------|- -x--|------------ <- z_OutUp, rt_OutUp
    //  |                               <------>
    //  |                       zLoPointed  ^  zHiPointed
    //  |                                   |
    //  |                                 dzMean
    //  |    ----------------- -y-------y- -----------------
    //  |    ----------------- --y-----y-- ----------------- <-|
    //  |                                                      | "InSeg"
    //  |    ----------------- ---y---y--- -----------------   | dz, drt, tl_axis is defined above
    //  |    ----------------- ----y-y---- ----------------- <-|
    //  |
    //  |
    //  |-----------------------------------------------------------------> z
    //
    //
    // We point it via the inner segment

    // direction estimate
    const float coshEta = dr3_InSeg / drt_InSeg;

    // Defining the error terms along the z-direction
    float dzErr = (zpitch_InLo + zpitch_OutUp) * (zpitch_InLo + zpitch_OutUp) * 2.f; //both sides contribute to direction uncertainty

    // Multiple scattering
    //FIXME (later) more realistic accounting of material effects is needed
    const float sdlThetaMulsF = 0.015f * sqrt(0.1f + 0.2 * (rt_OutUp - rt_InLo) / 50.f) * sqrt(r3_InLo / rt_InLo);
    const float sdlMuls = sdlThetaMulsF * 3.f / ptCut * 4.f; // will need a better guess than x4?
    dzErr += sdlMuls * sdlMuls * drt_OutUp_InLo * drt_OutUp_InLo / 3.f * coshEta * coshEta; //sloppy
    dzErr = sqrt(dzErr);

    // Constructing upper and lower bound
    const float dzMean = dz_InSeg / drt_InSeg * drt_OutUp_InLo;
    const float zWindow = dzErr / drt_InSeg * drt_OutUp_InLo + (zpitch_InLo + zpitch_OutUp); //FIXME for ptCut lower than ~0.8 need to add curv path correction
    const float zLoPointed = z_InLo + dzMean * (z_InLo > 0.f ? 1.f : dzDrtScale) - zWindow;
    const float zHiPointed = z_InLo + dzMean * (z_InLo < 0.f ? 1.f : dzDrtScale) + zWindow;

    //==========================================================================
    //
    // Cut #2: Pointed Z (Inner segment two MD points to outer segment inner MD)
    //
    //==========================================================================
    setRecoVars("zLoPointed", zLoPointed);
    setRecoVars("zHiPointed", zHiPointed);
    if (not (z_OutUp >= zLoPointed and z_OutUp <= zHiPointed))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #3 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " zLoPointed: " << zLoPointed <<  " z_OutUp: " << z_OutUp <<  " zHiPointed: " << zHiPointed <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::deltaZPointed);
    //--------------------------------------------------------------------------

    return true;

}

bool SDL::Triplet::passPointingConstraintBarrelBarrelEndcap(SDL::LogLevel logLevel)
{
    const float deltaZLum = 15.f;
    const float pixelPSZpitch = 0.15;
    const float strip2SZpitch = 5.0;
    const float zGeom =
        ((innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS ? pixelPSZpitch : strip2SZpitch)
        +
        ((outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS ? pixelPSZpitch : strip2SZpitch);
    const float rtIn = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float rtOut = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float zIn = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->z();
    const float zOut = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->z();
    const float rtOut_o_rtIn = rtOut / rtIn;
    const float kRinv1GeVf = (2.99792458e-3 * 3.8);
    const float k2Rinv1GeVf = kRinv1GeVf / 2.;
    const float ptCut = PTCUT;
    const float sinAlphaMax = 0.95;
    const float sdlSlope = std::asin(std::min(rtOut * k2Rinv1GeVf / ptCut, sinAlphaMax));
    const float dzDrtScale = tan(sdlSlope) / sdlSlope;//FIXME: need approximate value
    const float zLo = zIn + (zIn - deltaZLum) * (rtOut_o_rtIn - 1.f) * (zIn > 0.f ? 1.f : dzDrtScale) - zGeom; //slope-correction only on outer end

    // Cut #0: Preliminary
    // (Only here in endcap case)
    if (not (zIn * zOut > 0))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #0 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " zIn: " << zIn <<  std::endl;
            SDL::cout <<  " zOut: " << zOut <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::deltaZ);


    const float dLum = std::copysign(deltaZLum, zIn);
    const bool isOutSgInnerMDPS = (outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS;
    const float rtGeom1 = isOutSgInnerMDPS ? pixelPSZpitch : strip2SZpitch;//FIXME: make this chosen by configuration for lay11,12 full PS
    const float zGeom1 = std::copysign(zGeom, zIn); //used in B-E region
    const float rtLo = rtIn * (1.f + (zOut - zIn - zGeom1) / (zIn + zGeom1 + dLum) / dzDrtScale) - rtGeom1; //slope correction only on the lower end

    // Cut #1: rt condition
    if (not (rtOut >= rtLo))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #1 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " rtOut: " << rtOut <<  std::endl;
            SDL::cout <<  " rtLo: " << rtLo <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }

    float zInForHi = zIn - zGeom1 - dLum;
    if (zInForHi * zIn < 0)
        zInForHi = std::copysign(0.1f, zIn);
    const float rtHi = rtIn * (1.f + (zOut - zIn + zGeom1) / zInForHi) + rtGeom1;

    // Cut #2: rt condition
    if (not (rtOut >= rtLo and rtOut <= rtHi))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #2 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " rtOut: " << rtOut <<  std::endl;
            SDL::cout <<  " rtLo: " << rtLo <<  std::endl;
            SDL::cout <<  " rtHi: " << rtHi <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }

    const float rIn = sqrt(zIn*zIn + rtIn*rtIn);
    const float drtSDIn = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->rt() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float dzSDIn = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->z() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->z();
    const float dr3SDIn = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->r3() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->r3();
    const float coshEta = dr3SDIn / drtSDIn; //direction estimate
    const float dzOutInAbs = std::abs(zOut - zIn);
    const float multDzDr = dzOutInAbs * coshEta / (coshEta * coshEta - 1.f);
    const float zGeom1_another = pixelPSZpitch; // TODO-Q Why only one?
    const float kZ = (zOut - zIn) / dzSDIn;
    float drtErr = zGeom1_another * zGeom1_another * drtSDIn * drtSDIn / dzSDIn / dzSDIn * (1.f - 2.f * kZ + 2.f * kZ * kZ); //Notes:122316
    const float sdlThetaMulsF = 0.015f * sqrt(0.1f + 0.2 * (rtOut - rtIn) / 50.f) * sqrt(rIn / rtIn);
    const float sdlMuls = sdlThetaMulsF * 3.f / ptCut * 4.f; //will need a better guess than x4?
    drtErr += sdlMuls * sdlMuls * multDzDr * multDzDr / 3.f * coshEta * coshEta; //sloppy: relative muls is 1/3 of total muls
    drtErr = sqrt(drtErr);
    const float drtMean = drtSDIn * dzOutInAbs / std::abs(dzSDIn); //
    const float rtWindow = drtErr + rtGeom1;
    const float rtLo_another = rtIn + drtMean / dzDrtScale - rtWindow;
    const float rtHi_another = rtIn + drtMean + rtWindow;

    // Cut #3: rt-z pointed
    if (not (kZ >= 0 and rtOut >= rtLo and rtOut <= rtHi))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #3 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " kZ: " << kZ <<  std::endl;
            SDL::cout <<  " rtOut: " << rtOut <<  std::endl;
            SDL::cout <<  " rtLo: " << rtLo <<  std::endl;
            SDL::cout <<  " rtHi: " << rtHi <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::deltaZPointed);

    return true;
}

bool SDL::Triplet::passPointingConstraintBarrelEndcapEndcap(SDL::LogLevel logLevel)
{
    return passPointingConstraintBarrelBarrelEndcap(logLevel);
}

bool SDL::Triplet::passPointingConstraintEndcapEndcapEndcap(SDL::LogLevel logLevel)
{
    const float deltaZLum = 15.f;
    const float pixelPSZpitch = 0.15;
    const float strip2SZpitch = 5.0;
    const float zGeom =
        ((innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS ? pixelPSZpitch : strip2SZpitch)
        +
        ((outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS ? pixelPSZpitch : strip2SZpitch);
    const float rtIn = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float rtOut = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float zIn = innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->z();
    const float zOut = outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->z();
    const float rtOut_o_rtIn = rtOut / rtIn;
    const float kRinv1GeVf = (2.99792458e-3 * 3.8);
    const float k2Rinv1GeVf = kRinv1GeVf / 2.;
    const float ptCut = PTCUT;
    const float sinAlphaMax = 0.95;
    const float sdlSlope = std::asin(std::min(rtOut * k2Rinv1GeVf / ptCut, sinAlphaMax));
    const float dzDrtScale = tan(sdlSlope) / sdlSlope;//FIXME: need approximate value
    const float zLo = zIn + (zIn - deltaZLum) * (rtOut_o_rtIn - 1.f) * (zIn > 0.f ? 1.f : dzDrtScale) - zGeom; //slope-correction only on outer end

    // Cut #0: Preliminary
    // (Only here in endcap case)
    // https://github.com/slava77/cms-tkph2-ntuple/blob/superDoubletLinked-91X-noMock/doubletAnalysis.C#L3631-L3633
    if (not (zIn * zOut > 0))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #0 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " zIn: " << zIn <<  std::endl;
            SDL::cout <<  " zOut: " << zOut <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }


    const float dLum = std::copysign(deltaZLum, zIn);
    const bool isOutSgInnerMDPS = (outerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS;
    const bool isInSgInnerMDPS = (innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS;
    // https://github.com/slava77/cms-tkph2-ntuple/blob/superDoubletLinked-91X-noMock/doubletAnalysis.C#L3670-L3674
    // we're in mockMode == 3
    const float rtGeom = (isInSgInnerMDPS && isOutSgInnerMDPS ? 2.f * pixelPSZpitch
                 : (isInSgInnerMDPS || isOutSgInnerMDPS ) ? (pixelPSZpitch + strip2SZpitch)
                            : 2.f * strip2SZpitch);
    const float rtGeom1 = isOutSgInnerMDPS ? pixelPSZpitch : strip2SZpitch;//FIXME: make this chosen by configuration for lay11,12 full PS
    const float zGeom1 = std::copysign(zGeom, zIn); //used in B-E region
    const float dz = zOut - zIn;
    const float rtLo = rtIn * (1.f + dz / (zIn + dLum) / dzDrtScale) - rtGeom; //slope correction only on the lower end

    // Cut #1: rt condition
    // https://github.com/slava77/cms-tkph2-ntuple/blob/superDoubletLinked-91X-noMock/doubletAnalysis.C#L3679
    if (not (rtOut >= rtLo))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #1 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " rtOut: " << rtOut <<  std::endl;
            SDL::cout <<  " rtLo: " << rtLo <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }

    const float rtHi = rtIn * (1.f + dz / (zIn - dLum)) + rtGeom;

    // Cut #2: rt condition
    if (not (rtOut >= rtLo and rtOut <= rtHi))
    {
        if (logLevel >= SDL::Log_Debug3)
        {
            SDL::cout << "Failed Cut #2 in " << __FUNCTION__ << std::endl;
            SDL::cout <<  " rtOut: " << rtOut <<  std::endl;
            SDL::cout <<  " rtLo: " << rtLo <<  std::endl;
            SDL::cout <<  " rtHi: " << rtHi <<  std::endl;
        }
        passAlgo_ &= (0 << SDL::Default_TLAlgo);
        return false;
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::deltaZ);

    const bool isInSgOuterMDPS = (innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->getModule()).moduleType() == SDL::Module::PS;

    const float drOutIn = (rtOut - rtIn);
    const float drtSDIn = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->rt() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->rt();
    const float dzSDIn = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->z() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->z();
    const float dr3SDIn = innerSegmentPtr()->outerMiniDoubletPtr()->anchorHitPtr()->r3() - innerSegmentPtr()->innerMiniDoubletPtr()->anchorHitPtr()->r3();

    const float coshEta = dr3SDIn / drtSDIn; //direction estimate
    const float dzOutInAbs = std::abs(zOut - zIn);
    const float multDzDr = dzOutInAbs * coshEta / (coshEta * coshEta - 1.f);

    const float kZ = (zOut - zIn) / dzSDIn;

    const float sdlThetaMulsF = 0.015f * sqrt(0.1f + 0.2 * (rtOut - rtIn) / 50.f);
    const float sdlMuls = sdlThetaMulsF * 3.f / ptCut * 4.f; //will need a better guess than x4?

    float drtErr = pixelPSZpitch * pixelPSZpitch * 2.f / dzSDIn / dzSDIn * dzOutInAbs * dzOutInAbs; //both sides contribute to direction uncertainty
    drtErr += sdlMuls * sdlMuls * multDzDr * multDzDr / 3.f * coshEta * coshEta; //sloppy: relative muls is 1/3 of total muls
    drtErr = sqrt(drtErr);
    const float drtMean = drtSDIn * dzOutInAbs / std::abs(dzSDIn);
    const float rtWindow = drtErr + rtGeom; //
    const float rtLo_point = rtIn + drtMean / dzDrtScale - rtWindow;
    const float rtHi_point = rtIn + drtMean + rtWindow;

    // Cut #3: rt-z pointed
    // https://github.com/slava77/cms-tkph2-ntuple/blob/superDoubletLinked-91X-noMock/doubletAnalysis.C#L3765
    if (isInSgInnerMDPS and isInSgOuterMDPS) // If both PS then we can point
    {
        if (not (kZ >= 0 and rtOut >= rtLo_point and rtOut <= rtHi_point))
        {
            if (logLevel >= SDL::Log_Debug3)
            {
                SDL::cout << "Failed Cut #3 in " << __FUNCTION__ << std::endl;
                SDL::cout <<  " kZ: " << kZ <<  std::endl;
                SDL::cout <<  " rtOut: " << rtOut <<  std::endl;
                SDL::cout <<  " rtLo: " << rtLo <<  std::endl;
                SDL::cout <<  " rtHi: " << rtHi <<  std::endl;
            }
            passAlgo_ &= (0 << SDL::Default_TLAlgo);
            return false;
        }
    }
    // Flag the pass bit
    passBitsDefaultAlgo_ |= (1 << TripletSelection::deltaZPointed);
    return true;
}
