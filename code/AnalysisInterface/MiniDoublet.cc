# include "MiniDoublet.h"

SDL::MiniDoublet::MiniDoublet(float dz, float drt, float dphi, float dphichange, float dphinoshift, float dphichangenoshift, float dzCut, float drtCut, float miniCut, std::shared_ptr<SDL::Hit>& lowerHitPtr, std::shared_ptr<SDL::Hit>& upperHitPtr)
{
    dz_ = dz;
    drt_ = drt;
    dphi_ = dphi;
    dphichange_ = dphichange;
    dphinoshift_ = dphinoshift;
    dphichangenoshift_ = dphichangenoshift;
    lowerHitPtr_ = lowerHitPtr;
    upperHitPtr_ = upperHitPtr;
    dzCut_  = dzCut;
    drtCut_ = drtCut;
    drt_ = drt;
    miniCut_ = miniCut;
    setAnchorHit();
}

const std::shared_ptr<SDL::Hit>& SDL::MiniDoublet::lowerHitPtr() const
{
    return lowerHitPtr_;
}

const std::shared_ptr<SDL::Hit>& SDL::MiniDoublet::upperHitPtr() const
{
    return upperHitPtr_;
}

const std::shared_ptr<SDL::Hit>& SDL::MiniDoublet::anchorHitPtr() const
{
    return anchorHitPtr_;
}

const float& SDL::MiniDoublet::getDz() const
{
    return dz_;
}

const float& SDL::MiniDoublet::getDeltaPhi() const
{
    return dphi_;
}

const float& SDL::MiniDoublet::getDeltaPhiChange() const
{
    return dphichange_;
}

const float& SDL::MiniDoublet::getDeltaPhiNoShift() const
{
    return dphinoshift_;
}

const float& SDL::MiniDoublet::getDeltaPhiChangeNoShift() const
{
    return dphichangenoshift_;
}

const float& SDL::MiniDoublet::getDzCut() const
{
    return dzCut_;
}

const float& SDL::MiniDoublet::getDrtCut() const
{
    return drtCut_;
}

const float& SDL::MiniDoublet::getDrt() const
{
    return drt_;
}

const float& SDL::MiniDoublet::getMiniCut() const
{
    return miniCut_;
}

void SDL::MiniDoublet::setAnchorHit()
{
    const SDL::Module& lowerModule = lowerHitPtr()->getModule();

    // Assign anchor hit pointers based on their hit type
    if (lowerModule.moduleType() == SDL::Module::PS)
    {
        if (lowerModule.moduleLayerType() == SDL::Module::Pixel)
        {
            anchorHitPtr_ = lowerHitPtr();
        }
        else
        {
            anchorHitPtr_ = upperHitPtr();
        }
    }
    else
    {
        anchorHitPtr_ = lowerHitPtr();
    }
}

SDL::MiniDoublet::~MiniDoublet()
{
}
