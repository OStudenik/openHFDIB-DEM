/*---------------------------------------------------------------------------*\
                        _   _ ____________ ___________    ______ ______ _    _
                       | | | ||  ___|  _  \_   _| ___ \   |  _  \|  ___| \  / |
  ___  _ __   ___ _ __ | |_| || |_  | | | | | | | |_/ /   | | | || |_  |  \/  |
 / _ \| '_ \ / _ \ '_ \|  _  ||  _| | | | | | | | ___ \---| | | ||  _| | |\/| |
| (_) | |_) |  __/ | | | | | || |   | |/ / _| |_| |_/ /---| |/ / | |___| |  | |
 \___/| .__/ \___|_| |_\_| |_/\_|   |___/  \___/\____/    |___/  |_____|_|  |_|
      | |                     H ybrid F ictitious D omain - I mmersed B oundary
      |_|                                        and D iscrete E lement M ethod
-------------------------------------------------------------------------------
License

    openHFDIB-DEM is licensed under the GNU LESSER GENERAL PUBLIC LICENSE (LGPL).

    Everyone is permitted to copy and distribute verbatim copies of this license
    document, but changing it is not allowed.

    This version of the GNU Lesser General Public License incorporates the terms
    and conditions of version 3 of the GNU General Public License, supplemented
    by the additional permissions listed below.

    You should have received a copy of the GNU Lesser General Public License
    along with openHFDIB. If not, see <http://www.gnu.org/licenses/lgpl.html>.

InNamspace
    

Description
    Static class for storing contact zones for contact counterss

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "prtCounterBox.H"

using namespace Foam;

//---------------------------------------------------------------------------//
prtCounterBox::prtCounterBox
(
    const scalar& contactTime,
    const scalar& screeningTime,
    const scalar& particleDistance
)
:
contactTime_(contactTime),
screeningTime_(screeningTime),
activeTime_(0),
storedParticles_(0)
{}
//---------------------------------------------------------------------------//
void prtCounterBox::checkStoredParticles(const immersedBody& ib, const scalar& timeStep)
{
    currentParticles_.insert(ib.getBodyId());
    if(!storedParticles_.found(ib.getBodyId()))
    {
        storedParticles_.insert(ib.getBodyId(), autoPtr<presentParticle>(new presentParticle));
        storedParticles_[ib.getBodyId()]->presentTime += timeStep;
        storedParticles_[ib.getBodyId()]->id = ib.getBodyId();
    }
    else
    {
        storedParticles_[ib.getBodyId()]->presentTime += timeStep;
        if(storedParticles_[ib.getBodyId()]->presentTime > contactTime_)
        {
            storedParticles_[ib.getBodyId()]->compContact = false;
        }
        else
        {
            storedParticles_[ib.getBodyId()]->compContact = true;
            checkForContact_.append(ib.getBodyId());
        }
    }
}
//---------------------------------------------------------------------------//
void prtCounterBox::checkPresentParticles()
{
    for (auto stored :storedParticles_.toc())
    {
        if(!currentParticles_.found(stored))
        {
            storedParticles_.erase(stored);
        }        
    }
}
//---------------------------------------------------------------------------//
bool prtCounterBox::checkTimeCounter()
{
    if(activeTime_ > screeningTime_)
    {
        activeTime_ = 0;
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
void prtCounterBox::runPossibleContactScreening(PtrList<immersedBody>& immersedBodies)
{
    if(checkForContact_.size() > 2)// it has to be pair
    {
        //createPossibleContactPairs        
        DynamicList<Tuple2<label,label>> possibleContactPairs;
        forAll(checkForContact_,i)
        {
            forAll(checkForContact_,j)
            {
                if(i != j && i < j)
                {
                    possibleContactPairs.append(Tuple2<label,label>(checkForContact_[i],checkForContact_[j]));
                }
            }
        }
        for(auto cPair : possibleContactPairs)
        {
            if(checkPossibleContact1(immersedBodies[cPair.first()],immersedBodies[cPair.second()]))
            {
                recordedContactCount_.first()++;
            }
            if(checkPossibleContact2(immersedBodies[cPair.first()],immersedBodies[cPair.second()]))
            {
                recordedContactCount_.second()++;
            }
        }
    }
}
//---------------------------------------------------------------------------//
bool prtCounterBox::checkPossibleContact1
(
    immersedBody& cIb, 
    immersedBody& tIb
)
{

    scalar prtDist(mag(cIb.getGeomModel().getCoM() - tIb.getGeomModel().getCoM())-(cIb.getGeomModel().getDC()/2 + tIb.getGeomModel().getDC()/2));
    vector particleNormal(cIb.getGeomModel().getCoM() - tIb.getGeomModel().getCoM());
    particleNormal /= mag(particleNormal);
    scalar relVel(mag(cIb.getVel() - tIb.getVel()));
    scalar relVelAngle((cIb.getAxis() - tIb.getAxis()) & particleNormal);
    if(prtDist < particleDistance_ && relVelAngle < 0 && relVel > particleDistance_/screeningTime_)
    {
        return true;
    }
    return false;
}//---------------------------------------------------------------------------//
bool prtCounterBox::checkPossibleContact2
(
    immersedBody& cIb, 
    immersedBody& tIb
)

{
    scalar prtDist(mag(cIb.getGeomModel().getCoM() - tIb.getGeomModel().getCoM())-(cIb.getGeomModel().getDC()/2 + tIb.getGeomModel().getDC()/2));
    vector particleNormal(cIb.getGeomModel().getCoM() - tIb.getGeomModel().getCoM());
    particleNormal /= mag(particleNormal);
    scalar relVel(mag(cIb.getVel() - tIb.getVel()));
    scalar relVelAngle((cIb.getAxis() - tIb.getAxis()) & particleNormal);
    if(prtDist < particleDistance_ && relVelAngle < 0 && relVel > particleDistance_/screeningTime_)
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//