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
    const scalar& particleDistance,
    const vector& studyPlaneNormal
)
:
contactTime_(contactTime),
screeningTime_(screeningTime),
particleDistance_(particleDistance),
studyPlaneNormal_(studyPlaneNormal),
activeTime_(0),
storedParticles_(0)
{}
//---------------------------------------------------------------------------//
void prtCounterBox::checkStoredParticles(immersedBody& ib, const scalar& timeStep)
{
    currentParticles_.insert(ib.getBodyId());
    if(!storedParticles_.found(ib.getBodyId()))
    {
        storedParticles_.insert(ib.getBodyId(), autoPtr<presentParticle>(new presentParticle));
        storedParticles_[ib.getBodyId()]->presentTime += timeStep;
        storedParticles_[ib.getBodyId()]->id = ib.getBodyId();
        storedParticles_[ib.getBodyId()]->oldPos = ib.getGeomModel().getCoM();
    }
    else
    {
        storedParticles_[ib.getBodyId()]->presentTime += timeStep;

        if(storedParticles_[ib.getBodyId()]->presentTime > contactTime_)
        {
            storedParticles_[ib.getBodyId()]->compContact = true;
            checkForContact_.append(ib.getBodyId());
            if(checkScreeningTime() && storedParticles_[ib.getBodyId()]->presentTime > screeningTime_)
            {
                storedParticles_[ib.getBodyId()]->velocity = (ib.getGeomModel().getCoM() - storedParticles_[ib.getBodyId()]->oldPos)/screeningTime_;
                storedParticles_[ib.getBodyId()]->oldPos = ib.getGeomModel().getCoM();
                checkForContact2_.append(ib.getBodyId());
            }            
        }
        else
        {
            storedParticles_[ib.getBodyId()]->compContact = false;
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
    if(activeTime_ >= screeningTime_)
    {
        activeTime_ = 0;
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
void prtCounterBox::runPossibleContactScreening(PtrList<immersedBody>& immersedBodies)
{
    Info << "-- contact Filter -> Number of contact screening particles: " << checkForContact_.size() << endl;
    if(checkForContact_.size() > 2)// it has to be pair
    {
        //createPossibleContactPairs        
        DynamicList<Tuple2<label,label>> possibleContactPairs;
        forAll(checkForContact_,i)
        {
            forAll(checkForContact_,j)
            {
                if(i != j && checkForContact_[i] < checkForContact_[j])
                {
                    possibleContactPairs.append(Tuple2<label,label>(checkForContact_[i],checkForContact_[j]));
                }
            }
        }

        // Info << "-- contact Filter -> Number of possible contact pairs: " << possibleContactPairs.size() << endl;

        for(auto cPair : possibleContactPairs)
        {
            // Info << "-- contact Filter -> Checking possible contact pair: " << cPair.first() << " " << cPair.second() << endl;
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
void prtCounterBox::runPossibleContactScreening2(PtrList<immersedBody>& immersedBodies)
{

    Info << "-- contact Filter type2 -> Number of contact screening particles: " << checkForContact2_.size() << endl;
    if(checkForContact2_.size() > 2)// it has to be pair
    {
        //createPossibleContactPairs        
        DynamicList<Tuple2<label,label>> possibleContactPairs;
        forAll(checkForContact2_,i)
        {
            forAll(checkForContact2_,j)
            {
                if(i != j && checkForContact2_[i] < checkForContact2_[j])
                {
                    possibleContactPairs.append(Tuple2<label,label>(checkForContact2_[i],checkForContact2_[j]));
                }
            }
        }

        for(auto cPair : possibleContactPairs)
        {
            // Info << "-- contact Filter -> Checking possible contact pair: " << cPair.first() << " " << cPair.second() << endl;
            if(checkPossibleContact3(immersedBodies[cPair.first()],immersedBodies[cPair.second()]))
            {
                countedContacts_++;
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
    // scalar prtDistPlane(mag(projectToPlane(cIb.getGeomModel().getCoM()) - projectToPlane(tIb.getGeomModel().getCoM()))-(cIb.getGeomModel().getDC()/2 + tIb.getGeomModel().getDC()/2));
    vector particleNormal(cIb.getGeomModel().getCoM() - tIb.getGeomModel().getCoM());
    particleNormal /= mag(particleNormal);
    scalar relVel(mag(cIb.getVel() - tIb.getVel()));
    scalar relVelAngle((cIb.getAxis() - tIb.getAxis()) & particleNormal);
    
    // scalar relVelMag(mag(projectToPlane(cIb.getVel() - tIb.getVel())));
    // scalar relVelAnglePlane(projectToPlane((cIb.getVel() - tIb.getVel())) & projectToPlane(particleNormal));
    // bool condition1 = prtDistPlane < particleDistance_;
    // bool condition2 = relVelAnglePlane < 0;
    // bool condition3 = relVelMag > particleDistance_/screeningTime_;

    if(prtDist < particleDistance_ && relVelAngle < 0 && relVel > particleDistance_/screeningTime_)
    {
        // Info << "-- contact Filter 1 -> contact pair: " << cIb.getBodyId() << "-" << tIb.getBodyId() << endl;
        // Info << "-- contact Filter 1 -> prtDist: " << prtDist << " prtDistPlane " << prtDistPlane << " Status : " << condition1 << endl;
        // Info << "-- contact Filter 1 -> relVelAngle: " << relVelAngle << " relVelAnglePlane " << relVelAnglePlane << " Status : " << condition2 << endl;
        // Info << "-- contact Filter 1 -> relVelMag: " << relVel << " relVelMagPlane " << relVelMag  << " Status : " << condition3<< endl;
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
bool prtCounterBox::checkPossibleContact2
(
    immersedBody& cIb, 
    immersedBody& tIb
)
{
    scalar prtDist(mag(projectToPlane(cIb.getGeomModel().getCoM()) - projectToPlane(tIb.getGeomModel().getCoM()))-(cIb.getGeomModel().getDC()/2 + tIb.getGeomModel().getDC()/2));
    vector particleNormal(cIb.getGeomModel().getCoM() - tIb.getGeomModel().getCoM());
    particleNormal /= mag(particleNormal);
    vector relVel(cIb.getVel() - tIb.getVel());
    scalar relVelMag(mag(projectToPlane(relVel)));
    scalar relVelAngle(projectToPlane(relVel) & projectToPlane(particleNormal));
    // scalar relVelOutPlane(relVel & particleNormal);

    if(prtDist < particleDistance_ && relVelAngle < 0 && relVelMag > particleDistance_/screeningTime_)
    {
        return true;
    }
    return false;
}
//---------------------------------------------------------------------------//
// bool prtCounterBox::checkPossibleContact3
// (
//     immersedBody& cIb, 
//     immersedBody& tIb
// )
// {
//     scalar prtDist = mag(projectToPlane(cIb.getGeomModel().getCoM()) - projectToPlane(tIb.getGeomModel().getCoM()))-(cIb.getGeomModel().getDC()/2 + tIb.getGeomModel().getDC()/2);
//     vector particleNormal(projectToPlane(cIb.getGeomModel().getCoM()) - projectToPlane(tIb.getGeomModel().getCoM()));
//     particleNormal /= mag(particleNormal);
//     vector relVel(projectToPlane(storedParticles_[cIb.getBodyId()]->velocity) - projectToPlane(storedParticles_[tIb.getBodyId()]->velocity));
//     scalar relVelMag(mag(projectToPlane(relVel)));
//     scalar relVelAngle(relVel & particleNormal);

//     bool condition1 = prtDist < particleDistance_;
//     bool condition2 = relVelAngle < 0;
//     bool condition3 = relVelMag > particleDistance_/screeningTime_;

//     if(prtDist < particleDistance_ && relVelAngle < 0 && relVelMag > particleDistance_/screeningTime_)
//     {
//         // Info << "-- contact Filter 3 -> contact pair: " << cIb.getBodyId() << "-" << tIb.getBodyId() << endl;
//         // Info << "-- contact Filter 3 -> prtDist: " << prtDist << " Status : " << condition1 << endl;
//         // Info << "-- contact Filter 3 -> relVelAngle: " << relVelAngle << " Status : " << condition2 << endl;
//         // Info << "-- contact Filter 3 -> relVelMag: " << relVelMag  << " Status : " << condition3<< endl;
//         return true;
//     }
//     return false;
// }
//---------------------------------------------------------------------------//
bool prtCounterBox::checkPossibleContact3
(
      immersedBody& cIb, 
      immersedBody& tIb
)
{
    //in-plane projected CoMs
    vector cIbCoM(projectToPlane(cIb.getGeomModel().getCoM()));
    vector tIbCoM(projectToPlane(tIb.getGeomModel().getCoM()));

    //distance between in-plane projected particles (using full radii)
    scalar prtDist(
        mag(cIbCoM - tIbCoM) - 0.5*(cIb.getGeomModel().getDC() + tIb.getGeomModel().getDC())
    ); 

    //particle normal between in-plane projected particles 
    vector particleNormal(cIbCoM - tIbCoM); 
    particleNormal /= mag(particleNormal);

    //in-plane projected velocities
    // vector cIbVel(projectToPlane(cIb.getVel()));
    // vector tIbVel(projectToPlane(tIb.getVel()));
    vector cIbVel(projectToPlane(storedParticles_[cIb.getBodyId()]->velocity));
    vector tIbVel(projectToPlane(storedParticles_[tIb.getBodyId()]->velocity));
    // Note (MI): Ondra, here, we should plug-in the corrected particle velocity computation (based on the experiment) 

    //relative velocity (after projection)
    vector relVel(cIbVel-tIbVel);
    scalar relVelMag(mag(relVel));

    //angle between relative velocity and particle normal (I already work with projected data)
    scalar relVelAngle(relVel & particleNormal);
    relVelAngle /= mag(relVel);
 // Note (MI): particleNormal is a unit vector, relVel is not. This is just a normalization for the cases of extremely small particles  
            
      if(prtDist < particleDistance_ && relVelAngle < 0 && relVelMag > particleDistance_/screeningTime_)
      {
        return true;
      }
      return false;
}

