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

#ifndef prtCounterBox_H
#define prtCounterBox_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "presentParticle.H"
#include "immersedBody.H"

namespace Foam
{
    // Class for storing
    class prtCounterBox
    {
        private:

        typedef DynamicList<label>  DynamicLabelList;
        typedef DynamicList<scalar> DynamicScalarList;
        typedef DynamicList<vector> DynamicVectorList;
        typedef DynamicList<point>  DynamicPointList;
        
        scalar contactTime_= 0.0;
        scalar screeningTime_= 0.0;
        scalar activeTime_= 0.0;

        Tuple2<label,label> recordedContactCount_= Tuple2<label,label> (0,0);

        //screeningParameters
        scalar particleDistance_= 0.0;
        
        HashTable<autoPtr<presentParticle>,label,Hash<label>>  storedParticles_;
        labelHashSet currentParticles_;
        DynamicLabelList checkForContact_;


        
        public:
        // Constructors -------------------------------------------------
        prtCounterBox
        (
            const scalar& contactTime,
            const scalar& screeningTime,
            const scalar& particleDistance
        );
        // Destructors -----------------------------------------------------
        ~prtCounterBox(){}

        void checkStoredParticles
        (
            const immersedBody& ib,
            const scalar& timeStep
        );

        void checkStoredParticles();

        void clearData()
        {
            currentParticles_.clear();
            checkForContact_.clear();
        }

        label getParticleCount()
        {
            return storedParticles_.size();
        }

        label getCurrentContactReadyCount()
        {
            return checkForContact_.size();
        }

        void checkPresentParticles();

        void updateTimeCounter
        (
            const scalar& timeStep
        )
        {
            activeTime_ += timeStep;
        }

        scalar getActiveTime()
        {
            return activeTime_;
        }

        bool checkTimeCounter();

        void runPossibleContactScreening
        (
            PtrList<immersedBody>& immersedBodies
        );

        bool checkPossibleContact1
        (
            immersedBody& cIb,
            immersedBody& tIb
        );

        bool checkPossibleContact2
        (
            immersedBody& cIb,
            immersedBody& tIb
        );
        
        Tuple2<label,label> getContactCount()
        {
            return recordedContactCount_;
        }
    };
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace Foam
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#endif
// ************************************************************************* //
