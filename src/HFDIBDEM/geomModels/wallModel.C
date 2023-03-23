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
    Foam

Contributors
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*),
    Ondřej Studeník (2020-*)
\*---------------------------------------------------------------------------*/
#include "wallModel.H"

using namespace Foam;

//---------------------------------------------------------------------------//
wallModel::wallModel(List<string> contactPatches)
{
    forAll(contactPatches, i)
    {
        contactPlanes_.append(
            plane(
                wallPlaneInfo::getWallPlaneInfo()[wallName][1],
                wallPlaneInfo::getWallPlaneInfo()[wallName][0]
            )
        );
    }
}

wallModel::~wallModel()
{
}
//---------------------------------------------------------------------------//
volumeType wallModel::getPlaneVolumeType
(
    const plane& pl,
    const subVolume& sv
)
{
    tmp<pointField> points = sv.points();

    bool normalSide = false;
    bool flipSide = false;

    forAll(points, pI)
    {
        if (pl.sideOfPlane(points[pI])
            ==
            plane::side::NORMAL)
        {
            normalSide = true;
        }
        else
        {
            flipSide = true;
        }
    }

    if (normalSide && flipSide)
    {
        return volumeType::mixed;
    }
    else if (normalSide)
    {
        return volumeType::inside;
    }

    return volumeType::outside;
}
//---------------------------------------------------------------------------//
volumeType wallModel::getVolumeType
(
    subVolume& sv,
    bool cIb
)
{
    tmp<pointField> points = sv.points();
    List<volumeType> planesVT =
        List<volumeType>(points().size(), volumeType::unknown);

    forAll(contactPlanes_, cI)
    {
        planesVT[cI] = getPlaneVolumeType(contactPlanes_[cI], sv);

        if (planesVT[cI] == volumeType::inside)
        {
            return volumeType::inside;
        }
    }

    forAll(planesVT, pI)
    {
        if (planesVT[pI] == volumeType::mixed)
        {
            return volumeType::mixed;
        }
    }

    return volumeType::outside;
}
//---------------------------------------------------------------------------//
bool wallModel::limitFinalSubVolume
(
    const subVolume& sv,
    bool cIb,
    boundBox& limBBox
)
{
    limBBox = boundBox(sv.min(), sv.max());

    forAll(contactPlanes_, cI)
    {
        plane& planeI = contactPlanes_[cI];
        // Round normal to closest axis
        scalar magX = mag(planeI.normal().x());
        scalar magY = mag(planeI.normal().y());
        scalar magZ = mag(planeI.normal().z());

        if (magX > magY && magX > magZ)
        {
            if (sign(planeI.normal().x()) < 0)
            {
                limBBox.max().x() = planeI.refPoint().x();
                return true;
            }
            else
            {
                limBBox.min().x() = planeI.refPoint().x();
                return true;
            }
        }
        else if (magY > magX && magY > magZ)
        {
            if (sign(planeI.normal().y()) < 0)
            {
                limBBox.max().y() = planeI.refPoint().y();
                return true;
            }
            else
            {
                limBBox.min().y() = planeI.refPoint().y();
                return true;
            }
        }
        else if (magZ > magX && magZ > magY)
        {
            if (sign(planeI.normal().z()) < 0)
            {
                limBBox.max().z() = planeI.refPoint().z();
                return true;
            }
            else
            {
                limBBox.min().z() = planeI.refPoint().z();
                return true;
            }
        }
    }

    return false;
}
//---------------------------------------------------------------------------//
void wallModel::getClosestPointAndNormal
(
    const point& startPoint,
    const vector& span,
    point& closestPoint,
    vector& normal
)
{
    closestPoint = contactPlanes_[0].nearestPoint(startPoint);
    normal = contactPlanes_[0].normal();
}
//---------------------------------------------------------------------------//
