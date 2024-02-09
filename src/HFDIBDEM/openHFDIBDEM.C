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
    Federico Municchi (2016),
    Martin Isoz (2019-*), Martin Kotouč Šourek (2019-*)
\*---------------------------------------------------------------------------*/
#include "openHFDIBDEM.H"
#include "polyMesh.H"
#include "fvCFD.H"
#include "fvMatrices.H"
#include "geometricOneField.H"

#include "interpolationCellPoint.H"
#include "interpolationCell.H"

#include "scalarMatrices.H"
#include "OFstream.H"
#include <iostream>
#include "defineExternVars.H"
#include "parameters.H"

#define ORDER 2

using namespace Foam;
using namespace contactModel;

//---------------------------------------------------------------------------//
openHFDIBDEM::openHFDIBDEM(const Foam::fvMesh& mesh)
:
mesh_(mesh),
HFDIBDEMDict_
(
    IOobject
    (
        "HFDIBDEMDict",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
transportProperties_
(
    IOobject
    (
        "transportProperties",
        "constant",
        mesh_,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
),
bodyNames_(HFDIBDEMDict_.lookup("bodyNames")),
prtcInfoTable_(0),
stepDEM_(readScalar(HFDIBDEMDict_.lookup("stepDEM"))),
recordSimulation_(readBool(HFDIBDEMDict_.lookup("recordSimulation")))
{
    materialProperties::matProps_insert(
        "None",
        materialInfo("None", 1, 1, 1, 1, 1, 1)
    );

    if(HFDIBDEMDict_.found("recordFirstTimeStep"))
    {
        recordFirstTimeStep_ = readBool(HFDIBDEMDict_.lookup("recordFirstTimeStep"));
    }

    if(HFDIBDEMDict_.found("nSolidsInDomain"))
    {
        solverInfo::setNSolidsTreshnold(readLabel(HFDIBDEMDict_.lookup("nSolidsInDomain")));
    }
    
    dictionary demDic = HFDIBDEMDict_.subDict("DEM");
    dictionary materialsDic = demDic.subDict("materials");
    List<word> materialsNames = materialsDic.toc();
    forAll(materialsNames, matI)
    {
        dictionary matIDic = materialsDic.subDict(materialsNames[matI]);
        materialProperties::matProps_insert(
            materialsNames[matI],
            materialInfo(
                materialsNames[matI],
                readScalar(matIDic.lookup("Y")),
                readScalar(matIDic.lookup("nu")),
                readScalar(matIDic.lookup("gamma")),
                readScalar(matIDic.lookup("mu")),
                readScalar(matIDic.lookup("adhN")),
                readScalar(matIDic.lookup("eps"))
            )
        );
    }

    if(demDic.found("interfaceAdh"))
    {
        dictionary interfAdhDic = demDic.subDict("interfaceAdh");
        List<word> interNames = interfAdhDic.toc();
        forAll(interNames, interI)
        {
            dictionary interDicI = interfAdhDic.subDict(interNames[interI]);
            wordList interMat = interDicI.lookup("materials");
            string interKey;
            if(interMat[0] < interMat[1])
            {
                interKey += interMat[0];
                interKey += "-";
                interKey += interMat[1];
            }
            else
            {
                interKey += interMat[1];
                interKey += "-";
                interKey += interMat[0];
            }

            interAdhesion::interAdhesion_insert(
                interKey,
                readScalar(interDicI.lookup("value"))
            );
        }
    }

    if(demDic.found("LcCoeff"))
    {
        contactModelInfo::setLcCoeff(readScalar(demDic.lookup("LcCoeff")));
        // Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contactModelInfo::getLcCoeff() << endl;
    }
    else
    {
        contactModelInfo::setLcCoeff(4.0);
        // Info <<" -- Coefficient for characteristic Lenght Lc is set to : 4.0"<< endl;
    }

    if(demDic.found("betaCoeff"))
    {
        contactModelInfo::setBetaCoeff(readScalar(demDic.lookup("betaCoeff")));
        // Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contactModelInfo::getLcCoeff() << endl;
    }
    else
    {
        contactModelInfo::setBetaCoeff(5.0);
        // Info <<" -- Coefficient for characteristic Lenght Lc is set to : 4.0"<< endl;
    }
    
    if(demDic.found("useOldModel"))
    {
        contactModelInfo::setContactModel(readBool(demDic.lookup("useOldModel")));
        // Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contactModelInfo::getLcCoeff() << endl;
    }
    else
    {
        contactModelInfo::setContactModel(true);
        // Info <<" -- Coefficient for characteristic Lenght Lc is set to : 4.0"<< endl;
    }
    

    Info <<" -- Coefficient for characteristic Lenght Lc is set to : "<< contactModelInfo::getLcCoeff() << endl;
    Info <<" -- Coefficient for beta disipation term is set to  : "<< contactModelInfo::getBetaCoeff() << endl;
    if(contactModelInfo::getContactModel())
    {
        Info <<" -- Old dissipation term with gamma term is active : "<< endl;
    }
    else
    {
        Info <<" -- New dissipation term with Coefficient of resttution is active : "<< endl;
    }
    
    

    dictionary patchDic = demDic.subDict("collisionPatches");
    List<word> patchNames = patchDic.toc();
    forAll(patchNames, patchI)
    {
        word patchMaterial = patchDic.subDict(patchNames[patchI]).lookup("material");
        vector patchNVec = patchDic.subDict(patchNames[patchI]).lookup("nVec");
        vector planePoint = patchDic.subDict(patchNames[patchI]).lookup("planePoint");

        wallPlaneInfo::wallPlaneInfo_insert(
            patchNames[patchI],
            patchNVec,
            planePoint
        );

        wallMatInfo::wallMatInfo_insert(
            patchNames[patchI],
            materialProperties::getMatProps()[patchMaterial]
        );
    }

    if(demDic.found("cyclicPatches"))
    {
        Info << "CyclicPatches Found " << endl;
        dictionary cyclicPatchDic = demDic.subDict("cyclicPatches");
        List<word> cyclicPatchNames = cyclicPatchDic.toc();
        forAll(cyclicPatchNames, patchI)
        {
            vector patchNVec = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("nVec");
            vector planePoint = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("planePoint");
            word neighbourPatch = cyclicPatchDic.subDict(cyclicPatchNames[patchI]).lookup("neighbourPatch");

            cyclicPlaneInfo::insert(
                cyclicPatchNames[patchI],
                patchNVec,
                planePoint,
                neighbourPatch
            );
        }
        Info << "CyclicPatches  " <<  cyclicPatchNames <<endl;
    }
    
    if (HFDIBDEMDict_.found("geometricD"))
    {
        geometricD = HFDIBDEMDict_.lookup("geometricD");
    }
    else
    {
        geometricD = mesh_.geometricD();
    }

    forAll (geometricD, direction)
    {
        if (geometricD[direction] == -1)
        {
            case3D = false;
            emptyDir[direction] = 1;
            emptyDim = direction;
            break;
        }
    }

    if (HFDIBDEMDict_.isDict("virtualMesh"))
    {
        dictionary vMDic = HFDIBDEMDict_.subDict("virtualMesh");
        virtualMeshLevel::setVirtualMeshLevel(readScalar(vMDic.lookup("level")),readScalar(vMDic.lookup("charCellSize")));
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;

    }
    else
    {
        virtualMeshLevel::setVirtualMeshLevel(1,1);
        Info <<" -- VirtMesh Decomposition Level is set to        : "<< virtualMeshLevel::getVirtualMeshLevel() << endl;
        Info <<" -- VirtMesh charCellSize for boundary is set to  : "<< virtualMeshLevel::getCharCellSize() << endl;

    }

    recordOutDir_ = mesh_.time().rootPath() + "/" + mesh_.time().globalCaseName() + "/bodiesInfo";
}
//---------------------------------------------------------------------------//
openHFDIBDEM::~openHFDIBDEM()
{}
//---------------------------------------------------------------------------//
void openHFDIBDEM::initialize
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF,
    label recomputeM0,
    word runTime
)
{
    if(HFDIBDEMDict_.found("outputSetup"))
    {
        dictionary outputDic = HFDIBDEMDict_.subDict("outputSetup");
        bool basicOutput = readBool(outputDic.lookup("basic"));
        bool iBoutput = readBool(outputDic.lookup("iB"));
        bool DEMoutput = readBool(outputDic.lookup("DEM"));
        bool addModelOutput = readBool(outputDic.lookup("addModel"));
        bool parallelDEMOutput = readBool(outputDic.lookup("parallelDEM"));
        InfoH.setOutput(
            basicOutput,
            iBoutput,
            DEMoutput,
            addModelOutput,
            parallelDEMOutput
        );
    }

    // get data from HFDIBDEMDict
    //HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");
    preCalculateCellPoints();

    if(HFDIBDEMDict_.found("interpolationSchemes"))
    {
        HFDIBinterpDict_ = HFDIBDEMDict_.subDict("interpolationSchemes");

        if(HFDIBinterpDict_.found("method"))
        {
            word intMethod = HFDIBinterpDict_.lookup("method");

            if(intMethod == "leastSquares")
            {
                dictionary lsCoeffsDict
                    = HFDIBinterpDict_.subDict("leastSquaresCoeffs");
                ibInterp_.set(new leastSquaresInt(
                    mesh_,
                    readScalar(lsCoeffsDict.lookup("distFactor")),
                    readScalar(lsCoeffsDict.lookup("radiusFactor")),
                    readScalar(lsCoeffsDict.lookup("angleFactor")),
                    readScalar(lsCoeffsDict.lookup("maxCCRows"))
                ));
            }
            else if(intMethod == "line")
            {
                ibInterp_.set(new lineInt(HFDIBinterpDict_));
            }
        }
    }

    bool startTime0(runTime == "0");

    // initialize addModels
    addModels_.setSize(bodyNames_.size());
    immersedBodies_.setSize(0);                                         //on the fly creation
    refineF *= 0;
    recomputeM0_ = recomputeM0;

    if(!startTime0)
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            fileNameList entries(readDir(recordOutDir_,fileType::directory)); // OF version 8, For version 6 use fileName::DIRECTORY instead of fileType::directory
            scalar runTimeS(stod(runTime));
            forAll(entries,entry)
            {
                scalar dirTime(stod(entries[entry].name()));
                if(dirTime > runTimeS)
                {
                    word pathI(recordOutDir_ + "/" + entries[entry]);
                    rmDir(pathI);
                }
            }
        }

        restartSimulation(body, refineF, runTime);
    }
    else
    {
        if(!isDir(recordOutDir_))
            mkDir(recordOutDir_);
        else
        {
            rmDir(recordOutDir_);
            mkDir(recordOutDir_);
        }
    }

    #include "initializeAddModels.H"

    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);
        InfoH << basic_Info << "Creating immersed body based on: " << bodyName << endl;

        label maxAdditions(1000);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions and immersedBodies_.size() < solverInfo::getNSolidsTreshnold())
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));
            cAddition++;

            // initialize the immersed bodies
            if (addModels_[modelI].getBodyAdded())
            {
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                InfoH << addModel_Info << "Trying to set immersedBodies" << endl;
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );
                immersedBodies_[addIBPos].createImmersedBody(body,refineF);
                immersedBodies_[addIBPos].computeBodyCharPars();
                if (immersedBodies_[addIBPos].getStartSynced())
                {
                    immersedBodies_[addIBPos].initSyncWithFlow(U);
                }
                verletList_.addBodyToVList(immersedBodies_[addIBPos]);
                InfoH << addModel_Info << "Body based on: " << bodyName << " successfully added" << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info << "Body based on: "
                    << bodyName << " should have been added but was not "
                    << "(probably overlap with an already existing body)"
                    << endl;
            }
        }
    }

    verletList_.initialSorting();
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::createBodies(volScalarField& body,volScalarField& refineF)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].postContactUpdateBodyField(body,refineF);
        }
    }

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncCreateImmersedBody(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            immersedBodies_[bodyId].updateOldMovementVars();
        }
    }

    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].chceckBodyOp();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // create body or compute body-fluid coupling and estimate
            // potential contacts with walls
            immersedBodies_[bodyId].inContactWithStatic(false);

            immersedBodies_[bodyId].updateOldMovementVars();
            immersedBodies_[bodyId].printStats();
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::postUpdateBodies
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].clearIntpInfo();
            immersedBodies_[bodyId].postPimpleUpdateImmersedBody(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::recreateBodies
(
    volScalarField& body,
    volScalarField& refineF
)
{
    refineF *= 0;
    preCalculateCellPoints();
    forAll (addModels_,modelI)
    {
        addModels_[modelI].recreateBoundBox();
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].recreateBodyField(body,refineF);
        }
    }
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].syncCreateImmersedBody(body,refineF);
            immersedBodies_[bodyId].checkIfInDomain(body);
            if(immersedBodies_[bodyId].getrecomputeM0() > 0)
            {
                immersedBodies_[bodyId].computeBodyCharPars();
                immersedBodies_[bodyId].recomputedM0();
            }
            InfoH << iB_Info << "-- body "
                << immersedBodies_[bodyId].getBodyId() << " Re-created" << endl;
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::interpolateIB( volVectorField & V
                              ,volVectorField & Vs
                              ,volScalarField & body)
{
    if(ibInterp_.valid())
    {
        ibInterp_->resetInterpolator(V);
    }
    // reset imposed field
    Vs = V;

    // loop over all the immersed bodies
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            // update imposed field according to body
            immersedBodies_[bodyId].updateVectorField(Vs, V.name(),body);

            if(ibInterp_.valid())
            {
                ibInterp_->ibInterpolate
                (
                    immersedBodies_[bodyId].getIntpInfo(),
                    Vs,
                    immersedBodies_[bodyId].getUatIbPoints(),
                    mesh_
                );
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeBodiesInfo()
{
    if(!recordSimulation_)
        return;

    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());


    mkDir(curOutDir);
    mkDir(curOutDir +"/stlFiles");
    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }
    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;
    // Pout << "Processor "<< Pstream::myProcNo() << endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        // Pout <<"Processor "<< Pstream::myProcNo() << " writes Body " << bodyId << endl;
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateDEM(volScalarField& body,volScalarField& refineF)
{
    clockTime checkForCycleContactRun;
    if (cyclicPlaneInfo::getCyclicPlaneInfo().size() > 0)
    {
        forAll (immersedBodies_,bodyId)
        {
            if (!immersedBodies_[bodyId].getGeomModel().isCluster())
            {
                vector transVec = vector::zero;

                if (detectCyclicContact(
                    immersedBodies_[bodyId].getWallCntInfo(),
                    transVec
                ))
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));
                    std::shared_ptr<periodicBody> newPeriodicBody
                        = std::make_shared<periodicBody>(mesh_, thrSurf);

                    newPeriodicBody->setRhoS(immersedBodies_[bodyId].getGeomModel().getRhoS());
                    std::shared_ptr<geomModel> iBcopy(immersedBodies_[bodyId].getGeomModel().getCopy());
                    iBcopy->bodyMovePoints(transVec);
                    newPeriodicBody->addBodyToCluster(immersedBodies_[bodyId].getGeomModelPtr());
                    newPeriodicBody->addBodyToCluster(iBcopy);
                    immersedBodies_[bodyId].getGeomModelPtr() = newPeriodicBody;

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                    Info << "Periodic body created for body " << bodyId << endl;
                }
            }
            else
            {
                periodicBody& cBody = dynamic_cast<periodicBody&>(immersedBodies_[bodyId].getGeomModel());

                if(cBody.shouldBeUnclustered())
                {
                    verletList_.removeBodyFromVList(immersedBodies_[bodyId]);

                    immersedBodies_[bodyId].getGeomModelPtr() = cBody.getRemGeomModel();

                    verletList_.addBodyToVList(immersedBodies_[bodyId]);
                    Info << "Periodic body unclustered for body " << bodyId << endl;
                }
            }
        }
    }
    checkForCyclic_ += checkForCycleContactRun.timeIncrement();
    scalar deltaTime(mesh_.time().deltaT().value());
    scalar pos(0.0);
    scalar step(stepDEM_);
    // scalar timeStep(step*deltaTivim me);
    // list<Tuple2<label,contactType>>
    List<DynamicList<vector>> bodiesPositionList(Pstream::nProcs());
    List<DynamicList<pointField>> bodiesSTLPositionList(Pstream::nProcs());
    // Infos <<bodiesPositionList.size() << endl;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> syncOutForceKeyTable;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> contactResolvedKeyTable;
    HashTable <label,Tuple2<label, label>,Hash<Tuple2<label, label>>> contactSizesKeyTable;
    HashTable <label,label,Hash<label>> wallContactIBTable;
    bodiesPositionList[Pstream::myProcNo()].clear();
    while( pos < 1)
    {
        bodiesPositionList[Pstream::myProcNo()].clear();
        bodiesSTLPositionList[Pstream::myProcNo()].clear();
        // Info <<bodiesPositionList.size() << endl;
        // bodiesPositionList[Pstream::myProcNo()] = List<vector>(immersedBodies_.size());
        label possibleWallContacts(0);
        label resolvedWallContacts(0);
        label possiblePrtContacts(0);
        label resolvedPrtContacts(0);

        InfoH << DEM_Info << " Start DEM pos: " << pos
            << " DEM step: " << step << endl;

        InfoH << basic_Info << " DEM - CFD Time: "
            << mesh_.time().value() + deltaTime*pos << endl;
        clockTime updateMovementAndMovement;
        forAll (immersedBodies_,ib)
        {
            clockTime updateMovement;
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);
            updateMovementTime_ += updateMovement.timeIncrement();
        }
        label iter = 0;
        clockTime calcMoveIB;
        // Info << "cP1" << endl;
        forAll(immersedBodies_,iter)
        {
            if(Pstream::myProcNo() == 0 && immersedBodies_[iter].getGeomModel().getcType() == sphere)
            {
                immersedBodies_[iter].moveImmersedBody(deltaTime*step);
                bodiesPositionList[Pstream::myProcNo()].append(immersedBodies_[iter].getGeomModel().getBodyPosition());
            }

            if(Pstream::myProcNo() == 0 && immersedBodies_[iter].getGeomModel().getcType() != sphere && immersedBodies_[iter].getGeomModel().getcType() != cluster)
            {
                immersedBodies_[iter].moveImmersedBody(deltaTime*step);
                bodiesSTLPositionList[Pstream::myProcNo()].append(immersedBodies_[iter].getGeomModel().getSTLBodyPoints());
            }

        }
        // Info << "cP2" << endl;
        calcMoveBodyTime_ += calcMoveIB.timeIncrement();
        clockTime reduceTime;
        reduce(iter,maxOp<label>());
        // Info << "cP3" << endl;
        Pstream::gatherList(bodiesPositionList,0);
        Pstream::scatterList(bodiesPositionList,0);

        Pstream::gatherList(bodiesSTLPositionList,0);
        Pstream::scatterList(bodiesSTLPositionList,0);
        // Info << "cP4" << endl;
        moveIBTime_ += reduceTime.timeIncrement();
        clockTime moveIBSync;

        forAll (immersedBodies_,ib)
        {
            if(immersedBodies_[ib].getGeomModel().getcType() == sphere)
            {
                immersedBodies_[ib].getGeomModel().setBodyPosition(bodiesPositionList[0][ib]);
            }
            if(immersedBodies_[ib].getGeomModel().getcType() != sphere && immersedBodies_[iter].getGeomModel().getcType() != cluster)
            {
                immersedBodies_[ib].getGeomModel().setBodyPosition(bodiesSTLPositionList[0][ib]);
            }
        }
        moveBodySyncTime_ += moveIBSync.timeIncrement();

        bodiesPositionList[Pstream::myProcNo()].clear();
        bodiesSTLPositionList[Pstream::myProcNo()].clear();
        // Info << "cP5" << endl;
        updateMovementAndStuff_ += updateMovementAndMovement.timeIncrement();
        clockTime updateVerletListRun;
        verletList_.update(immersedBodies_);
        updateVerletList_ += updateVerletListRun.timeIncrement();
//OS Time effitiency Testing    
        clockTime wallContactStopWatch;
//OS Time effitiency Testing
        DynamicLabelList wallContactIB;
        DynamicList<wallSubContactInfo*> wallContactList;
        // Info << "CheckPoint #1" << endl;
        wallContactIBTable.clear();
        clockTime detectWallContactClock;
        forAll (immersedBodies_,bodyId)
        {
            immersedBody& cIb(immersedBodies_[bodyId]);
            if (cIb.getIsActive())
            {
                // set F_ and T_ to zero.
                cIb.resetContactForces();

                if(cIb.getbodyOperation() != 0)
                {
                    // detect wall contact
                    clockTime detectWallContactFunctionClock;
                    bool isContact(detectWallContact
                        (
                            mesh_,
                            cIb.getibContactClass(),
                            cIb.getWallCntInfo()
                    )); 
                    pureDetectWallContactTime_ += detectWallContactFunctionClock.timeIncrement();
                    if(isContact)
                    {
                        cIb.getibContactClass().setWallContact(true);
                        cIb.getibContactClass().inContactWithStatic(true);
                        wallContactIB.append(bodyId);
                        wallContactIBTable.insert(bodyId,wallContactIB.size()-1);
                        cIb.getWallCntInfo().registerSubContactList(wallContactList);
                    }
                }
            }
        }
        detectWallContactTime_ += detectWallContactClock.timeIncrement();
        // Info << "cP6" << endl;
        possibleWallContacts = wallContactIB.size();
        // Info << " -- wallContactIB.size() : "<< wallContactIB.size() << " wallContactList.size() : "<<wallContactList.size() << endl;
        List<bool> wallContactResolvedList(wallContactList.size(),false);
        clockTime wallContactEvalClock;
        if(wallContactIB.size() > 0)
        {
            label wallContactPerProc(ceil(double(wallContactList.size())/Pstream::nProcs()));
            // Info <<" wallContactPerProc : "<< wallContactPerProc << endl;
            if( wallContactList.size() <= Pstream::nProcs())
            {
                wallContactPerProc = 1;
            }
//OS Time effitiency Testing                
            clockTime wallContactParallelRun;
//OS Time effitiency Testing                
            for(int assignProc = Pstream::myProcNo()*wallContactPerProc; assignProc < min((Pstream::myProcNo()+1)*wallContactPerProc,wallContactList.size()); assignProc++)
            {
                resolvedWallContacts++;
                wallSubContactInfo* sCW = wallContactList[assignProc];
                immersedBody& cIb(immersedBodies_[sCW->getBodyId()]);
                bool resolved(solveWallContact(
                    mesh_,
                    cIb.getWallCntInfo(),
                    deltaTime*step,
                    *sCW
                    ));
                sCW->setResolvedContact(resolved);
                wallContactResolvedList[assignProc] += resolved;
            }
            // Info << "CheckPoint #3" << endl;
            //OS Time effitiency Testing                
            wallContactParallelTime_ += wallContactParallelRun.timeIncrement();
            clockTime wallContactSCRun;
            //OS Time effitiency Testing 
            // Info << "CheckPoint #4" << endl;
            reduce(wallContactResolvedList,sumOp<List<bool>>());

            List<vector> iBodyOutForceList(wallContactIB.size(),vector::zero);
            List<vector> iBodyOutTorqueList(wallContactIB.size(),vector::zero);

            forAll (wallContactList,iB)
            {
                wallSubContactInfo* sC = wallContactList[iB];
                label bodyId(sC->getBodyId());
                if(wallContactIBTable.found(bodyId))
                {
                    label cKey(wallContactIBTable[bodyId]);
                    if(wallContactResolvedList[iB])
                    {
                        iBodyOutForceList[cKey] += sC->getOutForce().F;
                        iBodyOutTorqueList[cKey] += sC->getOutForce().T;
                    }
                }
                
            }

            reduce(iBodyOutForceList,sumOp<List<vector>>());
            reduce(iBodyOutTorqueList,sumOp<List<vector>>());

            forAll (wallContactIB,iB)
            {
                immersedBody& cIb(immersedBodies_[wallContactIB[iB]]);
                forces cF;
                cF.F = iBodyOutForceList[iB];
                cF.T = iBodyOutTorqueList[iB];
                cIb.updateContactForces(cF);
                cIb.getWallCntInfo().clearOldContact();
            }
            
            //OS Time effitiency Testing                
            wallContactReduceTime_ += wallContactSCRun.timeIncrement();
            //OS Time effitiency Testing 
        }
        resolveWallContactTime_ += wallContactEvalClock.timeIncrement();
        // Info << "CheckPoint #5" << endl;
        wallContactList.clear();
        wallContactIB.clear();
        // Info << "cP7" << endl;
        wallContactTime_ += wallContactStopWatch.timeIncrement();
        clockTime reduceStatisticalDataWall;
        reduce(resolvedWallContacts,sumOp<label>());
        InfoH << basic_Info << " -- Possible Wall Contacts: " << possibleWallContacts
            << " Resolved Wall Contacts: " << resolvedWallContacts 
            << " contactPerProc : " << ceil(possibleWallContacts/Pstream::nProcs()) << endl;
        // Info << "CheckPoint #6" << endl;
        reduceStatisticalDataWallTime_ += reduceStatisticalDataWall.timeIncrement();
        
//OS Time effitiency Testing        
        clockTime particleContactStopWatch;
        clockTime particleContactSCRun;
//OS Time effitiency Testing
        // Pout <<" Survived 3 " << endl;
        // Info << "CPoint #2" << endl;
        DynamicList<prtSubContactInfo*> contactList;
        // check only pairs whose bounding boxes are intersected for the contact
        // Info << "CheckPoint ;#7" << endl;
        label vListSize(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);

            label cInd(cPair.first());
            bool cStatic(immersedBodies_[cInd].getbodyOperation() == 0);

            label tInd(cPair.second());
            bool tStatic(immersedBodies_[tInd].getbodyOperation() == 0);

            if((immersedBodies_[cInd].getIsActive() && immersedBodies_[tInd].getIsActive())
                &&
                !(cStatic && tStatic)

            )
            {
                if(cStatic)
                    immersedBodies_[tInd].inContactWithStatic(true);

                if(tStatic)
                    immersedBodies_[cInd].inContactWithStatic(true);

                prtContactInfo& prtcInfo(getPrtcInfo(
                    cPair)
                );
// 
                prtcInfo.clearData();
                getContacts(
                    mesh_,
                    prtcInfo
                );
                // Pout << " -- cpair :"<< cPair <<" - prtcInfo.getContactList().size() : " << prtcInfo.getContactListSize() << endl;
                clockTime prtcSyncInfo;
                // prtcInfo.syncContactList();
                syncParticleContactInfoTime_ += prtcSyncInfo.timeIncrement();
                prtcInfo.registerContactList(contactList);
            }
            vListSize++;
        }
        // Info << "CPoint #3" << endl;
//OS Time effitiency Testing
        // Info << "CheckPoint #8" << endl;
        prtContactReduceTime_ += particleContactSCRun.timeIncrement();
        clockTime particleContactParallelRun;
//OS Time effitiency Testing 
        // Info << "CPoint #4" << endl;
        possiblePrtContacts = contactList.size();      
        // Info << "CheckPoint #5.0  possiblePrtContacts" << possiblePrtContacts << endl;

        List<bool> contactResolved(contactList.size(),false);
        // Info << "CPoint #5.1 -- contactResolved.size() : " << contactResolved.size()<<endl;
        List<label> contactResolvedcKey(contactList.size(),0);
        List<label> contactResolvedtKey(contactList.size(),0);
        // Info << "CPoint #5.2 -- contactResolvedcKey.size() : " << contactResolvedcKey.size()<<endl;
        // Info << "CPoint #5.3 -- contactResolvedtKey.size() : " << contactResolvedtKey.size()<<endl;
        bool syncedData(true);
        reduce(syncedData, orOp<bool>());
        // Info << "cP8" << endl;        // Pout << " -- contactList.size() : " << contactList.size() << endl;
        // Info << "CPoint #6" << endl;
        if(contactList.size() > 0 )
        {
            label contactPerProc(ceil(double(contactList.size())/Pstream::nProcs()));
            if( contactList.size() <= Pstream::nProcs())
            {
                contactPerProc = 1;
            }

            for(int assignProc = Pstream::myProcNo()*contactPerProc; assignProc < min((Pstream::myProcNo()+1)*contactPerProc,contactList.size()); assignProc++)
            {
                prtSubContactInfo* sCI = contactList[assignProc];

                const Tuple2<label, label>& cPair = sCI->getCPair();

                contactResolvedcKey[assignProc] = cPair.first();
                contactResolvedtKey[assignProc] = cPair.second();

                ibContactClass& cClass(immersedBodies_[cPair.first()].getibContactClass());
                ibContactClass& tClass(immersedBodies_[cPair.second()].getibContactClass());

                if(detectPrtPrtContact(mesh_,cClass,tClass,*sCI))
                {
                    resolvedPrtContacts++;
                    prtContactInfo& prtcInfo(getPrtcInfo(cPair));
                    
                    bool resolved(solvePrtContact(mesh_, prtcInfo, *sCI, deltaTime*step));
                    sCI->setResolvedContact(resolved);
                    contactResolved[assignProc] += resolved;
                }
            }
        }
        clockTime syncSortTime3_1;
        // Info << "cP9" << endl;
        // Pout << contactResolved.size() << " contactResolved.size() " << endl;
        reduce(contactResolved,sumOp<List<bool>>());
        // Info << "cP9.01" << endl;
        reduce(contactResolvedcKey,sumOp<List<label>>());
        // Info << "cP9.02" << endl;
        reduce(contactResolvedtKey,sumOp<List<label>>());        
        // Info << "cP9.03" << endl;
        // Info << "CPoint #8" << endl;
        contactResolvedKeyTable.clear();
        // Info << "cP9.1" << endl;
        forAll(contactResolvedcKey,cKey)
        {   
            contactResolvedKeyTable.insert(Tuple2<label, label>(contactResolvedcKey[cKey],contactResolvedtKey[cKey]),cKey);
        }
        // Info << "cP9.2" << endl;
        syncResultTime_syncDataCost_ += syncSortTime3_1.timeIncrement();
//OS Time effitiency Testing

        prtContactParallelTime_ += particleContactParallelRun.timeIncrement();
        prtContactTime_ += particleContactStopWatch.timeIncrement();
        clockTime syncResultsRun;
//OS Time effitiency Testing
        // Info << "cP9.3" << endl;
        List<vector> cBodyOutForceList(vListSize,vector::zero);
        List<vector> cBodyOutTorqueList(vListSize,vector::zero);
        List<vector> tBodyOutForceList(vListSize,vector::zero);
        List<vector> tBodyOutTorqueList(vListSize,vector::zero);
        // Info << "cP9.4" << endl;
       syncOutForceKeyTable.clear();
        // Info << "CPoint #10" << endl;
       label nIter(0);
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            clockTime syncSortTime;
            // List<vector> outForce(4,vector::zero);
            // Info << "CheckPoint #9.1" << endl;
            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);
            // Info << "CheckPoint #9.2" << endl;
            // label nIters(0);
            prtContactInfo& prtcInfo(getPrtcInfo(cPair));
            
            syncResultTime_sortContactList_ += syncSortTime.timeIncrement();

            if(contactResolvedKeyTable.found(cPair))
            {
                label nSubContact(0);
                std::vector<std::shared_ptr<prtSubContactInfo>>& subCList
                    = prtcInfo.getPrtSCList();
                // Info << "cP9.5" << endl;
                for(auto sC : subCList)
                {
                    nSubContact++;
                    cBodyOutForceList[nIter] += sC->getOutForce().first().F;
                    cBodyOutTorqueList[nIter] += sC->getOutForce().first().T;
                    tBodyOutForceList[nIter] += sC->getOutForce().second().F;
                    tBodyOutTorqueList[nIter] += sC->getOutForce().second().T;
                }
                // Info << "cP9.6" << endl;
                
            }
            // bodiesOutForceList[Pstream::myProcNo()][nIter] = outForce;

            // outForce.clear();

            syncOutForceKeyTable.insert(cPair,nIter);
            nIter++;
        }
        clockTime syncSortTime3_8;
        // Info << "cP10" << endl;

        // Pstream::gatherList(bodiesOutForceList,0);
        // Pstream::scatterList(bodiesOutForceList,0);
        reduce(cBodyOutForceList,sumOp<List<vector>>());
        reduce(cBodyOutTorqueList,sumOp<List<vector>>());
        reduce(tBodyOutForceList,sumOp<List<vector>>());
        reduce(tBodyOutTorqueList,sumOp<List<vector>>());
        // Info << "CheckPoint #1" << endl;
        syncResultTime_syncDataCost_ += syncSortTime3_8.timeIncrement();
        label nvListIter(0);
        // Info << "CheckPoint #2" << endl;
        for (auto it = verletList_.begin(); it != verletList_.end(); ++it)
        {
            clockTime syncSortTime;

            const Tuple2<label, label> cPair = Tuple2<label, label>(it->first, it->second);
            label cInd(cPair.first());
            label tInd(cPair.second());

            syncResultTime_sortContactList_ += syncSortTime.timeIncrement();

            if(!contactResolvedKeyTable.found(cPair))
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
            else if(contactResolved[contactResolvedKeyTable[cPair]])
            {
                if(!syncOutForceKeyTable.found(cPair))
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " not found in syncOutForceKeyTable" << endl;
                    continue;
                }

                nvListIter = syncOutForceKeyTable[cPair];
                if(nvListIter > cBodyOutForceList.size())
                {
                    Pout <<" -- cPair  "<<cInd << " - "<<tInd << " nvListIter > bodiesOutForceList[Pstream::myProcNo()].size()" << endl;
                    continue;
                }
                                
                clockTime syncSortTime3;
                
                vector F1 = vector::zero;
                vector T1 = vector::zero;
                vector F2 = vector::zero;
                vector T2 = vector::zero;

                F1 += cBodyOutForceList[nvListIter];                    
                T1 += cBodyOutTorqueList[nvListIter];
                F2 += tBodyOutForceList[nvListIter];
                T2 += tBodyOutTorqueList[nvListIter];

                forces cF;
                cF.F = F1;
                cF.T = T1;
                forces tF;
                tF.F = F2;
                tF.T = T2;

                immersedBodies_[cInd].updateContactForces(cF);
                immersedBodies_[tInd].updateContactForces(tF);

                syncResultTime_syncDataCost_ +=syncSortTime3.timeIncrement();
            }
            else
            {
                if(prtcInfoTable_.found(cPair))
                {
                    prtcInfoTable_.erase(cPair);
                    continue;
                }
            }
            // Info << "cP11" << endl;

        
        }

        // Info << "CheckPoint #10" << endl;
        syncResultsTime_ += syncResultsRun.timeIncrement();
        clockTime reduceStatisticalDataPrt;
        reduce(resolvedPrtContacts,sumOp<label>());
        reduceStatisticalDataPrtTime_ += reduceStatisticalDataPrt.timeIncrement();
        
        clockTime DEMIntergrationRun;
        InfoH << basic_Info << " -- Possible Particle Contacts: " << possiblePrtContacts
            << " Resolved Particle Contacts: " << resolvedPrtContacts 
            << " contactPerProc : " << ceil(float(possiblePrtContacts)/Pstream::nProcs()) << endl;
        scalar maxCoNum = 0;
        label  bodyId = 0;
        // Info << "cP12" << endl;
        forAll (immersedBodies_,ib)
        {
            clockTime updateMovement;
            immersedBodies_[ib].updateMovement(deltaTime*step*0.5);
            updateMovementTime_ += updateMovement.timeIncrement();
            immersedBodies_[ib].printBodyInfo();
            clockTime computeCONum;
            if(!solverInfo::getOnlyDEM())
            {
                immersedBodies_[ib].computeBodyCoNumber();
                if (maxCoNum < immersedBodies_[ib].getCoNum())
                {
                    maxCoNum = immersedBodies_[ib].getCoNum();
                    bodyId = ib;
                }
            }
            computeCONumTime_ += computeCONum.timeIncrement();
        }
        // Info << "cP13" << endl;
        if(!solverInfo::getOnlyDEM())
        {
            InfoH << basic_Info << "Max CoNum = " << maxCoNum << " at body " << bodyId << endl;
        }

        //OS Time effitiency Testing            
        demItegrationTime_ += DEMIntergrationRun.timeIncrement(); 
        //OS Time effitiency Testing

        pos += step;
        
        if (pos + step + SMALL >= 1)
            step = 1 - pos;

    }
}
//---------------------------------------------------------------------------//
prtContactInfo& openHFDIBDEM::getPrtcInfo(Tuple2<label,label> cPair)
{
    if(!prtcInfoTable_.found(cPair))
    {
        prtcInfoTable_.insert(cPair, autoPtr<prtContactInfo>( new prtContactInfo(
            immersedBodies_[cPair.first()].getibContactClass(),
            immersedBodies_[cPair.first()].getContactVars(),
            immersedBodies_[cPair.second()].getibContactClass(),
            immersedBodies_[cPair.second()].getContactVars()
        )));
    }

    return prtcInfoTable_[cPair]();
}
//---------------------------------------------------------------------------//
// function to either add or remove bodies from the simulation
void openHFDIBDEM::addRemoveBodies
(
    volScalarField& body,
    volVectorField& U,
    volScalarField& refineF
)
{
    forAll (addModels_,modelI)
    {
        word bodyName(bodyNames_[modelI]);

        label maxAdditions(50);
        label cAddition(0);

        while (addModels_[modelI].shouldAddBody(body) and cAddition < maxAdditions)
        {
            InfoH << addModel_Info << "addModel invoked action, trying to add new body" << endl;
            std::shared_ptr<geomModel> bodyGeomModel(addModels_[modelI].addBody(body, immersedBodies_));

            cAddition++;

            if (addModels_[modelI].getBodyAdded())
            {
                InfoH << addModel_Info << "STL file correctly generated, registering the new body" << endl;

                // prepare pointer list for IBs (increase its size)
                label newIBSize(immersedBodies_.size()+1);
                label addIBPos(newIBSize - 1);
                immersedBodies_.setSize(newIBSize);

                // create the new body
                immersedBodies_.set
                (
                    addIBPos,
                    new immersedBody
                    (
                        bodyName,
                        mesh_,
                        HFDIBDEMDict_,
                        transportProperties_,
                        addIBPos,
                        recomputeM0_,
                        bodyGeomModel,
                        ibInterp_,
                        cellPoints_
                    )
                );

                // get reference for further processing
                immersedBody& nBody(immersedBodies_[addIBPos]);
                nBody.createImmersedBody(body,refineF);
                nBody.computeBodyCharPars();
                if (nBody.getStartSynced())
                {
                    nBody.initSyncWithFlow(U);
                }
                verletList_.addBodyToVList(nBody);

                InfoH << addModel_Info
                    << "new body included into the simulation" << endl;
                cAddition = 0;
            }
            else
            {
                InfoH << addModel_Info
                    << "new body should have been added but was not "
                    << "(probably overlap with an existing body)"
                    << endl;
            }
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::updateFSCoupling
(
    volScalarField& body,
    volVectorField& f
)
{
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            immersedBodies_[bodyId].pimpleUpdate(body,f);
        }
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::restartSimulation
(
    volScalarField& body,
    volScalarField& refineF,
    word runTime
)
{
    word timePath(recordOutDir_+"/"+runTime);
    fileNameList files(readDir(timePath));
    scalar thrSurf(readScalar(HFDIBDEMDict_.lookup("surfaceThreshold")));

    forAll(files,f)
    {
        IOdictionary bodyDict
        (
            IOobject
            (
                timePath + "/" + files[f],
                mesh_,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        );

        word bodyId(std::to_string(readLabel(bodyDict.lookup("bodyId"))));
        word bodyName(bodyDict.lookup("bodyName"));
        vector Vel(bodyDict.lookup("Vel"));
        scalar omega(readScalar(bodyDict.lookup("omega")));
        vector Axis(bodyDict.lookup("Axis"));
        bool isStatic(readBool(bodyDict.lookup("static")));
        label timeStepsInContWStatic(readLabel(bodyDict.lookup("timeStepsInContWStatic")));

        std::shared_ptr<geomModel> bodyGeomModel;
        word bodyGeom;
        // check if the immersedDict_ contains bodyGeom
        if (HFDIBDEMDict_.subDict(bodyName).found("bodyGeom"))
        {
            word input = HFDIBDEMDict_.subDict(bodyName).lookup("bodyGeom");
            bodyGeom = input;
            InfoH << iB_Info << "Found bodyGeom for "
                << bodyName << ", the body is: " << bodyGeom << endl;
        }
        else
        {
            bodyGeom = "convex";
            InfoH << iB_Info << "Did not find bodyGeom for "
                << bodyName << ", using bodyGeom: " << bodyGeom << endl;
        }

        if(bodyGeom == "convex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<convexBody>(mesh_,stlPath,thrSurf);
        }
        else if(bodyGeom == "nonConvex")
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath,thrSurf);
        }
        else if(bodyGeom == "sphere")
        {
            vector startPosition = bodyDict.subDict("sphere").lookup("position");
            scalar radius = readScalar(bodyDict.subDict("sphere").lookup("radius"));

            bodyGeomModel = std::make_shared<sphereBody>(mesh_,startPosition,radius,thrSurf);
        }
        else
        {
            word stlPath(timePath + "/stlFiles/"+bodyId+".stl");
            InfoH << iB_Info << "bodyGeom: " << bodyGeom
                << " not supported, using bodyGeom nonConvex" << endl;
            bodyGeom = "nonConvex";
            bodyGeomModel = std::make_shared<nonConvexBody>(mesh_,stlPath,thrSurf);
        }

        label newIBSize(immersedBodies_.size()+1);
        label addIBPos(newIBSize - 1);
        immersedBodies_.setSize(newIBSize);

        InfoH << iB_Info << "Restarting body: " << bodyId << " as "
            << addIBPos << " bodyName: " << bodyName << endl;
        immersedBodies_.set
        (
            addIBPos,
            new immersedBody
            (
                bodyName,
                mesh_,
                HFDIBDEMDict_,
                transportProperties_,
                addIBPos,
                recomputeM0_,
                bodyGeomModel,
                ibInterp_,
                cellPoints_
            )
        );

        immersedBodies_[addIBPos].createImmersedBody(body,refineF);
        immersedBodies_[addIBPos].computeBodyCharPars();
        immersedBodies_[addIBPos].setRestartSim(Vel,omega,Axis,isStatic,timeStepsInContWStatic);
        verletList_.addBodyToVList(immersedBodies_[addIBPos]);
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::preCalculateCellPoints()
{
    cellPoints_.clear();
    cellPoints_.setSize(mesh_.nCells());
    forAll(mesh_.C(), cellI)
    {
        cellPoints_[cellI] = mesh_.cellPoints()[cellI];
    }

    forAll (immersedBodies_,bodyId)
    {
        immersedBodies_[bodyId].getGeomModel().resetHashTable();
    }
}
//---------------------------------------------------------------------------//
void openHFDIBDEM::writeFirtsTimeBodiesInfo()
{
    word curOutDir(recordOutDir_ + "/" + mesh_.time().timeName());
    bool checkExistance(false);
    reduce(checkExistance,orOp<bool>());
    if(!recordSimulation_ || isDir(curOutDir))
        return;

    if(Pstream::myProcNo() == 0)
    {
        mkDir(curOutDir);
    }

    reduce(checkExistance,orOp<bool>());

    
    mkDir(curOutDir +"/stlFiles");
    DynamicLabelList activeIB;
    forAll (immersedBodies_,bodyId)
    {
        if (immersedBodies_[bodyId].getIsActive())
        {
            activeIB.append(bodyId);
        }
    }
    wordList bodyNames;
    scalar listZize(activeIB.size());
    label bodiesPerProc = ceil(listZize/Pstream::nProcs());
    InfoH << basic_Info << "Active IB listZize      : " << listZize<< endl;
    InfoH << basic_Info << "bodiesPerProc : " << bodiesPerProc<< endl;
    // Pout << "Processor "<< Pstream::myProcNo() << endl;

    for(int assignProc = Pstream::myProcNo()*bodiesPerProc; assignProc < min((Pstream::myProcNo()+1)*bodiesPerProc,activeIB.size()); assignProc++)
    {
        const label bodyId(activeIB[assignProc]);
        // Pout <<"Processor "<< Pstream::myProcNo() << " writes Body " << bodyId << endl;
        word path(curOutDir + "/body" + std::to_string(immersedBodies_[bodyId].getBodyId()) +".info");
        OFstream ofStream(path);
        IOobject outClass
            (
                path,
                mesh_,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            );
        IOdictionary outDict(outClass);

        outDict.writeHeader(ofStream);
        immersedBodies_[bodyId].recordBodyInfo(outDict,curOutDir);
        outDict.writeData(ofStream);
    }

}
//---------------------------------------------------------------------------//
void openHFDIBDEM::setSolverInfo()
{
    solverInfo::setOnlyDEM(true);
}
//---------------------------------------------------------------------------//
