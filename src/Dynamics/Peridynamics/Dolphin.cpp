#include "Dolphin.h"
#include "Topology/TetrahedronSet.h"
#include "Topology/TriangleSet.h"
#include "Topology/PointSet.h"
#include "Topology/MixSet.h"
#include "Topology/NeighborPointQuery.h"

#include "Mapping/PointSetToPointSet.h"
// #include "Mapping/CapsuleToMixSet.h"

#include "Peridynamics/Module/Peridynamics.h"
#include "Peridynamics/Module/ConstraintSolverModule.h"

#include "SharedFunc.h"



namespace dyno
{
    IMPLEMENT_TCLASS(Dolphin, TDataType)

    template<typename TDataType>
    Dolphin<TDataType>::Dolphin(std::string name)
        :ParticleSystem<TDataType>(name)
    {
        // Init
        auto mixSet = std::make_shared<MixSet<TDataType>>();
		this->stateTopology()->setDataPtr(mixSet);
        this->varHorizon()->setValue(0.0085);

        // Peridynamics
        auto peri = std::make_shared<Peridynamics<TDataType>>();
        {
            peri->varOutForce()->setValue(true);
            this->varTimeStep()->connect(peri->inTimeStep());
            this->statePosition()->connect(peri->inPosition());
            this->stateVelocity()->connect(peri->inVelocity());
            this->stateForce()->connect(peri->inForce());
            this->stateRestShape()->connect(peri->inRestShape());
            this->animationPipeline()->pushModule(peri);
        }

        // Generate rigid particle & constrain 
        auto gen = std::make_shared<RigidGenPoint<TDataType>>();
        {
            this->varRadius()->connect(gen->inRadius());
            this->statePosition()->connect(gen->inPosition());
            this->stateRigidPosition()->connect(gen->inRigidPosition());
            this->inBone()->connect(gen->inCapsule());
            this->inRotate()->connect(gen->inRotate());
            this->inTranslate()->connect(gen->inTranslate());
            this->animationPipeline()->pushModule(gen);
        }

        // Solver Constraint
        auto con = std::make_shared<ConstraintSolverModule<TDataType>>();
        {
            this->varTimeStep()->connect(con->inTimeStep());
            this->statePosition()->connect(con->inPosition());
            this->stateVelocity()->connect(con->inVelocity());
            this->stateRigidPosition()->connect(con->inRigidPosition());
            gen->outCapsuleId()->connect(con->inCapsuleId());
            gen->outConstraintRR()->connect(con->inConstraintRR());
            gen->outConstraintRE()->connect(con->inConstraintRE());
            peri->outElasticForce()->connect(con->inElasticForce());
            this->animationPipeline()->pushModule(con);
        }


        //TODO: Rigid System

		//Create a node for surface mesh rendering
        {
            m_surfaceNode = std::make_shared<Node>("Mesh");
            
            auto triSet = std::make_shared<TriangleSet<TDataType>>();
            m_surfaceNode->stateTopology()->setDataPtr(triSet);
            
            //Set the topology mapping from MixSet to TriangleSet
            auto surfaceMapping = this->template addTopologyMapping<PointSetToPointSet<TDataType>>("surface_mapping");
            auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->stateTopology()->getDataPtr());

            surfaceMapping->setFrom(ptSet);
            surfaceMapping->setTo(triSet);        
        }

        // Set the Topology mapping from Capsule(JointTree) to MixSet 
        // 1.Module.update()
        // 2.Mapping.apply()
        /*
        {
            auto capsuleMapping = this->template addTopologyMapping<CapsuleToMixSet<TDataType>>("joint_mapping");
            auto& jointMap = this->stateJointMap()->getDataPtr();

            capsuleMapping->setFrom(jointMap);
            capsuleMapping->setTo(mixSet);
            capsuleMapping->setCapsuleRadius(this->varRadius()->getData());
            
            this->stateVelocity()->connect(capsuleMapping->inVelocity());
            this->stateForce()->connect(capsuleMapping->inForce());
            this->varTimeStep()->connect(capsuleMapping->inTimeStep());
            // TODO: 修改用速度更新。
            capsuleMapping->outV0()->connect(this->outV0());
            capsuleMapping->outV1()->connect(this->outV1());
        }*/
    }

    template<typename TDataType>
    Dolphin<TDataType>::~Dolphin()
    {

    }
    
    template<typename TDataType>
    bool Dolphin<TDataType>::scale(Real s)
    {
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->stateTopology()->getDataPtr())->scale(s);
        TypeInfo::cast<MixSet<TDataType>>(this->stateTopology()->getDataPtr())->scale(s);
        // TODO
        // this->stateJointMap()->getData()[0]->scale(s); // Root

        return true;
    }

    template<typename TDataType>
    bool Dolphin<TDataType>::translate(Coord t)
    {
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->stateTopology()->getDataPtr())->translate(t);
        TypeInfo::cast<MixSet<TDataType>>(this->stateTopology()->getDataPtr())->translate(t);
        //TODO
        // this->stateJointMap()->getData()[0]->translate(t); // Root

        return true;
    }

    template<typename TDataType>
    void Dolphin<TDataType>::resetStates()
    {
        auto mixSet = TypeInfo::cast<MixSet<TDataType>>(this->stateTopology()->getDataPtr());
        auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->stateTopology()->getDataPtr());
		if (ptSet == nullptr) return;

		auto pts = ptSet->getPoints();
        
        // reset Position
		if (pts.size() > 0)
		{
			this->statePosition()->setElementCount(pts.size());
			this->stateVelocity()->setElementCount(pts.size());
			this->stateForce()->setElementCount(pts.size());

			this->statePosition()->getData().assign(pts);
			this->stateVelocity()->getDataPtr()->reset();
		}
		Node::resetStates();

        //Update Neighbor & RefPosition
		auto nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
 		this->varHorizon()->connect(nbrQuery->inRadius());
 		this->statePosition()->connect(nbrQuery->inPosition());
		nbrQuery->update();

		if (!this->statePosition()->isEmpty())
		{
			this->stateRestShape()->allocate();
			auto nbrPtr = this->stateRestShape()->getDataPtr();
			nbrPtr->resize(nbrQuery->outNeighborIds()->getData());
            
			constructRestShape(*nbrPtr, nbrQuery->outNeighborIds()->getData(), this->statePosition()->getData());

			this->stateReferencePosition()->allocate();
			this->stateReferencePosition()->getDataPtr()->assign(this->statePosition()->getData());

			this->stateNeighborIds()->allocate();
			this->stateNeighborIds()->getDataPtr()->assign(nbrQuery->outNeighborIds()->getData());
		}
    }

    template<typename TDataType>
    void Dolphin<TDataType>::updateTopology()
    {
        auto mixSet = TypeInfo::cast<MixSet<TDataType>>(this->stateTopology()->getDataPtr());
		auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->stateTopology()->getDataPtr());

        auto& pts = ptSet->getPoints();
        auto& curPos = this->statePosition()->getData();
        pts.assign(curPos);

        //if (int (this->varElapsedTime()->getData() * 240 ) % 2 == 0)
        {
            // TODO
            //Animation
            // for (auto joint : jointMap)
            // {
            //     joint->applyAnimationAll(this->varElapsedTime()->getData());
            //     // joint->applyAnimationAll(0.05);
            // }
        

            auto tMappings = this->getTopologyMappingList();
            for(auto iter = tMappings.begin(); iter != tMappings.end(); iter++){
                (*iter)->apply();
            }
        }

        //DEBUG
        curPos.assign(pts);
    }
    
	template<typename TDataType>
	std::shared_ptr<PointSetToPointSet<TDataType>> Dolphin<TDataType>::getTopologyMapping()
	{
		// auto mapping = this->template getModule<CapsuleToMixSet<TDataType>>("joint_mapping");
        auto mapping = this->template getModule<PointSetToPointSet<TDataType>>("surface_mapping");

		return mapping;
	}    

    template<typename TDataType>
    void Dolphin<TDataType>::loadMixFile(std::string filename)
    {
        TypeInfo::cast<MixSet<TDataType>>(this->stateTopology()->getDataPtr())->loadMixFile(filename);
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->stateTopology()->getDataPtr())->loadObjFile(filename.append("_surface.obj"));
    }

    DEFINE_CLASS(Dolphin);
}