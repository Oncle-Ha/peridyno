#include "Dolphin.h"
#include "Topology/TetrahedronSet.h"
#include "Topology/TriangleSet.h"
#include "Topology/PointSet.h"
#include "Topology/MixSet.h"
#include "Mapping/PointSetToPointSet.h"
// #include "Mapping/CapsuleToMixSet.h"
#include "Topology/NeighborPointQuery.h"
#include "Peridynamics/Peridynamics.h"
#include "SharedFunc.h"

namespace dyno
{
    IMPLEMENT_TCLASS(Dolphin, TDataType)

    template<typename TDataType>
    Dolphin<TDataType>::Dolphin(std::string name)
        :ParticleSystem<TDataType>(name)
    {

        auto mixSet = std::make_shared<MixSet<TDataType>>();
		this->currentTopology()->setDataPtr(mixSet);
        // this->currentPoints()->setDataPtr(mixSet->getPointSet());

        this->varHorizon()->setValue(0.0085);

		auto peri = std::make_shared<Peridynamics<TDataType>>();
		this->varTimeStep()->connect(peri->inTimeStep());
		this->statePosition()->connect(peri->inPosition());
		this->stateVelocity()->connect(peri->inVelocity());
		this->stateForce()->connect(peri->inForce());
		this->currentRestShape()->connect(peri->inRestShape());
		// this->animationPipeline()->pushModule(peri);// 暂时只控制点集

		//Create a node for surface mesh rendering
		// m_surfaceNode = this->template createAncestor<Node>("Mesh");
        m_surfaceNode = std::make_shared<Node>("Mesh");
        
		// auto triSet = m_surfaceNode->template setTopologyModule<TriangleSet<TDataType>>("surface_mesh");
        auto triSet = std::make_shared<TriangleSet<TDataType>>();
		m_surfaceNode->currentTopology()->setDataPtr(triSet);

		//Set the topology mapping from MixSet to TriangleSet
		// auto surfaceMapping = this->template addTopologyMapping<PointSetToPointSet<TDataType>>("surface_mapping");
        auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->currentTopology()->getDataPtr());

		// surfaceMapping->setFrom(ptSet);
		// surfaceMapping->setTo(triSet);        

        // Set the Topology mapping from Capsule(JointTree) to MixSet 
        // 1.Module.update()
        // 2.Mapping.apply()
        auto jointMapping = this->template addTopologyMapping<CapsuleToMixSet<TDataType>>("joint_mapping");
        jointMapping->setFrom(&m_jointMap);
        jointMapping->setTo(mixSet);
        // jointMapping->setCapsuleRadius(0.0125);
		jointMapping->setCapsuleRadius(0.0325);
        // jointMapping->setCapsuleRadius(0.055);
        // jointMapping->setCapsuleRadius(0.085);
        this->currentColor()->connect(jointMapping->outColor());
        // jointMapping->outColor()->connect(this->currentColor());

    }

    template<typename TDataType>
    Dolphin<TDataType>::~Dolphin()
    {

    }
    
    template<typename TDataType>
    bool Dolphin<TDataType>::scale(Real s)
    {
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->currentTopology()->getDataPtr())->scale(s);
        TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr())->scale(s);
        m_jointMap[0]->scale(s); // Root

        return true;
    }

    template<typename TDataType>
    bool Dolphin<TDataType>::translate(Coord t)
    {
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->currentTopology()->getDataPtr())->translate(t);
        TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr())->translate(t);
        m_jointMap[0]->translate(t); // Root

        return true;
    }

    template<typename TDataType>
    void Dolphin<TDataType>::resetStates()
    {
        auto mixSet = TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr());
        auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->currentTopology()->getDataPtr());
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
			this->currentRestShape()->allocate();
			auto nbrPtr = this->currentRestShape()->getDataPtr();
			nbrPtr->resize(nbrQuery->outNeighborIds()->getData());
            
			constructRestShape(*nbrPtr, nbrQuery->outNeighborIds()->getData(), this->statePosition()->getData());

			this->currentReferencePosition()->allocate();
			this->currentReferencePosition()->getDataPtr()->assign(this->statePosition()->getData());

			this->currentNeighborIds()->allocate();
			this->currentNeighborIds()->getDataPtr()->assign(nbrQuery->outNeighborIds()->getData());
		}
    }

    template<typename TDataType>
    void Dolphin<TDataType>::updateTopology()
    {
        auto mixSet = TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr());
		auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->currentTopology()->getDataPtr());

        auto& pts = ptSet->getPoints();
        auto& curPos = this->statePosition()->getData();
        pts.assign(curPos);

        //DEBUG
        for (auto joint : m_jointMap)
        {
            joint->applyAnimationAll(this->varElapsedTime()->getData());
            // joint->applyAnimationAll(0.05);
        }

        auto tMappings = this->getTopologyMappingList();
        for(auto iter = tMappings.begin(); iter != tMappings.end(); iter++){
            (*iter)->apply();
        }
        

        //DEBUG
        curPos.assign(pts);
    }
    
	template<typename TDataType>
	std::shared_ptr<CapsuleToMixSet<TDataType>> Dolphin<TDataType>::getTopologyMapping()
	{
		auto mapping = this->template getModule<CapsuleToMixSet<TDataType>>("joint_mapping");

		return mapping;
	}    

    template<typename TDataType>
    void Dolphin<TDataType>::loadMixFile(std::string filename)
    {
        TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr())->loadMixFile(filename);
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->currentTopology()->getDataPtr())->loadObjFile(filename.append("_surface.obj"));
    }

    DEFINE_CLASS(Dolphin);
}