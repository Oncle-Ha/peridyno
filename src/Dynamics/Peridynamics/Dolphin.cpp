#include "Dolphin.h"
#include "Topology/TetrahedronSet.h"
#include "Topology/TriangleSet.h"
#include "Topology/PointSet.h"
#include "Topology/MixSet.h"
#include "Mapping/PointSetToPointSet.h"
#include "Topology/NeighborPointQuery.h"
#include "Peridynamics/Peridynamics.h"
#include "SharedFunc.h"

namespace dyno
{
    IMPLEMENT_CLASS_1(Dolphin, TDataType)

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
		this->currentPosition()->connect(peri->inPosition());
		this->currentVelocity()->connect(peri->inVelocity());
		this->currentForce()->connect(peri->inForce());
		this->currentRestShape()->connect(peri->inRestShape());
		this->animationPipeline()->pushModule(peri);

		//Create a node for surface mesh rendering
		m_surfaceNode = this->template createAncestor<Node>("Mesh");

		auto triSet = m_surfaceNode->template setTopologyModule<TriangleSet<TDataType>>("surface_mesh");
		m_surfaceNode->currentTopology()->setDataPtr(triSet);

		//Set the topology mapping from MixSet to TriangleSet
		auto surfaceMapping = this->template addTopologyMapping<PointSetToPointSet<TDataType>>("surface_mapping");
		// auto mixSet = TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr());
        auto ptSet = TypeInfo::cast<PointSet<TDataType>>(this->currentTopology()->getDataPtr());

		surfaceMapping->setFrom(ptSet);
		surfaceMapping->setTo(triSet);        

        //Create a node for joint
        // m_jointNode = this->template createAncestor<Node>("Joint");
        // auto jointTree = m_jointNode->template setTopologyModule<JointTree<TDataType>>("joint_tree");
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

        return true;
    }

    template<typename TDataType>
    bool Dolphin<TDataType>::translate(Coord t)
    {
        TypeInfo::cast<TriangleSet<TDataType>>(m_surfaceNode->currentTopology()->getDataPtr())->translate(t);
        TypeInfo::cast<MixSet<TDataType>>(this->currentTopology()->getDataPtr())->translate(t);

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
			this->currentPosition()->setElementCount(pts.size());
			this->currentVelocity()->setElementCount(pts.size());
			this->currentForce()->setElementCount(pts.size());

			this->currentPosition()->getData().assign(pts);
			this->currentVelocity()->getDataPtr()->reset();
		}
		Node::resetStates();        

        //Update Neighbor & RefPosition
		auto nbrQuery = std::make_shared<NeighborPointQuery<TDataType>>();
 		this->varHorizon()->connect(nbrQuery->inRadius());
 		this->currentPosition()->connect(nbrQuery->inPosition());
		nbrQuery->update();

		if (!this->currentPosition()->isEmpty())
		{
			this->currentRestShape()->allocate();
			auto nbrPtr = this->currentRestShape()->getDataPtr();
			nbrPtr->resize(nbrQuery->outNeighborIds()->getData());
            
			constructRestShape(*nbrPtr, nbrQuery->outNeighborIds()->getData(), this->currentPosition()->getData());

			this->currentReferencePosition()->allocate();
			this->currentReferencePosition()->getDataPtr()->assign(this->currentPosition()->getData());

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
        pts.assign(this->currentPosition()->getData());

        auto tMappings = this->getTopologyMappingList();
        for(auto iter = tMappings.begin(); iter != tMappings.end(); iter++){
            (*iter)->apply();
        }

    }
    
	template<typename TDataType>
	std::shared_ptr<PointSetToPointSet<TDataType>> Dolphin<TDataType>::getTopologyMapping()
	{
		auto mapping = this->template getModule<PointSetToPointSet<TDataType>>("surface_mapping");

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