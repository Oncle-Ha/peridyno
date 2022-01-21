#pragma once
#include "ParticleSystem/ParticleSystem.h"
#include "NeighborData.h"
#include "Topology/PointSet.h"
#include "Topology/JointTree.h"
#include "Topology/Cluster.h"

namespace dyno
{
    template<typename> class PointSetToPointSet;

    /**
     * \class Dolphin
     * \brief 
     */
    template<typename TDataType>
	class Dolphin : public ParticleSystem<TDataType>
    {
        DECLARE_CLASS_1(Dolphin, TDataType)

    public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef TPair<TDataType> NPair;

        Dolphin(std::string name = "default");
		virtual ~Dolphin();

        void loadMixFile(std::string filename);
        
		bool translate(Coord t) override;
		bool scale(Real s) override;

        std::shared_ptr<Node> getSurfaceNode() {return m_surfaceNode;}
	    std::shared_ptr<PointSetToPointSet<TDataType>> getTopologyMapping();
        

	public:
		DEF_VAR(Real, Horizon, 0.01, "Horizon");

		DEF_EMPTY_CURRENT_ARRAY(ReferencePosition, Coord, DeviceType::GPU, "Reference position");

		DEF_EMPTY_CURRENT_ARRAYLIST(int, NeighborIds, DeviceType::GPU, "Storing the ids for neighboring particles");

		DEF_EMPTY_CURRENT_ARRAYLIST(NPair, RestShape, DeviceType::GPU, "Storing neighbors");

        // DEF_INSTANCE_STATE(PointSet<TDataType>, Points, ""); // TODO: remove this

        std::vector<std::shared_ptr<JointTree<TDataType>>> m_jointMap;
        std::vector<std::shared_ptr<Cluster<TDataType>>> m_clusters;

    protected:

        void resetStates() override; // set->current

        void updateTopology() override;// current->set

    private:
        std::shared_ptr<Node> m_surfaceNode;
	};
}