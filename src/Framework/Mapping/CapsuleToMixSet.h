#pragma once
#include "Module/TopologyMapping.h"
#include "Topology/PointSet.h"
#include "Topology/JointTree.h"
#include "Topology/MixSet.h"

namespace dyno
{
	template<typename TDataType> class PointSet;

	template<typename TDataType>
	class CapsuleToMixSet : public TopologyMapping
	{
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		CapsuleToMixSet();
		CapsuleToMixSet(std::shared_ptr<JointTree<TDataType>> from, std::shared_ptr<MixSet<TDataType>> to);
		~CapsuleToMixSet() override;

		void setCapsuleRadius(Real r) { m_radius = r; }

		void setFrom(std::shared_ptr<JointTree<TDataType>> from) { m_from = from; }
		void setTo(std::shared_ptr<MixSet<TDataType>> to) { m_to = to; }

		bool apply() override;

		void match(std::shared_ptr<JointTree<TDataType>> from, std::shared_ptr<MixSet<TDataType>> to);

	protected:
		bool initializeImpl() override;

	private:
		//Searching radius
		Real m_radius = 0.0125;

		DArrayList<int> mNeighborIds;

		std::shared_ptr<JointTree<TDataType>> m_from = nullptr;
		std::shared_ptr<MixSet<TDataType>> m_to = nullptr;

		std::shared_ptr<JointTree<TDataType>> m_initFrom = nullptr;
		std::shared_ptr<MixSet<TDataType>> m_initTo = nullptr;
	};
}