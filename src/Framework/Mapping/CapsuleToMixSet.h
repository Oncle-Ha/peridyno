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
		typedef std::vector<std::shared_ptr<JointTree<typename TDataType>>> JointList;

		CapsuleToMixSet();
		CapsuleToMixSet(JointList* from, std::shared_ptr<MixSet<TDataType>> to);
		~CapsuleToMixSet() override;

		void setCapsuleRadius(Real r) { m_radius = r; }

		void setFrom(JointList* from) { m_from = from; }
		void setTo(std::shared_ptr<MixSet<TDataType>> to) { m_to = to; }

		bool apply() override;

		// void match(JointList* from, std::shared_ptr<MixSet<TDataType>> to);
		void match();

	protected:
		bool initializeImpl() override;

	private:
		//Searching radius
		Real m_radius = 0.0125;

		DArrayList<int> m_pointClusters;
		DArrayList<int> m_tetClusters;
		DArrayList<int> m_triClusters;
		
		JointList* m_from = nullptr;
		// std::shared_ptr<JointTree<TDataType>> m_from = nullptr;
		std::shared_ptr<MixSet<TDataType>> m_to = nullptr;

		// JointList* m_initfrom = nullptr;
		// std::shared_ptr<JointTree<TDataType>> m_initFrom = nullptr;
		// std::shared_ptr<MixSet<TDataType>> m_initTo = nullptr;
	};
}