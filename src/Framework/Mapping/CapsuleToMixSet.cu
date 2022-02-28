#include "CapsuleToMixSet.h"

#include "Matrix/MatrixFunc.h"
#include "Topology/NeighborPointQuery.h"

#include <iostream>

namespace dyno
{
	template<typename TDataType>
	CapsuleToMixSet<TDataType>::CapsuleToMixSet()
		: TopologyMapping()
	{

	}

	template<typename TDataType>
	CapsuleToMixSet<TDataType>::CapsuleToMixSet(std::shared_ptr<JointTree<TDataType>> from, std::shared_ptr<MixSet<TDataType>> to)
	{
		m_from = from;s
		m_to = to;
	}

	template<typename TDataType>
	CapsuleToMixSet<TDataType>::~CapsuleToMixSet()
	{

	}


	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::initializeImpl()
	{
		match(m_from, m_to);
		return true;
	}

	//TODO: fix the problem
	
	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::apply()
	{
		uint pDim = cudaGridSize(m_to->getPoints().size(), BLOCK_SIZE);
        
        // Animation
	
		return true;
	}

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::match(std::shared_ptr<JointTree<TDataType>> from, std::shared_ptr<MixSet<TDataType>> to)
	{
		m_initFrom = std::make_shared<JointTree<TDataType>>();
		m_initTo = std::make_shared<MixSet<TDataType>>();

		m_initFrom->copyFrom(*from);
		m_initTo->copyFrom(*to);

        //Collision among Capsule Tet Tri

		// auto nbQuery = std::make_shared<NeighborPointQuery<TDataType>>();

		// nbQuery->inRadius()->setValue(m_radius);
		// nbQuery->inPosition()->allocate()->assign(m_initFrom->getPoints());
		// nbQuery->inOther()->allocate()->assign(m_initTo->getPoints());

		// nbQuery->update();

		// mNeighborIds.assign(nbQuery->outNeighborIds()->getData());
	}

	DEFINE_CLASS(CapsuleToMixSet);
}