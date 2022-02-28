#include "CapsuleToMixSet.h"

#include "Matrix/MatrixFunc.h"
#include "Topology/NeighborPointQueryJoint.h"

#include <iostream>

namespace dyno
{
	template<typename TDataType>
	CapsuleToMixSet<TDataType>::CapsuleToMixSet()
		: TopologyMapping()
	{

	}

	template<typename TDataType>
	CapsuleToMixSet<TDataType>::CapsuleToMixSet(JointList* from, std::shared_ptr<MixSet<TDataType>> to)
	{
		m_from = from;
		m_to = to;
	}

	template<typename TDataType>
	CapsuleToMixSet<TDataType>::~CapsuleToMixSet()
	{

	}


	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::initializeImpl()
	{
		match();
		return true;
	}

	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::apply()
	{
		uint pDim = cudaGridSize(m_to->getPoints().size(), BLOCK_SIZE);
        
        // Animation
	
		return true;
	}

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::match()
	{

        //Collision among Capsule Tet Tri

		auto nbQuery = std::make_shared<NeighborPointQueryJoint<TDataType>>();

		nbQuery->inRadius()->setValue(m_radius);
		nbQuery->inPosition()->allocate()->assign(m_to->getPoints());
		nbQuery->inJointSize()->setValue(m_from->size());
		
		for (auto joint : *m_from)
		{
			joint->getGlobalTransform();
			Vec4f tmp = joint->GlobalTransform * Vec4f(0, 0, 0, 1) ;
			joint->GlCoord = Coord(tmp[0] / tmp[3], tmp[1] / tmp[3], tmp[2] / tmp[3]);
		}

		std::vector<JCapsule>capsule_list;
		int id_joint = 0;
		for (auto joint : *m_from)
		{
			int id_cap = 0;
			for (auto joint_son : joint->children)
			{
				capsule_list.push_back(JCapsule{id_joint, id_cap, 
												joint->GlCoord, joint_son->GlCoord});
				++id_cap;
			}
			++id_joint;
		}

		nbQuery->inCapsule()->allocate()->assign(capsule_list);

		nbQuery->update();

		mClusters.assign(nbQuery->outCluster()->getData());
	}

	DEFINE_CLASS(CapsuleToMixSet);
}