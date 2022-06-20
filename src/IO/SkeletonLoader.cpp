#include "SkeletonLoader.h"

namespace dyno
{
	IMPLEMENT_TCLASS(SkeletonLoader, TDataType)

	template<typename TDataType>
	SkeletonLoader<TDataType>::SkeletonLoader()
		: Node()
	{

	}

	template<typename TDataType>
	SkeletonLoader<TDataType>::~SkeletonLoader()
	{
		
	}

	template<typename TDataType>
	void SkeletonLoader<TDataType>::resetStates()
	{
        //TODO
        // m_jointMap[0]->scale(this->varScale()->getData())
        // m_jointMap[0]->translate(this->varLocation()->getData())

        // init Bone
	}

	template<typename TDataType>
	void SkeletonLoader<TDataType>::updateTopology()
    {
        //TODO
        //Animation
        // for (auto joint : m_jointMap)
        // {
        //     joint->applyAnimationAll(this->varElapsedTime()->getData());
        //     // joint->applyAnimationAll(0.05);
        // }
        //

        //TODO
        // std::vector<JCapsule>capsule_list;
		// int id_joint = 0;
		// for (auto joint : *m_from)
		// {
		// 	int id_cap = 0;
		// 	for (auto joint_son : joint->children)
		// 	{
		// 		capsule_list.push_back(JCapsule{id_joint, id_cap, 
		// 										joint->GlCoord, joint_son->GlCoord});
		// 		++id_cap;
		// 		//DEBUG0
		// 		// float s = 1;
		// 		// printf("(%f,%f,%f) -> (%f,%f,%f)\n",
		// 		// 	joint->GlCoord[0] / s, joint->GlCoord[1] / s, joint->GlCoord[2] / s,
		// 		// 	joint_son->GlCoord[0] / s, joint_son->GlCoord[1] / s, joint_son->GlCoord[2] / s);
		// 	}
		// 	++id_joint;
		// }
    }

	DEFINE_CLASS(SkeletonLoader);
}