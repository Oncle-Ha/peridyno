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

	// 以重心为坐标原点的旋转、平移四元数转换
	template<typename TDataType>
	void SkeletonLoader<TDataType>::getCenterQuat(Coord v0, Coord v1, Quat<Real> &T, Quat<Real> &R)
	{
		Coord center = (v0 + v1) / 2.f;
		Coord dir = v1 - v0;
		float cos2 = dir[2];
		float cos1 = sqrtf((1 + cos2) / 2.0); 
		float sin1 = sqrtf((1 - cos2) / 2.0);
		Vec3f axis = Vec3f(-dir[1], dir[0], 0).normalize();
		Quat<Real> q(axis.x * sin1, axis.y * sin1, axis.z * sin1, cos1);
		Quat<Real> t(center[0], center[1], center[2], 0.f);

		R = q;
		T = t;
	}

	template<typename TDataType>
	bool SkeletonLoader<TDataType>::scale(Real s)
	{
		m_jointMap[0]->scale(s);
		
		return true;
	}
	
	template<typename TDataType>
	bool SkeletonLoader<TDataType>::translate(Coord t)
	{
		m_jointMap[0]->translate(t);
		
		return true;
	}

	template<typename TDataType>
	void SkeletonLoader<TDataType>::resetStates()
	{
		if (m_jointMap.empty())
		{
			printf("Load Skeleton failed.");
			return;
		}

        // init Bone
		{
			for (auto joint : m_jointMap)
			{
				joint->getGlobalQuat();
				joint->getGlobalCoord();
			}
			
			int id_joint = 0;
			int id_cap = 0;
			for (auto joint : m_jointMap)
			{
				for (auto joint_son : joint->children)
				{
					m_capLists.push_back(JCapsule{id_joint, id_cap, 
													joint->GlCoord, joint_son->GlCoord});
					Quat<Real> t, r;
					getCenterQuat(joint->GlCoord, joint_son->GlCoord, t, r);
					m_T.push_back(t);
					m_R.push_back(r);
					++id_cap;
				}
				++id_joint;
			}
			m_numCaps = id_cap;
			m_numjoints = id_joint;

			this->outBone()->allocate();
			this->outBone()->getData().resize(m_numCaps);
			this->outBone()->getData().assign(m_capLists);

			this->outTranslate()->allocate();
			this->outTranslate()->getData().resize(m_numCaps);
			this->outTranslate()->getData().assign(m_T);

			this->outRotate()->allocate();
			this->outRotate()->getData().resize(m_numCaps);
			this->outRotate()->getData().assign(m_R);
		}
	}

	template<typename TDataType>
	void SkeletonLoader<TDataType>::updateTopology()
	{
		if (m_jointMap.empty())
		{
			printf("Load Skeleton failed.");
			return;
		}		
        //Animation
        for (auto joint : m_jointMap)
        {
            joint->applyAnimationAll(this->varElapsedTime()->getData() / 100);
            // joint->applyAnimationAll(0.05);
        }
        
		for (auto joint : m_jointMap)
		{
			joint->getGlobalQuat();
			joint->getGlobalCoord();
		}

		{
			int index = 0;
			for (auto joint : m_jointMap)
			{
				for (auto joint_son : joint->children)
				{
					m_capLists[index].v0 = joint->GlCoord;
					m_capLists[index].v1 = joint_son->GlCoord;
					Quat<Real> t, r;
					getCenterQuat(joint->GlCoord, joint_son->GlCoord, t, r);
					m_T[index] = t;
					m_R[index] = r;
					index++;
				}
			}

			this->outBone()->getData().assign(m_capLists);
			this->outTranslate()->getData().assign(m_T);
			this->outRotate()->getData().assign(m_R);
		}
    }

	DEFINE_CLASS(SkeletonLoader);
}