#pragma once
#include "Node.h"

#include "Topology/JointTree.h"

namespace dyno
{
	/*!
	*	\class	SkeletonLoader
	*	\brief	Load a Skeleton 
	*/
	template<typename TDataType>
	class SkeletonLoader : public Node
	{
		DECLARE_TCLASS(SkeletonLoader, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

        typedef std::vector<std::shared_ptr<JointTree<typename TDataType>>> JointList;

		SkeletonLoader();
		virtual ~SkeletonLoader();

        void setJointMap(JointList* jointMap) { m_jointMap = jointMap; }
		// DEF_VAR(std::string, FileName, "", "");

	public:
        DEF_ARRAY_OUT(JCapsule, InitBone, DeviceType::GPU, "Init Bone Capsule");
        DEF_ARRAY_OUT(JCapsule, Bone, DeviceType::GPU, "current Bone Capsule");

        DEF_ARRAY_OUT(Coord, Velocity, DeviceType::GPU, "Bone Velocity");
        DEF_ARRAY_OUT(Coord, AngularVelocity, DeviceType::GPU, "Bone AngularVelocity");

        JointList* m_jointMap = nullptr;
        
	protected:
		void resetStates() override;
        void updateTopology() override;
	};
}