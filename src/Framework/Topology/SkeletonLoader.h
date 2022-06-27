/**
 * Copyright 2022 Xukun Luo
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
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

        void setJointMap(JointList &jointMap) { m_jointMap = jointMap; }
		void getCenterQuat(Coord v0, Coord v1, Quat<Real> &T, Quat<Real> &R);
		// DEF_VAR(std::string, FileName, "", "");

		bool translate(Coord t);
		bool scale(Real s);

	public:
        DEF_ARRAY_OUT(JCapsule, Bone, DeviceType::GPU, "Bone Capsule");

        // DEF_ARRAY_OUT(Coord, Velocity, DeviceType::GPU, "Bone Velocity");
        // DEF_ARRAY_OUT(Coord, AngularVelocity, DeviceType::GPU, "Bone AngularVelocity");
		DEF_ARRAY_OUT(Quat<Real>, Rotate, DeviceType::GPU, "Bone Rotate");
		DEF_ARRAY_OUT(Quat<Real>, Translate, DeviceType::GPU, "Bone Translate");

        JointList m_jointMap;

		std::vector<JCapsule> m_capLists;
		std::vector<Quat<Real>> m_T;
		std::vector<Quat<Real>> m_R;

        int m_numCaps = 0;
		int m_numjoints = 0;
	protected:
		void resetStates() override;
        void updateTopology() override;

	};
}