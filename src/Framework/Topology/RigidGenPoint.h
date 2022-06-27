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
#include "Module/ComputeModule.h"
#include "JointTree.h"

namespace dyno 
{
	// Generate Rigid's Particle 
	// initializeImpl:  init particle
	// compute: 		moving particle 
	
	template<typename TDataType>
	class RigidGenPoint : public ComputeModule
	{
		DECLARE_TCLASS(RigidGenPoint, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TopologyModule::Pair2 Pair2;
        typedef typename TopologyModule::Pair3 Pair3;
		
		typedef VectorND<Real, 3> Pair3f;


		RigidGenPoint();
		~RigidGenPoint() override;
		
		void compute() override;
        
		void genRigidPoint();
        
    protected:
        bool initializeImpl() override;

	private:
		// <CapsuleID, JointID, 顶点Id>
        DArray<Pair3> m_pointClusters;
		DArray<Quat<Real>> m_initQuatT;			
		DArray<Quat<Real>> m_initQuatR;
		
		Reduction<int> m_reduce;
		Scan m_scan;

		int m_numRigidPoints;
	public:
		/**
		* @brief Capsules radius
		*/
		DEF_VAR_IN(Real, Radius, "Capsules radius");

		/**
		* @brief The number of joints
		*/
		// DEF_VAR_IN(int, JointSize, "The number of joints");

		/**
		 * @brief A set of points to be required from.
		 */
		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "A set of points to be required from");

		/**
		 * @brief A set of Capsules from JointTree.
		 */
		DEF_ARRAY_IN(JCapsule, Capsule, DeviceType::GPU, "A set of Capsules from JointTree.");

		DEF_ARRAY_IN(Quat<Real>, Rotate, DeviceType::GPU, "Bone Rotate");
		DEF_ARRAY_IN(Quat<Real>, Translate, DeviceType::GPU, "Bone Translate");

		// DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "Bone Velocity");
        // DEF_ARRAY_IN(Coord, AngularVelocity, DeviceType::GPU, "Bone AngularVelocity");

		DEF_ARRAY_IN(Coord, RigidPosition, DeviceType::GPU, "Virtual Point Set.");


        /**
		 * @brief Capsule Id
		 */
		DEF_ARRAY_OUT(int, CapsuleId, DeviceType::GPU, "Rigid Particle's Bone Id");

		/**
		 * @brief Constraint Between rigid particle and elastic particle. <elastic particle ID, Joint ID>
		 */
		DEF_ARRAY_OUT(Pair2, ConstraintRE, DeviceType::GPU, "Constraint: rigid & elastic");

		/**
		 * @brief Distance Joint
		 */
		DEF_ARRAY_OUT(Pair2, ConstraintRR, DeviceType::GPU, "Constraint: rigid & rigid");
	};
}