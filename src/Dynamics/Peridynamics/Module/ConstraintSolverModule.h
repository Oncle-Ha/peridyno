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
#include "Module/ConstraintModule.h"
#include "Module/TopologyModule.h"
#include "ViewGPUData.h"

namespace dyno {

	/**
	  * @brief This is an implementation of elasticity based on projective peridynamics.
	  *		   For more details, please refer to[He et al. 2017] "Projective Peridynamics for Modeling Versatile Elastoplastic Materials"
	  */
	template<typename TDataType>
	class ConstraintSolverModule : public ConstraintModule
	{
		DECLARE_TCLASS(ConstraintSolverModule, TDataType)

	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TDataType::Matrix Matrix;
		typedef typename TopologyModule::Pair2 Pair2;
        typedef typename TopologyModule::Pair3 Pair3;

		ConstraintSolverModule();
		~ConstraintSolverModule() override;
		
		void constrain() override;
		void computeInitValue();
		void updatePosition();
		void updateVelocity();
		void solveConstraint();

	protected:
		bool initializeImpl() override;

	public:
		DEF_VAR_IN(Real, TimeStep, "");

		/**
		 * @brief Particle position
		 */
		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "Particle position");

		/**
		 * @brief Particle velocity
		 */
		DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "Particle Velocity");


		/**
		 * @brief Rigid Particle position
		 */
		DEF_ARRAY_IN(Coord, RigidPosition, DeviceType::GPU, "Rigid Particle position");

        /**
		 * @brief Capsule Id
		 */
		DEF_ARRAY_IN(int, CapsuleId, DeviceType::GPU, "Rigid Particle's Bone Id");

		/**
		 * @brief Constraint Between rigid particle and elastic particle.
		 */
		DEF_ARRAY_IN(Pair2, ConstraintRE, DeviceType::GPU, "Constraint: rigid & elastic");

		/**
		 * @brief Distance Joint
		 */
		DEF_ARRAY_IN(Pair2, ConstraintRR, DeviceType::GPU, "Constraint: rigid & rigid");

		/**
		 * @brief Elastic force
		 */
		DEF_ARRAY_IN(Coord, ElasticForce, DeviceType::GPU, "Particle's Elastic Force");

	public:
		/**
		 * @brief  parameters
		 */
		// DEF_VAR(Real, Mu, 0.001, "Lame parameters: mu");

		// DEF_VAR(Real, Lambda, 0.01, "Lame parameters: lambda");
		// DEF_VAR(Real, Mu, 0.0001, "Lame parameters: mu");

		// DEF_VAR(Real, Lambda, 0.001, "Lame parameters: lambda");		

		DEF_VAR(uint, IterationNumber, 50, "Iteration number");

	protected:
		int mN;
		int mC;
		int numEPos;
		int numRPos;

		// DArray<>


		DArray<Coord> mY;
		DArray<Coord> mX;
		DArray<Coord> mF;
		DArray<Coord> mV;

		// DArray<Matrix> mJ;

		
		DArray<Real> mInvK;
		DArray<Real> mA;
		DArray<Real> mAdy;

		//Jacobi Method
		DArray<Real> mWa;
		DArray<Real> mWb;
		DArray<Real> mWc;
		DArray<Real> mLambda;

		DArray<Pair2> mConstraint;
		DArray<Real> mPhi;
		DArray<Real> mDis;
		DArray<Real> mConstJ;
		DArray<Real> mLen; // = 0
	};
}
