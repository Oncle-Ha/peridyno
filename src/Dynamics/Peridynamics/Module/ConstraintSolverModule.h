/**
 * Copyright 2021 Xiaowei He
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
		typedef TPair<TDataType> NPair;

		ConstraintSolverModule();
		~ConstraintSolverModule() override;
		
		void constrain() override;

		void solveConstraint();

	protected:
		void preprocess() override;

	public:
		DEF_VAR_IN(Real, TimeStep, "");

		/**
		 * @brief Particle position
		 */
		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "Particle position");

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

		// DEF_VAR(uint, IterationNumber, 10, "Iteration number");

	protected:
		// DArray<>
	};
}