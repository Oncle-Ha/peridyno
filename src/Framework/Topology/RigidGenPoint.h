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
        
        
    protected:
        bool initializeImpl() override;

	private:
        DArray<Pair3> m_pointClusters;
		Reduction<int> m_reduce;
		Scan m_scan;

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

		DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "Bone Velocity");
        DEF_ARRAY_IN(Coord, AngularVelocity, DeviceType::GPU, "Bone AngularVelocity");
		/**
		 * @brief Virtual Point Set
		 */
		DEF_ARRAY_OUT(Coord, Position, DeviceType::GPU, "Return Virtual Point Set.");

        /**
		 * @brief Capsule Id
		 */
		DEF_ARRAY_OUT(int, CapsuleId, DeviceType::GPU, "Rigid Particle's Bone Id");

		/**
		 * @brief Constraint Between rigid particle and elastic particle.
		 */
		DEF_ARRAY_OUT(Pair2, ConstraintRE, DeviceType::GPU, "Constraint: rigid & elastic");

		/**
		 * @brief Distance Joint
		 */
		DEF_ARRAY_OUT(Pair2, ConstraintRR, DeviceType::GPU, "Constraint: rigid & rigid");
	};
}