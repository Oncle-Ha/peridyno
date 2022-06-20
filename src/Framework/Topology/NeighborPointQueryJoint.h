#pragma once
#include "Module/ComputeModule.h"
#include "JointTree.h"

namespace dyno 
{
	// return Point near Joint 
	template<typename TDataType>
	class NeighborPointQueryJoint : public ComputeModule
	{
		DECLARE_TCLASS(NeighborPointQueryJoint, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TopologyModule::Pair2 Pair2;
		typedef typename TopologyModule::Pair3 Pair3;
		
		typedef VectorND<Real, 4> Pair4f;


		NeighborPointQueryJoint();
		~NeighborPointQueryJoint() override;
		
		void compute() override;

	private:
		Reduction<int> m_reduce;
		Scan m_scan;

	public:
		DEF_VAR(uint, SizeLimit, 40, "Maximum number of neighbors");

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
		
		/**
		 * @brief <Joint, Point near Joint, Capsule>
		 */
		DEF_ARRAY_OUT(Pair3, PJPair, DeviceType::GPU, "Return <Joint, Point near Joint, Capsule>.");
	};
}