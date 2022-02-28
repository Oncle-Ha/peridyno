#pragma once
#include "Module/ComputeModule.h"
#include "JointTree.h"

namespace dyno 
{
	// return Number of Joint near Point
	template<typename TDataType>
	class NeighborPointQueryJoint : public ComputeModule
	{
		DECLARE_CLASS_1(NeighborPointQueryJoint, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		NeighborPointQueryJoint();
		~NeighborPointQueryJoint() override;
		
		void compute() override;

	private:
		void requestDynamicNeighborIds();

		void requestFixedSizeNeighborIds();

	public:
		// DEF_VAR(uint, SizeLimit, 0, "Maximum number of neighbors");

		/**
		* @brief Capsules radius
		*/
		DEF_VAR_IN(Real, Radius, "Capsules radius");

		/**
		 * @brief A set of points to be required from.
		 */
		DEF_ARRAY_IN(Coord, Position, DeviceType::GPU, "A set of points to be required from");

		/**
		 * @brief A set of Capsules from JointTree.
		 */
		DEF_ARRAY_IN(JointTree::JCapsule, Capsule, DeviceType::GPU, "A set of Capsules from JointTree.");
		
		/**
		 * @brief Joint ids
		 */
		DEF_ARRAY_OUT(int, JointIds, DeviceType::GPU, "Return Joint ids");
	};
}