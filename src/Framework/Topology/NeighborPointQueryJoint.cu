#include "NeighborPointQueryJoint.h"
#include "Topology/GridHash.h"

namespace dyno
{
	__constant__ int offset_nq[27][3] = { 
		0, 0, 0,
		0, 0, 1,
		0, 1, 0,
		1, 0, 0,
		0, 0, -1,
		0, -1, 0,
		-1, 0, 0,
		0, 1, 1,
		0, 1, -1,
		0, -1, 1,
		0, -1, -1,
		1, 0, 1,
		1, 0, -1,
		-1, 0, 1,
		-1, 0, -1,
		1, 1, 0,
		1, -1, 0,
		-1, 1, 0,
		-1, -1, 0,
		1, 1, 1,
		1, 1, -1,
		1, -1, 1,
		-1, 1, 1,
		1, -1, -1,
		-1, 1, -1,
		-1, -1, 1,
		-1, -1, -1
	};

	IMPLEMENT_CLASS_1(NeighborPointQueryJoint, TDataType)

	template<typename TDataType>
	NeighborPointQueryJoint<TDataType>::NeighborPointQueryJoint()
		: ComputeModule()
	{
		this->inOther()->tagOptional(true);
	}

	template<typename TDataType>
	NeighborPointQueryJoint<TDataType>::~NeighborPointQueryJoint()
	{
	}

	// 寻找关节所延伸胶囊体控制的顶点
	template<typename Coord, typename JCapsule,  typename TDataType>
	__global__ void K_ComputeNeighbor(
		DArray<Coord> position, 
		GridHash<TDataType> hash, 
		Real h,
		DArray<JCapsule> caps,
		DArrayList<int> potintIds,
		DArrayList<int> capIds)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= caps.size()) return;
		
		JCapsule cap = caps[pId];
		int3 vId0 = hash.getIndex3(cap.v0);
		int3 vId1 = hash.getIndex3(cap.v1);

		for 
	}

	template<typename TDataType>
	void NeighborPointQueryJoint<TDataType>::compute()
	{
		// Prepare inputs
		auto& points	= this->inPosition()->getData();
		auto& capsules  = this->inCapsule()->getData();
		auto h			= this->inRadius()->getData();

		// Prepare outputs
		if (this->outJointIds()->isEmpty())
			this->outJointIds()->allocate();

		auto& jointIds = this->outJointIds()->getData();

		uint numPt  = this->inPosition()->getDataPtr()->size();
		
		jointIds.resize(numPt);

		// Construct hash grid
		Reduction<Coord> reduce;
		Coord hiBound = reduce.maximum(points.begin(), points.size());
		Coord loBound = reduce.minimum(points.begin(), points.size());

		GridHash<TDataType> hashGrid;
		hashGrid.setSpace(h, loBound - Coord(h), hiBound + Coord(h));
		hashGrid.clear();
		hashGrid.construct(points);


		hashGrid.clear();
	}

	DEFINE_CLASS(NeighborPointQueryJoint);
}