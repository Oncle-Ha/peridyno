#include "RigidGenPoint.h"
#include <thrust/sort.h>
#include "Topology/NeighborPointQueryJoint.h"

namespace dyno
{

	IMPLEMENT_TCLASS(RigidGenPoint, TDataType)

	template<typename TDataType>
	RigidGenPoint<TDataType>::RigidGenPoint()
		: ComputeModule()
	{
	}

	template<typename TDataType>
	RigidGenPoint<TDataType>::~RigidGenPoint()
	{
	}

    template<typename Coord, typename Pair3>
	__global__ void RGP_GetVirtualCoord(
		DArray<Pair3> clusters,
		DArray<Coord> to,
		DArray<Coord> virTo,
		DArray<int> capsuleId)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		int p = clusters[pId][1];

        virTo[pId] = to[p];
		capsuleId[pId] = clusters[pId][2];
	}

    template<typename Pair3, typename Pair2>
	__global__ void RGP_GetConstraintRE(
		DArray<Pair3> clusters,
		DArray<Pair2> ConstraintRE)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		int p = clusters[pId][1];

		ConstraintRE[pId] = Pair2(pId, p);
	}

	template<typename TDataType>
	bool RigidGenPoint<TDataType>::initializeImpl()
	{
		this->outPosition()->allocate();
		this->outCapsuleId()->allocate();
		this->outConstraintRE()->allocate();
		this->outConstraintRR()->allocate();
        auto& virtualPosition = this->outPosition()->getData();
		auto& capsuleId = this->outCapsuleId()->getData();
		auto& constraintRE = this->outConstraintRE()->getData();
		auto& constraintRR = this->outConstraintRR()->getData();
		
        auto& pos = this->inPosition()->getData();
		auto nbQuery = std::make_shared<NeighborPointQueryJoint<TDataType>>();

		nbQuery->inRadius()->setValue(this->inRadius()->getValue());
		nbQuery->inPosition()->allocate()->assign(pos);
		nbQuery->inCapsule()->allocate()->assign(this->inCapsule()->getData());

		nbQuery->update();

		// 获取顶点所对应的几何体
		// 返回DArray<Pair<几何体ID, JointID>> 

		m_pointClusters.assign(nbQuery->outPJPair()->getData());
		int numPoint = pos.size();
		int numPair = m_pointClusters.size();

        DArray<int> count;

        virtualPosition.resize(numPair);
		capsuleId.resize(numPair);

		// TODO:sort by capsule ID

        cuExecute(numPair,
            RGP_GetVirtualCoord,
            m_pointClusters,
            pos,
            virtualPosition,
			capsuleId);
        cuSynchronize();

		constraintRE.resize(numPair);
        cuExecute(numPair,
            RGP_GetConstraintRE,
            m_pointClusters,
			capsuleId);
        cuSynchronize();		
		
		// constraintRR
		

		count.clear();

		return true;
	}

	template<typename TDataType>
	void RigidGenPoint<TDataType>::compute()
	{
		// update Position
	}

	DEFINE_CLASS(RigidGenPoint);
}