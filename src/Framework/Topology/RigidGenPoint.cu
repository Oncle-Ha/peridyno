#include "RigidGenPoint.h"
#include <thrust/sort.h>
#include "Topology/NeighborPointQueryJoint.h"
#include <thrust/sort.h>

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
		int p = clusters[pId][2];

        virTo[pId] = to[p];
		capsuleId[pId] = clusters[pId][0];
	}

    template<typename Pair3, typename Pair2>
	__global__ void RGP_GetConstraintRE(
		DArray<Pair3> clusters,
		DArray<Pair2> ConstraintRE)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;

		ConstraintRE[pId] = Pair2(clusters[pId][2], clusters[pId][1]);
	}

    template<typename Pair3>
	__global__ void RGP_AdjustOrder(
		DArray<Pair3> clusters)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		int p = clusters[pId][2];
		clusters[pId][2] = clusters[pId][0];
		clusters[pId][0] = p;
	}

	template<typename Coord, typename Real, template<typename Real> typename Quat>
	__global__ void RGP_ApplyTransformPointByQuat(
		DArray<int> capId,
		DArray<Quat<Real>> initT,
		DArray<Quat<Real>> initR,
		DArray<Quat<Real>> T,
		DArray<Quat<Real>> R,
		DArray<Coord> pos)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pos.size()) return;
		int p = pId; // point
		int j = capId[pId]; // Capsule
		// PT_d("point", point_i, pId);
		// PT_d("joint", couple[0], pId);
		Coord old_p = pos[p];
		Quat<Real> tmp_p(old_p[0], old_p[1], old_p[2], 0);
		
		tmp_p = initR[j].conjugate() * (tmp_p - initT[j]) * initR[j];
		tmp_p = R[j] * tmp_p * R[j].conjugate() + T[j];

		// PT_f("tmp3", tmp_p[3], pId);
		Coord new_p = Coord(tmp_p.x, tmp_p.y, tmp_p.z);
		// if (pId == 0)
			// printf("Vel Anim: [%f, %f, %f]", vel[p][0], vel[p][1], vel[p][2]);
		// vel[p] = (new_p - old_p) ;
		pos[p] = new_p;
	}	

	template<typename TDataType>
	void RigidGenPoint<TDataType>::genRigidPoint()
	{
		if(this->inRigidPosition()->isEmpty())
			this->inRigidPosition()->allocate();
		
		this->outCapsuleId()->allocate();
		this->outConstraintRE()->allocate();
		this->outConstraintRR()->allocate();
        auto& virtualPosition = this->inRigidPosition()->getData();
		auto& capsuleId = this->outCapsuleId()->getData();
		auto& constraintRE = this->outConstraintRE()->getData();
		auto& constraintRR = this->outConstraintRR()->getData();
		
		// Quat
		m_initQuatR.assign(this->inRotate()->getData());
		m_initQuatT.assign(this->inTranslate()->getData());
		
		// NeighborPointQueryJoint
        auto& pos = this->inPosition()->getData();
		auto nbQuery = std::make_shared<NeighborPointQueryJoint<TDataType>>();

		nbQuery->inRadius()->setValue(this->inRadius()->getValue());
		nbQuery->inPosition()->allocate()->assign(pos);
		nbQuery->inCapsule()->allocate()->assign(this->inCapsule()->getData());

		nbQuery->update();

		// 获取顶点所对应的几何体
		// 返回DArray<Pair<顶点Id, JointID, CapsuleID>> 
		DArray<Pair3> p_cluster;

		p_cluster.assign(nbQuery->outPJPair()->getData());
		int numPoint = pos.size();
		int numPair = p_cluster.size();
		m_numRigidPoints = numPair;

        DArray<int> count;

        virtualPosition.resize(numPair);
		capsuleId.resize(numPair);


		//sort by capsule ID
		cuExecute(numPair,
			RGP_AdjustOrder,
			p_cluster);
		cuSynchronize();

		thrust::sort(thrust::device, p_cluster.begin(), p_cluster.begin() + p_cluster.size());
		m_pointClusters.assign(p_cluster);

		// Get Rigid Coord
        cuExecute(numPair,
            RGP_GetVirtualCoord,
            m_pointClusters,
            pos,
            virtualPosition,
			capsuleId);
        cuSynchronize();


		// constraintRE
		constraintRE.resize(numPair);
        cuExecute(numPair,
            RGP_GetConstraintRE,
            m_pointClusters,
			constraintRE);
        cuSynchronize();		
		
		// constraintRR

		count.clear();
		p_cluster.clear();
	}

	template<typename TDataType>
	bool RigidGenPoint<TDataType>::initializeImpl()
	{
		genRigidPoint();
		return true;
	}

	template<typename TDataType>
	void RigidGenPoint<TDataType>::compute()
	{
		// update Position
		auto& quatR = this->inRotate()->getData();
		auto& quatT = this->inTranslate()->getData();

		auto& capsuleId = this->outCapsuleId()->getData();
		auto& rigidPos = this->inRigidPosition()->getData();

		cuExecute(m_numRigidPoints,
			RGP_ApplyTransformPointByQuat,
			capsuleId,
			m_initQuatT,
			m_initQuatR,
			quatT,
			quatR,
			rigidPos);
		cuSynchronize();

	}

	DEFINE_CLASS(RigidGenPoint);
}
