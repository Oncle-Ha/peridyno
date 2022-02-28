#include "NeighborPointQueryJoint.h"
#include "Topology/GridHash.h"

namespace dyno
{
	__constant__ int offset_nq2[27][3] = { 
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
		//this->inOther()->tagOptional(true);
	}

	template<typename TDataType>
	NeighborPointQueryJoint<TDataType>::~NeighborPointQueryJoint()
	{
	}

	// 寻找关节所延伸胶囊体控制的顶点
	template<typename Coord, typename JCapsule, typename TDataType>
	__global__ void K_ComputeNeighbor(
		DArray<Coord> position, 
		GridHash<TDataType> hash, 
		Real h,
		DArray<JCapsule> caps,
		DArrayList<int> pointIds, 
		DArrayList<int> capIds,
		DArrayList<Real> capdis)  
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= caps.size()) return;
		
		JCapsule cap = caps[pId];
		int3 vId0 = hash.getIndex3(cap.v0);
		int3 vId1 = hash.getIndex3(cap.v1);
		
		Coord m = cap.v0;
		Coord s = (cap.v1 - cap.v0);
		Real d = s.norm();

		// 遍历线段所覆盖的Grid
		int3 vId = vId0;
		while(true)
		{	
			int gId = hash.getIndex(vId.x, vId.y, vId.z);
			int totalNum = hash.getCounter(gId);
			for (int i = 0; i < totalNum; i++) {
				int nbId = hash.getParticleId(gId, i);
				Coord pos_i = position[nbId];
				Real d_v0 = (pos_i - cap.v0).norm();
				Real d_v1 = (pos_i - cap.v1).norm();
				Real d_line = (pos_i - m).dot(s) / d;
				Real min_d = min(d_v0, min(d_v1, d_line));
				if (min_d < h)
				{
					capdis[nbId].atomicInsert(min_d);
					capIds[nbId].atomicInsert(pId);
					pointIds[pId].atomicInsert(nbId);
				}
			}

			// 选取线段上最近Grid
			float next_t = -1;
			int3 next_c;
			for (int c = 1; c < 27; c++)
			{
				int3 cId;
				cId.x = vId.x + offset_nq2[c][0];
				cId.y = vId.y + offset_nq2[c][1];
				cId.z = vId.z + offset_nq2[c][2];
				if (cId.x >= 0 && cId.y >= 0 && cId.z >= 0) 
				{ 	
					// 线段与立方体求交
					Coord min_v = hash.getMin3(cId);
					Coord max_v = hash.getMax3(cId);
					Coord min_t = (min_v - m) / s;
					Coord max_t = (max_v - m) / s;
					float mint = max(min_t[0], max(min_t[1], min_t[2]));
					float maxt = min(max_t[0], min(max_t[1], max_t[2]));
					if (mint > 0 && mint < maxt && (next_t == -1 || mint < next_t))
					{
						next_t = mint;
						next_c = cId;
					}
				}
			}

			if (next_t < 0 || next_t > 1) break;
			
			vId = next_c;
		}
	}


	// 确定每个顶点所属关节
	template<typename Coord, typename JCapsule>
	__global__ void K_ComputeJoint(	
		DArray<Coord> position, 
		Real h,
		DArray<JCapsule> caps,
		DArrayList<int> pointIds,
		DArrayList<int> capIds,
		DArrayList<Real> capdis,
		DArrayList<int> cluster)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId > caps.size()) return;

		auto& ptIds = pointIds[pId];
		int size_i = ptIds.size();
		for (int i = 0; i < size_i; ++i)
		{
			int point = ptIds[i];
			auto& cIds = capIds[point];
			auto& cdis = capdis[point];
			int size_j = cIds.size();
			Real min_d = -1;
			int id_joint = -1;
			for (int j = 0; j < size_j; ++j)
			{
				JCapsule cap = caps[cIds[j]];
				if (min_d == -1 || cdis[j] < min_d) // 就近原则（暂时）
				{
					min_d = cdis[j];
					id_joint = cap.id_joint;
				}
			}
			if (id_joint != -1) cluster[id_joint].atomicInsert(point);
		}
	}

	template<typename TDataType>
	void NeighborPointQueryJoint<TDataType>::compute()
	{
		// Prepare inputs
		auto& points	= this->inPosition()->getData();
		auto& capsules  = this->inCapsule()->getData();
		auto h			= this->inRadius()->getData();

		// Prepare outputs
		if (this->outCluster()->isEmpty())
			this->outCluster()->allocate();

		auto& clusters = this->outCluster()->getData();

		uint numJt  = this->inJointSize()->getData();
		uint numPt  = this->inPosition()->getDataPtr()->size();
		uint numCp  = this->inCapsule()->getDataPtr()->size();

		clusters.resize(numJt);

		// Construct hash grid
		Reduction<Coord> reduce;
		Coord hiBound = reduce.maximum(points.begin(), points.size());
		Coord loBound = reduce.minimum(points.begin(), points.size());

		GridHash<TDataType> hashGrid;
		hashGrid.setSpace(h, loBound - Coord(h), hiBound + Coord(h));
		hashGrid.clear();
		hashGrid.construct(points);

		DArrayList<int> pointIds;
		DArrayList<int> capIds;
		DArrayList<Real> capdis;
		pointIds.resize(numCp);
		capIds.resize(numPt);
		capdis.resize(numPt); // TODO: 待优化空间

		cuExecute(numCp,
			K_ComputeNeighbor,
			points,
			hashGrid,
			h,
			capsules,
			pointIds,
			capIds,
			capdis);
		cuSynchronize();

		cuExecute(numCp,
			K_ComputeJoint,
			points,
			h,
			capsules,
			pointIds,
			capIds,
			capdis,
			clusters);
		cuSynchronize();

		pointIds.clear();
		capIds.clear();
		capdis.clear();
		hashGrid.clear();
	}

	DEFINE_CLASS(NeighborPointQueryJoint);
}