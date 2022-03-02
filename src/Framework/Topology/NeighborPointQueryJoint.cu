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

	// max_cap {O(Gird*Point)}
	// 寻找关节所延伸胶囊体控制的顶点
	template<typename Coord, typename JCapsule, typename TDataType>
	__global__ void K_ComputeNeighbor(
		DArray<Coord> position, 
		GridHash<TDataType> hash, 
		Real h,
		DArray<JCapsule> caps,
		DArrayList<int> pointIds,  // cap:[points..]
		DArrayList<int> capIds,	   // point:[caps..]
		DArrayList<Real> capdis)   // point:[dis..]
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
			//DEBUG

			int gId = hash.getIndex(vId.x, vId.y, vId.z);
			printf("Gird: %d, Pid: %d\n", gId, pId);
			if (gId == -1) break;

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
			//FIXME: 选取不正确
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
						//DEBUG
						//if (cId == vId) printf("Error cId == vId\n");
					}
				}
			}

			if (next_t < 0 || next_t > 1) break;
			
			vId = next_c;
		}
	}

	// max_cap {O(Gird*Point)}
	// 统计<顶点，关节>点对
	__global__ void K_CountPair(
		DArrayList<int> pointIds,
		DArray<int> count)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= count.size()) return;

		count[pId] = pointIds[pId].size();
		//if (count[pId] < 0) printf("Error <0\n");
		//if (count[pId] == 0) printf("Error =0 %d\n", pId);
	}

	// max_cap {O(Gird*Point)}
	// 存储<顶点，关节>点对
	template<typename JCapsule, typename Pair2>
	__global__ void K_SetPair(
		DArray<JCapsule> caps,
		DArrayList<int> pointIds,
		DArrayList<int> capIds,
		DArrayList<Real> capdis,
		DArray<int> count,
		DArray<Pair2> pairs)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= caps.size()) return;

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
			if (id_joint != -1)  // assert(id_joint != -1)
			{
				pairs[count[pId] + i] = Pair2(id_joint, point);
			}
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
		if (this->outPJPair()->isEmpty())
			this->outPJPair()->allocate();

		auto& pairs = this->outPJPair()->getData();

		// uint numJt  = this->inJointSize()->getData();
		uint numPt  = this->inPosition()->getDataPtr()->size();
		uint numCp  = this->inCapsule()->getDataPtr()->size();

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
		DArray<int> count;
		
		pointIds.resize(numCp);
		capIds.resize(numPt);
		capdis.resize(numPt); // TODO: 待优化空间
		count.resize(numCp);
		// FIXME 无输出
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

		// BUGXX
		int ppp = pointIds.elementSize();
		int ppp2 = capIds.elementSize();
		int ppp3 = capdis.elementSize();
		// std::cerr << pointIds.elementSize() << std::endl;

		cuExecute(numCp,
			K_CountPair,
			pointIds,
			count);
		cuSynchronize();
	
		int numPr = m_reduce.accumulate(count.begin(), count.size());
		m_scan.exclusive(count, true);
		pairs.resize(numPr);

		cuExecute(numCp,
			K_SetPair,
			capsules,
			pointIds,
			capIds,
			capdis,
			count,
			pairs);
		cuSynchronize();

		count.clear();
		pointIds.clear();
		capIds.clear();
		capdis.clear();
		hashGrid.clear();
	}

	DEFINE_CLASS(NeighborPointQueryJoint);
}