#include "NeighborPointQueryJoint.h"
#include "Topology/GridHash.h"
#include <thrust/sort.h>


#define PT_d(s, x, y) printf("%s: %d  pId: %d\n", s, x, y)
#define PT_f(s, x, y) printf("%s: %f  pId: %d\n", s, x, y)
#define PT_e(s, y) printf("[%s]  pId: %d\n", s, y)
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

	IMPLEMENT_TCLASS(NeighborPointQueryJoint, TDataType)

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
	// count<顶点, 距离, 关节>对
	template<typename Coord, typename JCapsule, typename TDataType>
	__global__ void K_CountNeighbor(
		DArray<Coord> position, 
		GridHash<TDataType> hash, 
		Real h,
		DArray<JCapsule> caps,
		DArray<int> count)  // cap:[<>..]
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= caps.size()) return;
		
		count[pId] = 0;
		JCapsule cap = caps[pId];

		// 削短胶囊体
		Coord tmp_s = (cap.v1 - cap.v0);
		Real tmp_d = tmp_s.norm();
		Real rate =  (h / tmp_d) * 0;
		if(rate > 0.5) rate = 0.5;
		cap.v0 += rate* tmp_s;
		cap.v1 -= rate * tmp_s;

		int3 vId0 = hash.getIndex3(cap.v0);
		int3 vId1 = hash.getIndex3(cap.v1);
		float eps = 1e-6;

		//DEBUG 
		// printf("[(%d,%d,%d)->(%d,%d,%d)] pId: %d\n", vId0.x, vId0.y, vId0.z, vId1.x, vId1.y, vId1.z, pId);

		Coord m = cap.v0;
		Coord s = (cap.v1 - cap.v0);
		Real d = s.norm();

		// 遍历线段所覆盖的Grid
		int3 vId = vId0;
		
		// 线段与立方体求交
		Coord min_v = hash.getMin3(vId);
		Coord max_v = hash.getMax3(vId);
		Coord time1 = (min_v - m) / s;
		Coord time2 = (max_v - m) / s;
		Coord max_t(max(time1[0],time2[0]), max(time1[1],time2[1]), max(time1[2],time2[2]));
		float next_t = min(max_t[0], min(max_t[1], max_t[2]));	//边界时间段	
		// PT_f("time", next_t, pId);
		while(true)
		{	
			int gId = hash.getIndex(vId.x, vId.y, vId.z);

			//DEBUG
			// PT_d("Gird", gId, pId);
			// printf("Gird:%d [%d,%d,%d] pId: %d\n", gId, vId.x, vId.y, vId.z, pId);

			if (gId == -1) break;

			//DEBUG 
			//PT_d("Num", totalNum, pId);
			
			for (int c = 0; c < 27; c++)
			{
				int3 cId;
				cId.x = vId.x + offset_nq2[c][0];
				cId.y = vId.y + offset_nq2[c][1];
				cId.z = vId.z + offset_nq2[c][2];
				if (cId.x >= 0 && cId.y >= 0 && cId.z >= 0) 
				{ 	
					int cNumId = hash.getIndex(cId.x, cId.y, cId.z);
					int totalNum = hash.getCounter(cNumId);
					for (int i = 0; i < totalNum; i++) {
						int nbId = hash.getParticleId(cNumId, i);
						Coord pos_i = position[nbId];
						
						Real d_v0 = (pos_i - cap.v0).norm();
						Real d_v1 = (pos_i - cap.v1).norm();
						Real d_t = (pos_i - m).dot(s);
						Real d_line = fabs(((pos_i - m).cross(s)).norm() / d);
						Real min_d = (d_t < d * d && d_t > 0)? d_line : min(d_v0, d_v1);

						if (min_d < h)
						{
							// PT_f("d_line", d_line, pId);
							count[pId] +=1;
						}
					}
				}
			}
			if (next_t > 1) break; //终点格子

			float tmp_t = -1;
			int3 next_c;	
			

			// 选取线段上最近Grid
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
					Coord time1 = (min_v - m) / s;
					Coord time2 = (max_v - m) / s;
					Coord min_t(min(time1[0],time2[0]), min(time1[1],time2[1]), min(time1[2],time2[2]));
					Coord max_t(max(time1[0],time2[0]), max(time1[1],time2[1]), max(time1[2],time2[2]));
					float mint = max(min_t[0], max(min_t[1], min_t[2]));
					float maxt = min(max_t[0], min(max_t[1], max_t[2]));

					// printf("GirdTest [%d,%d,%d] - (%f, %f) < (%f) pId: %d\n", cId.x, cId.y, cId.z, mint, maxt, tmp_t, pId);
					// printf("Checsk( %d )\n", (mint > next_t - eps));
					if (mint > 0 && mint < 1 && mint < maxt && (mint < tmp_t || tmp_t < 0) && mint > next_t - eps)
					{
						tmp_t = maxt;
						next_c = cId;
						//DEBUG
						//if (cId == vId) printf("Error cId == vId\n");
					}
				}
			}
			
			next_t = tmp_t;
			// PT_f("time", next_t, pId);
			if (next_t < 0) break;
			vId = next_c;
		}
		//DEBUG
		PT_d("Count", count[pId], pId);
	}

	// max_cap {O(Gird*Point)} TODO更新同Count一样
	// 寻找关节所延伸胶囊体控制的顶点
	template<typename Coord, typename JCapsule, typename TDataType, typename Pair4f>
	__global__ void K_ComputeNeighbor(
		DArray<Coord> position, 
		GridHash<TDataType> hash, 
		Real h,
		DArray<JCapsule> caps,
		DArray<Pair4f> capPairs,
		DArray<int> count)  // cap:[<>..]
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= caps.size()) return;
		
		JCapsule cap = caps[pId];
		// 削短胶囊体
		Coord tmp_s = (cap.v1 - cap.v0);
		Real tmp_d = tmp_s.norm();
		Real rate =  (h / tmp_d) * 0;
		if(rate > 0.5) rate = 0.5;
		cap.v0 += rate* tmp_s;
		cap.v1 -= rate * tmp_s;

		int3 vId0 = hash.getIndex3(cap.v0);
		int3 vId1 = hash.getIndex3(cap.v1);
		float eps = 1e-6;

		//DEBUG 
		// printf("[(%d,%d,%d)->(%d,%d,%d)] pId:%d\n", vId0.x, vId0.y, vId0.z, vId1.x, vId1.y, vId1.z, pId);

		Coord m = cap.v0;
		Coord s = (cap.v1 - cap.v0);
		Real d = s.norm();
		

		// 遍历线段所覆盖的Grid
		int3 vId = vId0;
		
		// 线段与立方体求交
		// 线段与立方体求交
		Coord min_v = hash.getMin3(vId);
		Coord max_v = hash.getMax3(vId);
		Coord time1 = (min_v - m) / s;
		Coord time2 = (max_v - m) / s;
		Coord max_t(max(time1[0],time2[0]), max(time1[1],time2[1]), max(time1[2],time2[2]));
		float next_t = min(max_t[0], min(max_t[1], max_t[2]));	//边界时间段
		int start = count[pId];
		int cnt = 0;
		while(true)
		{	
			int gId = hash.getIndex(vId.x, vId.y, vId.z);

			//DEBUG
			// PT_d("Gird", gId, pId);

			if (gId == -1) break;
			
			for (int c = 0; c < 27; c++)
			{
				int3 cId;
				cId.x = vId.x + offset_nq2[c][0];
				cId.y = vId.y + offset_nq2[c][1];
				cId.z = vId.z + offset_nq2[c][2];
				if (cId.x >= 0 && cId.y >= 0 && cId.z >= 0) 
				{ 	
					int cNumId = hash.getIndex(cId.x, cId.y, cId.z);
					int totalNum = hash.getCounter(cNumId);
					for (int i = 0; i < totalNum; i++) {
						int nbId = hash.getParticleId(cNumId, i);
						Coord pos_i = position[nbId];

						Real d_v0 = (pos_i - cap.v0).norm();
						Real d_v1 = (pos_i - cap.v1).norm();
						Real d_t = (pos_i - m).dot(s);
						Real d_line = fabs(((pos_i - m).cross(s)).norm() / d);
						Real min_d = (d_t < d * d && d_t > 0)? d_line : min(d_v0, d_v1);

						if (min_d < h)
						{
							// PT_f("MIN", min_d, pId);
							// printf("<id:%d dis:%f joint:%d>  pId:%d\n", nbId, min_d, cap.id_joint, pId);
							capPairs[cnt + start] = (Pair4f(nbId, -min_d, cap.id_joint, cap.id_cap));
							// PT_d("index", cnt + start, pId);
							cnt++;
						}
					}
				}
			}

			if (next_t > 1) break; //终点格子

			float tmp_t = -1;
			int3 next_c;	
			// 选取线段上最近Grid
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
					Coord time1 = (min_v - m) / s;
					Coord time2 = (max_v - m) / s;
					Coord min_t(min(time1[0],time2[0]), min(time1[1],time2[1]), min(time1[2],time2[2]));
					Coord max_t(max(time1[0],time2[0]), max(time1[1],time2[1]), max(time1[2],time2[2]));
					float mint = max(min_t[0], max(min_t[1], min_t[2]));
					float maxt = min(max_t[0], min(max_t[1], max_t[2]));
					if (mint > 0 && mint < 1 && mint < maxt && (mint < tmp_t || tmp_t < 0) && mint > next_t - eps)
					{
						tmp_t = maxt;
						next_c = cId;
						//DEBUG
						//if (cId == vId) printf("Error cId == vId\n");
					}
				}
			}

			next_t = tmp_t;
			// PT_f("time", next_t, pId);
			if (next_t < 0) break;
			vId = next_c;
		}
	}


	// max_cap 
	// 统计顶点数
	template<typename Pair4f>
	__global__ void K_CountPoint(
		DArray<Pair4f> capPairs,
		DArray<int> count)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= capPairs.size()) return;
		if( pId == capPairs.size() - 1 || int(capPairs[pId][0]) != int(capPairs[pId + 1][0]))
			count[pId] = 1;
		else 
		{
			// printf("[%d, %d] pId:%d\n", int(capPairs[pId][0]), int(capPairs[pId + 1][0]), pId);
			count[pId] = 0;
		}
		
	}

	// max_point
	// set out<关节, 顶点, 胶囊体>
	template<typename Pair4f, typename Pair3>
	__global__ void K_SetOutPair(
		DArray<Pair4f> capPairs,
		DArray<int> count,
		DArray<Pair3> outPairs)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= capPairs.size()) return;
		if( pId == capPairs.size() - 1 || int(capPairs[pId][0]) != int(capPairs[pId + 1][0]))
		{
			outPairs[count[pId]] = Pair3(capPairs[pId][2], capPairs[pId][0], capPairs[pId][3]);
			// printf("<joint:%d, point:%d>\n", outPairs[count[pId]][0], outPairs[count[pId]][1]);
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

		auto& outPairs = this->outPJPair()->getData();

		// uint numJt  = this->inJointSize()->getData();
		uint numPt  = this->inPosition()->getDataPtr()->size();
		uint numCp  = this->inCapsule()->getDataPtr()->size();
		uint sizeLimit = this->varSizeLimit()->getData();

		// Construct hash grid
		Reduction<Coord> reduce;
		Coord hiBound = reduce.maximum(points.begin(), points.size());
		Coord loBound = reduce.minimum(points.begin(), points.size());

		GridHash<TDataType> hashGrid;
		hashGrid.setSpace(h, loBound - Coord(h), hiBound + Coord(h));
		hashGrid.clear();
		hashGrid.construct(points);

		DArray<int> count;
		
		count.resize(numCp);
		cuExecute(numCp,
			K_CountNeighbor,
			points,
			hashGrid,
			h,
			capsules,
			count);
		cuSynchronize();

		int numPair = m_reduce.accumulate(count.begin(), count.size());
		m_scan.exclusive(count, true);

		//FixMe: numPair = 0?

		DArray<Pair4f>capJointPairs;
		capJointPairs.resize(numPair);
		cuExecute(numCp,
			K_ComputeNeighbor,
			points,
			hashGrid,
			h,
			capsules,
			capJointPairs,
			count);
		cuSynchronize();

		thrust::sort(thrust::device, capJointPairs.begin(), capJointPairs.begin() + capJointPairs.size());
		
		count.resize(numPair);
		cuExecute(numPair,
			K_CountPoint,
			capJointPairs,
			count);
		cuSynchronize();

		int numPoint = m_reduce.accumulate(count.begin(), count.size());
		m_scan.exclusive(count, true);

		assert(numPoint >= 0);

		outPairs.resize(numPoint);
		cuExecute(numPair,
			K_SetOutPair,
			capJointPairs,
			count,
			outPairs);
		cuSynchronize();

		capJointPairs.clear();
		count.clear();
		hashGrid.clear();
	}

	DEFINE_CLASS(NeighborPointQueryJoint);
}