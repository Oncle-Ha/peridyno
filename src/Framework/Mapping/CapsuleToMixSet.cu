#include "CapsuleToMixSet.h"

#include "Matrix/MatrixFunc.h"
#include "Topology/NeighborPointQueryJoint.h"

#include <iostream>

#define PT_d(s, x, y) printf("%s: %d  pId: %d\n", s, x, y)
#define PT_f(s, x, y) printf("%s: %f  pId: %d\n", s, x, y)
#define PT_e(s, y) printf("[%s]  pId: %d\n", s, y)

namespace dyno
{
	template<typename TDataType>
	CapsuleToMixSet<TDataType>::CapsuleToMixSet()
		: TopologyMapping()
	{

	}

	template<typename TDataType>
	CapsuleToMixSet<TDataType>::CapsuleToMixSet(JointList* from, std::shared_ptr<MixSet<TDataType>> to)
	{
		m_from = from;
		m_to = to;
	}

	template<typename TDataType>
	CapsuleToMixSet<TDataType>::~CapsuleToMixSet()
	{

	}


	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::initializeImpl()
	{
		match();
		return true;
	}

	template<typename Pair2, typename Mat4f, typename Triangle, typename Coord>
	__global__ void CM_ApplyTransformTri(
		DArray<Pair2> clusters,
		DArray<Mat4f> transform,
		DArray<Triangle> tris,
		DArray<Coord> to)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		auto& couple = clusters[pId];
		auto& tri = tris[couple[0]];
		for(int i = 0; i < 3; ++i)
		{
			Coord old_p = to[tri[i]];
			Vec4f tmp_p(old_p[0], old_p[1], old_p[2], 1);

			tmp_p = transform[couple[1]] * tmp_p;
			old_p[0] = tmp_p[0] / tmp_p[3];
			old_p[1] = tmp_p[1] / tmp_p[3];
			old_p[2] = tmp_p[2] / tmp_p[3];

			to[tri[i]] = old_p;
		}
	}

	template<typename Pair2, typename Mat4f, typename Tetrahedron, typename Coord>
	__global__ void CM_ApplyTransformTet(
		DArray<Pair2> clusters,
		DArray<Mat4f> transform,
		DArray<Tetrahedron> tets,
		DArray<Coord> to)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		auto& couple = clusters[pId];
		auto& tet = tets[couple[0]];
		for(int i = 0; i < 4; ++i)
		{
			Coord old_p = to[tet[i]];
			Vec4f tmp_p(old_p[0], old_p[1], old_p[2], 1);
			
			tmp_p = transform[couple[1]] * tmp_p;
			old_p[0] = tmp_p[0] / tmp_p[3];
			old_p[1] = tmp_p[1] / tmp_p[3];
			old_p[2] = tmp_p[2] / tmp_p[3];

			to[tet[i]] = old_p;
		}
	}	

	template<typename Pair2, typename Mat4f, typename Coord>
	__global__ void CM_ApplyTransformPoint(
		DArray<Pair2> clusters,
		DArray<Mat4f> transform,
		DArray<Coord> to)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		auto& couple = clusters[pId];
		int point_i = couple[1];
		
		// PT_d("point", point_i, pId);
		// PT_d("joint", couple[0], pId);
		Coord old_p = to[point_i];
		Vec4f tmp_p(old_p[0], old_p[1], old_p[2], 1);
		
		tmp_p = transform[couple[0]] * tmp_p;
		// PT_f("tmp3", tmp_p[3], pId);
		old_p[0] = tmp_p[0] / tmp_p[3];
		old_p[1] = tmp_p[1] / tmp_p[3];
		old_p[2] = tmp_p[2] / tmp_p[3];

		to[point_i] = old_p;
	}	

	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::apply()
	{
		uint pDim = cudaGridSize(m_to->getPoints().size(), BLOCK_SIZE);
        
		std::vector<Mat4f> p_transform;
		DArray<Mat4f> GlTransform;
		GlTransform.resize(m_from->size());

		for (int i = 0; i < m_from->size(); ++i)
		{
			auto& joint = (*m_from)[i];
			joint->getGlobalTransform();
			p_transform.push_back(joint->GlobalTransform * m_initInvTransform[i]);
		}
		GlTransform.assign(p_transform);

		// FIXME: 四面体和三角面片重合点会无法复原
		// TODO: 改为只控制顶点
        // Animation
		if(is_body)
		{
			cuExecute(m_triClusters.size(),
				CM_ApplyTransformTri,
				m_triClusters,
				GlTransform,
				m_to->getTriangles(),
				m_to->getAllPoints());
			cuSynchronize();
		
			cuExecute(m_tetClusters.size(),
				CM_ApplyTransformTet,
				m_tetClusters,
				GlTransform,
				m_to->getTetrahedrons(),
				m_to->getAllPoints());
			cuSynchronize();
		}
		else
		{
			cuExecute(m_pointClusters.size(),
				CM_ApplyTransformPoint,
				m_pointClusters,
				GlTransform,
				m_to->getAllPoints());
			cuSynchronize();
		}
		
		return true;
	}

	// 统计<几何体, 关节, 顶点>对数
	template<typename Pair2>
	__global__ void CM_CountPair3(
		DArray<Pair2> pairVerJoint,
		DArrayList<int> ver2X,
		DArray<int> count)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pairVerJoint.size()) return;
		count[pId] += ver2X[pairVerJoint[pId][1]].size();
	}

	// Set<几何体, 关节, 顶点>
	template<typename Pair2, typename Pair3>
	__global__ void CM_SetPair3(
		DArray<Pair2> pairVerJoint,
		DArrayList<int> ver2X,
		DArray<int> count,
		DArray<Pair3> pair3)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pairVerJoint.size()) return;
		
		int point = pairVerJoint[pId][1];
		int joint = pairVerJoint[pId][0];
		auto& list = ver2X[point];
		int size = list.size();
		for (int i = 0; i < size; ++i) {
			pair3[count[pId] + i] = Pair3(list[i], joint, point);
		}
	}

	// Count<几何体>
	template<typename Pair3>
	__global__ void CM_CountTopoList(
		DArray<Pair3> pair3,
		DArray<int> count)
	{
		int pId = threadIdx.x + (blockIdx.x + blockDim.x);
		if (pId >= pair3.size()) return;
		count[pId] =  (pId == pair3.size() - 1 || pair3[pId + 1][0] != pair3[pId][0]); // TODO: check m_scan
	}

	// Set 几何体:[<关节, 顶点>...]
	template<typename Pair3>
	__global__ void CM_SetTopoList(
		DArray<Pair3> pair3,
		DArray<int> count,
		DArrayList<Pair3> topoList)
	{
		int pId = threadIdx.x + (blockIdx.x + blockDim.x);
		if (pId >= pair3.size()) return;
		topoList[count[pId]].atomicInsert(pair3[pId]);
	}

	// 确定<几何体, 关节>
	template<typename Pair2, typename Pair3>
	__global__ void CM_ComputeJoint(
		DArrayList<Pair3> topolist,
		DArray<Pair2> pairJoint)
	{
		int pId = threadIdx.x + (blockIdx.x + blockDim.x);
		if (pId >= pairJoint.size()) return;
		auto& list = topolist[pId];

		int tmp_joint = list[0][1];
		int size_i = list.size();
		for (int i = 1; i < size_i; ++i)
		{
			if (list[i - 1][2] != list[i][2] && list[i - 1][1] == list[i][1])
				{
					tmp_joint = list[i][1];
					break;
				}
		}
		pairJoint[pId] = Pair2(list[0][0], tmp_joint);
	}

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::match()
	{
		// 获取胶囊体接触的顶点

		
		int for_cnt = 0;
		for (auto joint : *m_from)
		{
			++for_cnt;
			//FIXME : 坐标误差很大
			joint->getGlobalTransform();
			Vec4f tmp = joint->GlobalTransform * Vec4f(0, 0, 0, 1) ;
			joint->GlCoord = Coord(tmp[0] / tmp[3], tmp[1] / tmp[3], tmp[2] / tmp[3]);
			m_initInvTransform.push_back(joint->GlobalTransform.inverse());
			// std::cerr << for_cnt  <<" :  " <<joint->GlCoord[0] << ", " << joint->GlCoord[1] << ", " << joint->GlCoord[2] << "\n";
		}

		std::vector<JCapsule>capsule_list;
		int id_joint = 0;
		for (auto joint : *m_from)
		{
			int id_cap = 0;
			for (auto joint_son : joint->children)
			{
				capsule_list.push_back(JCapsule{id_joint, id_cap, 
												joint->GlCoord, joint_son->GlCoord});
				++id_cap;
				//DEBUG
				//float s = 1;
				//printf("(%f,%f,%f) -> (%f,%f,%f)\n",
				//	joint->GlCoord[0] / s, joint->GlCoord[1] / s, joint->GlCoord[2] / s,
				//	joint_son->GlCoord[0] / s, joint_son->GlCoord[1] / s, joint_son->GlCoord[2] / s);
			}
			++id_joint;
		}
 		auto nbQuery = std::make_shared<NeighborPointQueryJoint<TDataType>>();

		nbQuery->inRadius()->setValue(m_radius);
		nbQuery->inPosition()->allocate()->assign(m_to->getPoints());
		// nbQuery->inJointSize()->setValue(m_from->size());
		nbQuery->inCapsule()->allocate()->assign(capsule_list);

		nbQuery->update();

		// 获取顶点所对应的几何体
		// 返回DArray<Pair<几何体ID, JointID>> 
		DArray<Pair2> p_pairs2;

		p_pairs2.assign(nbQuery->outPJPair()->getData());

		DArray<int> count;
		
		auto fun_body = [=](DArrayList<int>& Ver2X, DArray<Pair2>& m_Clusters) mutable 
		{
			int pairVerJointSize = p_pairs2.size();
			count.resize(pairVerJointSize);

			cuExecute(pairVerJointSize,
				CM_CountPair3,
				p_pairs2,
				Ver2X,
				count);
			cuSynchronize();	

			int numPair3 = m_reduce.accumulate(count.begin(), count.size());
			m_scan.exclusive(count, true);
		
			DArray<Pair3> p_pairs3;
			p_pairs3.resize(numPair3);

			cuExecute(pairVerJointSize,
				CM_SetPair3,
				p_pairs2,
				Ver2X,
				count,
				p_pairs3);
			cuSynchronize();

			count.resize(numPair3);
			// sort?
			cuExecute(numPair3,
				CM_CountTopoList,
				p_pairs3,
				count);
			cuSynchronize();

			int numTopo = m_reduce.accumulate(count.begin(), count.size());			

			DArrayList<Pair3> p_pairsList;
			p_pairsList.resize(numTopo);

			cuExecute(numPair3,
				CM_SetTopoList,
				p_pairs3,
				count,
				p_pairsList);
			cuSynchronize();

			m_Clusters.resize(numTopo);

			cuExecute(numTopo,
				CM_ComputeJoint,
				p_pairsList,
				m_Clusters);
			cuSynchronize();
		};
		
		if (is_body && p_pairs2.size())
		{
			fun_body(m_to->getVer2Tri(), m_triClusters);
			fun_body(m_to->getVer2Tet(), m_tetClusters);
		}
		else
		{
			m_pointClusters.assign(p_pairs2);
		}
		


	}

	DEFINE_CLASS(CapsuleToMixSet);
}