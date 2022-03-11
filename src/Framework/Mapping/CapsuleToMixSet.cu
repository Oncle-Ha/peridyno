#include "CapsuleToMixSet.h"

#include "Matrix/MatrixFunc.h"
#include "Topology/NeighborPointQueryJoint.h"

#include <iostream>

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

	template<typename Pair2, typename Mat, typename Triangle, typename Coord>
	__global__ void CM_ApplyTransformTri(
		DArray<Pair2> clusters,
		DArray<Mat> transform,
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

	template<typename Pair2, typename Mat, typename Tetrahedron, typename Coord>
	__global__ void CM_ApplyTransformTet(
		DArray<Pair2> clusters,
		DArray<Mat> transform,
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

	template<typename Pair2, typename Mat, typename Coord>
	__global__ void CM_ApplyTransformPoint(
		DArray<Pair2> clusters,
		DArray<Mat> transform,
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


	template<typename Pair2, typename Coord, typename Real, template<typename Real> typename Quat>
	__global__ void CM_ApplyTransformPointByQuat(
		DArray<Pair2> clusters,
		DArray<Quat<Real>> initT,
		DArray<Quat<Real>> initR,
		DArray<Real> initS,
		DArray<Quat<Real>> T,
		DArray<Quat<Real>> R,
		DArray<Real> S,
		DArray<Coord> to)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		auto& couple = clusters[pId];
		int p = couple[1]; // point
		int j = couple[0]; // joint
		// PT_d("point", point_i, pId);
		// PT_d("joint", couple[0], pId);
		Coord old_p = to[p];
		Quat<Real> tmp_p(old_p[0], old_p[1], old_p[2], 0);
		
		tmp_p = (1. / initS[j]) * initR[j].conjugate() * (tmp_p - initT[j]) * initR[j];
		tmp_p = S[j] * R[j] * tmp_p * R[j].conjugate() + T[j];

		// PT_f("tmp3", tmp_p[3], pId);

		to[p] = Coord(tmp_p.x, tmp_p.y, tmp_p.z);
	}	

	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::apply()
	{
		// uint pDim = cudaGridSize(m_to->getPoints().size(), BLOCK_SIZE);
	
		
		// Matrix
		/*
		{
			std::vector<Mat> p_transform;
			std::vector<Mat> p_nextInvTransform;
			DArray<Mat> GlTransform;
			GlTransform.resize(m_from->size());
			for (int i = 0; i < m_from->size(); ++i)
			{
				auto& joint = (*m_from)[i];
				joint->getGlobalTransform();
				p_nextInvTransform.push_back(joint->GlobalTransform.inverse());
				p_transform.push_back(joint->GlobalTransform * m_initInvTransform[i]);
			}
			GlTransform.assign(p_transform);
		
			// Animation(Matrix)
			if(is_body)
			{
				// FIXME: 四面体和三角面片重合点会无法复原
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

			m_initInvTransform.assign(p_nextInvTransform.begin(), p_nextInvTransform.end());
		}
		*/
		
		// Quat
		{
			std::vector<Quat<Real>> p_T;
			std::vector<Quat<Real>> p_R;
			std::vector<Real> p_S;

			DArray<Quat<Real>> new_T;
			DArray<Quat<Real>> new_R;
			DArray<Real> new_S;
			
			new_T.resize(m_from->size());
			new_R.resize(m_from->size());
			new_S.resize(m_from->size());
			for (int i = 0; i < m_from->size(); ++i)
			{
				auto& joint = (*m_from)[i];
				joint->getGlobalQuat();
				p_T.push_back(joint->GlT);
				p_R.push_back(joint->GlR);
				p_S.push_back(joint->GlS);
			}
			new_T.assign(p_T);
			new_R.assign(p_R);
			new_S.assign(p_S);

			// Animation(Quat)
			if(is_body)
			{
				// // FIXME: 四面体和三角面片重合点会无法复原
				// cuExecute(m_triClusters.size(),
				// 	CM_ApplyTransformTri,
				// 	m_triClusters,
				// 	GlTransform,
				// 	m_to->getTriangles(),
				// 	m_to->getAllPoints());
				// cuSynchronize();
			
				// cuExecute(m_tetClusters.size(),
				// 	CM_ApplyTransformTet,
				// 	m_tetClusters,
				// 	GlTransform,
				// 	m_to->getTetrahedrons(),
				// 	m_to->getAllPoints());
				// cuSynchronize();
			}
			else
			{
				cuExecute(m_pointClusters.size(),
					CM_ApplyTransformPointByQuat,
					m_pointClusters,
					m_initQuatT,
					m_initQuatR,
					m_initS,
					new_T,
					new_R,
					new_S,
					m_to->getAllPoints());
				cuSynchronize();
			}
			m_initQuatT.assign(new_T);
			m_initQuatR.assign(new_R);
			m_initS.assign(new_S);
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
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
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
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= pair3.size()) return;
		topoList[count[pId]].atomicInsert(pair3[pId]);
	}

	// 确定<几何体, 关节>
	template<typename Pair2, typename Pair3>
	__global__ void CM_ComputeJoint(
		DArrayList<Pair3> topolist,
		DArray<Pair2> pairJoint)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
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

	template<typename Vec3f>
	__global__ void CM_UpdateColor1(
		DArray<Vec3f> colors)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= colors.size()) return;
		colors[pId] = Vec3f(0.f, 0.f, 0.f);
	}
	
	template<typename Vec3f, typename Pair2>
	__global__ void CM_UpdateColor2(
		DArray<Vec3f> colors,
		DArray<Pair2> clusters)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		colors[clusters[pId][1]] = Vec3f(0.f,1.f,0.f);
	}

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::match()
	{
		// 获取胶囊体接触的顶点

		m_initQuatT.resize(m_from->size());
		m_initQuatR.resize(m_from->size());
		m_initS.resize(m_from->size());

		std::vector<Quat<Real>> p_T;			
		std::vector<Quat<Real>> p_R;
		std::vector<Real> p_S;

		int for_cnt = 0;
		for (auto joint : *m_from)
		{
			++for_cnt;
			// Matrix
			joint->getGlobalTransform();
			m_initInvTransform.push_back(joint->GlobalTransform.inverse());// 求逆误差? 
			Coord cm = joint->getCoordByMatrix(Coord(0,0,0));	

			// Quat
			joint->getGlobalQuat();
			p_T.push_back(joint->GlT);
			p_R.push_back(joint->GlR);
			p_S.push_back(joint->GlS);
			Coord qm = joint->getCoordByQuat(Coord(0,0,0));

			joint->GlCoord = qm;
			// std::cerr << for_cnt  <<" (M) :  " <<cm[0] << ", " << cm[1] << ", " << cm[2] << "\n";
			//std::cerr << for_cnt  <<" (Q) :  " <<qm[0] << ", " << qm[1] << ", " << qm[2] << "\n";
			std::cerr << for_cnt  <<" :  " <<joint->GlCoord[0] << ", " << joint->GlCoord[1] << ", " << joint->GlCoord[2] << "\n";
		}
		
		//Quat
		m_initQuatT.assign(p_T);
		m_initQuatR.assign(p_R);
		m_initS.assign(p_S);

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
				// float s = 1;
				// printf("(%f,%f,%f) -> (%f,%f,%f)\n",
				// 	joint->GlCoord[0] / s, joint->GlCoord[1] / s, joint->GlCoord[2] / s,
				// 	joint_son->GlCoord[0] / s, joint_son->GlCoord[1] / s, joint_son->GlCoord[2] / s);
			}
			++id_joint;
		}
 		auto nbQuery = std::make_shared<NeighborPointQueryJoint<TDataType>>();

		nbQuery->inRadius()->setValue(m_radius);
		nbQuery->inPosition()->allocate()->assign(m_to->getPoints());
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
			int numPoint = m_to->getAllPoints().size();

			if (this->outColor()->isEmpty())
				this->outColor()->allocate();
			auto& out_color = this->outColor()->getData();

			out_color.resize(numPoint);
			cuExecute(numPoint,
				CM_UpdateColor1,
				out_color);
			cuSynchronize();

			printf("Num of Clusters:%d\n", m_pointClusters.size());
			cuExecute(m_pointClusters.size(),
				CM_UpdateColor2,
				out_color,
				m_pointClusters);
			cuSynchronize();			
		}
		


	}

	DEFINE_CLASS(CapsuleToMixSet);
}