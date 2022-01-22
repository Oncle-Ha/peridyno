#include "JointTreeToPointSet.h"

#include "Matrix/MatrixFunc.h"
#include "Topology/NeighborPointQuery.h"

#include <iostream>

namespace dyno
{
	template<typename TDataType>
	JointTreeToPointSet<TDataType>::JointTreeToPointSet()
		: TopologyMapping()
	{

	}

	template<typename TDataType>
	void JointTreeToPointSet<TDataType>::set(
			std::shared_ptr<PointSet<TDataType>> from,
			std::vector<std::shared_ptr<Cluster<TDataType>>>* clusters,
			std::vector<std::shared_ptr<JointTree<TDataType>>>* jointTree)
	{
		m_from = from;
		m_clusters = clusters;
		m_jointTree = jointTree;
	}

	template<typename TDataType>
	JointTreeToPointSet<TDataType>::~JointTreeToPointSet()
	{

	}



	// 对每个关节所控制的点集做动画更新
	template <typename Real, typename Coord>
	__global__ void UpdateAnimation(
		DArray<int> indices,
		DArray<Real> weights,
		Mat4f Mt,
		Mat4f Mtl,
		Mat4f GlobalTransform,
		DArray<Coord> points)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= indices.size()) return;	
		
		Coord old_p = points[indices[pId]];

		// TODO:Update Coord
		// old_p * Mt * Mtl^-1 * GlobalTransform
		
		points[indices[pId]] = old_p;
	}

	// 初始化每个点所受控制的关节关系
	template <typename Real>
	__global__ void InitMatch(
		DArray<int> indices,
		DArray<Real> weights,
		DArrayList<int> clusterIds,
		DArrayList<Real> clusterWeights,
		int clusterId) 
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= indices.size()) return;
		int pointId = indices[pId];
		clusterWeights[pointId].atomicInsert(weights[pId]);
		clusterIds[pointId].atomicInsert(clusterId);
	}

	template<typename TDataType>
	bool JointTreeToPointSet<TDataType>::initializeImpl()
	{
		match();
		return true;
	}

	template<typename TDataType>
	bool JointTreeToPointSet<TDataType>::apply()
	{
		for (auto joint : *m_jointTree)
		{
			joint->getGlobalTransform(); // DFS 遍历
		}

		std::shared_ptr<Cluster<TDataType>> v;
		for (int i = 0; i < this->m_clusters->size(); i++)
		{
			v = (*(this->m_clusters))[i];
			uint pDim = v->m_indices.size();
			cuExecute(pDim,
				UpdateAnimation,
				v->m_indices,
				v->m_weights,
				v->m_transform,
				v->m_transformLink,
				(*m_jointTree)[v->m_jointIndex]->GlobalTransform,
				// Mat4f(0.0),
				m_from->getPoints());
			cuSynchronize();
		}


		// for (auto v : *(this->m_clusters))
		// {
		// 	uint pDim = v->m_indices.size();
		// 	cuExecute(pDim,
		// 		UpdateAnimation,
		// 		v->m_indices,
		// 		v->m_weights,
		// 		v->m_transform,
		// 		v->m_transformLink,
		// 		(*m_jointTree)[v->m_jointIndex]->getGlobalTransform(),
		// 		this->m_from->getPoints());
		// 	cuSynchronize();
		// }
		return true;
	}



	template<typename TDataType>
	void JointTreeToPointSet<TDataType>::match()
	{
	}

	DEFINE_CLASS(JointTreeToPointSet);
}