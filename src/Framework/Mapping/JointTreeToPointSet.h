#pragma once
#include "Module/TopologyMapping.h"
#include "Topology/PointSet.h"
#include "Topology/JointTree.h"
#include "Topology/Cluster.h"

namespace dyno
{
	template<typename TDataType> class PointSet;
    template<typename TDataType> class JointTree;
    template<typename TDataType> class Cluster;

	template<typename TDataType>
	class JointTreeToPointSet : public TopologyMapping
	{
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;

		JointTreeToPointSet();

		~JointTreeToPointSet() override;

		bool apply() override;

		void match();

		void set(std::shared_ptr<PointSet<TDataType>> , 
			std::vector<std::shared_ptr<Cluster<TDataType>>>* ,
			std::vector<std::shared_ptr<JointTree<TDataType>>>*);

	// void UpdateAnimationFor(
	// 	DArray<int> indices,
	// 	DArray<Real> weights,
	// 	Mat4f Mt,
	// 	Mat4f Mtl,
	// 	Mat4f GlobalTransform,
	// 	int size2d,
	// 	int Size3d,
	// 	DArray<Coord> points);	
			
	protected:
		bool initializeImpl() override; 

	private:

        // PointSet -> Clusters -> Joint
		std::shared_ptr<PointSet<TDataType>> m_from = nullptr;
		int m_tetPointSize = 0; 
		int m_triPointSize = 0; 

		std::vector<std::shared_ptr<Cluster<TDataType>>>* m_clusters;

		std::vector<std::shared_ptr<JointTree<TDataType>>>* m_jointTree;// 顺序为DFS序
		
		// DArrayList<int> m_clusterIds;
		// DArrayList<Real> m_clusterWeights;
	};
}