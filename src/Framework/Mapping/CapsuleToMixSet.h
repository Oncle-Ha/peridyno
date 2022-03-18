#pragma once
#include "Module/TopologyMapping.h"
#include "Topology/PointSet.h"
#include "Topology/JointTree.h"
#include "Topology/MixSet.h"

namespace dyno
{
	template<typename TDataType> class PointSet;

	template<typename TDataType>
	class CapsuleToMixSet : public TopologyMapping
	{
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef typename TDataType::Matrix Mat3;
		typedef typename TopologyModule::Pair2 Pair2;
		typedef typename TopologyModule::Pair3 Pair3;
		typedef typename Mat4f Mat;

		typedef std::vector<std::shared_ptr<JointTree<typename TDataType>>> JointList;

		CapsuleToMixSet();
		CapsuleToMixSet(JointList* from, std::shared_ptr<MixSet<TDataType>> to);
		~CapsuleToMixSet() override;

		void setCapsuleRadius(Real r) { m_radius = r; }

		void setFrom(JointList* from) { m_from = from; }
		void setTo(std::shared_ptr<MixSet<TDataType>> to) { m_to = to; }

		bool apply() override;

		bool shapeMatch();

		void match();
		void animByQuat();
		void animByMatrix();

		void animByVir(); //更新虚粒子带来的弹力

		/**
		 * @brief Color of neighboring particles
		 */
		DEF_ARRAY_OUT(Vec3f, Color, DeviceType::GPU, "Return neighbor ids");

		/**
		 * @brief Velocity of neighboring particles
		 */
		DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "");

		/**
		 * @brief Force of neighboring particles
		 */
		DEF_ARRAY_IN(Coord, Force, DeviceType::GPU, "");

		DEF_VAR_IN(Real, TimeStep, "");
		
	protected:
		bool initializeImpl() override;

	private:
		Reduction<int> m_reduce;
		Scan m_scan;
		bool is_body = false;
		//Searching radius
		Real m_radius = 0.125;

	
		//映射关系
		DArray<Pair2> m_tetClusters; 
		DArray<Pair2> m_triClusters;
		DArray<Pair2> m_pointClusters;

		// shape match
		DArray<Coord> m_initCoord; //r_i

		JointList* m_from = nullptr;
		// std::shared_ptr<JointTree<TDataType>> m_from = nullptr;
		std::shared_ptr<MixSet<TDataType>> m_to = nullptr;
		DArray<Coord> m_virtualCoord; // 虚拟点坐标
	
		std::vector<Mat> m_initInvTransform; // X2 = M2 M1^-1 X1
		DArray<Quat<Real>> m_initQuatT;			
		DArray<Quat<Real>> m_initQuatR;
		DArray<Real> m_initS;
		
		// JointList* m_initfrom = nullptr;
		// std::shared_ptr<JointTree<TDataType>> m_initFrom = nullptr;
		// std::shared_ptr<MixSet<TDataType>> m_initTo = nullptr;
	};
}
