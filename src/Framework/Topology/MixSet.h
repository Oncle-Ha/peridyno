#pragma once
#include "TopologyConstants.h"
#include "Module/TopologyModule.h"

#include "TetrahedronSet.h"
#include "TriangleSet.h"
#include "PointSet.h"

namespace dyno
{
    enum NodeType
    {
        TwoD = 0,
        Joint,
        ThreeD,
    };   
    
    class FKey
	{
	public:
		DYN_FUNC FKey()
		{
			id[0] = EMPTY;
			id[1] = EMPTY;
			id[2] = EMPTY;
		}

		DYN_FUNC FKey(double v0, double v1, double v2)
		{
			id[0] = v0;
			id[1] = v1;
			id[2] = v2;

			swap(id[0], id[1]);
			swap(id[0], id[2]);
			swap(id[1], id[2]);
		}

		DYN_FUNC inline double operator[] (unsigned int i) { return id[i]; }
		DYN_FUNC inline double operator[] (unsigned int i) const { return id[i]; }

		DYN_FUNC inline bool operator>= (const FKey& other) const {
			if (id[0] >= other.id[0]) return true;
			if (id[0] == other.id[0] && id[1] >= other.id[1]) return true;
			if (id[0] == other.id[0] && id[1] == other.id[1] && id[2] >= other.id[2]) return true;

			return false;
		}

		DYN_FUNC inline bool operator> (const FKey& other) const {
			if (id[0] > other.id[0]) return true;
			if (id[0] == other.id[0] && id[1] > other.id[1]) return true;
			if (id[0] == other.id[0] && id[1] == other.id[1] && id[2] > other.id[2]) return true;

			return false;
		}

		DYN_FUNC inline bool operator<= (const FKey& other) const {
			if (id[0] <= other.id[0]) return true;
			if (id[0] == other.id[0] && id[1] <= other.id[1]) return true;
			if (id[0] == other.id[0] && id[1] == other.id[1] && id[2] <= other.id[2]) return true;

			return false;
		}

		DYN_FUNC inline bool operator< (const FKey& other) const {
			if (id[0] < other.id[0]) return true;
			if (id[0] == other.id[0] && id[1] < other.id[1]) return true;
			if (id[0] == other.id[0] && id[1] == other.id[1] && id[2] < other.id[2]) return true;

			return false;
		}

		DYN_FUNC inline bool operator== (const FKey& other) const {
			return id[0] == other.id[0] && id[1] == other.id[1] && id[2] == other.id[2];
		}

		DYN_FUNC inline bool operator!= (const FKey& other) const {
			return id[0] != other.id[0] || id[1] != other.id[1] || id[2] != other.id[2];
		}

	private:
		DYN_FUNC inline void swap(double& v0, double& v1)
		{
			double vt = v0;
			v0 = v0 < v1 ? v0 : v1;
			v1 = vt < v1 ? v1 : vt;
		}

		double id[3];
	};
    
	template<typename TDataType>
	class MixSet : public TopologyModule
	{
		DECLARE_CLASS_1(PointSet, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
        
		typedef typename TopologyModule::Triangle Triangle;
		typedef typename TopologyModule::Tetrahedron Tetrahedron;

		MixSet();
		~MixSet();
        // load by suffix
        void loadMixFile(std::string filename);
        
        // return coords
        DArray<Coord>&  get2DPoints() { return this->m_triSet->getPoints();}
        DArray<Coord>&  get3DPoints() { return this->m_tetSet->getPoints();}
        DArray<Coord>&  getAllPoints() { return this->m_ptSet->getPoints();}

		// return Set
		std::shared_ptr<PointSet<TDataType>>		getPointSet()		{ return this->m_ptSet; }
		std::shared_ptr<TriangleSet<TDataType>>		getTriangleSet()	{ return this->m_triSet;}
		std::shared_ptr<TetrahedronSet<TDataType>>	getTetrahedronSet() { return this->m_tetSet;}

		// Transform
		void scale(Real s);
		void scale(Coord s);
		void translate(Coord t);

        // sync coords
        void update2DPoints();
        void update3DPoints();
        void updateAllPoints();

        // get joint between 2D & 3D
        void getJointVer();

		void copyFrom(MixSet<TDataType> mixSet);

	protected:
        std::shared_ptr<TriangleSet<TDataType>> m_triSet;       // 2D
        std::shared_ptr<TetrahedronSet<TDataType>> m_tetSet;    // 3D 分别执行碰撞等

        std::shared_ptr<PointSet<TDataType>> m_ptSet; 		   // 2D&3D所有顶点用Peridynamic做弹性
        DArray<NodeType> m_nodeType;           // 顶点类型

        DArray<int> m_joints;    // 2D和3D交界点的对应点
	};
}