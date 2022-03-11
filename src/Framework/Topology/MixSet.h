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
	class MixSet : public PointSet<TDataType>
	{
		DECLARE_CLASS_1(PointSet, TDataType)
	public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
        
		typedef typename TopologyModule::Edge Edge;
		typedef typename TopologyModule::Triangle Triangle;
		typedef typename TopologyModule::Tetrahedron Tetrahedron;

		MixSet();
		~MixSet();
        // load by suffix
        void loadMixFile(std::string filename);
        
        // return coords
        DArray<Coord>&  getTriPoints() { m_triCoordsTmp.resize(m_triPointSize); m_triCoordsTmp.assign(m_coords, m_triPointSize, 0, 0); return m_triCoordsTmp;}
        DArray<Coord>&  getTetPoints() { m_tetCoordsTmp.resize(m_tetPointSize); m_tetCoordsTmp.assign(m_coords, m_tetPointSize, 0, m_triPointSize); return m_tetCoordsTmp;}
        DArray<Coord>&  getAllPoints() { return m_coords;}

		// return Set
		// std::shared_ptr<PointSet<TDataType>>		getPointSet()		{ return this->getPoints(); }
		// std::shared_ptr<TriangleSet<TDataType>>		getTriangleSet()	{ return this->m_triSet;}
		// std::shared_ptr<TetrahedronSet<TDataType>>	getTetrahedronSet() { return this->m_tetSet;}

		// return DArray Tri, Tet
		DArray<Triangle>& getTriangles() { return m_triangles; }
		DArray<Tetrahedron>& getTetrahedrons() { return m_tetrahedron; }

		// return DArrayList Ver2
		DArrayList<int>& getVer2Edge() {return m_ver2Edge;}
		DArrayList<int>& getVer2Tri() {return m_ver2Tri;}
		DArrayList<int>& getVer2Tet() {return m_ver2Tet;}

        // get joint between 2D & 3D
        void getJointVer();

		void copyFrom(MixSet<TDataType> mixSet);

		// Point
		void setTriPoints(std::vector<Coord>& pos);

		void setTetPoints(std::vector<Coord>& pos);

		// Edge
		void setEdges();

		// Tri
		void loadObjFile(std::string filename);

		void setTriangles(std::vector<Triangle>& triangles);

		void copyFrom(TriangleSet<TDataType>& triangleSet);

		// Tet
		void loadTetFile(std::string filename);

		void setTetrahedrons(std::vector<Tetrahedron>& tetrahedrons);

		void copyFrom(TetrahedronSet<TDataType> tetSet);

		// Ver2T
		void setVer2T();
	protected:
		// [		PointSet		] : m_coords, m_pointNeighbors, m_verType
		// [		EdgeSet			] : all edge
		// [	Triangle	]		  : non-body
		// 			[	Tetrahedron	] : body
		//			[ inter ]		  
		
		// PointSet (2D:m_triPointSize, 3D:m_tetPointSize)
        DArray<int> m_intersections;    // 重合点对

		DArray<NodeType> m_verType; 

		// DArrayList<TopoNumber> m_ver2Topo; 
		DArray<Coord> m_triCoordsTmp;
		DArray<Coord> m_tetCoordsTmp;
		
		// Edge
		DArray<Edge> m_edges;			
		DArray<Triangle> m_triangles; 		 
		DArray<Tetrahedron> m_tetrahedron;
		
		DArrayList<int> m_ver2Edge;
		DArrayList<int> m_ver2Tri;
		DArrayList<int> m_ver2Tet;

		int m_tetPointSize;
		int m_triPointSize;
	};
}