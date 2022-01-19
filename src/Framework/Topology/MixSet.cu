#include "MixSet.h"

#include <string>
#include <sstream>

#include <thrust/sort.h>

namespace dyno
{
    IMPLEMENT_CLASS_1(MixSet, TDataType)

    template<typename TDataType>
    MixSet<TDataType>::MixSet()
    {

    }
    template<typename TDataType>
    MixSet<TDataType>::~MixSet()
    {

    }

	template<typename TDataType>
	void MixSet<TDataType>::setTriPoints(std::vector<Coord>& pos)
	{
		//printf("%d\n", pos.size());
        m_triPointSize = pos.size();
		m_coords.resize(pos.size());
		m_coords.assign(pos);

		tagAsChanged();
	}

	template<typename TDataType>
	void MixSet<TDataType>::setTetPoints(std::vector<Coord>& pos)
	{
        assert(m_triPointSize > 0);

        DArray<Coord> tmp;
        tmp.resize(m_triPointSize);
        tmp.assign(m_coords);

        m_tetPointSize = pos.size();
		m_coords.resize(m_triPointSize + m_tetPointSize);
		m_coords.assign(tmp, m_triPointSize, 0, 0);
		//FIXME bug: assign(vector)
        // m_coords.assign(pos, m_tetPointSize, m_triPointSize, 0); 

        tmp.resize(m_tetPointSize);
        tmp.assign(pos);	
		m_coords.assign(tmp, m_tetPointSize, m_triPointSize, 0);

		tagAsChanged();
	}

	template<typename TDataType>
	void MixSet<TDataType>::setTriangles(std::vector<Triangle>& triangles)
	{
		m_triangles.resize(triangles.size());
		m_triangles.assign(triangles);
	}

	template<typename TDataType>
	void MixSet<TDataType>::setTetrahedrons(std::vector<Tetrahedron>& tetrahedrons)
	{
		m_tethedrons.resize(tetrahedrons.size());
		m_tethedrons.assign(tetrahedrons);
	}

    template<typename TDataType>
	void MixSet<TDataType>::loadObjFile(std::string filename)
	{
		if (filename.size() < 5 || filename.substr(filename.size() - 4) != std::string(".obj")) {
			std::cerr << "Error: Expected OBJ file with filename of the form <name>.obj.\n";
			exit(-1);
		}

		std::ifstream infile(filename);
		if (!infile) {
			std::cerr << "Failed to open. Terminating.\n";
			exit(-1);
		}

		int ignored_lines = 0;
		std::string line;
		std::vector<Coord> vertList;
		std::vector<Triangle> faceList;
		while (!infile.eof()) {
			std::getline(infile, line);

			//.obj files sometimes contain vertex normals indicated by "vn"
			if (line.substr(0, 1) == std::string("v") && line.substr(0, 2) != std::string("vn")) {
				std::stringstream data(line);
				char c;
				Coord point;
				data >> c >> point[0] >> point[1] >> point[2];
				vertList.push_back(point);
			}
			else if (line.substr(0, 1) == std::string("f")) {
				std::stringstream data(line);
				char c;
				int v0, v1, v2;
				data >> c >> v0 >> v1 >> v2;
				faceList.push_back(Triangle(v0 - 1, v1 - 1, v2 - 1));
			}
			else {
				++ignored_lines;
			}
		}
		infile.close();

		setTriPoints(vertList);
		setTriangles(faceList);
	}
            

    template<typename TDataType>
	void MixSet<TDataType>::loadTetFile(std::string filename)
	{
		std::string filename_node = filename;	filename_node.append(".node");
		std::string filename_ele = filename;	filename_ele.append(".ele");

		std::ifstream infile_node(filename_node);
		std::ifstream infile_ele(filename_ele);
		if (!infile_node || !infile_ele) {
			std::cerr << "Failed to open the tetrahedron file. Terminating.\n";
			exit(-1);
		}

		std::string line;
		std::getline(infile_node, line);
		std::stringstream ss_node(line);

		int node_num;
		ss_node >> node_num;
		std::vector<Coord> nodes;
		for (int i = 0; i < node_num; i++)
		{
			std::getline(infile_node, line);
			std::stringstream data(line);
			int id;
			Coord v;
			data >> id >> v[0] >> v[1] >> v[2];
			nodes.push_back(v);
		}

		
		std::getline(infile_ele, line);
		std::stringstream ss_ele(line);

		int ele_num;
		ss_ele >> ele_num;
		std::vector<Tetrahedron> tets;
		for (int i = 0; i < ele_num; i++)
		{
			std::getline(infile_ele, line);
			std::stringstream data(line);
			int id;
			Tetrahedron tet;
			data >> id >> tet[0] >> tet[1] >> tet[2] >> tet[3];
			tet[0] -= 1;
			tet[1] -= 1;
			tet[2] -= 1;
			tet[3] -= 1;
			tets.push_back(tet);
		}

		setTetPoints(nodes);
        setTetrahedrons(tets);
	}


    template<typename TDataType>
    void MixSet<TDataType>::loadMixFile(std::string filename)
    {
        std::string filename_2d = filename; filename_2d.append("_2d.obj");
        std::string filename_3d = filename; filename_3d.append("_3d");
        loadObjFile(filename_2d);
        loadTetFile(filename_3d);

        this->getJointVer();
    }

    template<typename Coord, typename FKey>
    __global__ void MS_SetupKeys(
        DArray<int> ids,
        DArray<FKey> keys,
        DArray<Coord> coords,
        DArray<NodeType> nodetype,
        DArray<int> joints,
        int size2d)
    {
        int tId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (tId >= keys.size()) return;
        ids[tId] = tId;
        keys[tId] = FKey(coords[tId][0], coords[tId][1], coords[tId][2]);
        if(tId < size2d) nodetype[tId] = NodeType::TwoD;
        else nodetype[tId] = NodeType::ThreeD;
        joints[tId] = -1;
    }

    template<typename FKey>
    __global__ void MS_SetupJoints(
        DArray<int> ids,
        DArray<FKey> keys,
		DArray<NodeType> nodetype,
        DArray<int> joints)
    {
        int tId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (tId >= keys.size()) return;
        if(tId != 0 && keys[tId] == keys[tId - 1])
        {
            joints[ids[tId - 1]] = ids[tId];
            joints[ids[tId]] = ids[tId - 1];
            nodetype[ids[tId - 1]] = NodeType::Joint;
            nodetype[ids[tId]] = NodeType::Joint;
        }
    }

    template<typename TDataType>
    void MixSet<TDataType>::getJointVer()
    {
    
        uint coordSize = m_coords.size();

        m_joints.resize(coordSize);
        m_verType.resize(coordSize);    

		DArray<int> coordIds;
        DArray<FKey> coordKeys;

        coordIds.resize(coordSize);
        coordKeys.resize(coordSize);

        cuExecute(coordSize,
			MS_SetupKeys,
			coordIds,
            coordKeys,
            m_coords,
            m_verType,
            m_joints,
            m_triPointSize);

        thrust::sort_by_key(thrust::device, coordKeys.begin(), coordKeys.begin() + coordKeys.size(), coordIds.begin());

        cuExecute(coordSize,
			MS_SetupJoints,
			coordIds,
            coordKeys,
			m_verType,
            m_joints);
    }

    template<typename TDataType>
    void MixSet<TDataType>::copyFrom(MixSet<TDataType> mixSet) 
    {
        if (m_coords.size() != mixSet.m_coords.size())
		{
			m_coords.resize(mixSet.m_coords.size());
		}
		m_coords.assign(mixSet.m_coords);

        m_joints.resize(mixSet.m_joints.size());
        m_joints.assign(mixSet.m_joints);

        m_verType.resize(mixSet.m_verType.size());
        m_verType.assign(mixSet.m_verType);     

        m_triangles.resize(mixSet.m_triangles.size());
        m_triangles.assign(mixSet.m_triangles);   

        m_tethedrons.resize(mixSet.m_tethedrons.size());
        m_tethedrons.assign(mixSet.m_tethedrons);   

        m_tetPointSize = mixSet.m_tetPointSize;
        m_triPointSize = mixSet.m_triPointSize;
    }

    DEFINE_CLASS(MixSet);
}