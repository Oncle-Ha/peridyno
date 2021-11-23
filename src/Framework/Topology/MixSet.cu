#include "MixSet.h"

#include <string>

#include <thrust/sort.h>

namespace dyno
{
    IMPLEMENT_CLASS_1(MixSet, TDataType)

    template<typename TDataType>
    MixSet<TDataType>::MixSet()
    {
        m_ptSet = std::make_shared<PointSet<TDataType>>();
        m_triSet = std::make_shared<TriangleSet<TDataType>>();
        m_tetSet = std::make_shared<TetrahedronSet<TDataType>>();
    }
    template<typename TDataType>
    MixSet<TDataType>::~MixSet()
    {

    }
        
    template<typename TDataType>
    void MixSet<TDataType>::loadMixFile(std::string filename)
    {
        std::string filename_2d = filename; filename_2d.append("_2d.obj");
        std::string filename_3d = filename; filename_3d.append("_3d");
        this->m_tetSet->loadTetFile(filename_3d);
        this->m_triSet->loadObjFile(filename_2d);

        auto& coords2d = this->m_triSet->getPoints();
        auto& coords3d = this->m_tetSet->getPoints();

        this->updateAllPoints();
        this->getJointVer();
    }

    template<typename TDataType>
    void MixSet<TDataType>::updateAllPoints()
    {
        auto& coords2d = this->m_triSet->getPoints();
        auto& coords3d = this->m_tetSet->getPoints();
        auto& coords = this->m_ptSet->getPoints();
        auto p_size = coords2d.size() + coords3d.size();
        m_joints.resize(p_size);
        m_nodeType.resize(p_size);
		coords.resize(p_size);
		coords.assign(coords2d, coords2d.size(), 0, 0);
		coords.assign(coords3d, coords3d.size(), coords2d.size(), 0);
    }
    
    template<typename TDataType>
    void MixSet<TDataType>::update2DPoints()
    {
        auto& coords = this->m_ptSet->getPoints();
        auto& coords2d = this->m_triSet->getPoints();
        coords2d.assign(coords, coords2d.size(), 0, 0);
    }

    template<typename TDataType>
    void MixSet<TDataType>::update3DPoints()
    {
        auto& coords = this->m_ptSet->getPoints();
        auto& coords2d = this->m_triSet->getPoints();
        auto& coords3d = this->m_tetSet->getPoints();
        coords3d.assign(coords, coords3d.size(), 0, coords2d.size());
    }

    template<typename TDataType>
    void MixSet<TDataType>::scale(Real s)
    {
        this->m_ptSet->scale(s);
        this->m_triSet->scale(s);
        this->m_tetSet->scale(s);
    }

    template<typename TDataType>
    void MixSet<TDataType>::scale(Coord s)
    {
        this->m_ptSet->scale(s);
        this->m_triSet->scale(s);
        this->m_tetSet->scale(s);
    }

    template<typename TDataType>
    void MixSet<TDataType>::translate(Coord t)
    {
        this->m_ptSet->translate(t);
        this->m_triSet->translate(t);
        this->m_tetSet->translate(t);
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
        auto& coords = this->m_ptSet->getPoints();
        uint coordSize = coords.size();

		DArray<int> coordIds;
        DArray<FKey> coordKeys;

        coordIds.resize(coordSize);
        coordKeys.resize(coordSize);

        cuExecute(coordSize,
			MS_SetupKeys,
			coordIds,
            coordKeys,
            coords,
            m_nodeType,
            m_joints,
            this->m_triSet->getPoints().size());

        thrust::sort_by_key(thrust::device, coordKeys.begin(), coordKeys.begin() + coordKeys.size(), coordIds.begin());

        cuExecute(coordSize,
			MS_SetupJoints,
			coordIds,
            coordKeys,
			m_nodeType,
            m_joints);
    }

    template<typename TDataType>
    void MixSet<TDataType>::copyFrom(MixSet<TDataType> mixSet) 
    {
        this->m_tetSet->copyFrom(*(mixSet.m_tetSet.get()));
        this->m_triSet->copyFrom(*(mixSet.m_triSet.get()));
        this->m_ptSet->copyFrom(*(mixSet.m_ptSet.get()));
    }

    DEFINE_CLASS(MixSet);
}