/**
 * Copyright 2022 Xukun Luo
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      https://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#pragma once
#include "ParticleSystem/ParticleSystem.h"
#include "NeighborData.h"
#include "Topology/PointSet.h"
#include "Topology/JointTree.h"
#include "Topology/Cluster.h"
#include "Topology/RigidGenPoint.h"
#include "Mapping/CapsuleToMixSet.h"

namespace dyno
{
    template<typename> class PointSetToPointSet;

    /**
     * \class Dolphin
     * \brief 
     */
    template<typename TDataType>
	class Dolphin : public ParticleSystem<TDataType>
    {
        DECLARE_TCLASS(Dolphin, TDataType)

    public:
		typedef typename TDataType::Real Real;
		typedef typename TDataType::Coord Coord;
		typedef TPair<TDataType> NPair;

        Dolphin(std::string name = "default");
		virtual ~Dolphin();

        void loadMixFile(std::string filename);
        
		bool translate(Coord t) override;
		bool scale(Real s) override;

        std::shared_ptr<Node> getSurfaceNode() {return m_surfaceNode;}
	    std::shared_ptr<PointSetToPointSet<TDataType>> getTopologyMapping();
        

	public:
		DEF_VAR(Real, Horizon, 0.01, "Horizon");
        DEF_VAR(Real, Radius, 0.0075, "Capsule Radius");

		DEF_ARRAY_STATE(Coord, ReferencePosition, DeviceType::GPU, "Reference position");

		DEF_ARRAYLIST_STATE(int, NeighborIds, DeviceType::GPU, "Storing the ids for neighboring particles");

		DEF_ARRAYLIST_STATE(NPair, RestShape, DeviceType::GPU, "Storing neighbors");

        //DEBUG
        // DEF_ARRAY_STATE(Vec3f, Color, DeviceType::GPU, "Color of point");
        // DEF_INSTANCE_STATE(PointSet<TDataType>, Points, ""); // TODO: remove this

        // TODO: 修改用速度更新。
		DEF_ARRAY_OUT(Coord, V0, DeviceType::GPU, "");
		DEF_ARRAY_OUT(Coord, V1, DeviceType::GPU, "");

        //  Particle Color
		DEF_ARRAY_OUT(Vec3f, Color, DeviceType::GPU, "");

        // Capsule info
        DEF_ARRAY_IN(JCapsule, Bone, DeviceType::GPU, "Bone Capsule");
		DEF_ARRAY_IN(Quat<Real>, Rotate, DeviceType::GPU, "Bone Rotate");
		DEF_ARRAY_IN(Quat<Real>, Translate, DeviceType::GPU, "Bone Translate");        
        // DEF_ARRAY_IN(Coord, Velocity, DeviceType::GPU, "Bone Velocity");
        // DEF_ARRAY_IN(Coord, AngularVelocity, DeviceType::GPU, "Bone AngularVelocity");
        
        DEF_ARRAY_STATE(Coord, RigidPosition, DeviceType::GPU, "Rigid Position");
        
        // std::vector<std::shared_ptr<JointTree<TDataType>>> m_jointMap;
        // std::vector<std::shared_ptr<Cluster<TDataType>>> m_clusters;

    protected:

        void resetStates() override;

        void updateTopology() override;

    private:
        std::shared_ptr<Node> m_surfaceNode;
        // DArray<JCapsule> m_bone;
	};
}