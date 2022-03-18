#include "CapsuleToMixSet.h"

#include "Matrix/MatrixFunc.h"
#include "Topology/NeighborPointQueryJoint.h"

#include <iostream>

namespace dyno
{
	const Vec3f default_color[] ={
		Vec3f(219, 68, 83) / 255.0,
		Vec3f(233, 87, 62) / 255.0,
		Vec3f(246, 187, 67) / 255.0,
		Vec3f(93, 156, 236) / 255.0,
		Vec3f(140, 192, 81) / 255.0,
		Vec3f(54, 188, 155) / 255.0,
		Vec3f(79, 192, 232) / 255.0,
		Vec3f(101, 109, 120) / 255.0,
		Vec3f(230, 233, 238) / 255.0,
		Vec3f(172, 146, 237) / 255.0,
		Vec3f(236, 135, 191) / 255.0,
		Vec3f(170, 178, 189) / 255.0};
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
		DArray<Coord> to,
		DArray<Coord> vel,
		Real dt)
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
		Coord new_p = Coord(tmp_p.x, tmp_p.y, tmp_p.z);
		//DEBUG
		// if (pId == 0)
			// printf("Vel Anim: [%f, %f, %f]", vel[p][0], vel[p][1], vel[p][2]);
		// vel[p] = (new_p - old_p) ;
		to[p] = new_p;
	}	

	template<typename Pair2, typename Coord, typename Real, template<typename Real> typename Quat>
	__global__ void CM_ApplyOriginTransformPointByQuat(
		DArray<Pair2> clusters,
		DArray<Quat<Real>> T,
		DArray<Quat<Real>> R,
		DArray<Real> S,
		DArray<Coord> to,
		DArray<Coord> newTo)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		auto& couple = clusters[pId];
		int p = pId; // point
		int j = couple[0]; // joint
		// PT_d("point", point_i, pId);
		// PT_d("joint", couple[0], pId);
		Coord old_p = to[p];
		Quat<Real> tmp_p(old_p[0], old_p[1], old_p[2], 0);
		tmp_p = S[j] * R[j] * tmp_p * R[j].conjugate() + T[j];
		// PT_f("tmp3", tmp_p[3], pId);
		Coord new_p = Coord(tmp_p.x, tmp_p.y, tmp_p.z);
		newTo[p] = new_p;
	}	

	template<typename Pair2, typename Coord, typename Real>
	__global__ void CM_UpdateForceByVir(
		DArray<Pair2> clusters,
		DArray<Coord> to,
		DArray<Coord> virTo,
		DArray<Coord> force,
		Real k)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		auto& couple = clusters[pId];
		int vir_p = pId; 	// vir_p
		int p = couple[1]; // point
		int j = couple[0]; // joint	
		force[p] = (virTo[vir_p] - to[p]) *k;
		if(pId == 0)
		{
			printf("virTo = (%f, %f, %f)\n", virTo[vir_p][0], virTo[vir_p][1], virTo[vir_p][2]);
			printf("To = (%f, %f, %f)\n", to[p][0], to[p][1], to[p][2]);
			printf("Force[%d] = (%f, %f, %f)\n", p, force[p][0], force[p][1], force[p][2]);
		}
	}

	// 计算重心坐标
	template<typename Coord, typename Pair2>
	__global__ void CM_CountClusterCoord(
		DArray<Coord> to,
		DArray<Pair2> clusters,
		DArray<int> sizeCluster,
		DArray<Coord> coordCluster)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		atomicAdd(&sizeCluster[clusters[pId][0]], 1);
		// Coord
		for (int i = 0; i < 3; ++i)
			atomicAdd(&coordCluster[clusters[pId][0]][i], to[clusters[pId][1]][i]);
	}

	// 计算A矩阵
	template<typename Coord, typename Pair2, typename Matrix>
	__global__ void CM_ComputeMatrixA(
		DArray<Coord> to,
		DArray<Pair2> clusters,
		DArray<int> sizeCluster,
		DArray<Coord> coordCluster,
		DArray<Coord> directDis,
		DArray<Matrix> matA1,
		DArray<Matrix> matA2)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		int point = clusters[pId][1];
		int joint = clusters[pId][0];
		Coord center = coordCluster[joint] / (1.0 * sizeCluster[joint]);
		//DEBUG
		if(pId == -1)
			printf("New center: [%f, %f, %f]\n", center[0], center[1], center[2]);

		Matrix a1 = Matrix(0);
		Matrix a2 = Matrix(0);
		Coord q, p;
		q = to[point] - center;
		p = directDis[pId];
			a1(0, 0) = q[0] * p[0]; a1(0, 1) = q[0] * p[1]; a1(0, 2) = q[0] * p[2];
			a1(1, 0) = q[1] * p[0]; a1(1, 1) = q[1] * p[1]; a1(1, 2) = q[1] * p[2];
			a1(2, 0) = q[2] * p[0]; a1(2, 1) = q[2] * p[1]; a1(2, 2) = q[2] * p[2];	
			
		q = directDis[pId];
		p = directDis[pId];
			a2(0, 0) = q[0] * p[0]; a2(0, 1) = q[0] * p[1]; a2(0, 2) = q[0] * p[2];
			a2(1, 0) = q[1] * p[0]; a2(1, 1) = q[1] * p[1]; a2(1, 2) = q[1] * p[2];
			a2(2, 0) = q[2] * p[0]; a2(2, 1) = q[2] * p[1]; a2(2, 2) = q[2] * p[2];		

		// Mat3
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
			{
				atomicAdd(&(matA1[joint](i, j)), a1(i, j));
				atomicAdd(&(matA2[joint](i, j)), a2(i, j));
			}
	}

	// R = Polar(A)
	template<typename Matrix>
	__global__ void CM_ComputeMatrixR(
		DArray<Matrix> matA1,
		DArray<Matrix> matA2)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= matA1.size()) return;
		Matrix A =  matA1[pId] * matA2[pId].inverse();
		Matrix R(0), U(0), D(0), V(0);
	
		
		polarDecomposition(A, R, U, D, V);
		matA1[pId] = R;
		//DEBUG 
		if (pId == -1)
		{
			printf("Matrix R:\n");
			for (int i = 0; i < 3; ++i)
				printf("[%f, %f, %f]\n", R(i, 0), R(i, 1), R(i, 2));	
			
			printf("Matrix A:\n");
			for (int i = 0; i < 3; ++i)
				printf("[%f, %f, %f]\n", A(i, 0), A(i, 1), A(i, 2));					
		}
			
	}

	// Shape Mathc更新点坐标
	template<typename Coord, typename Pair2, typename Matrix, typename NodeType>
	__global__ void CM_ShapeMatchPos(
		DArray<Coord> to,
		DArray<NodeType> nodetype,
		DArray<Coord> vel,
		DArray<Pair2> clusters,
		DArray<int> sizeCluster,
		DArray<Coord> coordCluster,
		DArray<Coord> directDis,
		DArray<Matrix> matR,
		Real dt)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		int point = clusters[pId][1];
		int joint = clusters[pId][0];

		if (nodetype[point] == NodeType::TwoD) return;

		Coord center = coordCluster[joint] / (1.0 * sizeCluster[joint]);
		if (joint >= matR.size())
			printf("Error\n");
		Coord new_to = center + matR[joint] * directDis[pId];
		//DEBUG
		if (pId == 0)
			printf("Vel match: [%f, %f, %f]", vel[point][0], vel[point][1], vel[point][2]);
		vel[point] += (new_to - to[point]) / dt;
		to[point] = new_to;

		//DEBUG
		if(pId == -1)
		{
			printf("Now center: [%f, %f, %f]\n", center[0], center[1], center[2]);
			printf("MatR:\n");
			for (int i = 0; i < 3; ++i)
				printf("[%f, %f, %f]\n", matR[joint](i, 0), matR[joint](i, 1), matR[joint](i, 2));
			printf("directDis: [%f, %f, %f]\n", directDis[pId][0], directDis[pId][1], directDis[pId][2]);
			printf("Now To: [%f, %f, %f]\n", to[point][0], to[point][1], to[point][2]);
		}
		
	}

	

	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::shapeMatch()
	{
		int numCluster = m_from->size();
		
		DArray<int> size_cluster;
		DArray<Coord> coord_cluster;

		int numPair = m_pointClusters.size();

		size_cluster.resize(numCluster);
		coord_cluster.resize(numCluster);
		size_cluster.reset();
		coord_cluster.reset();
		cuExecute(numPair,
			CM_CountClusterCoord,
			m_to->getPoints(),
			m_pointClusters,
			size_cluster,
			coord_cluster);
		cuSynchronize();

		DArray<Mat3> A1;
		DArray<Mat3> A2;
		A1.resize(numCluster);
		A2.resize(numCluster);
		A1.reset();
		A2.reset();
		cuExecute(numPair,
			CM_ComputeMatrixA,
			m_to->getPoints(),
			m_pointClusters,
			size_cluster,
			coord_cluster,
			m_initCoord,
			A1,
			A2);
		cuSynchronize();

		cuExecute(numCluster,
			CM_ComputeMatrixR,
			A1,
			A2);
		cuSynchronize();

		cuExecute(numPair,
			CM_ShapeMatchPos,
			m_to->getPoints(),
			m_to->getVerType(),
			this->inVelocity()->getData(),
			m_pointClusters,
			size_cluster,
			coord_cluster,
			m_initCoord,
			A1,
			this->inTimeStep()->getData());
		cuSynchronize();

		size_cluster.clear();
		coord_cluster.clear();
		A1.clear();
		A2.clear();
	}

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::animByMatrix()
	{
		// Matrix
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

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::animByQuat()
	{
		// Quat
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
				m_to->getAllPoints(),
				this->inVelocity()->getData(),
				this->inTimeStep()->getData());
			cuSynchronize();
		}
		m_initQuatT.assign(new_T);
		m_initQuatR.assign(new_R);
		m_initS.assign(new_S);

		new_T.clear();
		new_R.clear();
		new_S.clear();
	}

	template<typename TDataType>
	void CapsuleToMixSet<TDataType>::animByVir()
	{
		int numPair = m_pointClusters.size();
		DArray<Coord> p_newCoord;
		p_newCoord.resize(numPair);

		// Quat更新虚粒子坐标
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
			
			
			cuExecute(numPair,
				CM_ApplyOriginTransformPointByQuat,
				m_pointClusters,
				new_T,
				new_R,
				new_S,
				m_virtualCoord,
				p_newCoord);
			cuSynchronize();

			new_T.clear();
			new_R.clear();
			new_S.clear();		
		}

		// 虚粒子更新实粒子力
		{
			cuExecute(numPair,
				CM_UpdateForceByVir,
				m_pointClusters,
				m_to->getAllPoints(),
				p_newCoord,
				this->inForce()->getData(),
				400);
			cuSynchronize();
		}

		p_newCoord.clear();
	}

	template<typename TDataType>
	bool CapsuleToMixSet<TDataType>::apply()
	{
		// Shape match
		// shapeMatch();		
	
		// Matrix
		// animByMatrix();

		// Quat
		// animByQuat();
		
		// Vir
		animByVir();

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
		DArray<Pair2> clusters,
		DArray<Vec3f> default_color)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		// colors[clusters[pId][1]] = Vec3f(0.f, 1.f, 0.f);
		colors[clusters[pId][1]] = default_color[clusters[pId][0] % 12];
	}


	template<typename Coord, typename Pair2>
	__global__ void CM_CountPointCoord(
		DArray<Coord> initTo,
		DArray<Pair2> clusters,
		DArray<int> sizeCluster,
		DArray<Coord> coordCluster,
		DArray<Coord> directDis)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		Coord center = coordCluster[clusters[pId][0]] / (1.0 * sizeCluster[clusters[pId][0]]);
		//DEBUG
		if(pId == -1)
		{
			printf("center: [%f, %f, %f]\n", center[0], center[1], center[2]);
			printf("To: [%f, %f, %f]\n", initTo[clusters[pId][1]][0], initTo[clusters[pId][1]][1], initTo[clusters[pId][1]][2]);
		}
		directDis[pId] = initTo[clusters[pId][1]] - center;
	}


	template<typename Coord, typename Pair2, typename Real, template<typename Real> typename Quat>
	__global__ void CM_GetVirtualCoord(
		DArray<Pair2> clusters,
		DArray<Coord> to,
		DArray<Coord> virTo,
		DArray<Quat<Real>> initT,
		DArray<Quat<Real>> initR,
		DArray<Real> initS)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if (pId >= clusters.size()) return;
		int p = clusters[pId][1];
		int j = clusters[pId][0];

		Coord old_p = to[p];
		Quat<Real> tmp_p(old_p[0], old_p[1], old_p[2], 0);
		
		tmp_p = (1. / initS[j]) * initR[j].conjugate() * (tmp_p - initT[j]) * initR[j];
		Coord new_p = Coord(tmp_p.x, tmp_p.y, tmp_p.z);

		virTo[pId] = new_p;
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
		int numCluster = id_joint;

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
		
		// tet, tri
		/*
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
		*/

		m_pointClusters.assign(p_pairs2);
		int numPoint = m_to->getAllPoints().size();
		int numPair = m_pointClusters.size();

		// Set Color
		
		{
			std::vector<Vec3f> v_color(default_color, default_color + 12);
			DArray<Vec3f> d_color;
			d_color.resize(v_color.size());
			d_color.assign(v_color);

			if (this->outColor()->isEmpty())
				this->outColor()->allocate();
			auto& out_color = this->outColor()->getData();
			
			out_color.resize(numPoint);
			cuExecute(numPoint,
				CM_UpdateColor1,
				out_color);
			cuSynchronize();

			printf("Num of Clusters:%d\n", numPair);
			cuExecute(numPair,
				CM_UpdateColor2,
				out_color,
				m_pointClusters,
				d_color);
			cuSynchronize();
		}
		
		
		
		// Set Rigid Shape
		/*
		{
			DArray<Coord> barycentre;
			barycentre.resize(numCluster);
			count.resize(numCluster);
			
			// =0
			barycentre.reset();
			count.reset();
			
			// 获取每个joint的重心
			cuExecute(numPair,
				CM_CountClusterCoord,
				m_to->getAllPoints(), 
				m_pointClusters, 
				count, 
				barycentre);
			cuSynchronize();

			m_initCoord.resize(numPair);
			// 求解每个点的r_i
			cuExecute(numPair,
				CM_CountPointCoord,
				m_to->getAllPoints(), 
				m_pointClusters, 
				count, 
				barycentre,
				m_initCoord);
			
			barycentre.clear();
		}
		*/
		
		// Set m_virtualCoord
		{
			m_virtualCoord.resize(numPair);
			cuExecute(numPair,
				CM_GetVirtualCoord,
				m_pointClusters,
				m_to->getAllPoints(),
				m_virtualCoord,
				m_initQuatT,
				m_initQuatR,
				m_initS);
			cuSynchronize();
		}

		count.clear();
		p_pairs2.clear();

	}

	DEFINE_CLASS(CapsuleToMixSet);
}