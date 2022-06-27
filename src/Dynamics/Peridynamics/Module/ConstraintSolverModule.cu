#include "ConstraintSolverModule.h"
#include "Node.h"
#include "Matrix/MatrixFunc.h"
#include "ParticleSystem/Kernel.h"

namespace dyno
{
	IMPLEMENT_TCLASS(ConstraintSolverModule, TDataType)
	
	template<typename TDataType>
	ConstraintSolverModule<TDataType>::ConstraintSolverModule()
		: ConstraintModule()
	{
		// const
	}


	template<typename TDataType>
	ConstraintSolverModule<TDataType>::~ConstraintSolverModule()
	{
        // clear()
		mY.clear();
		mX.clear();
		mF.clear();
		mV.clear();

		mInvK.clear();
		mA.clear();

		//Jacobi Method
		mWa.clear();
		mWb.clear();
		mWc.clear();
		mLambda.clear();

		mConstraint.clear();
		mPhi.clear();
		mDis.clear();
		mConstJ.clear();
		mLen.clear(); // = 0
	}


	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::constrain()
	{
		this->computeInitValue();

		int itor = 0;
		uint maxIterNum = this->varIterationNumber()->getData();
		while (itor < maxIterNum) {
			this->solveConstraint();
			itor++;
		}

		this->updateVelocity();
		this->updatePosition();
	}

	template <typename Real, typename Coord>
	__global__ void CSM_CalParticleAB(
		DArray<Real> Wa,
		DArray<Real> Wb,
		DArray<Real> A,
		DArray<Coord> F,
		DArray<Coord> X)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= F.size()) return;
		int p0 = pId * 3;
		int p1 = p0 + 1;
		int p2 = p1 + 1;

		Wa[p0] = 1. / A[pId];
		Wa[p1] = Wa[p0];
		Wa[p2] = Wa[p0];
		Coord t = F[pId] + A[pId] * X[pId];
		Wb[p0] = t[0];
		Wb[p1] = t[1];
		Wb[p2] = t[2];
	}
	
	template <typename Coord, typename Real, typename Pair2>
	__global__ void CSM_CalConstraint(
		DArray<Real> Phi,
		DArray<Real> Dis,
		DArray<Real> ConstJ,
		DArray<Real> Len,
		DArray<Pair2> con,
		DArray<Coord> X)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= con.size()) return;
		int y1 = con[pId][0];
		int y2 = con[pId][1];
		Dis[pId] = (X[y1] - X[y2]).norm();
		Phi[pId] = Dis[pId] - Len[pId];
		ConstJ[pId] = 1. / sqrtf(Dis[pId]);
	}

	template <typename Coord, typename Real, typename Pair2>
	__global__ void CSM_CalLambdaAB(
		DArray<Real> Wa,
		DArray<Real> Wb,
		DArray<Real> invK,
		DArray<Real> Phi,
		DArray<Real> Dis,
		DArray<Coord> X,
		DArray<Pair2> con,
		int N,
		int C)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= C) return;
		int id = pId + N * 3;
		int y1 = con[pId][0];
		int y2 = con[pId][1];

		Wa[id] = 1. / invK[pId];
		Wb[id] = Phi[pId] - (X[y1].normSquared() + X[y2].normSquared()) / sqrtf(Dis[pId]);
	}
	
	template <typename Coord, typename Real, typename Pair2>
	__global__ void CSM_CalAllC(
		DArray<Pair2> con,
		DArray<Coord> X,
		DArray<Real> Wc,
		DArray<Real> ConstJ,
		DArray<Real> Dis,
		DArray<Real> Lambda,
		int N,
		int C)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= C) return;
		int id = pId + N * 3;
		int y1 = con[pId][0];
		int y2 = con[pId][1];

		// Lambda -> Y
		Coord t1 = - X[y1] * ConstJ[pId] * Lambda[pId];
		Coord t2 = - X[y2] * ConstJ[pId] * Lambda[pId];
		atomicAdd(&Wc[y1 * 3 + 0], t1[0]);
		atomicAdd(&Wc[y1 * 3 + 1], t1[1]);
		atomicAdd(&Wc[y1 * 3 + 2], t1[2]);

		atomicAdd(&Wc[y2 * 3 + 0], t2[0]);
		atomicAdd(&Wc[y2 * 3 + 1], t2[1]);
		atomicAdd(&Wc[y2 * 3 + 2], t2[2]);

		// Y -> Lambda
		Real l = - (X[y1].normSquared() + X[y2].normSquared()) * ConstJ[pId];
		atomicAdd(&Wc[id], l);
	}

		
	template <typename Coord, typename Real>
	__global__ void CSM_CalJacobi(
		DArray<Coord> X,
		DArray<Real> Lambda,
		DArray<Real> Wa,
		DArray<Real> Wb,
		DArray<Real> Wc,
		int N,
		int C)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= N + C) return;

		if(pId < N)
		{
			int p0 = pId * 3;
			int p1 = p0 + 1;
			int p2 = p1 + 1;
				
			Real t0 = Wa[p0] * (Wb[p0] - Wc[p0]);
			Real t1 = Wa[p1] * (Wb[p1] - Wc[p1]);
			Real t2 = Wa[p2] * (Wb[p2] - Wc[p2]);
			X[pId] = Coord(t0, t1, t2);
		}else
		{
			int p = pId + N * 3;
			int id = pId - N;
			Lambda[id] = Wa[p] * (Wb[p] - Wc[p]);
		}
	}

	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::solveConstraint()
	{
		// Init
		cuExecute(mC,
			CSM_CalConstraint,
			mPhi,
			mDis,
			mConstJ,
			mLen,
			mConstraint,
			mX);
		cuSynchronize();

		// Jacobi Method
		// x[3N+C] = a*(b-c)
		// [3N + C] :  a ,b

		cuExecute(mN,
			CSM_CalParticleAB,
			mWa,
			mWb,
			mA,
			mF,
			mX);
		cuSynchronize();

		cuExecute(mC,
			CSM_CalLambdaAB,
			mWa,
			mWb,
			mInvK,
			mPhi,
			mDis,
			mX,
			mConstraint,
			mN,
			mC);
		cuSynchronize();
			
		// for C : c  (y <-> \lambda)
		cuExecute(mC,
			CSM_CalAllC,
			mConstraint,
			mX,
			mWc,
			mConstJ,
			mDis,
			mLambda,
			mN,
			mC);
		cuSynchronize();

		// jacobi update pos
		cuExecute(mN + mC,
			CSM_CalJacobi,
			mX,
			mLambda,
			mWa,
			mWb,
			mWc,
			mN,
			mC);
		cuSynchronize();
	}
	template <typename Coord, typename Real, typename Pair2>
	__global__ void CSM_SetUpConstraint(
		DArray<Real> InvK,
		DArray<Pair2> InitconRE,
		DArray<Pair2> con,
		DArray<Real> Len,
		DArray<Coord> X,
		int numEPos,
		Real k)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= InvK.size()) return;

		InvK[pId] = 1. / k;
		con[pId] = Pair2(InitconRE[pId][0], pId + numEPos);
		Len[pId] = (X[con[pId][0]] - X[con[pId][1]]).norm();
	}

	template <typename Real>
	__global__ void CSM_SetUpLambda(
		DArray<Real> Lambda)
	{	
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= Lambda.size()) return;
		Lambda[pId] = 0;
	}

	template <typename Coord, typename Real>
	__global__ void CSM_SetUpParticle(
		DArray<Real> A,
		DArray<Coord> X,
		DArray<Coord> posE,
		DArray<Coord> velE,
		DArray<Coord> posR,
		int numE,
		int numR,
		Real dt,
		Real densityE,
		Real densityR)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= numE + numR) return;	
		if(pId < numE)	
		{
			A[pId] = densityE / (dt * dt);
			X[pId] = posE[pId] + velE[pId] * dt;
		}
		else
		{
			A[pId] = densityR / (dt * dt);
			X[pId] = posR[pId - numE];
		}
	}

	template <typename Coord, typename Real>
	__global__ void CSM_UpdateParticle(
		DArray<Coord> posE,
		DArray<Coord> velE,
		DArray<Coord> forceE,
		DArray<Coord> posR,
		DArray<Coord> Y,
		DArray<Coord> X,
		DArray<Coord> F,
		DArray<Coord> V,
		int numE,
		int numR,
		Real dt)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= numE + numR) return;

		if(pId < numE)
		{
			F[pId] = forceE[pId];
			X[pId] = posE[pId] + velE[pId] * dt;
			V[pId] = velE[pId];
		}
		else
		{
			F[pId] = Coord(0);
			X[pId] = posR[pId - numE];
			V[pId] = Coord(0);
		}
		Y[pId] = X[pId];
	}

	template <typename Coord, typename Real>
	__global__ void CSM_UpdateParticleVel(
		DArray<Coord> pos,
		DArray<Coord> cur_pos,
		DArray<Coord> vel,
		Real dt)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= pos.size()) return;

		vel[pId] = (cur_pos[pId] - pos[pId]) / dt;
		pos[pId] = cur_pos[pId];
	}

	template <typename Coord>
	__global__ void CSM_UpdateParticlePos(
		DArray<Coord> pos,
		DArray<Coord> posE,
		DArray<Coord> posR,
		int numE)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= pos.size()) return;

		if (pId < numE)
			posE[pId] = pos[pId];
		else
			posR[pId - numE] = pos[pId];
	}

	// 更新弹性粒子和刚体顶点坐标
	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::updatePosition()
	{
		cuExecute(mN,
			CSM_UpdateParticlePos,
			mY,
			this->inPosition()->getData(),
			this->inRigidPosition()->getData(),
			numEPos);
		cuSynchronize();
	}


	// 更新粒子速度
	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::updateVelocity()
	{
		cuExecute(mN,
			CSM_UpdateParticleVel,
			mY,
			mX,
			this->inVelocity()->getData(),
			this->inTimeStep()->getData());
		cuSynchronize();

		// TODO: Rigid Capsule Vel
	}

	// 计算迭代变量的初值
	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::computeInitValue()
	{
		cuExecute(mC,
			CSM_SetUpLambda,
			mLambda);
		cuSynchronize();

		cuExecute(mN,
			CSM_UpdateParticle,
			this->inPosition()->getData(),
			this->inVelocity()->getData(),
			this->inElasticForce()->getData(),
			this->inRigidPosition()->getData(),
			mY,
			mX,
			mF,
			mV,
			numEPos,
			numRPos,
			this->inTimeStep()->getData());
		cuSynchronize();
	}

	template<typename TDataType>
	bool ConstraintSolverModule<TDataType>::initializeImpl()
	{
		
		numEPos = this->inPosition()->getElementCount();
		numRPos = this->inRigidPosition()->getElementCount();
		
		mN = numEPos + numRPos;

		int numERCon = this->inConstraintRE()->getElementCount();
		// int numRRCon = this->inConstraintRR()->getElementCount();
		int numRRCon= 0;

		mC = numERCon + numRRCon;

		mY.resize(mN);
		mX.resize(mN);
		mF.resize(mN);
		mV.resize(mN);
		mA.resize(mN);

		cuExecute(mN,
			CSM_SetUpParticle,
			mA,
			mX,
			this->inPosition()->getData(),
			this->inVelocity()->getData(),
			this->inRigidPosition()->getData(),
			numEPos,
			numRPos,
			this->inTimeStep()->getData(),
			Real(1.0),
			Real(1.0));
		cuSynchronize();

		mInvK.resize(mC);
		mConstraint.resize(mC);
		mLen.resize(mC);
		mDis.resize(mC);
		mPhi.resize(mC);
		mConstJ.resize(mC);

		cuExecute(mC,
			CSM_SetUpConstraint,
			mInvK,
			this->inConstraintRE()->getData(),
			mConstraint,
			mLen,
			mX,
			numEPos,
			Real(1.0));
		cuSynchronize();

		mWa.resize(3 * mN + mC);
		mWb.resize(3 * mN + mC);
		mWc.resize(3 * mN + mC);
		mLambda.resize(mC);
		
		printf("Elastic Particle: %d\n", numEPos);
		printf("Ridig Particle: %d\n", numRPos);
		printf("Constraint: %d\n", mC);

		return true;
	}

	DEFINE_CLASS(ConstraintSolverModule);
}