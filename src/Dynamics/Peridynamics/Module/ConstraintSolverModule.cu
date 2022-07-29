#include "ConstraintSolverModule.h"
#include "Node.h"
#include "Matrix/MatrixFunc.h"
#include "ParticleSystem/Kernel.h"
#define CSM_BIAS 0
#define CSM_EPS 1e-6
#define CSM_BOUND 4
#define CSM_BIASCOORD Coord(-1e-5, -1e-5, -1e-5)
#define iszero(x, num) if (fabs(x) < CSM_EPS) printf("div zero in line %d\n", num)
#define isover(x, num) if (fabs(x) > CSM_BOUND) printf("over bounding in line %d\n", num)
#define Fline(num)  printf("overflow in line %d\n", num)

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
		// Jacobi Method
		// x[3N+C] = a*(b-c)

		// DEBUG
		auto& ef = this->inVelocity()->getData();
		ViewGPUData<TDataType> p_view[3];
		p_view[0].resize(mN);
		p_view[0].viewIm(ef, 0);
		p_view[1].resize(mN);
		p_view[1].viewIm(ef, 1);
		p_view[2].resize(mN);
		p_view[2].viewIm(ef, 2);

		this->computeInitValue();// [3N + C] :  a ,b

		int itor = 0;
		uint maxIterNum = this->varIterationNumber()->getData();
		while (itor < maxIterNum) {
			this->solveConstraint();// c=\sum a x
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

		Wa[p0] = 1. / A[pId]; iszero(A[pId], 76);
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
		// DArray<Real> ConstJ,
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
		// ConstJ[pId] = 1. / sqrtf(Dis[pId]); iszero(Dis[pId], 101);
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
		Wb[id] = Phi[pId] - Dis[pId] * sqrtf(Dis[pId]);
		// iszero(Dis[pId], 123);
	}
	
	template <typename Coord, typename Real, typename Pair2>
	__global__ void CSM_CalAllC(
		DArray<Pair2> con,
		DArray<Coord> X,
		DArray<Coord> oldX,
		DArray<Real> Wc,
		// DArray<Real> ConstJ,
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

		if (fabs(Dis[pId]) < CSM_EPS ) return;	// 0 * 1 / 0 = 0 当为原长时，特判为不做贡献（近似方法，本应该去考虑矩阵求解的细节）

		// Lambda -> Y
		Real sd = 1 / sqrtf(Dis[pId]); iszero(Dis[pId], 153);
		Real t = Lambda[pId] * sd;
		Coord t1 = - (oldX[y1] - oldX[y2]) * t;
		Coord t2 = - (oldX[y2] - oldX[y1]) * t;
		atomicAdd(&Wc[y1 * 3 + 0], t1[0]);
		atomicAdd(&Wc[y1 * 3 + 1], t1[1]);
		atomicAdd(&Wc[y1 * 3 + 2], t1[2]);

		atomicAdd(&Wc[y2 * 3 + 0], t2[0]);
		atomicAdd(&Wc[y2 * 3 + 1], t2[1]);
		atomicAdd(&Wc[y2 * 3 + 2], t2[2]);

		// Y -> Lambda
		Real l = sd * (X[y1].dot(oldX[y1] - oldX[y2]) + X[y2].dot(oldX[y2] - oldX[y1]));
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
			isover(t0, pId); isover(t1, pId); isover(t2, pId);
			if (pId == 2 || pId == 3)  //7
			{
				// printf("X[%d]: (%.3f %.3f %.3f)\n", pId, t0, t1, t2);
			}
		}else
		{
			int p = pId + N * 2;
			int id = pId - N;
			Lambda[id] = Wa[p] * (Wb[p] - Wc[p]);
		}
	}

	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::solveConstraint()
	{
		// for C : c  (y <-> \lambda)
		cuExecute(mC,
			CSM_CalAllC,
			mConstraint,
			mX,
			mY,
			mWc,
			// mConstJ,
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

		InvK[pId] = 1. / k;	iszero(k, 270);
		con[pId] = Pair2(InitconRE[pId][0], pId + numEPos);
		Len[pId] = (X[con[pId][0]] - X[con[pId][1]]).norm() + CSM_BIAS;
	}
	template <typename Real>
	__global__ void CSM_SetUpDynamicA(
		DArray<Real> A,
		DArray<Real> Ady,
		Real dt)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= A.size()) return;		
		Ady[pId] = A[pId] / (dt * dt); iszero(dt, 282);
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
			// A[pId] = densityE / (dt * dt);	iszero(dt, 301);
			A[pId] = densityE;
			X[pId] = posE[pId] + velE[pId] * dt;
			if (pId > posE.size()) Fline(315);
			if (pId > velE.size()) Fline(316);
		}
		else
		{
			// A[pId] = densityR / (dt * dt); iszero(dt, 306);
			A[pId] = densityR;
			X[pId] = posR[pId - numE];
			if (pId - numE > posR.size()) Fline(323);
		}
		isover(X[pId][0], 327); isover(X[pId][1], 327); isover(X[pId][2], 327);
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
			if (pId > posE.size()) Fline(347);
			if (pId > velE.size()) Fline(347);			
		}
		else
		{
			F[pId] = Coord(0);
			X[pId] = posR[pId - numE];
			// X[pId] += CSM_BIASCOORD;
			V[pId] = Coord(0);
			if (pId - numE > posR.size()) Fline(357);
		}
		Y[pId] = X[pId];
		isover(X[pId][0],363); isover(X[pId][1], 363); isover(X[pId][2], 363);
	}

	template <typename Coord, typename Real>
	__global__ void CSM_UpdateParticleVel(
		DArray<Coord> pos,
		DArray<Coord> cur_pos,
		DArray<Coord> vel,
		Real dt,
		int numE)
	{
		int pId = threadIdx.x + (blockIdx.x * blockDim.x);
		if(pId >= pos.size()) return;
		if (pId < numE)
			vel[pId] = (cur_pos[pId] - pos[pId]) / dt;  
		iszero(dt, 353);
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
		{
			// pos[pId] -= CSM_BIASCOORD;
			posR[pId - numE] = pos[pId];
		}
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
		// DEBUG
		auto& vel = this->inVelocity()->getData();
		ViewGPUData<TDataType> p_view[3];
		p_view[0].resize(mN);
		p_view[2].resize(mN);
		p_view[1].resize(mN);
		p_view[0].viewIm(vel, 0);
		p_view[1].viewIm(vel, 1);
		p_view[2].viewIm(vel, 2);

		Real dt = this->inTimeStep()->getData();
		cuExecute(mN,
			CSM_UpdateParticleVel,
			mY,
			mX,
			this->inVelocity()->getData(),
			dt,
			// 0);
			numEPos);
		cuSynchronize();

		// DEBUG
		p_view[0].viewIm(vel, 0);
		p_view[1].viewIm(vel, 1);
		p_view[2].viewIm(vel, 2);

		// TODO: Rigid Capsule Vel
	}

	// 计算迭代变量的初值
	template<typename TDataType>
	void ConstraintSolverModule<TDataType>::computeInitValue()
	{
		Real dt = this->inTimeStep()->getData();

		cuExecute(mN,
			CSM_SetUpDynamicA,
			mA,
			mAdy,
			dt);
		cuSynchronize();

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
			dt);
		cuSynchronize();

		// Init
		cuExecute(mC,
			CSM_CalConstraint,
			mPhi,
			mDis,
			// mConstJ,
			mLen,
			mConstraint,
			mX);
		cuSynchronize();	


		cuExecute(mN,
			CSM_CalParticleAB,
			mWa,
			mWb,
			mAdy,
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
	}

	template<typename TDataType>
	bool ConstraintSolverModule<TDataType>::initializeImpl()
	{
		Real dt = this->inTimeStep()->getData();	// dt = 0.0

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
		mAdy.resize(mN);

		cuExecute(mN,
			CSM_SetUpParticle,
			mA,
			mX,
			this->inPosition()->getData(),
			this->inVelocity()->getData(),
			this->inRigidPosition()->getData(),
			numEPos,
			numRPos,
			dt,
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