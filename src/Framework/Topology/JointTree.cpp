#pragma once
#include "JointTree.h"

#define cos_angle(angle) cos(double(angle * 3.1415926f / 180.0))
#define sin_angle(angle) sin(double(angle * 3.1415926f / 180.0))

namespace dyno
{
    IMPLEMENT_CLASS_1(JointTree, TDataType)
    
    template<typename TDataType>
    JointTree<TDataType>::JointTree()
    {
        id = -1;
        RotationActive = true;
        PreRotation = Coord(0);
        PreTranslation = Coord(0);
        PreScaling = Coord(1); 

        LclTranslation = Coord(0);
        LclRotation = Coord(0);
        LclScaling = Coord(1);
		
        tmp = Coord(0, 0, 0);

        GlT = Quat<Real>(0, 0, 0, 0);
        GlR = Quat<Real>(0, 0, 0, 0);
        GlS = Quat<Real>(0, 0, 0, 0);
    }
    
    template<typename TDataType>
    JointTree<TDataType>::~JointTree()
    {

    }


	template<typename TDataType>
    Mat4f JointTree<TDataType>::getTransform(Coord & T, Coord& R, Coord& S)
    {
		// TODO:改为四元数
        Mat4f translation = Mat4f(
            1, 0, 0, T[0],
            0, 1, 0, T[1],
            0, 0, 1, T[2],
            0, 0, 0, 1);
        
        // R[X,Y,Z] -> [Z,X,Y]轴

        double X = R[0];
        Mat4f rotation_x = Mat4f(
            1, 0, 0, 0,
            0, cos_angle(X), -sin_angle(X), 0,
            0, sin_angle(X), cos_angle(X), 0,
            0, 0, 0, 1);

		double Y = R[1];
        Mat4f rotation_y = Mat4f(
            cos_angle(Y), 0, sin_angle(Y), 0,
            0, 1, 0, 0,
            -sin_angle(Y), 0, cos_angle(Y), 0,
            0, 0, 0, 1);

		double Z = R[2];
        Mat4f rotation_z = Mat4f(
            cos_angle(Z), -sin_angle(Z), 0, 0,
            sin_angle(Z), cos_angle(Z), 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1);

        Mat4f scaling= Mat4f(
            S[0], 0, 0, 0,
            0, S[1], 0, 0,
            0, 0, S[2], 0,
            0, 0, 0, 1);

        return  translation * scaling * rotation_x * rotation_y * rotation_z;
    }

    // 遍历关节层次时，顺便更新
	template<typename TDataType>
    void JointTree<TDataType>::getGlobalTransform()
    {
        // 注意顺序
        this->GlobalTransform = getTransform(this->LclTranslation, this->LclRotation, this->LclScaling);
        if(this->parent != nullptr)
            this->GlobalTransform = this->parent->GlobalTransform * this->GlobalTransform;
        else
        {
            //DEBUG 简单平移
            this->GlobalTransform = getTransform(this->PreTranslation, this->PreRotation, this->PreScaling) * this->GlobalTransform;
            
            // T
            Coord tmp_t(0,0,0);
            for (int i = 0; i < 3; ++i)
                tmp_t[i] += dist(generator)  * 0.0005;
            // this->PreTranslation += tmp_t;
            
            // R
            this->PreRotation += tmp;

            tmp = Coord(0);

        }
            
    }

    template<typename TDataType>
    void JointTree<TDataType>::getQuat(Coord &T, Coord &R, Coord &S)
    {

        // Scaling
        Quat<Real> s(S[0], S[1], S[2], 0.f);
        // Translation
        Quat<Real> t(T[0], T[1], T[2], 0.f);
        // Rotate X
        Quat<Real> q_x(R[0], Vec3f(1.f, 0.f, 0.f));
        Quat<Real> q_y(R[1], Vec3f(0.f, 1.f, 0.f));
        Quat<Real> q_z(R[2], Vec3f(0.f, 0.f, 1.f));
        
        GlT = GlT + GlS * GlR * t * GlR.conjugate();
        GlS = GlS * s;
        GlR = GlR * (q_x * q_y * q_z);

        // TODO:
        // GlR.normalize(); 
    }

    // 遍历关节层次时，顺便更新
	template<typename TDataType>
    void JointTree<TDataType>::getGlobalTransform2()
    {
        if(this->parent != nullptr)
        {
            this->GlT = this->parent->GlT;
            this->GlS = this->parent->GlS;
            this->GlR = this->parent->GlR;
            getQuat(this->LclTranslation, this->LclRotation, this->LclScaling);
        }else
        {
            this->GlT = Quat<Real>(0.f, 0.f, 0.f, 0.f);
            this->GlS = Quat<Real>(1.f, 1.f, 1.f, 0.f);
            this->GlR = Quat<Real>(0.f, 0.f, 0.f, 1.f);
            getQuat(this->PreTranslation, this->PreRotation, this->PreScaling);
            getQuat(this->LclTranslation, this->LclRotation, this->LclScaling);
        }
    }

    template<typename TDataType>
    void JointTree<TDataType>::copyFrom(JointTree<TDataType> jointTree)
    {
        this->id = jointTree.id;
        this->PreRotation = jointTree.PreRotation;
        this->LclTranslation = jointTree.LclTranslation;
        this->LclRotation = jointTree.LclRotation;
        this->LclScaling = jointTree.LclScaling;
        this->GlCoord = jointTree.GlCoord;
        this->GlobalTransform = jointTree.GlobalTransform;
        this->RotationActive = jointTree.RotationActive;
        this->children.assign(jointTree.children.begin(), jointTree.children.end());
        this->parent = jointTree.parent;
    }

    template<typename TDataType>
    void JointTree<TDataType>::scale(Real s)
    {
		PreScaling *= s;
    }

    template<typename TDataType>
    void JointTree<TDataType>::translate(Coord t)
    {
        PreTranslation += t;
    }
    
#ifdef PRECISION_FLOAT
	template class JointTree<DataType3f>;
#else
	template class JointTree<DataType3d>;
#endif    
}