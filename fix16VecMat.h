#ifndef _FIT16_VEC_MAT_H_
#define _FIT16_VEC_MAT_H_

#ifdef __cplusplus
extern "C"
{
#endif

#include "fix16.h"
#include "fix_config.h"

typedef struct {
    fix16_t x;
    fix16_t y;
} v2d;
typedef struct {
	fix16_t x;
	fix16_t y;
	fix16_t z;
} v3d;
typedef struct {
	fix16_t x;
	fix16_t y;
	fix16_t z;
    fix16_t w;
} v4d;

typedef struct {
    fix16_t a; // Real part
    fix16_t b; // i
    fix16_t c; // j
    fix16_t d; // k
} qf16;

typedef struct{
    fix16_t m00;
    fix16_t m01;
    fix16_t m02;
    fix16_t m10;
    fix16_t m11;
    fix16_t m12;
    fix16_t m20;
    fix16_t m21;
    fix16_t m22;
}mat3_t;

typedef struct{
    fix16_t m00;
    fix16_t m01;
    fix16_t m02;
    fix16_t m03;
    fix16_t m10;
    fix16_t m11;
    fix16_t m12;
    fix16_t m13;
    fix16_t m20;
    fix16_t m21;
    fix16_t m22;
    fix16_t m23;
    fix16_t m30;
    fix16_t m31;
    fix16_t m32;
    fix16_t m33;
}mat4_t;
// Basic arithmetic
extern void v2d_add(v2d *dest, const v2d *a, const v2d *b);
extern void v2d_sub(v2d *dest, const v2d *a, const v2d *b);
extern void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b);
extern void v2d_div_s(v2d *dest, const v2d *a, fix16_t b);
extern void v3d_add(v3d *dest, const v3d *a, const v3d *b);
extern void v3d_sub(v3d *dest, const v3d *a, const v3d *b);
extern void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b);
extern void v3d_div_s(v3d *dest, const v3d *a, fix16_t b);
// Norm
extern fix16_t v2d_norm(const v2d *a);
extern void v2d_normalize(v2d *dest, const v2d *a);
extern fix16_t v3d_norm(const v3d *a);
extern void v3d_normalize(v3d *dest, const v3d *a);
// Dot product
extern fix16_t v2d_dot(const v2d *a, const v2d *b);
extern fix16_t v3d_dot(const v3d *a, const v3d *b);
// Rotation (positive direction = counter-clockwise, angle in radians)
extern void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle);
// Cross product
extern void v3d_cross(v3d *dest, const v3d *a, const v3d *b);
// Conjugate of quaternion
extern void qf16_conj(qf16 *dest, const qf16 *q);
// Multiply two quaternions, dest = q * r.
extern void qf16_mul(qf16 *dest, const qf16 *q, const qf16 *r);
// Add two quaternions, dest = q + r
extern void qf16_add(qf16 *dest, const qf16 *q, const qf16 *r);
// Multiply quaternion by scalar
extern void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s);
// Divide quaternion by scalar
extern void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s);
// Dot product of two quaternions
extern fix16_t qf16_dot(const qf16 *q, const qf16 *r);
// Quaternion norm
extern fix16_t qf16_norm(const qf16 *q);
// Normalize quaternion
extern void qf16_normalize(qf16 *dest, const qf16 *q);
// Quaternion power (exponentation)
extern void qf16_pow(qf16 *dest, const qf16 *q, fix16_t power);
// Weighted average of two quaternions
// Think of it as q = w * q1 + (1 - w) * q2, but the internal algorithm considers attitudes.
extern void qf16_avg(qf16 *dest, const qf16 *q1, const qf16 *q2, fix16_t weight);
// Unit quaternion from axis and angle.
// Axis should have unit length and angle in radians.
extern void qf16_from_axis_angle(qf16 *dest, const v3d *axis, fix16_t angle);
// Unit quaternion to rotation matrix
extern void qf16_to_matrix(mat3_t *pMat, const qf16 *q);
// Rotate vector using quaternion
extern void qf16_rotate(v3d *dest, const qf16 *q, const v3d *v);




extern void InitUnitMat4(mat4_t* pMat);
extern void Mat4XMat4(mat4_t *pMat1,mat4_t *pMat2, mat4_t *pMat3);

extern void CreateO2WMat(mat3_t* pRMat, v3d* pTranslation, v3d* pScale, mat4_t* pResult);
extern void GenerateW2CMat(mat4_t* pResult, mat3_t* pSource, vect3_t* pVect) ;
extern void MakeClipMatrix_perspective( fix16_t near_plane, fix16_t far_plane,
    fix16_t focalLength, fix16_t aspectRatio, mat4_t* pMat) ;
#ifdef __cplusplus
}
#endif

#endif
