#include "fit16VecMat.h"
static fix16_t scale_value(fix16_t value, int_fast8_t scale);
#ifdef __GNUC__
// Count leading zeros, using processor-specific instruction if available.
#define clz(x) __builtin_clzl(x)
#else
static uint8_t clz(uint32_t x)
{
  uint8_t result = 0;
  if (x == 0) return 32;
  while (!(x & 0xF0000000)) { result += 4; x <<= 4; }
  while (!(x & 0x80000000)) { result += 1; x <<= 1; }
  return result;
}
#endif

static fix16_t scale_value(fix16_t value, int_fast8_t scale)
{
    if (scale > 0)
    {
        fix16_t temp = value << scale;
        if (temp >> scale != value)
            return fix16_overflow;
        else
            return temp;
    }
    else if (scale < 0)
    {
        return value >> -scale;
    }
    else
    {
        return value;
    }
}

// Basic arithmetic
void v2d_add(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_add(a->x, b->x);
    dest->y = fix16_add(a->y, b->y);
}

void v2d_sub(v2d *dest, const v2d *a, const v2d *b)
{
    dest->x = fix16_sub(a->x, b->x);
    dest->y = fix16_sub(a->y, b->y);
}

void v2d_mul_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_mul(a->x, b);
    dest->y = fix16_mul(a->y, b);
}

void v2d_div_s(v2d *dest, const v2d *a, fix16_t b)
{
    dest->x = fix16_div(a->x, b);
    dest->y = fix16_div(a->y, b);
}

void v3d_add(v3d *dest, const v3d *a, const v3d *b)
{
        dest->x = fix16_add(a->x, b->x);
        dest->y = fix16_add(a->y, b->y);
        dest->z = fix16_add(a->z, b->z);
}

void v3d_sub(v3d *dest, const v3d *a, const v3d *b)
{
        dest->x = fix16_sub(a->x, b->x);
        dest->y = fix16_sub(a->y, b->y);
        dest->z = fix16_sub(a->z, b->z);
}

void v3d_mul_s(v3d *dest, const v3d *a, fix16_t b)
{
        dest->x = fix16_mul(a->x, b);
        dest->y = fix16_mul(a->y, b);
        dest->z = fix16_mul(a->z, b);
}

void v3d_div_s(v3d *dest, const v3d *a, fix16_t b)
{
        dest->x = fix16_div(a->x, b);
        dest->y = fix16_div(a->y, b);
        dest->z = fix16_div(a->z, b);
}

// Norm
fix16_t v2d_norm(const v2d *a){

    fix16_t max = 0;
    
    // Calculate inclusive OR of all values to find out the maximum.
    max |= fix16_abs(a->x);
    max |= fix16_abs(a->y);
    // To avoid overflows, the values before squaring can be max 128.0,
    // i.e. v & 0xFF800000 must be 0. Also, to avoid overflow in sum,
    // we need additional log2(n) bits of space.
    int_fast8_t scale = clz(max) - 10;
    fix16_t val = scale_value(a->x, scale);
    fix16_t sum = fix16_mul(val, val);

    val = scale_value(a->y, scale);
    sum += fix16_mul(val, val);

#ifndef FIXMATH_NO_OVERFLOW
    if (sum == fix16_overflow)
        return sum;
#endif
    fix16_t result = fix16_sqrt(sum);
    return scale_value(result, -scale);
}

void v2d_normalize(v2d *dest, const v2d *a)
{
    v2d_div_s(dest, a, v2d_norm(a));
}

// Norm
fix16_t v3d_norm(const v3d *a)
{
    
    fix16_t max = 0;
    
    // Calculate inclusive OR of all values to find out the maximum.
    max |= fix16_abs(a->x);
    max |= fix16_abs(a->y);
    max |= fix16_abs(a->z);
    // To avoid overflows, the values before squaring can be max 128.0,
    // i.e. v & 0xFF800000 must be 0. Also, to avoid overflow in sum,
    // we need additional log2(n) bits of space.
    int_fast8_t scale = clz(max) - 10;
    fix16_t val = scale_value(a->x, scale);
    fix16_t sum = fix16_mul(val, val);
    val = scale_value(a->y, scale);
    sum += fix16_mul(val, val);
    val = scale_value(a->z, scale);
    sum += fix16_mul(val, val);
#ifndef FIXMATH_NO_OVERFLOW
    if (sum == fix16_overflow)
        return sum;
#endif
    fix16_t result = fix16_sqrt(sum);
    return scale_value(result, -scale);
}

void v3d_normalize(v3d *dest, const v3d *a)
{
    v3d_div_s(dest, a, v3d_norm(a));
}

// Dot product
fix16_t v2d_dot(const v2d *a, const v2d *b)
{
    return fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y));
}

fix16_t v3d_dot(const v3d *a, const v3d *b)
{
    return fix16_add(fix16_add(fix16_mul(a->x, b->x), fix16_mul(a->y, b->y)),fix16_mul(a->z, b->z)) ;
}

// Rotation (positive direction = counter-clockwise, angle in radians)
void v2d_rotate(v2d *dest, const v2d *a, fix16_t angle)
{
    fix16_t c = fix16_cos(angle);
    fix16_t s = fix16_sin(angle);
    
    dest->x = fix16_add(fix16_mul(c, a->x), fix16_mul(-s, a->y));
    dest->y = fix16_add(fix16_mul(s, a->x), fix16_mul(c, a->y));
}

// Cross product
void v3d_cross(v3d *dest, const v3d *a, const v3d *b)
{
    fix16_t x,y,z;
    
    x = fix16_sub(fix16_mul(a->y, b->z), fix16_mul(a->z, b->y));
    y = fix16_sub(fix16_mul(a->z, b->x), fix16_mul(a->x, b->z));
    z = fix16_sub(fix16_mul(a->x, b->y), fix16_mul(a->y, b->x));
    dest->x = x;
    dest->y = y;
    dest->z = z;
}
//quaternion--------------------------------------------------------------------------------------------------------
static inline void qf16_from_v3d(qf16 *q, const v3d *v, fix16_t a);
static inline void qf16_to_v3d(v3d *v, const qf16 *q);
static void fa16Copy(qf16 *pSource,qf16 *pTarget);
static void fa16_unalias(qf16 *dest, qf16 **a, qf16 **b, qf16 *tmp);


static void fa16Copy(qf16 *pSource,qf16 *pTarget){
    pTarget->a = pSource->a;
    pTarget->b = pSource->b;
    pTarget->c = pSource->c;
    pTarget->d = pSource->d;
}

static void fa16_unalias(qf16 *dest, qf16 **a, qf16 **b, qf16 *tmp){
    if (dest == *a)
    {
        fa16Copy((*a),tmp);
        *a = tmp;
        
        if (dest == *b)
            *b = tmp;
    }
    else if (dest == *b)
    {
        fa16Copy((*b),tmp);
        *b = tmp;
    }
}

static inline void qf16_from_v3d(qf16 *q, const v3d *v, fix16_t a)
{
    q->a = a;
    q->b = v->x;
    q->c = v->y;
    q->d = v->z;
}

static inline void qf16_to_v3d(v3d *v, const qf16 *q)
{
    v->x = q->b;
    v->y = q->c;
    v->z = q->d;
}


// Conjugate of quaternion
void qf16_conj(qf16 *dest, const qf16 *q)
{
    dest->a = q->a;
    dest->b = - q->b;
    dest->c = - q->c;
    dest->d = - q->d;
}

// Multiply two quaternions, dest = a * b.
void qf16_mul(qf16 *dest, const qf16 *q, const qf16 *r)
{
    fix16_t a,b,c,d;
    
    a = fix16_mul(q->a, r->a) - fix16_mul(q->b, r->b) - fix16_mul(q->c, r->c) - fix16_mul(q->d, r->d);
    b = fix16_mul(q->a, r->b) + fix16_mul(q->b, r->a) + fix16_mul(q->c, r->d) - fix16_mul(q->d, r->c);
    c = fix16_mul(q->a, r->c) - fix16_mul(q->b, r->d) + fix16_mul(q->c, r->a) + fix16_mul(q->d, r->b);
    d = fix16_mul(q->a, r->d) + fix16_mul(q->b, r->c) - fix16_mul(q->c, r->b) + fix16_mul(q->d, r->a);
    dest->a = a;
    dest->b = b;
    dest->c = c;
    dest->d = d;
}

void qf16_add(qf16 *dest, const qf16 *q, const qf16 *r)
{
    dest->a = q->a + r->a;
    dest->b = q->b + r->b;
    dest->c = q->c + r->c;
    dest->d = q->d + r->d;
}

// Multiply quaternion by scalar
void qf16_mul_s(qf16 *dest, const qf16 *q, fix16_t s)
{
    dest->a = fix16_mul(q->a, s);
    dest->b = fix16_mul(q->b, s);
    dest->c = fix16_mul(q->c, s);
    dest->d = fix16_mul(q->d, s);
}

// Divide quaternion by scalar
void qf16_div_s(qf16 *dest, const qf16 *q, fix16_t s)
{
    dest->a = fix16_div(q->a, s);
    dest->b = fix16_div(q->b, s);
    dest->c = fix16_div(q->c, s);
    dest->d = fix16_div(q->d, s);
}

fix16_t qf16_dot(const qf16 *q, const qf16 *r)
{
    return fix16_mul(q->a,r->a)+fix16_mul(q->b,r->b)+fix16_mul(q->c,r->c)+fix16_mul(q->d,r->d); 
}

// Quaternion norm
fix16_t qf16_norm(const qf16 *q)
{    
    fix16_t max = 0;
    // Calculate inclusive OR of all values to find out the maximum.
    max |= fix16_abs(a->a);
    max |= fix16_abs(a->b);
    max |= fix16_abs(a->c);
    max |= fix16_abs(a->d);
    // To avoid overflows, the values before squaring can be max 128.0,
    // i.e. v & 0xFF800000 must be 0. Also, to avoid overflow in sum,
    // we need additional log2(n) bits of space.
    int_fast8_t scale = clz(max) - 10;
    fix16_t val = scale_value(a->a, scale);
    fix16_t sum = fix16_mul(val, val);
    val = scale_value(a->b, scale);
    sum += fix16_mul(val, val);
    val = scale_value(a->c, scale);
    sum += fix16_mul(val, val);
    val = scale_value(a->d, scale);
    sum += fix16_mul(val, val);
#ifndef FIXMATH_NO_OVERFLOW
    if (sum == fix16_overflow)
        return sum;
#endif
    fix16_t result = fix16_sqrt(sum);
    return scale_value(result, -scale);
}

// Normalize quaternion
void qf16_normalize(qf16 *dest, const qf16 *q)
{
    qf16_div_s(dest, q, qf16_norm(q));
}

// Quaternion power
void qf16_pow(qf16 *dest, const qf16 *q, fix16_t power)
{
    fix16_t old_half_angle = fix16_acos(q->a);
    fix16_t new_half_angle = fix16_mul(old_half_angle, power);
    fix16_t multiplier = 0;
    
    if (old_half_angle > 10) // Guard against almost-zero divider
    {
        multiplier = fix16_div(fix16_sin(new_half_angle),
                               fix16_sin(old_half_angle));
    }
    
    dest->a = fix16_cos(new_half_angle);
    dest->b = fix16_mul(q->b, multiplier);
    dest->c = fix16_mul(q->c, multiplier);
    dest->d = fix16_mul(q->d, multiplier);
}

// Weighted average
// See http://www.acsu.buffalo.edu/~johnc/ave_sfm07.pdf
void qf16_avg(qf16 *dest, const qf16 *q1, const qf16 *q2, fix16_t weight)
{
    // z = sqrt((w1 - w2)^2 + 4 w1 w2 (q1' q2)^2
    // <=>
    // z = sqrt((2 w1 - 1)^2 + 4 w1 (1 - w1) (q1' q2)^2)
    fix16_t dot = qf16_dot(q1, q2);
    fix16_t temp = 2 * weight - fix16_one;
    fix16_t z = fix16_mul(temp,temp)
            + fix16_mul(4 * weight, fix16_mul((fix16_one - weight), fix16_mul(dot,dot)));
    z = fix16_sqrt(z);
    
    // q = 2 * w1 * (q1' q2) q1 + (w2 - w1 + z) q2
    // <=>
    // q = 2 * w1 * (q1' q2) q1 + (1 - 2 * w1 + z) q2
    qf16 tmp1;
    qf16_mul_s(&tmp1, q1, fix16_mul(2 * weight, dot));
    
    qf16 tmp2;
    qf16_mul_s(&tmp2, q2, fix16_one - 2 * weight + z);
    
    qf16_add(dest, &tmp1, &tmp2);
    qf16_normalize(dest, dest);
}

void qf16_from_axis_angle(qf16 *dest, const v3d *axis, fix16_t angle){
    angle /= 2;
    fix16_t scale = fix16_sin(angle);
    
    dest->a = fix16_cos(angle);
    dest->b = fix16_mul(axis->x, scale);
    dest->c = fix16_mul(axis->y, scale);
    dest->d = fix16_mul(axis->z, scale);
}

// Unit quaternion to rotation matrix
void qf16_to_matrix(mat3_t *pMat, const qf16 *q){
    fix16_t x = pQuat->b; fix16_t y = pQuat->c; fix16_t z = pQuat->d; fix16_t w = pQuat->a;
    pMat->m00 = fix16_one - 2 * (fix16_mul(y , y) + fix16_mul(z , z));
    pMat->m01 = 2 * (fix16_mul(x , y) - fix16_mul(w , z));
    pMat->m02 = 2 * (fix16_mul(x , z) + fix16_mul(w , y));
    pMat->m10 = 2 * (fix16_mul(x , y) + fix16_mul(w , z));
    pMat->m11 = fix16_one - 2 * (fix16_mul(x , x) + fix16_mul(z , z));
    pMat->m12 = 2 * (fix16_mul(y , z) - fix16_mul(w , x));
    pMat->m20 = 2 * (fix16_mul(x , z) - fix16_mul(w , y));
    pMat->m21 = 2 * (fix16_mul(y , z) + fix16_mul(w , x));
    pMat->m22 = fix16_one - 2 * (fix16_mul(x , x) + fix16_mul(y , y));
}

void qf16_rotate(v3d *dest, const qf16 *q, const v3d *v)
{
    qf16 vector, q_conj;
    
    qf16_from_v3d(&vector, v, 0);
    qf16_conj(&q_conj, q);
    
    qf16_mul(&vector, q, &vector);
    qf16_mul(&vector, &vector, &q_conj);
    
    qf16_to_v3d(dest, &vector);
}





//matrix------------------------------------------------------------------------------------------------------------

void InitUnitMat4(mat4_t* pMat) {
#define M(x,y) (pMat)->m##x##y
    M(0, 0) = fix16_one; M(1, 0) = 0; M(2, 0) = 0; M(3, 0) = 0;
    M(0, 1) = 0; M(1, 1) = fix16_one; M(2, 1) = 0; M(3, 1) = 0;
    M(0, 2) = 0; M(1, 2) = 0; M(2, 2) = fix16_one; M(3, 2) = 0;
    M(0, 3) = 0; M(1, 3) = 0; M(2, 3) = 0; M(3, 3) = fix16_one;
#undef M
}

void Mat4XMat4(mat4_t *pMat1,mat4_t *pMat2, mat4_t *pMat3){
    fix16_t t0, t1, t2, t3;
    fix16_t s0, s1, s2, s3;
#define M(x,y) (pMat1)->m##x##y
#define N(x,y) (pMat2)->m##x##y
#define O(x,y) (pMat3)->m##x##y
    s0 = M(0, 0); s1 = M(1, 0); s2 = M(2, 0); s3 = M(3, 0);
    t0 = fix16_mul(s0 , N(0, 0)) + fix16_mul(s1 , N(0, 1)) + fix16_mul(s2 , N(0, 2)) + fix16_mul(s3 , N(0, 3));
    t1 = fix16_mul(s0 , N(1, 0)) + fix16_mul(s1 , N(1, 1)) + fix16_mul(s2 , N(1, 2)) + fix16_mul(s3 , N(1, 3));
    t2 = fix16_mul(s0 , N(2, 0)) + fix16_mul(s1 , N(2, 1)) + fix16_mul(s2 , N(2, 2)) + fix16_mul(s3 , N(2, 3));
    t3 = fix16_mul(s0 , N(3, 0)) + fix16_mul(s1 , N(3, 1)) + fix16_mul(s2 , N(3, 2)) + fix16_mul(s3 , N(3, 3));
    O(0, 0) = t0; O(1, 0) = t1; O(2, 0) = t2; O(3, 0) = t3;
    s0 = M(0, 1); s1 = M(1, 1); s2 = M(2, 1); s3 = M(3, 1);
    t0 = fix16_mul(s0 , N(0, 0)) + fix16_mul(s1 , N(0, 1)) + fix16_mul(s2 , N(0, 2)) + fix16_mul(s3 , N(0, 3));
    t1 = fix16_mul(s0 , N(1, 0)) + fix16_mul(s1 , N(1, 1)) + fix16_mul(s2 , N(1, 2)) + fix16_mul(s3 , N(1, 3));
    t2 = fix16_mul(s0 , N(2, 0)) + fix16_mul(s1 , N(2, 1)) + fix16_mul(s2 , N(2, 2)) + fix16_mul(s3 , N(2, 3));
    t3 = fix16_mul(s0 , N(3, 0)) + fix16_mul(s1 , N(3, 1)) + fix16_mul(s2 , N(3, 2)) + fix16_mul(s3 , N(3, 3));
    O(0, 1) = t0; O(1, 1) = t1; O(2, 1) = t2; O(3, 1) = t3;
    s0 = M(0, 2); s1 = M(1, 2); s2 = M(2, 2); s3 = M(3, 2);
    t0 = fix16_mul(s0 , N(0, 0)) + fix16_mul(s1 , N(0, 1)) + fix16_mul(s2 , N(0, 2)) + fix16_mul(s3 , N(0, 3));
    t1 = fix16_mul(s0 , N(1, 0)) + fix16_mul(s1 , N(1, 1)) + fix16_mul(s2 , N(1, 2)) + fix16_mul(s3 , N(1, 3));
    t2 = fix16_mul(s0 , N(2, 0)) + fix16_mul(s1 , N(2, 1)) + fix16_mul(s2 , N(2, 2)) + fix16_mul(s3 , N(2, 3));
    t3 = fix16_mul(s0 , N(3, 0)) + fix16_mul(s1 , N(3, 1)) + fix16_mul(s2 , N(3, 2)) + fix16_mul(s3 , N(3, 3));
    O(0, 2) = t0; O(1, 2) = t1; O(2, 2) = t2; O(3, 2) = t3;
    s0 = M(0, 3); s1 = M(1, 3); s2 = M(2, 3); s3 = M(3, 3);
    t0 = fix16_mul(s0 , N(0, 0)) + fix16_mul(s1 , N(0, 1)) + fix16_mul(s2 , N(0, 2)) + fix16_mul(s3 , N(0, 3));
    t1 = fix16_mul(s0 , N(1, 0)) + fix16_mul(s1 , N(1, 1)) + fix16_mul(s2 , N(1, 2)) + fix16_mul(s3 , N(1, 3));
    t2 = fix16_mul(s0 , N(2, 0)) + fix16_mul(s1 , N(2, 1)) + fix16_mul(s2 , N(2, 2)) + fix16_mul(s3 , N(2, 3));
    t3 = fix16_mul(s0 , N(3, 0)) + fix16_mul(s1 , N(3, 1)) + fix16_mul(s2 , N(3, 2)) + fix16_mul(s3 , N(3, 3));
    O(0, 3) = t0; O(1, 3) = t1; O(2, 3) = t2; O(3, 3) = t3;
#undef M
#undef N
#undef O

}

void CreateO2WMat(mat3_t* pRMat, v3d* pTranslation, v3d* pScale, mat4_t* pResult) {
    fix16_t sx = pScale->x; fix16_t sy = pScale->y; fix16_t sz = pScale->z;
    pResult->m03 = pTranslation->x; pResult->m13 = pTranslation->y; pResult->m23 = pTranslation->z;
    pResult->m00 = fix16_mul(pRMat->m00 , sx); pResult->m01 = fix16_mul(pRMat->m01 , sy) pResult->m02 = fix16_mul(pRMat->m02 , sz);
    pResult->m10 = fix16_mul(pRMat->m10 , sx); pResult->m11 = fix16_mul(pRMat->m11 , sy) pResult->m12 = fix16_mul(pRMat->m12 , sz);
    pResult->m20 = fix16_mul(pRMat->m20 , sx); pResult->m21 = fix16_mul(pRMat->m21 , sy) pResult->m22 = fix16_mul(pRMat->m22 , sz);
    pResult->m30 = 0; pResult->m31 = 0; pResult->m32 = 0; pResult->m33 = fix16_one;
}

void GenerateW2CMat(mat4_t* pResult, mat3_t* pSource, vect3_t* pVect) {
    fix16_t x = pVect->x;fix16_t y = pVect->y;fix16_t z = pVect->z;
    pResult->m00 = pSource->m00;
    pResult->m01 = pSource->m01;
    pResult->m02 = pSource->m02;
    pResult->m10 = pSource->m10;
    pResult->m11 = pSource->m11;
    pResult->m12 = pSource->m12;  
    pResult->m20 = pSource->m20;
    pResult->m21 = pSource->m21;
    pResult->m22 = pSource->m22; 
    pResult->m30 = 0; pResult->m31 = 0; pResult->m32 = 0; pResult->m33 = fix16_one;
    pResult->m03 =  fix16_mul(x , pResult->m00) +  fix16_mul(y , pResult->m01) +  fix16_mul(z , pResult->m02);
    pResult->m13 =  fix16_mul(x , pResult->m10) +  fix16_mul(y , pResult->m11) +  fix16_mul(z , pResult->m12);
    pResult->m23 =  fix16_mul(x , pResult->m20) +  fix16_mul(y , pResult->m21) +  fix16_mul(z , pResult->m22);

}

void MakeClipMatrix_perspective( fix16_t near_plane, fix16_t far_plane,
    fix16_t focalLength, fix16_t aspectRatio, mat4_t* pMat) {
#define M(x,y) (pMat)->m##x##y
    M(0, 0) = focalLength; M(1, 0) = 0;   M(2, 0) = 0;   M(3, 0) = 0;
    M(0, 1) = 0;   M(1, 1) = fix16_mul(focalLength , aspectRatio); M(2, 1) = 0;   M(3, 1) = 0;
    M(0, 2) = 0;   M(1, 2) = 0;   M(2, 2) = fix16_div(far_plane , (far_plane - near_plane)); M(3, 2) = fix16_one;
    M(0, 3) = 0;   M(1, 3) = 0;   M(2, 3) = -fix16_div(fix16_mul(near_plane ,far_plane) , (far_plane - near_plane));   M(3, 3) = 0;
#undef M
}
