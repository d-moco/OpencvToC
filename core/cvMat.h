#ifndef CVMAT_H
#define CVMAT_H

#include <stdint.h>
#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif

#define true        1
#define false       0

#define CV_CN_MAX     0x200//512
#define CV_CN_SHIFT   3
#define CV_DEPTH_MAX  (1 << CV_CN_SHIFT)

#define CV_8U   0
#define CV_8S   1
#define CV_16U  2
#define CV_16S  3
#define CV_32S  4
#define CV_32F  5
#define CV_64F  6
#define CV_16F  7

    //CV_MAT_DEPTH 将二进制位高于 CV_CN_SHIFT 的部分全部置 0，确保 depth 的范围在 [0, 7]
#define CV_MAT_DEPTH_MASK       (CV_DEPTH_MAX - 1)
#define CV_MAT_DEPTH(flags)     ((flags) & CV_MAT_DEPTH_MASK)

//CV_MAKETYPE 后半段将通道cn减1，越过 CV_MAT_DEPTH 放置
#define CV_MAKETYPE(depth,cn) (CV_MAT_DEPTH(depth) + (((cn)-1) << CV_CN_SHIFT))

//CV_MAT_CN 取回通道 cn 的值
#define CV_MAT_CN_MASK          ((CV_CN_MAX - 1) << CV_CN_SHIFT)
#define CV_MAT_CN(flags)        ((((flags) & CV_MAT_CN_MASK) >> CV_CN_SHIFT) + 1)

/* Size of each channel item,
   0x124489 = 1000 0100 0100 0010 0010 0001 0001 ~ array of sizeof(arr_type_elem) */
#define CV_ELEM_SIZE1(type) \
    ((((sizeof(size_t)<<28)|0x8442211) >> CV_MAT_DEPTH(type)*4) & 15)

   /* 0x3a50 = 11 10 10 01 01 00 00 ~ array of log2(sizeof(arr_type_elem)) */
#define CV_ELEM_SIZE(type) \
    (CV_MAT_CN(type) << ((((sizeof(size_t)/4+1)*16384|0x3a50) >> CV_MAT_DEPTH(type)*2) & 3))

#define CV_8UC1 CV_MAKETYPE(CV_8U,1)
#define CV_8UC2 CV_MAKETYPE(CV_8U,2)
#define CV_8UC3 CV_MAKETYPE(CV_8U,3)
#define CV_8UC4 CV_MAKETYPE(CV_8U,4)
#define CV_8UC(n) CV_MAKETYPE(CV_8U,(n))

#define CV_8SC1 CV_MAKETYPE(CV_8S,1)
#define CV_8SC2 CV_MAKETYPE(CV_8S,2)
#define CV_8SC3 CV_MAKETYPE(CV_8S,3)
#define CV_8SC4 CV_MAKETYPE(CV_8S,4)
#define CV_8SC(n) CV_MAKETYPE(CV_8S,(n))

#define CV_16UC1 CV_MAKETYPE(CV_16U,1)
#define CV_16UC2 CV_MAKETYPE(CV_16U,2)
#define CV_16UC3 CV_MAKETYPE(CV_16U,3)
#define CV_16UC4 CV_MAKETYPE(CV_16U,4)
#define CV_16UC(n) CV_MAKETYPE(CV_16U,(n))

#define CV_16SC1 CV_MAKETYPE(CV_16S,1)
#define CV_16SC2 CV_MAKETYPE(CV_16S,2)
#define CV_16SC3 CV_MAKETYPE(CV_16S,3)
#define CV_16SC4 CV_MAKETYPE(CV_16S,4)
#define CV_16SC(n) CV_MAKETYPE(CV_16S,(n))

#define CV_32SC1 CV_MAKETYPE(CV_32S,1)
#define CV_32SC2 CV_MAKETYPE(CV_32S,2)
#define CV_32SC3 CV_MAKETYPE(CV_32S,3)
#define CV_32SC4 CV_MAKETYPE(CV_32S,4)
#define CV_32SC(n) CV_MAKETYPE(CV_32S,(n))

#define CV_32FC1 CV_MAKETYPE(CV_32F,1)
#define CV_32FC2 CV_MAKETYPE(CV_32F,2)
#define CV_32FC3 CV_MAKETYPE(CV_32F,3)
#define CV_32FC4 CV_MAKETYPE(CV_32F,4)
#define CV_32FC(n) CV_MAKETYPE(CV_32F,(n))

#define CV_64FC1 CV_MAKETYPE(CV_64F,1)
#define CV_64FC2 CV_MAKETYPE(CV_64F,2)
#define CV_64FC3 CV_MAKETYPE(CV_64F,3)
#define CV_64FC4 CV_MAKETYPE(CV_64F,4)
#define CV_64FC(n) CV_MAKETYPE(CV_64F,(n))

#define CV_16FC1 CV_MAKETYPE(CV_16F,1)
#define CV_16FC2 CV_MAKETYPE(CV_16F,2)
#define CV_16FC3 CV_MAKETYPE(CV_16F,3)
#define CV_16FC4 CV_MAKETYPE(CV_16F,4)
#define CV_16FC(n) CV_MAKETYPE(CV_16F,(n))

#define Sf( y, x ) ((float*)(srcdata + y*srcstep))[x]
#define Sd( y, x ) ((double*)(srcdata + y*srcstep))[x]
#define Df( y, x ) ((float*)(dstdata + y*dststep))[x]
#define Dd( y, x ) ((double*)(dstdata + y*dststep))[x]


    /****************************************************************************************\
*                                 Determinant of the matrix                              *
\****************************************************************************************/

#define det2(m)   ((double)m(0,0)*m(1,1) - (double)m(0,1)*m(1,0))
#define det3(m)   (m(0,0)*((double)m(1,1)*m(2,2) - (double)m(1,2)*m(2,1)) -  \
                   m(0,1)*((double)m(1,0)*m(2,2) - (double)m(1,2)*m(2,0)) +  \
                   m(0,2)*((double)m(1,0)*m(2,1) - (double)m(1,1)*m(2,0)))
// matrix decomposition types
    enum {
        DECOMP_LU = 0,
        DECOMP_SVD = 1,
        DECOMP_EIG = 2,
        DECOMP_CHOLESKY = 3,
        DECOMP_QR = 4,
        DECOMP_NORMAL = 16
    };

    typedef struct _mat
    {
        int type;//CV_8UC3
        int step;
        int* refcount;
        unsigned char* ptr;
        short* s;
        int* i;
        float* fl;
        double* db;

        int rows;
        int cols;
        int channels;
        int Empty;  // 新增的成员变量，表示数组是否为空
    }CvMat;
   
    typedef struct _size
    {
        int width;
        int height;
    }CvSize;

    typedef struct {
        CvMat* data;
        size_t size;
        size_t capacity;
       
    } UMatVector;
    typedef struct {
        float* data;
        int size;
    } AutoBuffer;
    CvMat createMat(int rows, int cols, int type);
    void createMat_s(CvMat* mat, int rows, int cols, int type);
    void zerosMat(CvMat *datya);
    void freeMat(CvMat *_mat);
    CvMat getCvMat(CvMat* mat);

    void cvmat_vector_init(UMatVector* vec, size_t capacity);
    void cvmat_vector_push_back(UMatVector* vec, CvMat value);
    void cvmat_vector_free(UMatVector* vec);

    void saveImage(CvMat* data, const char* zOut);

    CvSize cvGetMatSize(const CvMat* mat);
    double cvGetArea(CvSize cs);
    void cvMatDev(const CvMat*src ,CvMat* dst,  float dev);
    void cvMatMul(const CvMat* src, double mul);
    void cvconvertTo(const CvMat* src, CvMat* dst,int type);
    void cvcopyTo(const CvMat* src,CvMat* dst);
    void cvAppendMatAdd(CvMat* src, int row, int col, float value);
    float cvGetMatValue(CvMat* src,int row,int col);
    CvMat cvInvert(CvMat* src,int method);
    void cvSetIdentity(CvMat* src);
    void cvDiag(CvMat* src,const int d);
    
#ifdef __cplusplus
}
#endif /* end of __cplusplus */





#endif