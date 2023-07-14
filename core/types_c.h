#ifndef TYPES_C_H
#define TYPES_C_H

#ifdef __cplusplus
extern "C" {
#endif
#ifndef NULL
#define NULL ((void*)0)
#endif
    // 定义CvTermCriteria结构体
    typedef struct CvTermCriteria {
        int type;     /**< 选择优化算法的类型 */
        int max_iter; /**< 最大迭代次数 */
        double epsilon; /**< 精度条件 */
    } CvTermCriteria;
    // 定义CvPoint2D32f结构体1
    typedef struct CvPoint2D32f {
        float x; /**< X坐标 */
        float y; /**< Y坐标 */
    } CvPoint2D32f;


    enum NormTypes {
        /**
        \f[
        norm =  \forkthree
        {\|\texttt{src1}\|_{L_{\infty}} =  \max _I | \texttt{src1} (I)|}{if  \(\texttt{normType} = \texttt{NORM_INF}\) }
        {\|\texttt{src1}-\texttt{src2}\|_{L_{\infty}} =  \max _I | \texttt{src1} (I) -  \texttt{src2} (I)|}{if  \(\texttt{normType} = \texttt{NORM_INF}\) }
        {\frac{\|\texttt{src1}-\texttt{src2}\|_{L_{\infty}}    }{\|\texttt{src2}\|_{L_{\infty}} }}{if  \(\texttt{normType} = \texttt{NORM_RELATIVE | NORM_INF}\) }
        \f]
        */
        NORM_INF = 1,
        /**
        \f[
        norm =  \forkthree
        {\| \texttt{src1} \| _{L_1} =  \sum _I | \texttt{src1} (I)|}{if  \(\texttt{normType} = \texttt{NORM_L1}\)}
        { \| \texttt{src1} - \texttt{src2} \| _{L_1} =  \sum _I | \texttt{src1} (I) -  \texttt{src2} (I)|}{if  \(\texttt{normType} = \texttt{NORM_L1}\) }
        { \frac{\|\texttt{src1}-\texttt{src2}\|_{L_1} }{\|\texttt{src2}\|_{L_1}} }{if  \(\texttt{normType} = \texttt{NORM_RELATIVE | NORM_L1}\) }
        \f]*/
        NORM_L1 = 2,
        /**
        \f[
        norm =  \forkthree
        { \| \texttt{src1} \| _{L_2} =  \sqrt{\sum_I \texttt{src1}(I)^2} }{if  \(\texttt{normType} = \texttt{NORM_L2}\) }
        { \| \texttt{src1} - \texttt{src2} \| _{L_2} =  \sqrt{\sum_I (\texttt{src1}(I) - \texttt{src2}(I))^2} }{if  \(\texttt{normType} = \texttt{NORM_L2}\) }
        { \frac{\|\texttt{src1}-\texttt{src2}\|_{L_2} }{\|\texttt{src2}\|_{L_2}} }{if  \(\texttt{normType} = \texttt{NORM_RELATIVE | NORM_L2}\) }
        \f]
        */
        NORM_L2 = 4,
        /**
        \f[
        norm =  \forkthree
        { \| \texttt{src1} \| _{L_2} ^{2} = \sum_I \texttt{src1}(I)^2} {if  \(\texttt{normType} = \texttt{NORM_L2SQR}\)}
        { \| \texttt{src1} - \texttt{src2} \| _{L_2} ^{2} =  \sum_I (\texttt{src1}(I) - \texttt{src2}(I))^2 }{if  \(\texttt{normType} = \texttt{NORM_L2SQR}\) }
        { \left(\frac{\|\texttt{src1}-\texttt{src2}\|_{L_2} }{\|\texttt{src2}\|_{L_2}}\right)^2 }{if  \(\texttt{normType} = \texttt{NORM_RELATIVE | NORM_L2SQR}\) }
        \f]
        */
        NORM_L2SQR = 5,
        /**
        In the case of one input array, calculates the Hamming distance of the array from zero,
        In the case of two input arrays, calculates the Hamming distance between the arrays.
        */
        NORM_HAMMING = 6,
        /**
        Similar to NORM_HAMMING, but in the calculation, each two bits of the input sequence will
        be added and treated as a single bit to be used in the same calculation as NORM_HAMMING.
        */
        NORM_HAMMING2 = 7,
        NORM_TYPE_MASK = 7, //!< bit-mask which can be used to separate norm type from norm flags
        NORM_RELATIVE = 8, //!< flag
        NORM_MINMAX = 32 //!< flag
    };

    //! type of the threshold operation
    enum {
        THRESH_BINARY = 0, // value = value > threshold ? max_value : 0
        THRESH_BINARY_INV = 1, // value = value > threshold ? 0 : max_value
        THRESH_TRUNC = 2, // value = value > threshold ? threshold : value
        THRESH_TOZERO = 3, // value = value > threshold ? value : 0
        THRESH_TOZERO_INV = 4, // value = value > threshold ? 0 : value
        THRESH_MASK = 7,
        THRESH_OTSU = 8  // use Otsu algorithm to choose the optimal threshold value
    };

   /* typedef struct rect rect;
    struct rect {
        short xmin;
        short ymin;
        short xmax;
        short ymax;
        short classNum;
    };

    typedef struct box_rect box_rect;
    struct box_rect
    {
        int xmin;
        int ymin;
        int width;
        int height;
    };*/
  

    typedef struct rect rect;
    struct rect {
        short xmin;
        short ymin;
        short xmax;
        short ymax;
        short classNum;
    };

    typedef struct box_rect box_rect;
    struct box_rect
    {
        int xmin;
        int ymin;
        int width;
        int height;
    };

#ifdef __cplusplus
}
#endif

#endif  // TYPES_C_H

