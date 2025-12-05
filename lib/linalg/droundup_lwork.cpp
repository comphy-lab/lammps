#ifdef __cplusplus
extern "C" {
#endif
#include "lmp_f2c.h"
doublereal droundup_lwork__(integer *lwork)
{
    const doublereal eps = 2.2204460492503131E-16;
    doublereal ret_val;
    ret_val = (doublereal)(*lwork);
    if ((integer)ret_val < *lwork) {
        ret_val *= (1.000000000000000 + eps);
    }
    return ret_val;
}
#ifdef __cplusplus
}
#endif
