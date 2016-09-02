/*
最小二乘法C++实现
输入 ： x
输出： 预测的y  
*/
#include "stdafx.h"
#include "lls.h"

void LeastSquare::init(const float *x, const float *y, int len)
{
        double t1=0, t2=0, t3=0, t4=0;
        for(int i=0; i<len; ++i)
        {
            t1 += x[i]*x[i];
            t2 += x[i];
            t3 += x[i]*y[i];
            t4 += y[i];
        }
        a = (t3*len - t2*t4) / (t1*len - t2*t2);  // 求得β1 
        b = (t1*t4 - t2*t3) / (t1*len - t2*t2);        // 求得β2
}

float LeastSquare::getY(const float x)
{
        return (a*x + b);
}