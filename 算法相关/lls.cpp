/*
��С���˷�C++ʵ��
���� �� x
����� Ԥ���y  
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
        a = (t3*len - t2*t4) / (t1*len - t2*t2);  // ��æ�1 
        b = (t1*t4 - t2*t3) / (t1*len - t2*t2);        // ��æ�2
}

float LeastSquare::getY(const float x)
{
        return (a*x + b);
}