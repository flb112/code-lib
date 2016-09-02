#ifndef _LEAST_SQUARE_
#define _LEAST_SQUARE_
class LeastSquare{
    double a,b;
public:
    LeastSquare(){};

    void init(const float *x, const float *y, int len);
    float getY(const float x);

};

#endif //_LEAST_SQUARE_