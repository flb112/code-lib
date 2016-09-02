#ifndef _BANDPASS_H_
#define _BANDPASS_H_

#define FREQ 6400	//采用频率6.4K
typedef struct{
	float re;
	float im;
}complex;


class bandpass{
public:
	int freq;//采样频率
	float wpl;	//通带下限截止频率
	float wph;	//通带上限截止频率
	float wsl;	//下阻带截止频率
	float wsh;	//上阻带截止频率
	float ap;	//通带最大衰减(dB)通常为3dB
	float as;	//阻带最大衰减(dB)通常为40dB

public:
	bandpass(void);
	float filter(float val);
private:
	complex csub(complex c1,complex c2);
	complex cadd(complex c1,complex c2);
	complex cmul(complex c1,complex c2);
	complex cdiv(complex c1,complex c2);
	complex rdiv(complex c1,float r);
	complex rmul(complex c1,float r);
};



#endif

	
	
