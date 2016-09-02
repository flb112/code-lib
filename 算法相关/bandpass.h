#ifndef _BANDPASS_H_
#define _BANDPASS_H_

#define FREQ 6400	//����Ƶ��6.4K
typedef struct{
	float re;
	float im;
}complex;


class bandpass{
public:
	int freq;//����Ƶ��
	float wpl;	//ͨ�����޽�ֹƵ��
	float wph;	//ͨ�����޽�ֹƵ��
	float wsl;	//�������ֹƵ��
	float wsh;	//�������ֹƵ��
	float ap;	//ͨ�����˥��(dB)ͨ��Ϊ3dB
	float as;	//������˥��(dB)ͨ��Ϊ40dB

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

	
	
