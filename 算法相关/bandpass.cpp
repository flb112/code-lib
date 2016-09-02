#include "bandpass.h"
#include<math.h>

float PI=3.141592654;

//函数声明

bandpass::bandpass(void)
{
	freq = FREQ;
	wpl = 300;	//通带下限截止频率
	wph = 1500;	//通带上限截止频率
	wsl = 500;	//下阻带截止频率
	wsh = 1000;	//上阻带截止频率
	ap = 0;	//通带最大衰减(dB)通常为3dB
	as = 0;	//阻带最大衰减(dB)通常为40dB
}
//主函数
float bandpass::filter(float val)
{
	//定义块
	complex z,s,p,eye,Gkp[1024];
	int Nk,m,k;//采样频率
	float C=0,N=0,Wsl=0,Wsh=0,Wpl=0,Wph=0,Ws=0,mone=0,Ntemp=0;//参照书上定义命名
	float ysl=0,ysh=0,ypl=0,yph=0,ys=0;
	for(m=0;m<1024;m++)
	{
		Gkp[m].re=1;
		Gkp[m].im=0;
	}

	//处理块，，，数字滤波器的指标转化为模拟滤波器的指标
//	N=ws;
//	C=2*PI;
	wpl=2*PI*wpl/freq;
	wph=2*PI*wph/freq;
	wsl=2*PI*wsl/freq;
	wsh=2*PI*wsh/freq;
	Wpl=tan(wpl/2);
	Wph=tan(wph/2);
	Wsl=tan(wsl/2);
	Wsh=tan(wsh/2);
	mone=Wph-Wpl;
	//归一化处理
	ysl=Wsl/mone;
	ysh=Wsh/mone;
	ypl=Wpl/mone;
	yph=Wph/mone;
	ys=(pow(ysh,2)-ypl*yph)/ysh;
	C=sqrt(pow(10,ap/10)-1);
	Ntemp=sqrt((pow(10,as/10)-1)/(pow(10,ap/10)-1));
	N=log(Ntemp)/log(ys);
	Nk=int(N)+1;
	for(m=1;m<=freq;m++)
	{
		z.re=cos(m*2*PI/freq);
		z.im=sin(m*2*PI/freq);
		eye.re=1;
		eye.im=0;
		s=cdiv(csub(z,eye),cadd(z,eye));
		p=cdiv(cadd(cmul(s,s),rmul(eye,Wph*Wpl)),rmul(s,mone));
		Ntemp=Nk%2;
		if(Ntemp==0)
		{
			for(k=1;k<=Nk/2;k++)
			{
			Gkp[m-1]=cdiv(Gkp[m-1],cadd(csub(cmul(p,p),rmul(p,2*cos((2*k+Nk-1)*PI/2/Nk))),eye));
			}
		}
		else
		{
			for(k=1;k<=(Nk-1)/2;k++)
			{
				Gkp[m-1]=cdiv(Gkp[m-1],cadd(csub(cmul(p,p),rmul(p,2*cos((2*k+Nk-1)*PI/2/Nk))),eye));
			}
			Gkp[m-1]=cdiv(Gkp[m-1],cadd(p,eye));
		}
		if(m<100)
		{
			cout<<Gkp[m-1].re<<'+'<<Gkp[m-1].im<<'i'<<endl;	
		}
	}
}



//定义复数减法
complex csub(complex c1,complex c2)
{
	complex c;
	c.re=c1.re-c2.re;
	c.im=c1.im-c2.im;
	return c;
}
//定义复数加法
complex cadd(complex c1,complex c2)
{
	complex c;
	c.re=c1.re+c2.re;
	c.im=c1.im+c2.im;
	return c;
}
//定义复数乘法
complex cmul(complex c1,complex c2)
{
	complex c;
	c.re=c1.re*c2.re-c1.im*c2.im;
	c.im=c1.re*c2.im+c1.im*c2.re;
	return c;
}
//定义复数除法
complex cdiv(complex c1,complex c2)
{
	complex c;
	float temp;
	temp=c2.re*c2.re+c2.im*c2.im;
	c.re=c1.re*c2.re+c1.im*c2.im;
	c.im=-c1.re*c2.im+c1.im*c2.re;
	c.re=c.re/temp;
	c.im=c.im/temp;
	return c;
}
complex rdiv(complex c1,float r)
{
	complex c;
	c.re=c1.re/r;
	c.im=c1.im/r;
	return c;
}
complex rmul(complex c1,float r)
{
	complex c;
	c.re=c1.re*r;
	c.im=c1.im*r;
	return c;
}

	
	
