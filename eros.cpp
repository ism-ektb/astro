// EROS.cpp : main project file.

#include "stdafx.h"
#include "iostream"
#include "fstream"
#define _USE_MATH_DEFINES
#include <math.h>
#include <string>
#include <ctime>
#include <Windows.h>
#pragma comment(linker, "/STACK:10000000")
using namespace System;
using namespace System::IO;
using namespace std;
using namespace DeLibrary;

double XSTART[3], VSTART[3], XEND[3], VEND[3], F[30], X[30], V[30], W[13], U[13];
double C[78], D[78], R[78], XI[78], H[14], Julian_Date[1000], eph_t[1000], poz[6];
double ul[10000], del[10000], al[10000], abss[10000], planet_centric[10000];
double *BT, *BE, *F1, *FJ, *Y, *Z, *B, TDIF, DIR, PW, TP, SR, W1, TM, SM, SS, Second, TS, TF, TSTART, TEND, time_abb, rad, alpha, delta;
double T_user_start, T_user_finish, user_step, T_user, a_Kepler, e_Kepler, i_Kepler, knot_Kepler, arg_Kepler, M_Kepler, poz1[6];
double h_alpha_eph, m_alpha_eph, s_alpha_eph, d_delta_eph, m_delta_eph, s_delta_eph;
double parameter_G, abs_mag, magnitude;
double X_abber[30], V_abber[30];
double X_in[30], V_in[30];
double sigma;
double alpha_dom, delta_dom;
double magnitude_vis;
double X_boul, Y_boul, Z_boul, VX_boul, VY_boul, VZ_boul;
struct Asteroid
{
	double a;
	double q;
	double p;
	double e;
	double i;
	double w;
	double node;
};

int nclass, nor, ni, ll, nv ,NOR, NCLASS, NV, LL, N_CLASS, NI, NF, KD, KF, KE;
int NPQ, NSF, NES, NCL, NPER, NCOUNT, NL, NS, Year, Month, Day, Hour, Minute, eph_count, aster_count;
int NW[14], MC[12];
int count_iter;
int s_c;

float xl, XL;

string str_spisok;
string path_out;

char * binPathForHeader = "In\\header.bin"; 
char * binPathForData = "In\\data.bin";

double massiv_mass[14]={6023600, 408523.71, 332946.050895, 3098708, 1047.34860, 3497.8980, 22902.980,
		                19412.240, 135200000, 27068700.387534, 1, 0.21276595744680851063829787234043e+10,
	                    1e+10, 0.76923076923076923076923076923077e+10};          

int force_var[16]={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  //               0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15
double TAST=2455500.5;
double massiv_Ceres[6]={1.155539633848469760, -2.335075520271295040, -1.335376831103344320, 
                     0.008978098952189530, 0.003812198314654010, -0.000032100530148940};

double massiv_Pallas[6]={-0.423783444796719840, -3.071749139849307520, 0.638279965102497280, 
                     0.008480892389052840, -0.002872960453692960, -0.000022944424181820};

double massiv_Vesta[6]={-1.553049265833779520, -1.474955583536473280, -0.384432998870605760, 
                     0.008702756162718960, -0.007371566342941650, -0.004073518782575940};

double massiv_Test[6]={-1.876167736180090880, 1.483497519729738560, 0.354186848392725840, 
					-0.009214306680071010, -0.006941186207782940, -0.000936185761031680};

const double AE=149597870.691;
const double mu=0.0002959122082855911025;
const double Earth_Garmonica=0.00108263;
const double Sun_Garmonica=1.9E-7;
const double Jupiter_Garmonica=1.4736E-2;
const double Radius_Earth=6378.140;
const double Radius_Sun=695508;
const double Radius_Jupiter=71492;
const double i_Sun=7.25*M_PI/180;
const double cos_i_Sun=cos(i_Sun);
const double sin_i_Sun=sin(i_Sun);
const double Sun_oblateness=mu*Sun_Garmonica*(Radius_Sun/AE)*(Radius_Sun/AE);
const double eps=(23*60*60+26*60+21.448)*M_PI/648000;
const double cc=299792.458/AE*86400;
const double mu1=1.0027379093;
const double relative_aux=mu/(cc*cc);
const double Earth_sec_garmonica=mu/massiv_mass[2]*Earth_Garmonica*((Radius_Earth*Radius_Earth)/(AE*AE));
const double Jupiter_sec_garmonica=mu/massiv_mass[4]*Jupiter_Garmonica*((Radius_Jupiter*Radius_Jupiter)/(AE*AE));

double Date_JD (int Year_, int Month_, int Day_, int Hour_, int Minute_, double Second_)
{
	double b=0,d, data;
	int m0, y0;
	long int y, m;
	m=Month_;
	y=Year_;
	y0=y;
	m0=m;
	d=Day_;
	
	if(m<3)
	{
		m0=m0+12;
		y0=y0-1;
	}

	if ((y0>1582)||(y0=1582)&&(m0>10)||(y0=1582)&&(m0=10)&&(d>=15.0))
	{
		b=2-int(y0/100)+int(int(y0/100)/4);
		data=int(365.25*y0)+int((30.6001)*(m0+1))+d+b+1720994.5;
	}

    data=data+(Hour_/24.)+(Minute_/1440.)+(Second_/86400.);
	return (data);
}
//вычисляет по календарной дате юлианскую дату
void Date_Grig (double Grig_time)
{
	int ich, icm, ice;
	double ch, f4; long int a1, a5, b0, c0, d0, e3, z0;
	z0=int (Grig_time+0.5);
	f4=Grig_time+0.5-z0;
	if(z0<2299161)
	{
		a5=z0;
	}
	else 
	{
		a1=int((z0-1867216.25)/36524.25);
		a5=z0+1+a1-int(a1/4.0);
	}

	b0=a5+1524;
	c0=int((b0-122.1)/365.25);
	d0=int(365.25*c0);
	e3=int((b0-d0)/30.6001);
	ch=b0-d0-int(30.6001*e3)+f4;
	ich=int(ch);
	
	if(e3>13.5)
	{
		icm=e3-13;
	}
	else
	{
		icm=e3-1;
		ice=c0-4715;
	}

	if(icm>2.5)
	{
		ice=ice-1;
	}
	int h, m;
	int s;
	h=int((ch-ich)*24);
	m=int((((ch-ich)*24)-h)*60);
	s=int(((((ch-ich)*24)-h)*60-m)*60);

	Year=ice; Month=icm; Day=ich; Hour=h; Minute=m; Second=s;
} 
//вычисляет по юлианской дате календарную дату
void eph_open()
{	 
	ifstream fine;
    eph_count=0; 
	int Year_, Month_, Day_;
	int Hour_=0, Minute_=0; 
	double Second_=0;
    fine.open("In\\dtime.txt");

	while(!fine.eof())
	{
		fine>>Year_>>Month_>>Day_>>eph_t[eph_count];
		Julian_Date[eph_count]=Date_JD(Year_, Month_, Day_, Hour_, Minute_, Second_);
		eph_count++;
	}
	
	eph_count=eph_count-1;
	fine.close();
} 
//Открывает список эфемеридных поправок
double eph_time(double epher)
{
	eph_open();
	int i=0;
	double eph_time;
	double J1, J2, J0;
	double deph1, deph2, deph0;
	J0=epher;
	J2=Julian_Date[eph_count];
	if(J0>J2)
	{
		deph0=eph_t[eph_count];
	}
	else
	{	
		J2=Julian_Date[0];
		if(J0<J2)
		{
			deph0=eph_t[0];
		}
		else
		{
			while(J0>J2)
			{
				i++;
				J2=Julian_Date[i];
			}
			J2=Julian_Date[i];
			J1=Julian_Date[i-1];
			deph1=eph_t[i-1];
			deph2=eph_t[i];
			deph0=deph1+(deph2-deph1)/(J2-J1)*(J0-J1);
		}
	}
	eph_time=J0+deph0/(60*60*24);
	return(eph_time);
}
 //переводит юлианскую дату к эфемеридной шкале
double star_time(double time_star)
{
	double tt;
	double s_time;
	tt=(time_star-2451545.0)/36525.;
	s_time=(24110.54841+8640184.812866*tt+0.0931048*tt*tt-0.00000621*tt*tt*tt);
	return(s_time);
}
//переводит юлианскую дату в звездное время
string read_boul(string name_ast)
{   
	ifstream fin_cat;
	fin_cat.open("Irina_in\\astorb.dat");
	if (!fin_cat) exit(500);
	int m=-1, num=7, j;
	string str_cat;

	while(m!=0) 
	{   getline(fin_cat,str_cat);
		if(str_cat=="")
		{cout<<name_ast<<endl<<"error";
		break;}
		int r=0; 
		for(j=num;j<name_ast.length()+num;j++)
		{
			if(name_ast.at(j-num)==str_cat.at(j))
			{
				r++;}}
		int qq=str_cat.at(num+name_ast.length());
		if(r==name_ast.length() && qq==32)
		{m=0; }
		r=0;}
	fin_cat.close();
	return str_cat;
}
//открывает каталог Боуэлла
double boul_abs(string str_cat) 
{   double abs=0; 
	abs=abs+0.1*(str_cat.at(45)-48); 
	if(str_cat.at(46)!=32)
	{abs=abs+0.01*(str_cat.at(46)-48);}
	abs=abs+(str_cat.at(43)-48);
	if(str_cat.at(42)!=32)
	{abs=abs+10.*(str_cat.at(42)-48);}
	return abs;}
//звездная величина из боуэлла
double boul_osc(string str_cat)
{   int y_start, m_start, d_start;
	y_start=0;
	int k=3;
	for(int i=106;i<=110;i++)
	{
		y_start+=(str_cat.at(i)-48)*pow(10,k);
		k--;
	}
	m_start=(str_cat.at(110)-48)*10+str_cat.at(111)-48;
	d_start=(str_cat.at(112)-48)*10+str_cat.at(113)-48;
	return Date_JD(y_start, m_start, d_start, 0, 0, 0);
}
//момент времени из боуэлла
double boul_axis(string str_cat)
{   double axis=0;
	int r=-1, j;
	for(j=173;j<181;j++)
	{axis=axis+pow(10.,r)*(str_cat.at(j)-48);
	 r--;}
	axis=axis+(str_cat.at(171)-48);
	return axis;}
//большая полуось из боуэлла
double boul_i(string str_cat)
{   double i=0;
	int r=-1, j=0;
	for(j=151;j<=156;j++)
	{i=i+pow(10.,r)*(str_cat.at(j)-48);
	 r--;}
	i=i+(str_cat.at(149)-48);
	if(str_cat.at(148)!=32)
	{i=i+10.*(str_cat.at(148)-48);}
	return i;}
//наклонение из боуэлла
double boul_e(string str_cat)
{   double e=0;
	int r=-1, j;
	for(j=160;j<=167;j++)
	{e=e+pow(10.,r)*(str_cat.at(j)-48);
	 r--;}
	e=e+int(str_cat.at(158)-48);
	return e;}
//эксцентриситет из боуэлла
double boul_arg(string str_cat)
{   double arg=0;
	int r=-1, j;
	for(j=130;j<=135;j++)
	{arg=arg+pow(10.,r)*(str_cat.at(j)-48);
	 r--;}
	arg=arg+(str_cat.at(128)-48);
	if(str_cat.at(127)!=32)
	{arg=arg+10.*(str_cat.at(127)-48);}
	if(str_cat.at(126)!=32)
	{arg=arg+100.*(str_cat.at(126)-48);}
	return arg;}
//аргумент перицентра из боуэлла
double boul_anomaly(string str_cat)
{   double anomaly=0;
	int r=-1, j;
	for(j=119;j<=124;j++)
	{anomaly=anomaly+pow(10.,r)*(str_cat.at(j)-48);
	 r--;}
	anomaly=anomaly+(str_cat.at(117)-48);
	if(str_cat.at(116)!=32)
	{anomaly=anomaly+10.*(str_cat.at(116)-48);}
	if(str_cat.at(115)!=32)
	{anomaly=anomaly+100.*(str_cat.at(115)-48);}
	return anomaly;}
//средняя аномалия из боуэлла
double boul_knot(string str_cat)
{   double knot=0;
	int r=-1, i;
	for(i=141;i<=146;i++)
	{knot=knot+pow(10.,r)*(str_cat.at(i)-48);
	 r--;}
	knot=knot+(str_cat.at(139)-48);
	if(str_cat.at(137)!=32)
	{knot=knot+100.*(str_cat.at(137)-48);}
	if(str_cat.at(138)!=32)
	{knot=knot+10.*(str_cat.at(138)-48);}
	return knot;}
//долгота восходящего узла из боуэлла
double boul_g(string str_cat)
{
double g=0;
int r=-1,
i;
for(i=51;i<=52;i++)
{
g=g+pow(10.,r)*(str_cat.at(i)-48);
r--;
}
g=g+(str_cat.at(49)-48);
return g;
}
//параметр G из боуэлла
string read_observatory(string name_num)
{   ifstream fin_spisok;
	fin_spisok.open("In\\obser.dat");
	int m=-1, j;
	while(m!=0)
	{getline(fin_spisok,str_spisok);
	if(str_spisok=="")
		{cout<<name_num<<endl<<"error";
		break;}
		int r=0; 
		for(j=0;j<3;j++)
		{if(name_num.at(j)==str_spisok.at(j))
			{r++; }}
		if(r==3)
		{m=0; }
		r=0; }
	fin_spisok.close();
	return str_spisok;}
//Открывает каталог обсерваторий
double names_lon(string str_spisok)
{   double lon=0;
	if(str_spisok.at(4)!=32)
	{lon=lon+100*(str_spisok.at(4)-48);}
	if(str_spisok.at(5)!=32)
	{lon=lon+10*(str_spisok.at(5)-48);}
	lon=lon+str_spisok.at(6)-48;
	int j, r=-1;
	for(j=8;j<=12;j++)
	{if(str_spisok.at(j)!=32)
		{lon=lon+pow(10.,r)*(str_spisok.at(j)-48);}
		r--;}
	return lon;}
//Считывает первую координату из каталога Обс
double names_alt(string str_spisok)
{   double alt=0;
	alt=alt+10*(str_spisok.at(14)-48);
	int j, r=-1;
	for(j=16;j<=21;j++)
	{if(str_spisok.at(j)!=32)
		{alt=alt+pow(10.,r)*(str_spisok.at(j)-48);}
		r--;}
	return alt;}
//Считывает вторую координату из каталога Обс
double names_hren(string str_spisok)
{   double hren=0;
	hren=hren+10*(str_spisok.at(24)-48);
	int j, r=-1;
	for(j=26;j<=31;j++)
	{if(str_spisok.at(j)!=32)
		{hren=hren+pow(10.,r)*(str_spisok.at(j)-48);}
		r--;}
	if(str_spisok.at(23)==45)
	{hren=hren*(-1);}
	return hren;}
//Считывает третью координату из каталога Обс
void createfondbin(){
	char * savePathForHeader = "In\\header.bin"; 
	char * savePathForData = "In\\data.bin";
	char * pathtoAsciiheader = strcat( DeAcsessor::PathToAssembly(), "\\read405\\header.405" );
	char *asciiData[31] = {
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1600.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1620.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1640.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1660.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1680.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1700.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1720.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1740.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1760.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1780.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1800.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1820.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1840.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1860.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1880.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1900.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1920.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1940.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1960.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp1980.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2000.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2020.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2040.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2060.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2080.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2100.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2120.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2140.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2160.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2180.405" ),
		strcat(DeAcsessor::PathToAssembly(), "\\read405\\ascp2200.405" )};
		DeAcsessor^ acs = gcnew DeAcsessor(pathtoAsciiheader,31,&asciiData[0],savePathForData,savePathForHeader);}
//создает бинарные файлы каталогов
const int N=100;
void read_obs(ifstream&fin,double data_ul[],double del[],double al[],double abs[],string obs[],int&n,double s_n[],double s_x[],double s_y[],double s_z[])
{
	string str;
	double ul;
	int i, qq=0, j=0, ii, k;
	for(i=0;i<5000;i++)
	{
		obs[i]="   ";
	}
	for(i=0;i<1000;i++)
	{
		s_n[i]=.5;
	}
	i=0;
	n=0;
	while(!fin.eof())
	{
		getline(fin,str);
	
		if ( str.length() > 0 )
		{

		if(str.at(14)!=82 && str.at(14)!=114 && qq!=1)
		{
			if(str.at(14)==83 || str.at(14)==115)
			{
				qq=1;
			}
		int y=(str.at(15)-48)*1000+(str.at(16)-48)*100+(str.at(17)-48)*10+(str.at(18)-48),
			m=(str.at(20)-48)*10+(str.at(21)-48),
			d=(str.at(23)-48)*10+(str.at(24)-48);
		double q=((str.at(26)-48)*.1+(str.at(27)-48)*.01+(str.at(28)-48)*.001+(str.at(29)-48)*.0001+(str.at(30)-48)*.00001);
		if(str.at(31)!=32)
		{
			q=q+(str.at(31)-48)*.000001;
		}
		int h=q*24,
			min=(q*24-h)*60;
		double s=((q*24-h)*60-min)*60;
		ul=Date_JD(y,m,d,h,min,s);
		if(i==0 || ul!=data_ul[i-1])
		{
			if(qq==1)
			{
				s_n[j]=i;
			}
			
				data_ul[i]=ul;
			
			
			q=(str.at(45)-48)*10+(str.at(46)-48)+((str.at(48)-48)*10+(str.at(49)-48))/60.+((str.at(51)-48)*10+(str.at(52)-48)+(str.at(54)-48)*.1)/3600.;
			if(str.at(55)!=32)
			{
				q=q+(str.at(55)-48)*.01/3600;
			}
			q=q*M_PI/180;
			if(str.at(44)==45)
			{
				q=q*(-1);
			}

			
				del[i]=q;
			
			
			q=(str.at(32)-48)*10+(str.at(33)-48)+((str.at(35)-48)*10+(str.at(36)-48))/60.+((str.at(38)-48)*10+(str.at(39)-48)+(str.at(41)-48)*.1+(str.at(42)-48)*.01)/3600.;
			if(str.at(43)!=32)
			{
				q=q+(str.at(43)-48)*.001/3600;
			}
			q=q*15*M_PI/180;

			
				al[i]=q;
			
			
			if(str.at(65)==32)
			{
				
					abs[i]=NULL;
				
			}
			else
			{
				q=(str.at(65)-48)*10+(str.at(66)-48)+(str.at(68)-48)*.1;
				if(str.at(69)!=32)
				{
					q=q+(str.at(69)-48)*.01;
				}
				
					abs[i]=q;
				
			}
			
			
				obs[i].at(0)=str.at(77);
				obs[i].at(1)=str.at(78);
				obs[i].at(2)=str.at(79);
				i++;
				n++;			
		}
		}
		if(qq==1)
		{
			getline(fin,str);
			k=3;
			s_x[j]=0;
			for(ii=36;ii<=44;ii++)
			{
				if(str.at(ii)!=46 && str.at(ii)!=32)
				{
					s_x[j]+=(str.at(ii)-48)*pow(10,k);
					k--;
				}
				if(str.at(ii)==32)
				{
					k--;
				}
			}
			s_x[j]/=AE;
			if(str.at(34)==45)
			{
				s_x[j]*=-1;
			}
			k=3;
			s_y[j]=0;
			for(ii=48;ii<=56;ii++)
			{
				if(str.at(ii)!=00 && str.at(ii)!=32)
				{
					s_y[j]+=(str.at(ii)-48)*pow(10,k);
					k--;
				}
				if(str.at(ii)==32)
				{
					k--;
				}
			}
			s_y[j]/=AE;
			if(str.at(46)==45)
			{
				s_y[j]*=-1;
			}
			k=3;
			s_z[j]=0;
			for(ii=60;ii<=68;ii++)
			{
				if(str.at(ii)!=46 && str.at(ii)!=32)
				{
					s_z[j]+=(str.at(ii)-48)*pow(10,k);
					k--;
				}
				if(str.at(ii)==32)
				{
					k--;
				}
			}
			s_z[j]/=AE;
			if(str.at(58)==45)
			{
				s_z[j]*=-1;
			}
			qq=0;
			j++;
		}
	    }	
	}
	
	for(j=i;j<5000;j++)
	{
		data_ul[j]=0;
		del[j]=0;
		al[j]=0;
		abs[j]=0;
	}
	
}
//работает с файлом наблюдений	
void gelio_to_planet (double *X, double *V, double *mass_coor)
{int i;
 for(i=0; i<3; i++)
 {planet_centric[i]=X[i]-mass_coor[i];
  planet_centric[i+3]=V[i]-mass_coor[i+3];}}
//переводит из гелиоцентрической в планетоцентрическую
void eklekv(double ekl[], double ekv[])
{	
	ekv[0]=ekl[0];
	ekv[3]=ekl[3];
	ekv[1]=ekl[1]*cos(eps)-ekl[2]*sin(eps);
	ekv[2]=ekl[1]*sin(eps)+ekl[2]*cos(eps);
	ekv[4]=ekl[4]*cos(eps)-ekl[5]*sin(eps);
	ekv[5]=ekl[4]*sin(eps)+ekl[5]*cos(eps);
}
//переводит координаты из эклиптической системы в экваториальную
void ekvekl(double ekl[],double ekv[])
{
	ekl[0]=ekv[0];
	ekl[3]=ekv[3];
	ekl[1]=ekv[1]*cos(eps)+ekv[2]*sin(eps);
	ekl[2]=-ekv[1]*sin(eps)+ekv[2]*cos(eps);
	ekl[4]=ekv[4]*cos(eps)+ekv[5]*sin(eps);
	ekl[5]=-ekv[4]*sin(eps)+ekv[5]*cos(eps);
 }
//переводит координаты из экваториальной системы в эклиптическую
void h_a(string name_num, double jd0, double alfa, double del, int&h_grad, int&h_min, double&h_sec, int&a_grad, int&a_min, double&a_sec)
{double hh;
string str_spisok=read_observatory(name_num);
double t0=floor(jd0+0.5)-0.5;
double t=(star_time(t0));
t=t+mu1*(jd0-t0)*24*60*60;
while(t>86400)
{t=t-86400;}
while(t<0)
{t=t+86400;}
t=t/3600.;
hh=names_lon(str_spisok);
t=t+hh/15;
t=t-alfa;
t=t/12.*M_PI;
double hren=names_hren(str_spisok);
double alt=names_alt(str_spisok);
double tan_fi=names_hren(str_spisok)/names_alt(str_spisok),
       rad_fi=atan(tan_fi);
if(rad_fi>M_PI/2)
{rad_fi=rad_fi-M_PI;}
else
{if(rad_fi<M_PI/(-2))
{rad_fi=rad_fi+M_PI;}}
del=del/180*M_PI;
double cos_z=sin(rad_fi)*sin(del)+cos(rad_fi)*cos(del)*cos(t);
	double tan_a=(cos(del)*sin(t))/((-cos(rad_fi)*sin(del)+sin(rad_fi)*cos(del)*cos(t)));
    double arc_a=atan(tan_a);
	double sin_a=cos(del)*sin(t);
	double cos_a=-cos(rad_fi)*sin(del)+sin(rad_fi)*cos(del)*cos(t);
	if(cos_a<0)
		{arc_a=arc_a+M_PI;}

double rad_z=acos(cos_z);
if(rad_z>M_PI)
{rad_z=2*M_PI-rad_z;}
else
{if(rad_z<0)
{rad_z=abs(rad_z);}}
double grad_h=90-rad_z*180/M_PI;
h_grad=int(grad_h);
h_min=int((grad_h-h_grad)*60);
h_sec=((grad_h-h_grad)*60-h_min)*60;
double grad_a=arc_a*180/M_PI;
if(grad_a<0)
{grad_a=grad_a+360;}
if(grad_a>360)
{grad_a=grad_a-360;}
a_grad=int(grad_a);
a_min=int((grad_a-a_grad)*60);
a_sec=((grad_a-a_grad)*60-a_min)*60;
}
//перевод вторую экваториальную систему в горизонтальную
void XYZ_boul(string name_asteroid)
{ string str;
  str=read_boul(name_asteroid);


  double axis_boul, i_boul, e_boul, arg_boul, anomaly_boul, knot_boul;
  axis_boul=boul_axis(str);
  i_boul=boul_i(str);
  e_boul=boul_e(str);
  arg_boul=boul_arg(str);
  anomaly_boul=boul_anomaly(str);
  knot_boul=boul_knot(str);
  i_boul=i_boul/180.*M_PI;
  arg_boul=arg_boul/180.*M_PI;
  anomaly_boul=anomaly_boul/180.*M_PI;
  knot_boul=knot_boul/180.*M_PI;

  double E_boul=0, E_boul1=0;
  E_boul=anomaly_boul;
  while(abs(E_boul-E_boul1)>10e-15)
  {E_boul1=E_boul;
	  E_boul=anomaly_boul+e_boul*sin(E_boul1);}
  double cv_boul, cu_boul, sv_boul, su_boul, r_boul;
  sv_boul=(sqrt(1-e_boul*e_boul)*sin(E_boul))/(1-e_boul*cos(E_boul));
  cv_boul=(cos(E_boul)-e_boul)/(1-e_boul*cos(E_boul));
  su_boul=sv_boul*cos(arg_boul)+cv_boul*sin(arg_boul);
  cu_boul=cv_boul*cos(arg_boul)-sv_boul*sin(arg_boul);
  r_boul=axis_boul*(1-e_boul*cos(E_boul));
  X_boul=r_boul*(cu_boul*cos(knot_boul)-su_boul*sin(knot_boul)*cos(i_boul));
  Y_boul=r_boul*(cu_boul*sin(knot_boul)+su_boul*cos(knot_boul)*cos(i_boul));
  Z_boul=r_boul*su_boul*sin(i_boul);
  double p_boul;
  p_boul=axis_boul*(1-e_boul*e_boul);
  double VR_boul, VN_boul;
  VR_boul=sqrt(mu/p_boul)*e_boul*sv_boul;
  VN_boul=sqrt(mu/p_boul)*(1+e_boul*cv_boul);
  VX_boul=(X_boul/r_boul)*VR_boul+(-su_boul*cos(knot_boul)-cu_boul*sin(knot_boul)*cos(i_boul))*VN_boul;
  VY_boul=(Y_boul/r_boul)*VR_boul+(-su_boul*sin(knot_boul)+cu_boul*cos(knot_boul)*cos(i_boul))*VN_boul;
  VZ_boul=(Z_boul/r_boul)*VR_boul+cu_boul*sin(i_boul)*VN_boul;}
//переводит из Кеплеровых элементов в прямоугольные
void force(double *X, double *V, double TS,double *F, DeAcsessor^ acs)
	{
	double r[4]; int i, j;
    aster_count=1+force_var[11]*3;
    for(i=0; i<aster_count; i++)
	{r[i]=sqrt(X[i*3]*X[i*3]+X[i*3+1]*X[i*3+1]+X[i*3+2]*X[i*3+2]);}

	
	double T_=TS;
	double mass_coor[10][6];
	double massiv_aster[4][4];
	for(j=0;j<10; j++)
		{   
	        acs->GetPlanetPoz(T_, j, true, poz);
			for(i=0; i<6; i++)
	    {mass_coor[j][i]=poz[i];}}        

	double massiv_r[10];
	     for(i=0; i<10; i++)
		 {massiv_r[i]=sqrt(mass_coor[i][0]*mass_coor[i][0]+mass_coor[i][1]*mass_coor[i][1]+mass_coor[i][2]*mass_coor[i][2]);}
		 for(j=0; j<aster_count; j++)
		 {for(i=0; i<aster_count; i++)
		 {massiv_aster[j][i]=sqrt((X[j*3]-X[i*3])*(X[j*3]-X[i*3])+(X[j*3+1]-X[i*3+1])*(X[j*3+1]-X[i*3+1])+(X[j*3+2]-X[i*3+2])*(X[j*3+2]-X[i*3+2]));}}
	double massiv_vector[4][10];
	for(j=0; j<aster_count; j++)    
	{for(i=0; i<10; i++)
	{massiv_vector[j][i]=sqrt((mass_coor[i][0]-X[3*j])*(mass_coor[i][0]-X[3*j])+(mass_coor[i][1]-X[3*j+1])*(mass_coor[i][1]-X[3*j+1])+(mass_coor[i][2]-X[3*j+2])*(mass_coor[i][2]-X[3*j+2]));}}
	int k;
	for(k=0; k<aster_count; k++)
	{for (i=0; i<3; i++)
	{F[k*3+i]= (-1*mu*X[i+k*3])/(r[k]*r[k]*r[k]);
	      for(j=0; j<10; j++)
{F[k*3+i]=F[k*3+i]+force_var[j]*mu*(mass_coor[j][i]-X[k*3+i])/(massiv_vector[k][j]*massiv_vector[k][j]*massiv_vector[k][j])/massiv_mass[j];
 F[k*3+i]=F[k*3+i]-force_var[j]*mu*mass_coor[j][i]/(massiv_r[j]*massiv_r[j]*massiv_r[j])/massiv_mass[j];}
	     for(j=1;j<aster_count;j++)
			 if(j!=k)
			 {
		  F[k*3+i]=F[k*3+i]+force_var[11]*mu*(X[j*3+i]-X[k*3+i])/(massiv_aster[k][j]*massiv_aster[k][j]*massiv_aster[k][j])/massiv_mass[10+j];
		  F[k*3+i]=F[k*3+i]-(force_var[11]*mu*X[j*3+i]/(r[j]*r[j]*r[j]))/massiv_mass[10+j];}
	}
	}
	if(force_var[12]==1)
	{double planet_c[3], sin_phi;
	gelio_to_planet(X, V, mass_coor[2]);
	for(i=0; i<3; i++)
	{planet_c[i]=planet_centric[i];}
	sin_phi=planet_c[2]/massiv_vector[0][2];
	double Garmonica_aux;
	double p2=(3*sin_phi*sin_phi-1)/2;
	Garmonica_aux=Earth_sec_garmonica/pow(massiv_vector[0][2], 5.)*(3*p2+3*sin_phi*sin_phi);
	F[0]=F[0]+planet_c[0]*Garmonica_aux;
	F[1]=F[1]+planet_c[1]*Garmonica_aux;
	F[2]=F[2]+planet_c[2]*Garmonica_aux-Earth_sec_garmonica*3*sin_phi/pow(massiv_vector[0][2], 4.);}

	double Garmonica_aux2[3];
	double Garmonica_aux3[3];

	if(force_var[13]==1)
	{
	 ekvekl(Garmonica_aux2, X);
	 double omega_Sun;
	  omega_Sun= (73.667 + (TS-2396758.5)/365.2425*0.01396)/180*M_PI;
	  double cc[3];
	cc[0]=sin_i_Sun*sin(omega_Sun);
	cc[1]=-sin_i_Sun*cos(omega_Sun);
	cc[2]=cos_i_Sun;
	double z, s;
	  z=cc[0]*Garmonica_aux2[0]+cc[1]*Garmonica_aux2[1]+cc[2]*Garmonica_aux2[2];
      s=3*z-r[0]*r[0];
	  double pe[3];

	  for(int i=0; i<3; i++)
	  { pe[i]=(Sun_oblateness/(r[0]*r[0]*r[0]))*((3*cc[i]*z-Garmonica_aux2[i])*(r[0]*r[0])-2*Garmonica_aux2[i]*s)/(r[0]*r[0]*r[0]);}
	    
	  eklekv(pe, Garmonica_aux3);

	 F[0]=F[0]+Garmonica_aux3[0];
	 F[1]=F[1]+Garmonica_aux3[1];
	 F[2]=F[2]+Garmonica_aux3[2];}

	if(force_var[14]==1)
	{double rel_vel, rel_vel2;
	 rel_vel=X[0]*V[0]+X[1]*V[1]+X[2]*V[2];
	 rel_vel2=V[0]*V[0]+V[1]*V[1]+V[2]*V[2];
	 for(i=0; i<3; i++)
	 {F[i]=F[i]+relative_aux*((4*mu/r[0]-rel_vel2)*X[i]+4*rel_vel*V[i])/(r[0]*r[0]*r[0]);}}
 
	if(force_var[15]==1)
		
	{double Jup_X[3], Jup_V[3];
	 Jup_X[0]=X[0];
	 Jup_X[1]=X[1];
	 Jup_X[2]=X[2];
	 Jup_V[0]=V[0];
	 Jup_V[1]=V[1];
	 Jup_V[2]=V[2];
	 gelio_to_planet(Jup_X, Jup_V, mass_coor[4]);
	  Jup_X[0]=planet_centric[0];
	  Jup_X[1]=planet_centric[1];
	  Jup_X[2]=planet_centric[2];

	  double T=(TS-2451545)/36525;
  double aux_X1=cos(0.00015707963267948965*T);
  double aux_Y1=cos(0.00005235987755982989*T);
  double aux_X2=sin(0.00015707963267948965*T); 
  double aux_Y2=sin(0.00005235987755982989*T);

  double aux_A[3];
  aux_A[0]=aux_X2*(-0.43041922177685454*aux_Y1-0.9019874904580493*aux_Y2)+aux_X1*(-0.014654512120466917*aux_Y1-0.030710028601556798*aux_Y2);
  aux_A[1]=aux_X1*(-0.43041922177685454*aux_Y1-0.9019874904580493*aux_Y2)+aux_X2*( 0.014654512120466917*aux_Y2+0.030710028601556798*aux_Y2);
  aux_A[2]=0.9025101322420253*aux_Y1-0.430668621100356*aux_Y2;
  double  aux_Z=aux_A[0]*Jup_X[0]+aux_A[1]*Jup_X[1]+aux_A[2]*Jup_X[2];


  double Jup_to_Ast=massiv_vector[0][4]*massiv_vector[0][4];
  double Jup_to_Ast1=Jup_to_Ast*Jup_to_Ast*Jup_to_Ast;
  double useless=3*aux_Z-Jup_to_Ast;
  double Jup_F[3];
   for (int Jup_count=0; Jup_count<3; Jup_count++)
          {  Jup_F[Jup_count]=Jupiter_sec_garmonica*((3*aux_A[Jup_count]*aux_Z-Jup_X[Jup_count])*Jup_to_Ast-2*Jup_X[Jup_count]*useless)/Jup_to_Ast1;
             F[Jup_count]=F[Jup_count]+Jup_F[Jup_count];}
	}

	if(nv>12)
	{double GMA[3][3];
	 int trig;
	 double summa;
	 int ii, jj, kk, kkk, jjj;
	 for(ii=0; ii<3; ii++)
		{for(jj=0; jj<3; jj++)
		{	if(ii==jj)
				{trig=1;}
			else
				{trig=0;}
	        summa=0;
	        for(kk=0; kk<10; kk++)
			{summa+=(mu/massiv_mass[kk])/(massiv_vector[0][kk]*massiv_vector[0][kk]*massiv_vector[0][kk])*(3*(mass_coor[kk][ii]-X[ii])*(mass_coor[kk][jj]-X[jj])/(massiv_vector[0][kk]*massiv_vector[0][kk])-trig);}
			GMA[ii][jj]=mu/(r[0]*r[0]*r[0])*(3*X[ii]*X[jj]/(r[0]*r[0])-trig)+summa;}}
	 jjj=3;
	 for(kkk=3; kkk<9; kkk++)
	 {F[jjj+3*force_var[11]]=GMA[0][0]*X[kkk+3*force_var[11]]+GMA[0][1]*X[kkk+6+3*force_var[11]]+GMA[0][2]*X[kkk+12+3*force_var[11]];
	  F[jjj+6+3*force_var[11]]=GMA[1][0]*X[kkk+3*force_var[11]]+GMA[1][1]*X[kkk+6+3*force_var[11]]+GMA[1][2]*X[kkk+12+3*force_var[11]];
	  F[jjj+12+3*force_var[11]]=GMA[2][0]*X[kkk+3*force_var[11]]+GMA[2][1]*X[kkk+6+3*force_var[11]]+GMA[2][2]*X[kkk+12+3*force_var[11]];
	  jjj++;}}
	}
//правые части
int radamaker(double *X, double  *V, double TF, DeAcsessor^ acs)
{int j, k, l, m, n, j2, la, jdm;
 double *be_g_ptr, *be_ptr, *be_i_ptr, *b_g_ptr, *b_ptr, *b_i_ptr, *bt_g_ptr, *bt_ptr, *w_ptr, *u_ptr, *r_ptr, *c_ptr;
 double s, q, res, hsum, temp, val, bdouble, t, tval, t2, eps = 1.0e-10L;

while(1)
 {for (k=0, be_g_ptr=BE+NV*(KE-1), b_g_ptr=B+NV*(KE-1); k < NV; k++, be_g_ptr++, b_g_ptr++)
  {*be_g_ptr = *b_g_ptr/W[KE-1];
   for (j=0, be_ptr = BE+k, b_ptr=B+k, w_ptr=W; j < KD; j++, be_ptr+=NV, b_ptr+=NV, w_ptr++)
   {*be_ptr = *b_ptr/(*w_ptr);
    for (l=j+1, b_i_ptr=b_ptr+NV; l<KE; l++, b_i_ptr+=NV)
    {n = NW[l] + j; *be_ptr += *b_i_ptr * D[n];}}}

t = TP; tval = fabs(t); t2 = pow(t, (double) N_CLASS);

for (m=0; m<NL; m++)
  {for (j=1, j2=1; j<KF; j++)
   {la = NW[j-1] - 1;
    jdm = j - 1;
    s = H[j];
    q = pow(s, (double) (N_CLASS - 1));
    if (NPQ)
    {for (k=0, b_g_ptr=B+(KE-1)*NV; k<NV; k++, b_g_ptr++)
     {res = *b_g_ptr;
      for (l=0, b_ptr=b_g_ptr-NV; l<KD; l++, b_ptr-=NV)
      {res = *b_ptr + s*res;}
      Y[k] = X[k] + q*(t*V[k] + t2*s*(F1[k]*W1 + s*res));}}
    else
    {for (k=0, b_g_ptr=B+(KE-1)*NV; k<NV; k++, b_g_ptr++)
     {res = *b_g_ptr;
      temp = res * U[KE-1]; 
      for (l=0, b_ptr=b_g_ptr-NV, u_ptr=U+KE-2; l<KD; l++, b_ptr-=NV, u_ptr--)
      {res = *b_ptr + s*res;
       temp = *b_ptr * (*u_ptr) + s*temp;}
      Y[k] = X[k] + q*(t*V[k] + t2*s*(F1[k]*W1 + s*res));
      Z[k] = V[k] + s*t*(F1[k] + s*temp);}}
    force(Y, Z, TM + s*t, FJ, acs); NF++;
    if (j2)
    {j2 = 0;
     for (k=0, be_ptr=BE,b_ptr=B; k<NV; k++, be_ptr++, b_ptr++)
     {temp = *be_ptr;
      res = (FJ[k] - F1[k])/s;
      *be_ptr = res;
      *b_ptr += (res - temp) * W[0];}}
    else
    {for (k=0, be_ptr=BE+(j-1)*NV, b_ptr=B+(j-1)*NV; k<NV; k++, be_ptr++, b_ptr++)
     {temp = *be_ptr;
      res = (FJ[k] - F1[k])/s;
      for (l=0, r_ptr=R+la+1, be_i_ptr=BE+k; l<jdm; l++, r_ptr++, be_i_ptr+=NV)
      {res = (res - *be_i_ptr)*(*r_ptr);}
      *be_ptr = res; temp = res - temp; *b_ptr += temp * W[j-1];
      for (l=0, c_ptr=C+la+1, b_i_ptr=B+k; l<jdm; l++, c_ptr++, b_i_ptr+=NV)
      {*b_i_ptr += *c_ptr * temp;}}}}

if (m < NI-1) continue;

   hsum = 0.0L;
   val = pow(tval, (double) (-KE));
   for (k=0, b_ptr=B+(KE-1)*NV; k<NV; k++, b_ptr++)
   {bdouble = *b_ptr;
    hsum += bdouble*bdouble;}
   hsum = val*sqrt(hsum);
   if (NSF) continue; if (fabs(hsum-SM) <= 0.01L * hsum) break;
   SM = hsum;}
  if (NSF == 0)
  {if (hsum != 0.0L) TP = pow(SS/hsum, PW) * DIR;

   if (NES)
    TP = XL;
   else
   {if (TP/t <= 1.0L)
    {TP *= 0.8L; NCOUNT++;
     if (NCOUNT > 10) return 0; else return 1;}}
   NSF = 1;}

   for (k=0, b_g_ptr=B+(KE-1)*NV; k<NV; k++, b_g_ptr++)
  {res = *b_g_ptr;
   for (l=0, b_ptr=B+k; l<KD; l++, b_ptr+=NV)
   {res += *b_ptr;}
   X[k] += V[k]*t + t2*(F1[k]*W1 + res);
   if (NCL) continue;
   res = *b_g_ptr * U[KE-1];
   for (l=0, b_ptr=B+k; l<KD; l++, b_ptr+=NV)
   {res += *b_ptr * U[l];}
   V[k] += t*(F1[k] + res);}

  TM += t;
  NS++;
  if (NPER) return 0;

   force( X, V, TM, F1, acs); NF++;
  if (NES)
   TP = XL;
  else
  {if (hsum != 0.0L) TP = pow(SS/hsum, PW) * DIR;
   if (TP/t > SR) TP = SR*t;}
   if (DIR*(TM+TP) >= DIR*TF - eps)
  {TP = TF - TM; NPER = 1;}
  q = TP/t;
  for (k=0, b_g_ptr=B, bt_g_ptr=BT; k<NV; k++, b_g_ptr++, bt_g_ptr++)
  {res = 1.0L;
   for (j=0, b_ptr=b_g_ptr, bt_ptr=bt_g_ptr; j<KE; j++, b_ptr+=NV, bt_ptr+=NV)
   {if (NS>1) *bt_ptr = *b_ptr - *bt_ptr;
    if (j<KE-1)
    {m = MC[j] - 1;
     for (l=j+1, b_i_ptr=b_ptr+NV; l<KE; l++, b_i_ptr +=NV, m++)
     {*b_ptr += XI[m] * (*b_i_ptr);}}
    res *= q;
    temp = res * (*b_ptr);
    *b_ptr = temp + (*bt_ptr);
    *bt_ptr = temp;}}
    NL = NI;}
	}
void rada27 (double *X, double *V, double TI, double TF, DeAcsessor^ acs)
{
double HH[48] = {
      0.212340538239152E+00, 0.590533135559265E+00, 0.911412040487296E+00,
      0.0985350857988264E+00,0.3045357266463639E+00,0.5620251897526139E+00,
      0.8019865821263918E+00,0.9601901429485313E+00,0.0562625605369221E+00,
      0.1802406917368924E+00,0.3526247171131696E+00,0.5471536263305554E+00,
      0.7342101772154105E+00,0.8853209468390957E+00,0.9775206135612875E+00,
      0.0362578128832095E+00,0.1180789787899987E+00,0.2371769848149604E+00,
      0.3818827653047059E+00,0.5380295989189891E+00,0.6903324200723622E+00,
      0.8238833438370047E+00,0.9256126102908040E+00,0.9855875903511235E+00,
      0.0252736203975203E+00,0.0830416134474051E+00,0.1691751003771814E+00,
      0.2777967151090320E+00,0.4015027202328608E+00,0.5318623869104160E+00,
      0.6599918420853348E+00,0.7771593929561621E+00,0.8753807748555569E+00,
      0.9479645488728194E+00,0.9899817195383196E+00,0.0186103650109879E+00,
      0.0614755408992690E+00,0.1263051786933106E+00,0.2098429717265625E+00,
      0.3078989982803983E+00,0.4155560359786595E+00,0.5274156139958823E+00,
      0.6378686027177611E+00,0.7413764592942375E+00,0.8327489886084423E+00,
      0.9074047753009974E+00,0.9616018612603216E+00,0.9926353489739107E+00};

int mc[12] =     {1,13,24,34,43,51,58,64,69,73,76,78};
int nw[14] =     {0,0,1,3,6,10,15,21,28,36,45,55,66,78};
int NXI[78] =	 {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 3, 6, 10, 15, 21, 28, 36, 45,
		          55, 66, 78, 4, 10, 20, 35, 56, 84, 120, 165, 220, 286, 5, 15, 35, 70, 
				  126, 210, 330, 495, 715, 6, 21, 56, 126, 252, 462, 792, 1287, 7, 28, 84, 
				  210, 462, 924, 1716, 8, 36, 120, 330, 792, 1716, 9, 45, 165, 495, 1287, 
				  10, 55, 220, 715, 11, 66, 286, 12, 78, 13};

int kd2, la, lc, j, l, k, m, n, jd, jdm, lb, ld, le, off;
double *mem;

NV = nv; NCLASS = nclass; NOR = nor; LL = ll; XL = xl; NI = ni;

for(j=0; j<14; j++) {NW[j] = nw[j]; H[j] = 0.0L;}
for(j=0; j<13; j++) U[j] = W[j] = 0.0L;
for(j=0; j<12; j++) MC[j] = mc[j];
for(j=0; j<78; j++) XI[j] = C[j] = D[j] = R[j] = 0.0L;

delete F1;
F1 = new double [43*NV]; FJ = F1+NV; Y = FJ + NV; Z = Y + NV;
B  = Z + NV; BE = B + 13*NV; BT = BE + 13*NV;

KD  = (NOR - 3)/2; kd2 = KD/2; KE  = KD + 1;
KF  = KD + 2; PW  = 1.0L/((double) (KD + 3));
NCL = NPQ = NES = 0;
if (NCLASS == 1) NCL = 1; if (NCLASS < 2)  NPQ = 1;

if (NV == 1) SR = 1.2L; else SR = 1.5L;
if (LL < 0)	  NES = 1;
N_CLASS = abs(NCLASS);
for(n=1, la=kd2*kd2-1; n<KF; n++, la++)
{H[n] = HH[la];
 W[n-1] = 1.0L/((double) ((n+1)*(1 + (n+1)*(N_CLASS-1))));
 U[n-1] = n+2;}

W1 = 1.0L/((double) N_CLASS);

for(j=0; j<KD; j++)
 {m = MC[j] - 1;
  for(l=j+1; l<KE; l++, m++)
  {XI[m] = (double) (NXI[m]) * W[j]/W[l];}}

C[0] = -H[1] * W[0]; D[0] =  H[1] / W[1]; R[0] =  1.0L /(H[2] - H[1]);

for(k=2, la=0, lc=0; k<KE; k++)
 {lb = la; la = lc + 1; lc = NW[k+1] - 1;
  jd = lc - la - 1;
  C[la] = -H[k] * C[lb];
  C[lc] = (C[la-1]/W[jd] - H[k]) * W[jd+1];
  D[la] = H[1] * D[lb] * W[k-1]/W[k];
  D[lc] = (D[la-1]*W[k-1] + H[k]) / W[k];
  R[la] = 1.0L/(H[k+1] - H[1]);
  R[lc] = 1.0L/(H[k+1] - H[k]);
 if (k != 2)
  {for(l=3; l<=k; l++)
   {ld = la + l - 2;
    le = lb + l - 3;
    jdm= ld - la - 1;
    C[ld] = W[jdm+1] * C[le] / W[jdm] - H[k] * C[le+1];
    D[ld] = (D[le] + H[l-1] * D[le+1]) * W[k-1]/W[k];
    R[ld] = 1.0L/(H[k+1] - H[l-1]);}}}

double *v_ptr, *b_ptr, *b_g_ptr, *bt_ptr, *bt_g_ptr, *be_ptr, *be_g_ptr;

NSF = NPER = 0; TDIF = TF - TI; DIR = TDIF/fabs(TDIF);
SS = pow(10.0, -LL); NL = NI + 30;

if (NES) TP = XL = fabs(XL)*DIR;
else TP = (((double) NOR) /11.0L) * pow(0.5L,0.4L* (double) LL) * DIR/2.0L;
if (TP/TDIF > 0.5L) TP = 0.5L*TDIF;

for(k=0, v_ptr=V, bt_g_ptr=BT, b_g_ptr=B, be_g_ptr=BE; k<NV; k++, v_ptr++, bt_g_ptr++, b_g_ptr++, be_g_ptr++)
 {if (NCL) *v_ptr = 0.0L;
  for(l=0, bt_ptr=bt_g_ptr, b_ptr=b_g_ptr, be_ptr=be_g_ptr; l<KE; l++, b_ptr+=NV, bt_ptr+=NV, be_ptr+=NV)
  {*b_ptr = *bt_ptr = 0.0L;}}

NF = 0; NCOUNT = 0;

do
  {NS = 0; TM = TI; SM = 10000.0L;
   force(X, V, TM, F1, acs); NF++;} 
  while(radamaker(X, V, TF, acs));}
void abberation(double *X, double *V, double JD0, DeAcsessor^ acs)
{double Earth_coor, radius_vremia, radius_vremia_=0, time_abb;
int count=0;
 acs->GetPlanetPoz(JD0, 2, true, poz);
 radius_vremia=sqrt((X[0]-poz[0])*(X[0]-poz[0])+(X[1]-poz[1])*(X[1]-poz[1])+(X[2]-poz[2])*(X[2]-poz[2]))/cc;
 poz1[0]=poz[0];
 poz1[1]=poz[1];
 poz1[2]=poz[2];
 poz1[3]=poz[3];
 poz1[4]=poz[4];
 poz1[5]=poz[5];
 double xpred[30], vpred[30];
 int cis=0;
 for(cis=0; cis<nv; cis++)
 {xpred[cis]=X[cis];
  vpred[cis]=V[cis];}

 while(abs(radius_vremia-radius_vremia_)>1e-14 && count<20)
 {
 time_abb=JD0-radius_vremia;
  for(cis=0; cis<nv; cis++)
 {X[cis]=xpred[cis];
  V[cis]=vpred[cis];}

 rada27(X, V, JD0, time_abb, acs);
  radius_vremia_=radius_vremia;
  radius_vremia=sqrt((X[0]-poz1[0])*(X[0]-poz1[0])+(X[1]-poz1[1])*(X[1]-poz1[1])+(X[2]-poz1[2])*(X[2]-poz1[2]))/cc;
  count=count+1;}
 time_abb=JD0;}
//производит учет абберации
void flat_to_sphere(double *X)
{       rad=sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]);
		delta=asin(X[2]/rad);
		alpha=atan(X[1]/X[0]);
		if(sin(alpha)*X[1]<0)
		{alpha=alpha+M_PI;}}
//переводит прямоугольные координаты в сферические
void flat_to_Kepler(double *X, double *V)
{       double m=2.959122082855911025e-4,
			   r=sqrt(X[0]*X[0]+X[1]*X[1]+X[2]*X[2]),
			   v=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2]),
			   h=v*v/2-m/r,
			   c1=X[1]*V[2]-X[2]*V[1],
			   c2=X[2]*V[0]-X[0]*V[2],
			   c3=X[0]*V[1]-X[1]*V[0],
			   l1=-m*X[0]/r+(V[1]*c3-V[2]*c2),
			   l2=-m*X[1]/r+(V[2]*c1-V[0]*c3),
			   l3=-m*X[2]/r+(V[0]*c2-V[1]*c1),
			   c=sqrt(c1*c1+c2*c2+c3*c3),
			   l=sqrt(l1*l1+l2*l2+l3*l3);
        a_Kepler=-m/(2*h);
		e_Kepler=l/m;
		i_Kepler=acos(c3/c);
		knot_Kepler=asin(c1/(c*sin(i_Kepler)));
		arg_Kepler=asin(l3/(l*sin(i_Kepler)));
		double u=asin(X[2]/(r*i_Kepler)),
			   sine=(sqrt(1-e_Kepler*e_Kepler)*sin(u-arg_Kepler))/(1+e_Kepler*cos(u-knot_Kepler)),
			   cose=(cos(u-knot_Kepler)+e_Kepler)/(1+e_Kepler*cos(u-knot_Kepler)),
			   E=u-knot_Kepler+atan((sine*cos(u-knot_Kepler)-cose*sin(u-knot_Kepler))/(cose*cos(u-cos(u-knot_Kepler+sine*sin(u-knot_Kepler)))));
		M_Kepler=E-e_Kepler*sine;}
//переводит прямоугольные координаты в Кеплеровы элементы
void inversion(double A[6][6], int N)
{
    double temp;
    double **E = new double *[N];
    for (int i = 0; i < N; i++)
    E[i] = new double [N];
 
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {E[i][j] = 0.0; 
            if (i == j)
                E[i][j] = 1.0;}
 
    for (int k = 0; k < N; k++)
    {temp = A[k][k]; 
        for (int j = 0; j < N; j++)
        {A[k][j] /= temp;
         E[k][j] /= temp;}
 
        for (int i = k + 1; i < N; i++)
        {temp = A[i][k];
 for (int j = 0; j < N; j++)
            {A[i][j] -= A[k][j] * temp;
             E[i][j] -= E[k][j] * temp;}}}
 
    for (int k = N - 1; k > 0; k--)
    {for (int i = k - 1; i >= 0; i--)
     {temp = A[i][k]; 
            for (int j = 0; j < N; j++)
            {A[i][j] -= A[k][j] * temp;
             E[i][j] -= E[k][j] * temp;}}}
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            A[i][j] = E[i][j];
    for (int i = 0; i < N; i++)
        delete [] E[i]; 
    delete [] E;
}
 // обращение матрицы
void ephemerida(double *X, double *V, double ephem, double JD0, string num_obs, DeAcsessor^ acs)
{
double moment_of_time=0, T0=0, ST, coor_obs[3], spi, gamma, alpha1, delta1;
int cis=0;	    
for(cis=0; cis<nv; cis++)
{X_in[cis]=X[cis];
 V_in[cis]=V[cis];}
abberation(X, V, ephem, acs);
for(cis=0; cis<nv; cis++)
{X_abber[cis]=X[cis];
 V_abber[cis]=V[cis];}
 gelio_to_planet(X, V, poz1);
 flat_to_sphere(planet_centric);
 moment_of_time=JD0;
 T0=floor(moment_of_time+0.5)-0.5;
 ST=star_time(T0)+mu1*(moment_of_time-T0)*86400;
 while(ST<0)
 {ST=ST+86400;}
 while(ST>86400)
 {ST=ST-86400;}

 str_spisok=read_observatory(num_obs);

 
coor_obs[0]=names_lon(str_spisok);
coor_obs[0]=coor_obs[0]/180*M_PI;
coor_obs[1]=names_alt(str_spisok);
coor_obs[2]=names_hren(str_spisok);

 ST=ST/86400*2*M_PI+coor_obs[0];
 spi=(Radius_Earth/AE)/rad;
 if(num_obs!="500")
 {alpha1=alpha-atan(coor_obs[1]*spi*sin(ST-alpha)/(cos(delta)-coor_obs[1]*spi*cos(ST-alpha)));
 gamma=atan((coor_obs[2]/coor_obs[1])*cos((alpha-alpha1)/2)*1/cos((ST-(alpha+alpha1)/2)));
 delta1=delta+atan(coor_obs[2]*spi/sin(gamma)*sin(delta-gamma)/(1-coor_obs[2]*spi/sin(gamma)*cos(delta-gamma)));}
 else
 {alpha1=alpha;
  delta1=delta;}
 alpha1=alpha1/M_PI*180;
 delta1=delta1/M_PI*180;
 alpha1=alpha1/15.;
 if(alpha1<0)
 {alpha1=alpha1+24;}
 if(alpha1>24)
 {alpha1=alpha1-24;}
 h_alpha_eph=int(alpha1);
 m_alpha_eph=int((alpha1-h_alpha_eph)*60);
 s_alpha_eph=((alpha1-h_alpha_eph)*60-m_alpha_eph)*60;
 d_delta_eph=int(delta1);
 m_delta_eph=int((delta1-d_delta_eph)*60);
 s_delta_eph=((delta1-d_delta_eph)*60-m_delta_eph)*60;}
//расчет эфемерид
void ephemerida(double *X, double *V, double ephem, double JD0, string num_obs, double s_n[],double s_x[],double s_y[],double s_z[], int number_obs, DeAcsessor^ acs)
{
	
double moment_of_time=0, T0=0, ST, coor_obs[3], spi, gamma, alpha1, delta1;
int cis=0;
for(cis=0; cis<nv; cis++)
{X_in[cis]=X[cis];
 V_in[cis]=V[cis];}
 abberation(X, V, ephem, acs);
for(cis=0; cis<nv; cis++)
{X_abber[cis]=X[cis];
 V_abber[cis]=V[cis];}
 gelio_to_planet(X, V, poz1);

  
 if(number_obs==s_n[s_c])
 {planet_centric[0]=planet_centric[0]-s_x[s_c];
  planet_centric[1]=planet_centric[1]-s_y[s_c];
  planet_centric[2]=planet_centric[2]-s_z[s_c];}


 flat_to_sphere(planet_centric);
 moment_of_time=JD0;
 T0=floor(moment_of_time+0.5)-0.5;
 ST=star_time(T0)+mu1*(moment_of_time-T0)*86400;
 while(ST<0)
 {ST=ST+86400;}
 while(ST>86400)
 {ST=ST-86400;}


 if(number_obs!=s_n[s_c])
 {str_spisok=read_observatory(num_obs);
 coor_obs[0]=names_lon(str_spisok);
 coor_obs[0]=coor_obs[0]/180*M_PI;
 coor_obs[1]=names_alt(str_spisok);
 coor_obs[2]=names_hren(str_spisok);
 ST=ST/86400*2*M_PI+coor_obs[0];
 spi=(Radius_Earth/AE)/rad;
 if(num_obs!="500")
 {alpha1=alpha-atan(coor_obs[1]*spi*sin(ST-alpha)/(cos(delta)-coor_obs[1]*spi*cos(ST-alpha)));
 gamma=atan((coor_obs[2]/coor_obs[1])*cos((alpha-alpha1)/2)*1/cos((ST-(alpha+alpha1)/2)));
 delta1=delta+atan(coor_obs[2]*spi/sin(gamma)*sin(delta-gamma)/(1-coor_obs[2]*spi/sin(gamma)*cos(delta-gamma)));}
 else
 {alpha1=alpha;
  delta1=delta;}
 }
 else
 {s_c++;
 alpha1=alpha;
 delta1=delta;}

 alpha1=alpha1/M_PI*180;
 delta1=delta1/M_PI*180;
 alpha1=alpha1/15.;
 if(alpha1<0)
 {alpha1=alpha1+24;}
 if(alpha1>24)
 {alpha1=alpha1-24;}
 h_alpha_eph=int(alpha1);
 m_alpha_eph=int((alpha1-h_alpha_eph)*60);
 s_alpha_eph=((alpha1-h_alpha_eph)*60-m_alpha_eph)*60;
 d_delta_eph=int(delta1);
 m_delta_eph=int((delta1-d_delta_eph)*60);
 s_delta_eph=((delta1-d_delta_eph)*60-m_delta_eph)*60;
 }
//расчет эфемерид для улучшения
double average_moment (double data_moment [], int number)
{double target;
 double sum=0;
 for(int i=0; i<number; i++)
 {sum+=data_moment[i];}
 target=sum/number;
 return target;
}
//вычисление среднего момента времени
void improve_elements(double in_X[], double in_V[], double time, double data_ul[], double del[], double al[], string obs[], int count, double s_n[],double s_x[],double s_y[],double s_z[], DeAcsessor^ acs)
{double epher[5000];
 double aaa[30000][6];
 double aa1[10000][6];
 double bb[10000];
 int i=0;
 double vector_error[6];
 double in_time;
 double aa2[6][6];
 double demph=1;
 double compare_sigma=1000000;

 ofstream fout12, fout13, fout15, fout16;
 fout12.open("result-residual.txt");
 fout13.open("result-epher.txt");
 fout15.open("covar.txt");
 fout16.open("search.txt");
 fout13.precision(12);
 fout16.precision(16);

 for(i=0; i<count; i++)
 {epher[i]=eph_time(data_ul[i]);
  fout13<<epher[i]<<endl;}
 in_time=average_moment(epher, count);
 in_time=int(in_time);
 in_time=in_time+0.5;
  X[0]=in_X[0];
  X[1]=in_X[1];
  X[2]=in_X[2];
  V[0]=in_V[0];
  V[1]=in_V[1];
  V[2]=in_V[2];
 nv=3+9*force_var[11];

  rada27(X, V, time, in_time, acs);

  in_X[0]=X[0];
  in_X[1]=X[1];
  in_X[2]=X[2];
  in_V[0]=V[0];
  in_V[1]=V[1];
  in_V[2]=V[2];
  time=in_time;

 vector_error[0]=100;
 vector_error[1]=100;
 vector_error[2]=100;

	 count_iter=0;

 while(sqrt(vector_error[0]*vector_error[0]+vector_error[1]*vector_error[1]+vector_error[2]*vector_error[2])>1e-10&&count_iter<40)
 {   s_c=0;
	 count_iter++;
	 fout12<<"Number of iteration "<<count_iter<<endl<<endl;
	 cout<<count_iter<<endl;
	 time=in_time;
	 X[0]=in_X[0];
     X[1]=in_X[1];
	 X[2]=in_X[2];
	 V[0]=in_V[0];
	 V[1]=in_V[1];
	 V[2]=in_V[2];
	  nv=3+9*force_var[11]+18;

	 if(force_var[11]==0)
	{X[3]=1;  X[4]=0;   X[5]=0;  X[6]=0;  X[7]=0; X[8]=0; 
     X[9]=0;  X[10]=1;  X[11]=0; X[12]=0; X[13]=0; X[14]=0;
	 X[15]=0; X[16]=0;  X[17]=1; X[18]=0; X[19]=0; X[20]=0;

	 V[3]=0;  V[4]=0;   V[5]=0;  V[6]=1;  V[7]=0; V[8]=0; 
     V[9]=0;  V[10]=0;  V[11]=0; V[12]=0; V[13]=1; V[14]=0;
	 V[15]=0; V[16]=0;  V[17]=0; V[18]=0; V[19]=0; V[20]=1;}

 else
      {X[12]=1;  X[13]=0;  X[14]=0;  X[15]=0; X[16]=0; X[17]=0; 
       X[18]=0;  X[19]=1;  X[20]=0;  X[21]=0; X[22]=0; X[23]=0;
	   X[24]=0;  X[25]=0;  X[26]=1;  X[27]=0; X[28]=0; X[29]=0;

	   V[12]=0;  V[13]=0;  V[14]=0;  V[15]=1; V[16]=0; V[17]=0; 
       V[18]=0;  V[19]=0;  V[20]=0;  V[21]=0; V[22]=1; V[23]=0;
	   V[24]=0;  V[25]=0;  V[26]=0;  V[27]=0; V[28]=0; V[29]=1;}

  double summa=0;
  double X_aux[30], V_aux[30];
  int sat_count=0;

  int iii, jjj, kkk;

  for(int opr=0; opr<30; opr++)
  {X_aux[opr]=X[opr];
  V_aux[opr]=V[opr];}

  for(kkk=0; kkk<count; kkk++)
  {
  for(int opr=0; opr<30; opr++)
  {X[opr]=X_aux[opr];
  V[opr]=V_aux[opr];}
  rada27(X, V, time, epher[kkk], acs);
  for(int opr=0; opr<30; opr++)
  {X_aux[opr]=X[opr];
  V_aux[opr]=V[opr];}
  time=epher[kkk];
  
 ephemerida(X, V, epher[kkk], data_ul[kkk], obs[kkk], s_n, s_x, s_y, s_z, kkk, acs);
 double alpha_im=((h_alpha_eph+m_alpha_eph/60.+s_alpha_eph/3600.)*15)/180*M_PI;
 double delta_im=(d_delta_eph+m_delta_eph/60.+s_delta_eph/3600.)/180*M_PI;

 for(iii=0; iii<3; iii++)
 {for(jjj=0; jjj<6; jjj++)
    {aaa[iii+6*kkk][jjj]=X_abber[3+6*iii+jjj+3*force_var[11]];}}
 for(iii=0; iii<3; iii++)
 {for(jjj=0; jjj<6; jjj++)
    {aaa[3+iii+6*kkk][jjj]=V_abber[3+6*iii+jjj+3*force_var[11]];}}
 for(iii=0; iii<6; iii++)
 {aa1[2*kkk][iii]=(cos(alpha_im)*aaa[1+6*kkk][iii]-sin(alpha_im)*aaa[6*kkk][iii])/rad;
 aa1[2*kkk+1][iii]=(cos(delta_im)*aaa[2+6*kkk][iii]-sin(delta_im)*(cos(alpha_im)*aaa[6*kkk][iii]+sin(alpha_im)*aaa[1+6*kkk][iii]))/rad;}

 bb[2*kkk]=(al[kkk]-alpha_im)*cos(del[kkk]);
 bb[1+2*kkk]=(del[kkk]-delta_im);
 fout12<<kkk+1<<" "<<bb[2*kkk]*206265<<" "<<bb[1+2*kkk]*206265<<endl;
 cout<<kkk+1<<" "<<bb[2*kkk]*206265<<" "<<bb[1+2*kkk]*206265<<endl;
 summa+=(bb[2*kkk]*180/M_PI)*3600*(bb[2*kkk]*180/M_PI)*3600;
 summa+=(bb[1+2*kkk]*180/M_PI)*3600*(bb[1+2*kkk]*180/M_PI)*3600;
 fout12<<summa<<endl;}
 

 sigma=sqrt(summa/(count*2-6));

 if(sigma>100&&sigma>compare_sigma)
 {demph=demph/10.;}

 compare_sigma=sigma;
 fout12<<endl<<"sigma = "<<sigma<<endl<<endl;
 cout<<sigma<<endl;
 int N=6;

 for(iii=0; iii<6; iii++)
 {for(jjj=0; jjj<6; jjj++)
 {aa2[iii][jjj]=0;
  for(kkk=0; kkk<count*2; kkk++)
  {aa2[iii][jjj]+=aa1[kkk][iii]*aa1[kkk][jjj];}}}
  inversion(aa2, 6);
  double aa4[6][10000];
  for(iii=0; iii<6; iii++)
  {for(jjj=0; jjj<count*2; jjj++)
  {aa4[iii][jjj]=0;
   for(kkk=0; kkk<6; kkk++)
   {aa4[iii][jjj]+=aa2[iii][kkk]*aa1[jjj][kkk];}}}
  for(iii=0; iii<6; iii++)
  {vector_error[iii]=0;
   for(jjj=0; jjj<count*2; jjj++)
   {vector_error[iii]+=aa4[iii][jjj]*bb[jjj];}
    vector_error[iii]=vector_error[iii]*demph;}

in_X[0]=in_X[0]+vector_error[0];
in_X[1]=in_X[1]+vector_error[1];
in_X[2]=in_X[2]+vector_error[2];
in_V[0]=in_V[0]+vector_error[3];
in_V[1]=in_V[1]+vector_error[4];
in_V[2]=in_V[2]+vector_error[5];
fout12.precision(15);
fout12<<"Improved coordinate"<<endl<<endl;
fout12<<in_X[0]<<endl;
fout12<<in_X[1]<<endl;
fout12<<in_X[2]<<endl;
fout12<<in_V[0]<<endl;
fout12<<in_V[1]<<endl;
fout12<<in_X[2]<<endl<<endl;
fout12.precision(6);
double ekv_aux[6];
ekv_aux[0]=in_X[0];
ekv_aux[1]=in_X[1];
ekv_aux[2]=in_X[2];
ekv_aux[3]=in_V[0];
ekv_aux[4]=in_V[1];
ekv_aux[5]=in_V[2];
double ekl_aux[6];
ekvekl(ekl_aux, ekv_aux);
double iin_X[3], iin_V[3];
iin_X[0]=ekl_aux[0];
iin_X[1]=ekl_aux[1];
iin_X[2]=ekl_aux[2];
iin_V[0]=ekl_aux[3];
iin_V[1]=ekl_aux[4];
iin_V[2]=ekl_aux[5];
flat_to_Kepler(iin_X, iin_V);
fout12<<"Elements"<<endl<<endl;
fout12<<a_Kepler<<endl;
fout12<<e_Kepler<<endl;
fout12<<i_Kepler*180/M_PI<<endl;
fout12<<arg_Kepler*180/M_PI<<endl;
fout12<<knot_Kepler*180/M_PI<<endl;
fout12<<M_Kepler*180/M_PI<<endl<<endl; 
fout12<<"Damping factors"<<endl<<endl;
fout12<<demph<<endl<<endl;}

fout15.precision(10);
for(i=0; i<6; i++)
{for(int j=0; j<6; j++)
{aa2[i][j]=aa2[i][j]*(sigma*M_PI/(180.*3600.))*(sigma*M_PI/(180.*3600.));}
fout15<<aa2[i][0]<<" "<<aa2[i][1]<<" "<<aa2[i][2]<<" "<<aa2[i][3]<<" "<<aa2[i][4]<<" "<<aa2[i][5]<<endl;
}

 fout16<<in_time<<endl<<in_X[0]<<" "<<in_X[1]<<" "<<in_X[2]<<" "<<in_V[0]<<" "<<in_V[1]<<" "<<in_V[2]<<" "<<endl;
 for(i=0; i<6; i++)
{fout16<<aa2[i][0]<<" "<<aa2[i][1]<<" "<<aa2[i][2]<<" "<<aa2[i][3]<<" "<<aa2[i][4]<<" "<<aa2[i][5]<<endl;}

 fout12.close(); 
 fout13.close();
 fout15.close();
 fout16.close();
}
//улучшение элементов орбит
void improve_elements(double in_X[], double in_V[], double time, double data_ul[], double del[], double al[], string obs[], int count, double s_n[],double s_x[],double s_y[],double s_z[], double sigma_im, DeAcsessor^ acs)
{
	double epher[5000];
	double aaa[30000][6];
	double aa1[10000][6];
	double bb[10000];

	int i=0;
	double vector_error[6];
	double in_time;
	double aa2[6][6];
	double demph=1;
	double compare_sigma=1000000;
	int mass_rejection[5000];
	for(i=0; i<5000; i++)
	{
		mass_rejection[i]=1;
	}

	ofstream fout12, fout13, fout15, fout16;
	fout12.open(path_out + "result-residual.txt");
	fout13.open(path_out + "result-epher.txt");
	fout15.open(path_out + "covar.txt");
	fout16.open(path_out + "search.txt");
	fout13.precision(12);
	fout16.precision(16);

	for(i=0; i<count; i++)
	{
		epher[i]=eph_time(data_ul[i]);
		fout13<<epher[i]<<endl;
	}

	nv=3+9*force_var[11];
	vector_error[0]=100;
	vector_error[1]=100;
	vector_error[2]=100;

	count_iter=0;
	in_time=time;

	while(sqrt(vector_error[0]*vector_error[0]+vector_error[1]*vector_error[1]+vector_error[2]*vector_error[2])>1e-10&&count_iter<40)
	{   
		s_c=0;
		count_iter++;
		fout12<<endl<<"Number of iteration "<<count_iter<<endl<<endl;
		cout<<count_iter<<endl;
		time=in_time;

		X[0]=in_X[0];
		X[1]=in_X[1];
		X[2]=in_X[2];

		V[0]=in_V[0];
		V[1]=in_V[1];
		V[2]=in_V[2];

		nv=3+9*force_var[11]+18;

		if(force_var[11]==0)
		{
			X[3]=1;  X[4]=0;   X[5]=0;  X[6]=0;  X[7]=0; X[8]=0; 
			X[9]=0;  X[10]=1;  X[11]=0; X[12]=0; X[13]=0; X[14]=0;
			X[15]=0; X[16]=0;  X[17]=1; X[18]=0; X[19]=0; X[20]=0;

			V[3]=0;  V[4]=0;   V[5]=0;  V[6]=1;  V[7]=0; V[8]=0; 
			V[9]=0;  V[10]=0;  V[11]=0; V[12]=0; V[13]=1; V[14]=0;
			V[15]=0; V[16]=0;  V[17]=0; V[18]=0; V[19]=0; V[20]=1;
		}
		else
		{
			X[12]=1;  X[13]=0;  X[14]=0;  X[15]=0; X[16]=0; X[17]=0; 
			X[18]=0;  X[19]=1;  X[20]=0;  X[21]=0; X[22]=0; X[23]=0;
			X[24]=0;  X[25]=0;  X[26]=1;  X[27]=0; X[28]=0; X[29]=0;

			V[12]=0;  V[13]=0;  V[14]=0;  V[15]=1; V[16]=0; V[17]=0; 
			V[18]=0;  V[19]=0;  V[20]=0;  V[21]=0; V[22]=1; V[23]=0;
			V[24]=0;  V[25]=0;  V[26]=0;  V[27]=0; V[28]=0; V[29]=1;
		}

		double summa=0;
		double X_aux[30], V_aux[30];
		int sat_count=0;

		int iii, jjj, kkk;

		for(int opr=0; opr<30; opr++)
		{
			X_aux[opr]=X[opr];
			V_aux[opr]=V[opr];
		}
  
		for(kkk=0; kkk<count; kkk++)
		{   
			if(mass_rejection[kkk]==0&&kkk==s_n[s_c])
			{
				s_c++;
			}

			if(mass_rejection[kkk]!=0)
			{
				for(int opr=0; opr<30; opr++)
				{
					X[opr]=X_aux[opr];
					V[opr]=V_aux[opr];
				}

				rada27(X, V, time, epher[kkk], acs);

				for(int opr=0; opr<30; opr++)
				{
					X_aux[opr]=X[opr];
					V_aux[opr]=V[opr];
				}

				time=epher[kkk];

				ephemerida(X, V, epher[kkk], data_ul[kkk], obs[kkk], s_n, s_x, s_y, s_z, kkk, acs);

				double alpha_im=((h_alpha_eph+m_alpha_eph/60.+s_alpha_eph/3600.)*15)/180*M_PI;
				double delta_im=(d_delta_eph+m_delta_eph/60.+s_delta_eph/3600.)/180*M_PI;

				for(iii=0; iii<3; iii++)
				{
					for(jjj=0; jjj<6; jjj++)
					{
						aaa[iii+6*kkk][jjj]=X_abber[3+6*iii+jjj+3*force_var[11]];
					}
				}

				for(iii=0; iii<3; iii++)
				{
					for(jjj=0; jjj<6; jjj++)
					{
						aaa[3+iii+6*kkk][jjj]=V_abber[3+6*iii+jjj+3*force_var[11]];
					}
				}

				for(iii=0; iii<6; iii++)
				{
					aa1[2*kkk][iii]=(cos(alpha_im)*aaa[1+6*kkk][iii]-sin(alpha_im)*aaa[6*kkk][iii])/rad;
					aa1[2*kkk+1][iii]=(cos(delta_im)*aaa[2+6*kkk][iii]-sin(delta_im)*(cos(alpha_im)*aaa[6*kkk][iii]+sin(alpha_im)*aaa[1+6*kkk][iii]))/rad;
				}

				bb[2*kkk]=(al[kkk]-alpha_im)*cos(del[kkk]);
				bb[1+2*kkk]=(del[kkk]-delta_im);

				if(abs((bb[2*kkk]/M_PI)*180*3600)>3*sigma_im||abs((bb[1+2*kkk]/M_PI)*180*3600)>3*sigma)
				{
					mass_rejection[kkk]=0;
				}

				fout12<<kkk+1<<" "<<bb[2*kkk]*206265<<" "<<bb[1+2*kkk]*206265<<endl;
				cout<<kkk+1<<" "<<bb[2*kkk]*206265<<" "<<bb[1+2*kkk]*206265<<endl;
				summa+=(bb[2*kkk]*180/M_PI)*3600*(bb[2*kkk]*180/M_PI)*3600;
				summa+=(bb[1+2*kkk]*180/M_PI)*3600*(bb[1+2*kkk]*180/M_PI)*3600;
			}
		}

		double num_rejection=0;
		for(i=0; i<count; i++)
		{
			if(mass_rejection[i]==0)
			{
				num_rejection++;
			}
		}

		sigma=sqrt(summa/((count-num_rejection)*2-6));

		if(sigma>100&&sigma>compare_sigma)
		{
			demph=demph/10.;
		}

		compare_sigma=sigma;
		fout12<<endl<<"sigma = "<<sigma<<endl<<endl;
		cout<<sigma<<endl;

		int N=6;

		for(iii=0; iii<6; iii++)
		{
			for(jjj=0; jjj<6; jjj++)
			{
				aa2[iii][jjj]=0;
				for(kkk=0; kkk<count*2; kkk++)
				{
					aa2[iii][jjj]+=aa1[kkk][iii]*aa1[kkk][jjj]*mass_rejection[int(kkk/2)];
				}
			}
		}

		inversion(aa2, 6);
		double aa4[6][10000];
		for(iii=0; iii<6; iii++)
		{
			for(jjj=0; jjj<count*2; jjj++)
			{
				aa4[iii][jjj]=0;
				for(kkk=0; kkk<6; kkk++)
				{
					aa4[iii][jjj]+=aa2[iii][kkk]*aa1[jjj][kkk]*mass_rejection[int(jjj/2)];
				}
			}
		}

		for(iii=0; iii<6; iii++)
		{
			vector_error[iii]=0;
			for(jjj=0; jjj<count*2; jjj++)
			{
				vector_error[iii]+=aa4[iii][jjj]*bb[jjj]*mass_rejection[int(jjj/2)];
			}

			vector_error[iii]=vector_error[iii]*demph;
		}

		in_X[0]=in_X[0]+vector_error[0];
		in_X[1]=in_X[1]+vector_error[1];
		in_X[2]=in_X[2]+vector_error[2];
		
		in_V[0]=in_V[0]+vector_error[3];
		in_V[1]=in_V[1]+vector_error[4];
		in_V[2]=in_V[2]+vector_error[5];

		fout12.precision(15);
		fout12<<"Improved coordinate"<<endl<<endl;
		fout12<<in_X[0]<<endl;
		fout12<<in_X[1]<<endl;
		fout12<<in_X[2]<<endl;
		fout12<<in_V[0]<<endl;
		fout12<<in_V[1]<<endl;
		fout12<<in_X[2]<<endl<<endl;
		fout12.precision(6);
		double ekv_aux[6];
		ekv_aux[0]=in_X[0];
		ekv_aux[1]=in_X[1];
		ekv_aux[2]=in_X[2];
		ekv_aux[3]=in_V[0];
		ekv_aux[4]=in_V[1];
		ekv_aux[5]=in_V[2];
		double ekl_aux[6];
		ekvekl(ekl_aux, ekv_aux);
		double iin_X[3], iin_V[3];
		iin_X[0]=ekl_aux[0];
		iin_X[1]=ekl_aux[1];
		iin_X[2]=ekl_aux[2];
		iin_V[0]=ekl_aux[3];
		iin_V[1]=ekl_aux[4];
		iin_V[2]=ekl_aux[5];
		flat_to_Kepler(iin_X, iin_V);
		fout12<<"Elements"<<endl<<endl;
		fout12<<a_Kepler<<endl;
		fout12<<e_Kepler<<endl;
		fout12<<i_Kepler*180/M_PI<<endl;
		fout12<<arg_Kepler*180/M_PI<<endl;
		fout12<<knot_Kepler*180/M_PI<<endl;
		fout12<<M_Kepler*180/M_PI<<endl<<endl; 
		fout12<<"Damping factors"<<endl<<endl;
		fout12<<demph<<endl<<endl;
		fout12<<"Rejected observations"<<endl<<endl;

		for(i=0; i<count; i++)
		{
			if(mass_rejection[i]==0)
			{
				fout12<<i<<" ";
			}
		}
	}
				
	fout15.precision(10);
	for(i=0; i<6; i++)
	{
		for(int j=0; j<6; j++)
		{
			aa2[i][j]=aa2[i][j]*(sigma*M_PI/(180.*3600.))*(sigma*M_PI/(180.*3600.));
		}

		fout15<<aa2[i][0]<<" "<<aa2[i][1]<<" "<<aa2[i][2]<<" "<<aa2[i][3]<<" "<<aa2[i][4]<<" "<<aa2[i][5]<<endl;
	}
		
	fout16<<in_time<<endl<<in_X[0]<<" "<<in_X[1]<<" "<<in_X[2]<<" "<<in_V[0]<<" "<<in_V[1]<<" "<<in_V[2]<<" "<<endl;

	for(i=0; i<6; i++)
	{
		fout16<<aa2[i][0]<<" "<<aa2[i][1]<<" "<<aa2[i][2]<<" "<<aa2[i][3]<<" "<<aa2[i][4]<<" "<<aa2[i][5]<<endl;
	}

	fout12.close(); 
	fout13.close();
	fout15.close();
	fout16.close();
}
//улучшение элементов орбит c отбраковкой
double magnitude_observ(double XM[3], double abs_magnitude, double GM, double epher_time_mag, DeAcsessor^ acs)
{
	double magnitude_observ;
	double mass_Earth[3], planet_centric[3], mass_Sun[3];
	
	acs->GetPlanetPoz(epher_time_mag, 2, true, mass_Earth);
	for(int i=0; i<3; i++)
	{
		planet_centric[i]=XM[i]-mass_Earth[i];
	}
	for(int i=0; i<3; i++)
	{
		mass_Sun[i]=mass_Earth[i]*(-1);
	}

	double distance_to_Earth, distance_Earth_to_Sun, distance_to_Sun;
	distance_to_Earth=sqrt(planet_centric[0]*planet_centric[0]+planet_centric[1]*planet_centric[1]+planet_centric[2]*planet_centric[2]);
	distance_to_Sun=sqrt(XM[0]*XM[0]+XM[1]*XM[1]+XM[2]*XM[2]);
	distance_Earth_to_Sun=sqrt(mass_Earth[0]*mass_Earth[0]+mass_Earth[1]*mass_Earth[1]+mass_Earth[2]*mass_Earth[2]);

	double var_p;
	var_p=(distance_Earth_to_Sun+distance_to_Earth+distance_to_Sun)/2;
	
	double var_SS;
	var_SS=sqrt(var_p*(var_p-distance_to_Sun)*(var_p-distance_Earth_to_Sun)*(var_p-distance_to_Earth));
	
	double var_tb;
	var_tb=var_SS/var_p/(var_p-distance_Earth_to_Sun);
	
	double var_betta;
	var_betta=2*atan(var_tb)*180/M_PI;
	
	double var_f[2];
	const double varmas_a[2]={3.33, 1.87}, varmas_b[2]={0.63, 1.22};
	int i=0;
	for(i=0; i<2; i++)
	{
		var_f[i]=exp(-varmas_a[i]*pow(var_tb, varmas_b[i]));
	}
	
	double var_a1, var_b1;
	var_a1=distance_to_Sun*distance_to_Earth;
	var_b1=(1-GM)*var_f[0]+GM*var_f[1];

	if(var_a1!=0 && var_b1!=0)
	{
		magnitude_observ=abs_magnitude+5.*log10(var_a1)-2.5*log10(var_b1);
	}

	return magnitude_observ;
}
//вычисление видимой звездной величины
void domain (double magnitude_d, double G_d, double domain_time, string num_obs, DeAcsessor^ acs)
{
	double covar_d[6][6], improve_time, coor_d[30], vel_d[30];
	ifstream fin;

	int i, j;
	for(i=0; i<30; i++)
	{
		coor_d[i]=0;
	}

	for(i=0; i<30; i++)
	{
		vel_d[i]=0;
	}

	fin.open(path_out + "search.txt");
	fin>>improve_time>>coor_d[0]>>coor_d[1]>>coor_d[2]>>vel_d[0]>>vel_d[1]>>vel_d[2];
	for(i=0; i<6; i++)
	{
		for(j=0; j<6; j++)
		{
			fin>>covar_d[i][j];
		}
	}
	fin.close();

	if(force_var[11]==0)
	{
		coor_d[3]=1;  coor_d[4]=0;   coor_d[5]=0;  coor_d[6]=0;  coor_d[7]=0; coor_d[8]=0; 
		coor_d[9]=0;  coor_d[10]=1;  coor_d[11]=0; coor_d[12]=0; coor_d[13]=0; coor_d[14]=0;
		coor_d[15]=0; coor_d[16]=0;  coor_d[17]=1; coor_d[18]=0; coor_d[19]=0; coor_d[20]=0;

		vel_d[3]=0;  vel_d[4]=0;   vel_d[5]=0;  vel_d[6]=1;  vel_d[7]=0; vel_d[8]=0; 
		vel_d[9]=0;  vel_d[10]=0;  vel_d[11]=0; vel_d[12]=0; vel_d[13]=1; vel_d[14]=0;
		vel_d[15]=0; vel_d[16]=0;  vel_d[17]=0; vel_d[18]=0; vel_d[19]=0; vel_d[20]=1;
	}
	else
	{
		coor_d[12]=1;  coor_d[13]=0;  coor_d[14]=0;  coor_d[15]=0; coor_d[16]=0; coor_d[17]=0; 
		coor_d[18]=0;  coor_d[19]=1;  coor_d[20]=0;  coor_d[21]=0; coor_d[22]=0; coor_d[23]=0;
		coor_d[24]=0;  coor_d[25]=0;  coor_d[26]=1;  coor_d[27]=0; coor_d[28]=0; coor_d[29]=0;

		vel_d[12]=0;  vel_d[13]=0;  vel_d[14]=0;  vel_d[15]=1; vel_d[16]=0; vel_d[17]=0; 
		vel_d[18]=0;  vel_d[19]=0;  vel_d[20]=0;  vel_d[21]=0; vel_d[22]=1; vel_d[23]=0;
		vel_d[24]=0;  vel_d[25]=0;  vel_d[26]=0;  vel_d[27]=0; vel_d[28]=0; vel_d[29]=1;
	}

	double ephertime;
	ephertime=eph_time(domain_time);

	rada27(coor_d, vel_d, improve_time, domain_time, acs); 

	ephemerida(coor_d, vel_d, ephertime, domain_time, num_obs, acs);
	magnitude_vis=magnitude_observ(coor_d, magnitude_d, G_d, ephertime, acs);
	double DA[6], DD[6];
	double alpha_rad_d=(h_alpha_eph+m_alpha_eph/60.+s_alpha_eph/3600.)*M_PI/12.;
	double delta_rad_d=(d_delta_eph+m_delta_eph/60.+s_alpha_eph/3600.)*M_PI/180.;
	if(force_var[11]==0)
	{
		for(int count=0; count<6; count++)
		{
			DA[count]=-coor_d[9+count]*sin(alpha_rad_d)+coor_d[15+count]*cos(alpha_rad_d);
			DD[count]=-coor_d[9+count]*sin(delta_rad_d)*cos(alpha_rad_d)-coor_d[15+count]*sin(alpha_rad_d)*sin(delta_rad_d)+vel_d[3+count]*cos(delta_rad_d);
			DA[count]=DA[count]/(rad*cos(delta_rad_d));
			DD[count]=DD[count]/rad;
		}
		double C11=0, C12=0, C22=0;

		for(i=0; i<6; i++)
		{
			for(j=0; j<6; j++)
			{
				C11=C11+DA[i]*covar_d[i][j]*DA[j];
				C12=C12+DA[i]*covar_d[i][j]*DD[j];
				C22=C22+DD[i]*covar_d[i][j]*DD[j];
			}
		}

		alpha_dom=sqrt(C11)*(12*60*60)/M_PI;
		delta_dom=sqrt(C22)*(180*60*60)/M_PI;
	}   
	else
	{
		for(int count=0; count<6; count++)
		{
			DA[count]=-coor_d[8+count+9]*sin(alpha_rad_d)+coor_d[14+count+9]*cos(alpha_rad_d);
			DD[count]=-coor_d[8+count+9]*sin(delta_rad_d)*cos(alpha_rad_d)-coor_d[14+count+9]*sin(alpha_rad_d)*sin(delta_rad_d)+vel_d[3+count+9]*cos(delta_rad_d);
			DA[count]=DA[count]/(rad*cos(delta_rad_d));
			DD[count]=DD[count]/rad;
		}

		double C11=0, C12=0, C22=0;
		for(i=0; i<6; i++)
		{
			for(j=0; i<6; j++)
			{
				C11=C11+DA[i]*covar_d[i][j]*DA[j];
				C12=C12+DA[i]*covar_d[i][j]*DD[j];
				C22=C22+DD[i]*covar_d[i][j]*DD[j];
			}
		}

		alpha_dom=sqrt(C11)*206265/15.;
		delta_dom=sqrt(C22)*206265;
	}
}
//вычисление доверительной области линейным методом
double to_rad(double d)
{
	return d*M_PI / 180;
}
//в радианы
double metric(Asteroid ast1, Asteroid ast2)
{
	double m, cosi, p1, p2, cosp;
	cosi = cos(to_rad(ast1.i))*cos(to_rad(ast2.i)) + sin(to_rad(ast1.i))*sin(to_rad(ast2.i))*cos(to_rad(ast1.node - ast2.node));
	p1 = (cos(to_rad(ast1.w))*cos(to_rad(ast2.w)) + cos(to_rad(ast1.i))*cos(to_rad(ast2.i))*sin(to_rad(ast1.w))*sin(to_rad(ast2.w)))*cos(to_rad(ast1.node - ast2.node));
	p2 = (cos(to_rad(ast2.i))*cos(to_rad(ast1.w))*sin(to_rad(ast2.w)) - cos(to_rad(ast1.i))*sin(to_rad(ast1.w))*cos(to_rad(ast2.w)))*sin(to_rad(ast1.node - ast2.node));
	cosp = sin(to_rad(ast1.i))*sin(to_rad(ast2.i))*sin(to_rad(ast1.w))*sin(to_rad(ast2.w)) + p1 + p2;
	m = (1 + pow(ast1.e, 2))*ast1.p + (1 + pow(ast2.e, 2))*ast2.p - 2 * sqrt(ast1.p*ast2.p)*(cosi + ast1.e*ast2.e*cosp);
	return m;
}
//вычисление метрики

/*
int main(void)
{

	//порядок интегратора
	nor=19;

//точность интегратора
	ll=11;

//настройка модели сил
	force_var[0]=1; //Меркурий
	force_var[1]=1; //Венера
	force_var[2]=1; //Земля
	force_var[3]=1; //Марс
	force_var[4]=1; //Юпитер
	force_var[5]=1; //Сатурн
	force_var[6]=1; //Уран
	force_var[7]=1; //Нептун
	force_var[8]=1; //Плутон
	force_var[9]=1; //Луна
	force_var[10]=1; //Солнце
	force_var[11]=0; //Церера, Паллада, Веста
	force_var[12]=1; //Сжатие Земли
	force_var[13]=1; //Сжатие Солнца
	force_var[14]=1; //Релятивистские эффекты
	force_var[15]=1; //Сжатие Юпитера

	nclass=2, ni=2, nv=3+9*force_var[11], xl=1.0e-9;

	
double testX[3], testY[3];
double F123[3];
testX[0]=0.092236854809664520;
testX[1]=-0.182559066416997630;
testX[2]=-0.031559776760198310;
testY[0]=0.047308464384152910;
testY[1]=-0.010961810814000040;
testY[2]= 0.016389372724113570;
double time=2457100.5;
DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
force(testX, testY, time, F123, acs);

return(0);}
*/

/*
int main()
{
	//порядок интегратора
	nor = 15;

	//точность интегратора
	ll = 11;

	//настройка модели сил
	force_var[0] = 1; //Меркурий
	force_var[1] = 1; //Венера
	force_var[2] = 1; //Земля
	force_var[3] = 1; //Марс
	force_var[4] = 1; //Юпитер
	force_var[5] = 1; //Сатурн
	force_var[6] = 1; //Уран
	force_var[7] = 1; //Нептун
	force_var[8] = 1; //Плутон
	force_var[9] = 1; //Луна
	force_var[10] = 1; //Солнце
	force_var[11] = 0; //Церера, Паллада, Веста
	force_var[12] = 1; //Сжатие Земли
	force_var[13] = 1; //Сжатие Солнца
	force_var[14] = 1; //Релятивистские эффекты
	force_var[15] = 0; //Сжатие Юпитера

	nclass = 2, ni = 2, nv = 3 + 9 * force_var[11], xl = 1.0e-9;

	string str_obs = "";
	//string str_date = "";
	//double date;
	double date = Date_JD(2017, 6, 1, 8, 4, 2);
	date = eph_time(date);
	ifstream fin_obs; //fin_date;
	fin_obs.open("obs.txt");
	//fin_date.open("date.txt");
	//ofstream fout2;
	ifstream fin168;
	fin168.open("search.txt");
	//fout2.precision(8);
	//fout2.open("domain.txt");
	//fout2<<"Name\t"<<"Data\t"<<"Magnitude\t"<<"Alpha\t"<<"Domain\t"<<"Delta\t"<<"Domain\t"<<"Iteration\t"<<"Sigma\t"<<endl;
	while (!fin_obs.eof())
	{
		getline(fin_obs, str_obs);
		//getline(fin_date, str_date);
		string ast = str_obs;
		str_obs += ".txt";

		//date=stod(str_date);

		ifstream fin5;
		fin5.open(str_obs);
		if (!fin5) exit(5);
		ofstream fout;

		fout.precision(18);
		fout.open("result-improve.txt");
		double ulian[8000], declination[8000], right[8000], magni[8000];
		string observatory[8000];
		int number = 0;
		double sat_num[5000], sat_x[5000], sat_y[5000], sat_z[5000];

		DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
		read_obs(fin5, ulian, declination, right, magni, observatory, number, sat_num, sat_x, sat_y, sat_z);
		fin5.close();
		string boul_m;
		boul_m = read_boul(ast);
		double time_m;
		//a_m, e_m, i_m, node_m, arg_m, mean_m, mag_m, g_m;
		time_m = boul_osc(boul_m);
		//mag_m = boul_abs(boul_m);
		//g_m = boul_g(boul_m);
		XYZ_boul(ast);
		double ekl_m[6], ekv[6];
		ekl_m[0] = X_boul;
		ekl_m[1] = Y_boul;
		ekl_m[2] = Z_boul;
		ekl_m[3] = VX_boul;
		ekl_m[4] = VY_boul;
		ekl_m[5] = VZ_boul;
		eklekv(ekl_m, ekv);
		double X_m[3], V_m[3];
		X_m[0] = ekv[0];
		X_m[1] = ekv[1];
		X_m[2] = ekv[2];
		V_m[0] = ekv[3];
		V_m[1] = ekv[4];
		V_m[2] = ekv[5];

		improve_elements(X_m, V_m, time_m, ulian, declination, right, observatory, number, sat_num, sat_x, sat_y, sat_z, acs);
		fin168 >> time_m >> X_m[0] >> X_m[1] >> X_m[2] >> V_m[0] >> V_m[1] >> V_m[2];
		fout << X_in[0] << " " << X_in[1] << " " << X_in[2] << " " << V_in[0] << " " << V_in[1] << " " << V_in[2] << " " << sigma;
		//improve_elements(X_m, V_m, time_m, ulian, declination, right, observatory, number, sat_num, sat_x, sat_y, sat_z, sigma, acs);
		domain2(ast);
		fout.close();
		fin5.close();
	}
	return 0;
}
*/
// мейн с улучшением для постройки областей, этот старый!!

/*
int main(array<System::String ^> ^args)
{
   //порядок интегратора
  nor=19;

  //точность интегратора
  ll=14;

  //настройка модели сил
  force_var[0]=1; //Меркурий
  force_var[1]=1; //Венера
  force_var[2]=1; //Земля
  force_var[3]=1; //Марс
  force_var[4]=1; //Юпитер
  force_var[5]=1; //Сатурн
  force_var[6]=1; //Уран
  force_var[7]=1; //Нептун
  force_var[8]=1; //Плутон
  force_var[9]=1; //Луна
  force_var[10]=1; //Солнце
  force_var[11]=0; //Церера, Паллада, Веста
  force_var[12]=1; //Сжатие Земли
  force_var[13]=1; //Сжатие Солнца
  force_var[14]=1; //Релятивистские эффекты

  nclass=2, ni=2, nv=3+9*force_var[11], xl=1.0e-9;

	string str_obs="";
	ifstream fin_obs;
	fin_obs.open("obs.txt");
	while(!fin_obs.eof()) {

	getline(fin_obs,str_obs);
	string ast=str_obs;
	str_obs+=".eph";

	ofstream fout5;
	fout5.open(str_obs);


	DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);

	string boul_m;
	boul_m=read_boul(ast);
	double time_m, a_m, e_m, i_m, node_m, arg_m, mean_m, mag_m, g_m;
	time_m=boul_osc(boul_m);
	mag_m=boul_abs(boul_m);
	g_m=boul_g(boul_m);
	XYZ_boul(ast);
	double ekl_m[6], ekv[6];
	ekl_m[0]=X_boul;
	ekl_m[1]=Y_boul;
	ekl_m[2]=Z_boul;
	ekl_m[3]=VX_boul;
	ekl_m[4]=VY_boul;
	ekl_m[5]=VZ_boul;
	eklekv(ekl_m, ekv);
	double X_m[3], V_m[3];
	X_m[0]=ekv[0];
	X_m[1]=ekv[1];
	X_m[2]=ekv[2];
	V_m[0]=ekv[3];
	V_m[1]=ekv[4];
	V_m[2]=ekv[5];


	double test_data=Date_JD(2015, 5, 1, 0, 0, 0);
	double test_data_2=Date_JD(2015, 5, 2, 0, 0, 0);

	double ephemer_time;
	ephemer_time=eph_time(test_data);
	double ephemer_time_2;
	ephemer_time_2=eph_time(test_data_2);
	rada27(X_m, V_m, time_m, ephemer_time, acs);
	double step_eph=10./(60.*24.);
	int d_height_eph, m_height_eph, d_asim_eph, m_asim_eph;
	double s_height_eph, s_asim_eph;
	double rogue;
	double mag_eph;
	double alpha_eph_test, delta_eph_test;
	
	for(rogue=test_data; rogue<=test_data_2; rogue=rogue+step_eph)
	{ephemerida(X_m, V_m, eph_time(rogue), rogue, "168", acs);
	 mag_eph=magnitude_observ(X_m, mag_m, g_m, eph_time(rogue), acs);
	 alpha_eph_test=h_alpha_eph+m_alpha_eph/60.+s_alpha_eph/3600.;
	 delta_eph_test=d_delta_eph+m_delta_eph/60.+s_delta_eph/3600.;
	 h_a("168", rogue, alpha_eph_test, delta_eph_test, d_height_eph, m_height_eph, s_height_eph, d_asim_eph, m_asim_eph, s_asim_eph);
	 Date_Grig(rogue);
	 fout5<<Year<<"\t"<<Month<<"\t"<<Day<<"\t"<<Hour<<"\t"<<Minute<<"\t"<<Second<<"\t"<<h_alpha_eph<<"\t"<<m_alpha_eph<<"\t"<<s_alpha_eph<<"\t"<<d_delta_eph<<"\t"<<abs(m_delta_eph)<<"\t"<<abs(s_delta_eph)<<"\t"<<mag_eph<<"\t"<<d_height_eph<<"\t"<<abs(m_height_eph)<<"\t"<<abs(s_height_eph)<<"\t"<<d_asim_eph<<"\t"<<m_asim_eph<<"\t"<<s_asim_eph<<endl;
	 rada27(X_m, V_m, eph_time(rogue), eph_time(rogue+step_eph), acs);}	
	fout5.close();
	}
	fin_obs.close();
	return 0;
}
*/
//мейн для эфемерид

/*int main (array<System::String ^> ^args)
{
   //порядок интегратора
  nor=19;

  //точность интегратора
  ll=11;

  //настройка модели сил
  force_var[0]=1; //Меркурий
  force_var[1]=1; //Венера
  force_var[2]=1; //Земля
  force_var[3]=1; //Марс
  force_var[4]=1; //Юпитер
  force_var[5]=1; //Сатурн
  force_var[6]=1; //Уран
  force_var[7]=1; //Нептун
  force_var[8]=1; //Плутон
  force_var[9]=1; //Луна
  force_var[10]=1; //Солнце
  force_var[11]=0; //Церера, Паллада, Веста
  force_var[12]=1; //Сжатие Земли
  force_var[13]=1; //Сжатие Солнца
  force_var[14]=1; //Релятивистские эффекты

  nclass=2, ni=2, nv=3+9*force_var[11], xl=1.0e-9;

    double Tupoi;
    string str_obs="";
	ifstream fin_obs;
	fin_obs.open("obs.txt");
	ofstream fout2;
	fout2.precision(8);
	fout2.open("opposition.txt");
	double scal, modX, modpoz, cosangle, angle, mag;
	double asterX[3], asterV[3];
	bool nikitka=false;
	int count=0;
	DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
	while(!fin_obs.eof()) {
	count++;
	Tupoi=Date_JD(2015, 10, 2, 0, 0 , 0);
	getline(fin_obs,str_obs);
	string ast=str_obs;
	XYZ_boul(ast);
	string aux=read_boul(ast);
	double T_boul=boul_osc(aux);
	double ekl_m[6], ekv[6];
	ekl_m[0]=X_boul;
	ekl_m[1]=Y_boul;
	ekl_m[2]=Z_boul;
	ekl_m[3]=VX_boul;
	ekl_m[4]=VY_boul;
	ekl_m[5]=VZ_boul;
	eklekv(ekl_m, ekv);
	double X_m[3], V_m[3];
	X_m[0]=ekv[0];
	X_m[1]=ekv[1];
	X_m[2]=ekv[2];
	V_m[0]=ekv[3];
	V_m[1]=ekv[4];
	V_m[2]=ekv[5];
	rada27(X_m, V_m, T_boul, Tupoi, acs);
	acs->GetPlanetPoz(Tupoi, 2, true, poz);
	gelio_to_planet(X_m, V_m, poz);
	poz[0]=poz[0]*-1;
	poz[1]=poz[1]*-1;
	poz[2]=poz[2]*-1;

	scal=planet_centric[0]*poz[0]+planet_centric[1]*poz[1]+planet_centric[2]*poz[2];
	modX=sqrt(planet_centric[0]*planet_centric[0]+planet_centric[1]*planet_centric[1]+planet_centric[2]*planet_centric[2]);
	modpoz=sqrt(poz[0]*poz[0]+poz[1]*poz[1]+poz[2]*poz[2]);
	cosangle=scal/(modX*modpoz);
	angle=acos(cosangle)*180/M_PI;
	mag=magnitude_observ(X_m, boul_abs(aux), boul_g(aux), Tupoi, acs);
	if((angle>40&&angle<320)&&mag<23)
	{nikitka=true;}
	while(nikitka==false)

	{Tupoi++;
	 cout<<Tupoi<<" "<<angle<<" "<<mag<<endl;
	 acs->GetPlanetPoz(Tupoi, 2, true, poz);
	 rada27(X_m, V_m, Tupoi-1, Tupoi, acs);
	 gelio_to_planet(X_m, V_m, poz);
	 poz[0]=poz[0]*-1;
	 poz[1]=poz[1]*-1;
	 poz[2]=poz[2]*-1;
	 scal=planet_centric[0]*poz[0]+planet_centric[1]*poz[1]+planet_centric[2]*poz[2];
	 modX=sqrt(planet_centric[0]*planet_centric[0]+planet_centric[1]*planet_centric[1]+planet_centric[2]*planet_centric[2]);
	 modpoz=sqrt(poz[0]*poz[0]+poz[1]*poz[1]+poz[2]*poz[2]);
	 cosangle=scal/(modX*modpoz);
	 angle=acos(cosangle)*180/M_PI;
	 mag=magnitude_observ(X_m, boul_abs(aux), boul_g(aux), Tupoi, acs);
	 if((angle>60&&angle<300)&&mag<22||Tupoi>2492500)
	 {nikitka=true;}
	}
	nikitka=false;
	Date_Grig(Tupoi);
	if(Tupoi>=2492500)
	{fout2<<-4000<<" "<<-4000<<" "<<-4000<<" "<<ast<<endl;
	 cout<<count<<" "<<ast<<" "<<"out of interval"<<endl;}
	else
	{fout2<<Tupoi<<" "<<angle<<" "<<mag<<" "<<ast<<endl;
	cout<<count<<" "<<ast<<" "<<"done"<<endl;}

	}
}

//мейн для оппозиций

/*
	int main(array<System::String ^> ^args)
{
   //порядок интегратора
  nor=19;

  //точность интегратора
  ll=11;

  //настройка модели сил
  force_var[0]=1; //Меркурий
  force_var[1]=1; //Венера
  force_var[2]=1; //Земля
  force_var[3]=1; //Марс
  force_var[4]=1; //Юпитер
  force_var[5]=1; //Сатурн
  force_var[6]=1; //Уран
  force_var[7]=1; //Нептун
  force_var[8]=1; //Плутон
  force_var[9]=1; //Луна
  force_var[10]=1; //Солнце
  force_var[11]=0; //Церера, Паллада, Веста
  force_var[12]=1; //Сжатие Земли
  force_var[13]=1; //Сжатие Солнца
  force_var[14]=1; //Релятивистские эффекты

  nclass=2, ni=2, nv=3+9*force_var[11], xl=1.0e-9;

  	string str_obs="";
	ifstream fin_obs;
	fin_obs.open("obs.txt");
	ofstream fout2;
	fout2.open("vidimost.txt");
	fout2.precision(8);
	double ulian[5000], declination[5000], right[5000], magni[5000];
	string observatory[5000]; 
	int number=0;
	double sat_num[1000], sat_x[1000], sat_y[1000], sat_z[1000];
	DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
	int i;
	double poz1[3], poz2[3];
	double aplha, delta;
	double angle;
	double rad_alpha1, rad_alpha2, rad_delta1, rad_delta2, alpha_diff;

	while(!fin_obs.eof()) {

	getline(fin_obs,str_obs);
	string ast=str_obs;
	str_obs+=".obs";

	ifstream fin5;
    fin5.open(str_obs);

	
	read_obs(fin5, ulian, declination, right, magni, observatory, number, sat_num, sat_x, sat_y, sat_z);
	for(i=0; i<number; i++)
	{cout<<i<<endl;
	 if(observatory[i]!="500")
	 {fout2<<eph_time(ulian[i])<<" "<<magni[i]<<" ";
	 acs->GetPlanetPoz(eph_time(ulian[i]), 2, true, poz);
	 poz1[0]=poz[0]*(-1);
	 poz1[1]=poz[1]*(-1);
	 poz1[2]=poz[2]*(-1);
	 poz2[0]=poz[3]*(-1);
	 poz2[1]=poz[4]*(-1);
	 poz2[2]=poz[5]*(-1);
	 ephemerida(poz1, poz2, eph_time(ulian[i]), ulian[i], observatory[i], acs);
	 alpha=(h_alpha_eph+m_alpha_eph/60.+s_alpha_eph/3600.)*15;
	 delta=d_delta_eph+m_delta_eph/60.+s_delta_eph/3600.;
	 rad_delta1=declination[i];
	 rad_delta2=delta/180*M_PI;
	 rad_alpha1=right[i];
	 rad_alpha2=alpha/180*M_PI;
	 alpha_diff=rad_alpha1-rad_alpha2;
	 if(alpha_diff<0)
	 {alpha_diff=alpha_diff+2*M_PI;}
	 if(alpha_diff>2*M_PI)
	 {alpha_diff=alpha_diff-2*M_PI;}
	 angle=acos(cos(M_PI/2.-rad_delta1)*cos(M_PI/2.-rad_delta2)+sin(M_PI/2.-rad_delta1)*sin(M_PI/2.-rad_delta2)*cos(alpha_diff));
	 angle=angle*180/M_PI;
	 fout2<<angle<<" "<<ast<<endl;}
	}
	 fin5.close(); 
	 cout<<ast<<" done"<<endl;


	}
	 
	fin_obs.close();
	fout2.close();
}


*/

//мейн для видимости
/*int main()
{
	//порядок интегратора
	nor = 19;

	//точность интегратора
	ll = 11;

	//настройка модели сил
	force_var[0] = 1; //Меркурий
	force_var[1] = 1; //Венера
	force_var[2] = 1; //Земля
	force_var[3] = 1; //Марс
	force_var[4] = 1; //Юпитер
	force_var[5] = 1; //Сатурн
	force_var[6] = 1; //Уран
	force_var[7] = 1; //Нептун
	force_var[8] = 1; //Плутон
	force_var[9] = 1; //Луна
	force_var[10] = 1; //Солнце
	force_var[11] = 0; //Церера, Паллада, Веста
	force_var[12] = 1; //Сжатие Земли
	force_var[13] = 1; //Сжатие Солнца
	force_var[14] = 1; //Релятивистские эффекты

	string s1 = "2000 SH89";
	string str_cat1;
	str_cat1 = read_boul(s1);
	//string str_cat1;
	//str_cat1=read_boul(s1);
	//double t1 = boul_osc(str_cat1);
	//double a1 = boul_axis(str_cat1);
	//double i1 = boul_i(str_cat1);
	//double e1 = boul_e(str_cat1);
	//double w1 = boul_arg(str_cat1);
	//double knot1 = boul_knot(str_cat1);
	//double A1 = boul_anomaly(str_cat1);
	

	XYZ_boul(s1);
	X_in[0] = X_boul;
	X_in[1] = Y_boul;
	X_in[2] = Z_boul;
	V_in[0] = VX_boul;
	V_in[1] = VY_boul;
	V_in[2] = VZ_boul;
	double t1 = boul_osc(str_cat1);
	ifstream fin1;
	fin1.open("165230.txt");
	if (!fin1) exit(1);
	int n=0;
	string obs[8000];
	double s_n[8000], s_x[8000], s_y[8000], s_z[8000];
	DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
	read_obs(fin1, ul, del, al, abss, obs, n, s_n, s_x, s_y, s_z);
	improve_elements(X_in, V_in, t1, ul, del, al, obs, n, s_n, s_x, s_y, s_z, acs);
	double mag1 = boul_abs(str_cat1);
	double g1 = boul_g(str_cat1);
	double date = Date_JD(2017, 3, 14, 0, 0, 0);
	domain(mag1, g1,date,"500", acs);
	return 0;
}*/
//Ирина

int main()
{
	//порядок интегратора
	nor = 19;

	//точность интегратора
	ll = 11;

	//настройка модели сил
	force_var[0] = 1; //Меркурий
	force_var[1] = 1; //Венера
	force_var[2] = 1; //Земля

	force_var[3] = 1; //Марс
	force_var[4] = 1; //Юпитер
	force_var[5] = 1; //Сатурн
	force_var[6] = 1; //Уран
	force_var[7] = 1; //Нептун
	force_var[8] = 1; //Плутон
	force_var[9] = 1; //Луна
	force_var[10] = 1; //Солнце
	force_var[11] = 0; //Церера, Паллада, Веста
	force_var[12] = 1; //Сжатие Земли
	force_var[13] = 1; //Сжатие Солнца
	force_var[14] = 1; //Релятивистские эффекты
	force_var[15] = 0; //Сжатие Юпитера
	nclass = 2, ni = 2, nv = 3 + 9 * force_var[11], xl = 1.0e-9;

	string str_obs = "";
	string ast;
	double X[3], V[3];
	double date_boul, date = Date_JD(2017, 6, 1, 8, 4, 2);
	ifstream fin_obs; //fin_date;
	fin_obs.open("Irina_in\\obs.txt");
	if (!fin_obs) exit(3);
	ofstream fout4;
	fout4.open("Irina_out\\nominal\\obs.kepler.t.txt");
	if (!fout4) exit(4);
	while (!fin_obs.eof())
	{
		getline(fin_obs, str_obs);
		cout << str_obs;
		ast=read_boul(str_obs);
		XYZ_boul(str_obs);
		X[0] = X_boul;
		X[1] = Y_boul;
		X[2] = Z_boul;
		V[0] = VX_boul;
		V[1] = VY_boul;
		V[2] = VZ_boul;
		cout<< X[0] << " " << X[1] << " " << X[2] << " " << V[0] << " " << V[1] << " " << V[2] << endl;
		date_boul=boul_osc(ast);
		cout << date<<endl;
		cout << date_boul<<endl;
		DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
		rada27(X, V, date_boul, date, acs);
		//fout4 << X[0] << " " << X[1] << " " << X[2] << " " << V[0] << " " << V[1] << " " << V[2] << endl;
		flat_to_Kepler(X, V);
		fout4.setf(ios::fixed);
		fout4.precision(16);
		fout4<< a_Kepler << " " << e_Kepler << " " << i_Kepler << " " << arg_Kepler << " " << knot_Kepler << endl;

	}
	fout4.close();
	fin_obs.close();
	Asteroid ast1, ast2;
	fin_obs.open("Irina_out\\nominal\\obs.kepler.t.txt");
	fin_obs>> ast1.a >> ast1.e >> ast1.i >> ast1.w >> ast1.node;
	ast1.q = ast1.a*(1 - ast1.e);
	ast1.p = ast1.a*(1 - ast1.e*ast1.e);
	fin_obs >> ast2.a >> ast2.e >> ast2.i >> ast2.w >> ast2.node;
	ast2.q = ast2.a*(1 - ast2.e);
	ast2.p = ast2.a*(1 - ast2.e*ast2.e);
	double metr = metric(ast1, ast2);
	double r = sqrt(metr);
	fout4.open("Irina_out\\nominal\\nominal_metric.txt");
	fout4.setf(ios::fixed);
	fout4.precision(16);
	fout4 << r << endl;
	fout4 << r*149597870.7;
	fin_obs.close();
	fout4.close();

}

//метрика номинальная орбита
/*
int main()
{
	int i = 1;
	while (1) {

		cout << i << endl;
		//порядок интегратора
		nor = 15;

		//точность интегратора
		ll = 11;

		//настройка модели сил
		force_var[0] = 1; //Меркурий
		force_var[1] = 1; //Венера
		force_var[2] = 1; //Земля
		force_var[3] = 1; //Марс
		force_var[4] = 1; //Юпитер
		force_var[5] = 1; //Сатурн
		force_var[6] = 1; //Уран
		force_var[7] = 1; //Нептун
		force_var[8] = 1; //Плутон
		force_var[9] = 1; //Луна
		force_var[10] = 1; //Солнце
		force_var[11] = 0; //Церера, Паллада, Веста
		force_var[12] = 1; //Сжатие Земли
		force_var[13] = 1; //Сжатие Солнца
		force_var[14] = 1; //Релятивистские эффекты
		force_var[15] = 0; //Сжатие Юпитера

		nclass = 2, ni = 2, nv = 3 + 9 * force_var[11], xl = 1.0e-9;


		ifstream fin_obs, fin_date;
		fin_obs.open("Irina_in\\obs.txt");
		fin_date.open("Irina_in\\time.txt");
		string str_obs = "";
		string str_date = "";
		double date;
		getline(fin_date, str_date);
		string t = str_date;
		cout << t << endl;
		date = stod(str_date);
		DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
		cout << "step 1" << endl;
		while (!fin_obs.eof())
		{
			getline(fin_obs, str_obs);
			string ast = str_obs;
			double timeimp;
			str_obs = ast + "_xyz_timp.txt";

			ifstream fin;
			string path_to_xyz_timp = "Irina_out\\main1\\";
			path_to_xyz_timp.append(str_obs);
			fin.open(path_to_xyz_timp);
			if (!fin) exit(1);
			double X[3], V[3];
			str_obs = ast + t + "_xyz_t.txt";
			ofstream fout1;
			string path_to_xyz_t = "Irina_out\\main2\\";
			path_to_xyz_t.append(str_obs);
			fout1.open(path_to_xyz_t);
			str_obs = ast + t + "_kepler_t.txt";
			ofstream fout2;
			string path_to_kepler_t = "Irina_out\\main2\\";
			path_to_kepler_t.append(str_obs);
			fout2.open(path_to_kepler_t);
			fin >> timeimp;
			cout << "step 2" << endl;
			while (!fin.eof())
			{
				fin >> X[0] >> X[1] >> X[2] >> V[0] >> V[1] >> V[2];
				rada27(X, V, timeimp, date, acs);
				fout1.precision(16);
				fout1 << X[0] << " " << X[1] << " " << X[2] << " " << V[0] << " " << V[1] << " " << V[2] << endl;
				flat_to_Kepler(X, V);
				fout2.precision(16);
				fout2 << a_Kepler << " " << e_Kepler << " " << i_Kepler << " " << arg_Kepler << " " << knot_Kepler << endl;

			}
			cout << "step 3" << endl;
			fin.close();
		}
		cout << "step 4" << endl;
		fin_obs.close();
		Asteroid ast1, ast2;
		int i = 1, a, end;
		double metr, r;
		string str1, str2;
		fin_obs.open("Irina_in\\obs.txt");
		while (!fin_obs.eof())
		{
			getline(fin_obs, str_obs);
			str1 = str_obs;
			cout << str1 << endl;
			getline(fin_obs, str_obs);
			str2 = str_obs;
			cout << str2 << endl;

			ifstream fin11;
			string path_to_kepler_t_1 = "Irina_out\\main2\\";
			path_to_kepler_t_1.append(str1);
			path_to_kepler_t_1.append(t);
			path_to_kepler_t_1.append("_kepler_t.txt");
			fin11.open(path_to_kepler_t_1);
			if (!fin11) exit(11);
			
			ifstream fin12;
			string path_to_kepler_t_2 = "Irina_out\\main2\\";
			path_to_kepler_t_2.append(str2);
			path_to_kepler_t_2.append(t);
			path_to_kepler_t_2.append("_kepler_t.txt");
			fin12.open(path_to_kepler_t_2);
			if (!fin12) exit(12);
			fin12.seekg(0, ios::end);
			end = fin12.tellg();
			ofstream foutm;
			string file = "Irina_out\\metric\\" + str1 + "_" + str2 + "_" + t + ".txt";
			foutm.open(file);
			if (!foutm) exit(8);

			while (!fin11.eof())
			{
				fin12.seekg(0, ios::beg);
				int s = fin12.tellg();
				fin11 >> ast1.a >> ast1.e >> ast1.i >> ast1.w >> ast1.node;
				ast1.q = ast1.a*(1 - ast1.e);
				ast1.p = ast1.a*(1 - ast1.e*ast1.e);
				do
				{
					s = fin12.tellg();
					fin12 >> ast2.a >> ast2.e >> ast2.i >> ast2.w >> ast2.node;
					ast2.q = ast2.a*(1 - ast2.e);
					ast2.p = ast2.a*(1 - ast2.e*ast2.e);
					metr = metric(ast1, ast2);
					r = sqrt(metr);
					foutm.setf(ios::fixed);
					foutm.precision(7);
					foutm << ast1.a << " " << ast1.e << " " << ast1.i << " " << ast1.w << " " << ast1.node << " ";
					foutm << ast2.a << " " << ast2.e << " " << ast2.i << " " << ast2.w << " " << ast2.node << " ";
					foutm << r << endl;
					a = fin12.tellg();
				} while (abs(a - end) > 18 * 5 + 1);
			}
			fin12.close();
			fin11.close();
			fin_date.seekg(0, ios::beg);
			i++;
		}
		cout << "step 5" << endl;
		fin_obs.close();
		fin_date.close();

		fin_date.open("Irina_in\\time.txt");
		double dateJD[1000];
		int n = 0;
		while (!fin_date.eof())
		{
			getline(fin_date, str_date);
			dateJD[n] = atof(str_date.c_str());
			n++;
		}

		fin_date.close();
		ofstream fout_date;
		fout_date.open("Irina_in\\time.txt", ios::trunc);
		fout_date.setf(ios::fixed);
		fout_date.precision(1);
		for (int m = 1; m <= n; m++)
		{
			fout_date << dateJD[m] << endl;
		}
		cout << "step end" << endl;
		i++;
		 
	}
	system("pause");

	return 0;
}
*/
//вычисление метрики, второй. +перезапись файла времени

/*
int main()
{
		//порядок интегратора
		nor = 15;

		//точность интегратора
		ll = 11;

		//настройка модели сил
		force_var[0] = 1; //Меркурий
		force_var[1] = 1; //Венера
		force_var[2] = 1; //Земля
		force_var[3] = 1; //Марс
		force_var[4] = 1; //Юпитер
		force_var[5] = 1; //Сатурн
		force_var[6] = 1; //Уран
		force_var[7] = 1; //Нептун
		force_var[8] = 1; //Плутон
		force_var[9] = 1; //Луна
		force_var[10] = 1; //Солнце
		force_var[11] = 0; //Церера, Паллада, Веста
		force_var[12] = 1; //Сжатие Земли
		force_var[13] = 1; //Сжатие Солнца
		force_var[14] = 1; //Релятивистские эффекты
		force_var[15] = 0; //Сжатие Юпитера

		nclass = 2, ni = 2, nv = 3 + 9 * force_var[11], xl = 1.0e-9;


		string str_obs = "";
	string str_date = "";
	double date;
	//double date = Date_JD(2017, 6, 1, 8, 4, 2);
	//date = eph_time(date);
	ifstream fin_obs, fin_date;
	fin_obs.open("Irina_in\\obs.txt");
	if (!fin_obs) exit(1);
	//fin_date.open("Irina_in\\date.txt");
	//if (!fin_date) exit(2);
	//ofstream fout2;
	ifstream fin168;
	fin168.open("search.txt");
	//fout2.precision(8);
	//fout2.open("domain.txt");
	//fout2<<"Name\t"<<"Data\t"<<"Magnitude\t"<<"Alpha\t"<<"Domain\t"<<"Delta\t"<<"Domain\t"<<"Iteration\t"<<"Sigma\t"<<endl;
	while (!fin_obs.eof())
	{
		getline(fin_obs, str_obs);
		//while (!fin_date.eof())
		//{
		//getline(fin_date, str_date);
		string ast = str_obs;
		//string t = str_date;
		str_obs += ".txt";

		//date = stod(str_date);

		ifstream fin5;
		string path_to_ast = "Irina_in\\";
		path_to_ast.append(str_obs);
		fin5.open(path_to_ast);
		if (!fin5) exit(5);
		ofstream fout;

		fout.precision(18);
		fout.open("Irina_out\\result-improve.txt");
		double ulian[8000], declination[8000], right[8000], magni[8000];
		string observatory[8000];
		int number = 0;
		double sat_num[5000], sat_x[5000], sat_y[5000], sat_z[5000];
		cout << "step 1" << endl;
		DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
		cout << "step 2" << endl;
		read_obs(fin5, ulian, declination, right, magni, observatory, number, sat_num, sat_x, sat_y, sat_z);
		cout << "step 3" << endl;
		fin5.close();
		string boul_m;
		boul_m = read_boul(ast);
		cout << "step 4" << endl;
		double time_m;
		//a_m, e_m, i_m, node_m, arg_m, mean_m, mag_m, g_m;
		time_m = boul_osc(boul_m);
		//mag_m = boul_abs(boul_m);
		//g_m = boul_g(boul_m);
		XYZ_boul(ast);
		double ekl_m[6], ekv[6];
		ekl_m[0] = X_boul;
		ekl_m[1] = Y_boul;
		ekl_m[2] = Z_boul;
		ekl_m[3] = VX_boul;
		ekl_m[4] = VY_boul;
		ekl_m[5] = VZ_boul;
		eklekv(ekl_m, ekv);
		double X_m[3], V_m[3];
		X_m[0] = ekv[0];
		X_m[1] = ekv[1];
		X_m[2] = ekv[2];
		V_m[0] = ekv[3];
		V_m[1] = ekv[4];
		V_m[2] = ekv[5];
		cout << "step 5" << endl;
		improve_elements(X_m, V_m, time_m, ulian, declination, right, observatory, number, sat_num, sat_x, sat_y, sat_z, acs);
		cout << "step 6" << endl;
		fin168 >> time_m >> X_m[0] >> X_m[1] >> X_m[2] >> V_m[0] >> V_m[1] >> V_m[2];
		fout << X_in[0] << " " << X_in[1] << " " << X_in[2] << " " << V_in[0] << " " << V_in[1] << " " << V_in[2] << " " << sigma;
		//improve_elements(X_m, V_m, time_m, ulian, declination, right, observatory, number, sat_num, sat_x, sat_y, sat_z, sigma, acs);
		domain2(ast);
		cout << "step 7" << endl;
		fout.close();
		fin5.close();

		
		//}
		//fin_date.seekg(0, ios::beg);
	}
	system("pause");
}
*/
//первый, постройка областей
/*
int main()
{
	//порядок интегратора
	nor = 15;

	//точность интегратора
	ll = 11;

	//настройка модели сил
	force_var[0] = 1; //Меркурий
	force_var[1] = 1; //Венера
	force_var[2] = 1; //Земля
	force_var[3] = 1; //Марс
	force_var[4] = 1; //Юпитер
	force_var[5] = 1; //Сатурн
	force_var[6] = 1; //Уран
	force_var[7] = 1; //Нептун
	force_var[8] = 1; //Плутон
	force_var[9] = 1; //Луна
	force_var[10] = 1; //Солнце
	force_var[11] = 0; //Церера, Паллада, Веста
	force_var[12] = 1; //Сжатие Земли
	force_var[13] = 1; //Сжатие Солнца
	force_var[14] = 1; //Релятивистские эффекты
	force_var[15] = 0; //Сжатие Юпитера

	nclass = 2, ni = 2, nv = 3 + 9 * force_var[11], xl = 1.0e-9;


	ifstream fin_obs, fin_date;
	fin_obs.open("obs.txt");
	fin_date.open("time.txt");
	string str_obs = "";
	string str_date = "";
	double date;
	DeAcsessor^ acs = gcnew DeAcsessor(binPathForHeader, binPathForData);
	while (!fin_obs.eof())
	{
		getline(fin_obs, str_obs);
		string ast = str_obs;
		double timeimp;
		str_obs = ast + ".xyz.timp.txt";
		while (!fin_date.eof())
		{
			getline(fin_date, str_date);
			string t = str_date;
			date = stod(str_date);
			ifstream fin;
			fin.open(str_obs);
			if (!fin) exit(111);
			double X[3], V[3];
			str_obs = ast + t + ".xyz.t.txt";
			ofstream fout1;
			fout1.open(str_obs);
			str_obs = ast + t + ".kepler.t.txt";
			ofstream fout2;
			fout2.open(str_obs);
			fin >> timeimp;
			while (!fin.eof())
			{
				fin >> X[0] >> X[1] >> X[2] >> V[0] >> V[1] >> V[2];
				rada27(X, V, timeimp, date, acs);
				fout1.precision(16);
				fout1 << X[0] << " " << X[1] << " " << X[2] << " " << V[0] << " " << V[1] << " " << V[2] << endl;
				flat_to_Kepler(X, V);
				fout2.precision(16);
				fout2 << a_Kepler << " " << e_Kepler << " " << i_Kepler << " " << arg_Kepler << " " << knot_Kepler << endl;

			}
			fin.close();
		}
		fin_date.seekg(0, ios::beg);
	}
	fin_obs.close();
	Asteroid ast1, ast2;
	int i = 1, a, end;
	double metr, r;
	string str1, str2;
	fin_obs.open("obs.txt");
	while (!fin_obs.eof())
	{
		getline(fin_obs, str_obs);
		str1 = str_obs;
		cout << str1 << endl;
		getline(fin_obs, str_obs);
		str2 = str_obs;
		cout << str2;
		while (!fin_date.eof())
		{
			getline(fin_date, str_date);
			str1 += str_date += ".kepler.t.txt";
			ifstream fin11;
			fin11.open(str1);
			if (!fin11) exit(11);
			i++;
			//while (!fin_obs.eof())
			//{
			str2 += str_date += ".kepler.t.txt";
			ifstream fin12;
			fin12.open(str2);
			if (!fin12) exit(12);
			fin12.seekg(0, ios::end);
			end = fin12.tellg();
			ofstream foutm;
			string file = str1 + "**" + str2 + ".txt";
			foutm.open(file);
			if (!foutm) exit(8);

			while (!fin11.eof())
			{
				fin12.seekg(0);
				int s = fin12.tellg();
				fin11 >> ast1.a >> ast1.e >> ast1.i >> ast1.w >> ast1.node;
				ast1.q = ast1.a*(1 - ast1.e);
				ast1.p = ast1.a*(1 - ast1.e*ast1.e);
				do
				{
					s = fin12.tellg();
					fin12 >> ast2.a >> ast2.e >> ast2.i >> ast2.w >> ast2.node;
					ast2.q = ast2.a*(1 - ast2.e);
					ast2.p = ast2.a*(1 - ast2.e*ast2.e);
					metr = metric(ast1, ast2);
					r = sqrt(metr);
					foutm.setf(ios::fixed);
					foutm.precision(7);
					foutm << ast1.a << " " << ast1.e << " " << ast1.i << " " << ast1.w << " " << ast1.node << " ";
					foutm << ast2.a << " " << ast2.e << " " << ast2.i << " " << ast2.w << " " << ast2.node << " ";
					foutm << r << endl;
					a = fin12.tellg();
				} while (abs(a - end) > 18 * 6);
			}
			fin12.close();
			//}
			fin11.close();
		}
		fin_date.seekg(0, ios::beg);
	}
	fin_obs.close();
	system("pause");
	return 0;
}
*/

//второй, работает сразу с файлом где времени много, не работает