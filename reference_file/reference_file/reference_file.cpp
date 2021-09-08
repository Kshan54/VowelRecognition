// dc_shift.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include<iostream>
#include<fstream>
#include <sstream>
#include <math.h> 
#include<vector>
#include <stdlib.h>
#include <string>
using namespace std;
vector<long double> steadyf(vector<long double> &v){
	long double max=0;
	long double sum=0;
	long long int max_ind=0;
	for(long long int i=0;i<v.size();i++)
	{
		if((i+1)%320==0)    
		{
			if(sum>max)
			{
				max=sum;
				max_ind=i;
			}
			sum=0;
		}
		else
		sum+=v[i]*v[i];
	}
	vector<long double> vst;
	long long int start_ind=max_ind-(320)*3;
	for(int i=start_ind;i<start_ind+(320)*5;i++)
	{
		vst.push_back(v[i]);
	}
	return vst;
}
void cep(vector<long double> &cepstral,vector<long double> &lpca)
{
	for(int i=1;i<=12;i++)
	{
		long double sum=0;
		for(int j=1;j<=i-1;j++)
		{
			long double k1=j;
			long double k2=i;
			sum+=(k1/k2)*(cepstral[j]*((lpca[i-j-1])));
		}
		sum+=((lpca[i-1]));
		cepstral.push_back(sum);
	}
}
void levinson_durbin(vector<long double> rval,vector<vector<long double>> &lpc_a)
{
	vector<long double> eval;
	vector<long double> kval;
	vector<long double> fval;
	long double arr[13][13];
	eval.push_back(rval[0]);
	kval.push_back(0);
	for(int i=1;i<13;i++)
	{
		double kk_init=0;
		for(int j=1;j<=i-1;j++)
		{
			kk_init+=arr[i-1][j]*rval[i-j];
		}
		kk_init=(rval[i]-kk_init)/eval[i-1];	
		kval.push_back(kk_init);
		arr[i][i]=kval[i];
		for(int j=1;j<=i-1;j++)
		{
			arr[i][j]=arr[i-1][j]-kval[i]*arr[i-1][i-j];
		}
		double g=(1-kval[i]*kval[i])*eval[i-1];
		eval.push_back(g);
	}
	for(int i=1;i<13;i++)
	{
		fval.push_back(arr[12][i]);
	}
	lpc_a.push_back(fval);
}
void normalization(vector<long double> &v1 ){
	long double max_amp=0.0;
	long double nmax_amp=0.0;
	for(int i=0;i<v1.size();i++){
		if(v1[i]>max_amp){
			max_amp=v1[i];
		}
		if(v1[i]<nmax_amp)
		{
			nmax_amp=v1[i];
		}
	}
	nmax_amp=-nmax_amp;
	if(nmax_amp>max_amp)
	{
		max_amp=nmax_amp;
	}
	long double ratio= 10000.0/max_amp;
	for(int i=0;i<v1.size();i++){
		v1[i]=v1[i]*ratio;
	}
}
vector<long double> hamming(vector<long double> v1,int start_index,long double threshold){
	vector<long double> rval;
	//Applying Hamming window
	long long int ind=0;
	for(int i=start_index;i<start_index+320;i++){
		v1[i]=v1[i]*(0.54-0.46*cos((2*3.141*ind)/319));
		ind++;
	}
	int flag=1;
	for(int i=0;i<13;i++)
	{
	    long double sum=0;
		for(int j=start_index;j<(start_index+320)-i;j++){
			sum+=v1[j]*v1[j+i];
		}
		if(i==0 && flag==1)
		{
			if(sum>threshold)
			{
				flag=0;
			}
			else
			{
				break;
			}
		}
		rval.push_back(sum);
	}
	return rval;

}
int _tmain(int argc, _TCHAR* argv[])
{
	int t=10;
	int count=1;
	vector< vector<long double>> cepstral;
	while(t--)
	{
		ifstream File;
		string s1 ="204101051_a_";
		stringstream ss;
		ss << count;
		string str = ss.str();
		s1+=str;
		string s3=".txt";
		s1+=s3;
		File.open(s1);
		cout<<s1<<endl;
		vector<long double> v1;
		long double k;
		while(File >> k){
			v1.push_back(k);
		}
		File.close();
		//cout<<"v1 size() "<<v1.size()<<endl;
		for(long long int i=0;i<v1.size();i++)
		{
			v1[i]=v1[i]+0.0220778;
		}
		normalization(v1);
		v1=steadyf(v1);
		vector<vector<long double>> v_r;
		for(int i=0;i<v1.size();i+=320)
		{
			if(i+320<=v1.size())
			{
				// threshold 61312
				vector<long double> result=hamming(v1,i,61312);
				if(result.size()!=0)
				{
					v_r.push_back(result);
				}
			}
		}
		vector<vector<long double>> lpc_a;

		for(int i=0;i<v_r.size();i++)
		{
			levinson_durbin(v_r[i],lpc_a);
		}

		
		for(int i=0;i<lpc_a.size();i++)
		{
			long double nsigma=0;
			nsigma=v_r[i][0]*v_r[i][0];
			nsigma=log(nsigma);
			vector<long double> ceps;
			ceps.push_back(nsigma);
			for(int j=0;j<lpc_a[i].size();j++)
			{
				cep(ceps,lpc_a[i]);
			}
			cepstral.push_back(ceps);
			ceps.resize(0);
		}
		count++;
	}
	for(int i=0;i<cepstral.size();i++)
	{
		for(int j=1;j<=12;j++)
		{
			long double k1=j;
			cepstral[i][j]=(1+6*sin((3.141*k1)/12))*cepstral[i][j];
		}
	}
	vector<long double> fcepstral;
	vector<vector<long double>> fi;
	
	for(int i=0;i<5;i++)
	{
		fcepstral.resize(12,0);
		for(int j=i;j<50;j+=5)
		{
			for(int k=1;k<=12;k++)
			{
				long double sum=cepstral[j][k];
				fcepstral[k-1]+=sum;
			}
		}
		for(int i=0;i<fcepstral.size();i++)
		{
			fcepstral[i]=fcepstral[i]/10.0;
		}
		fi.push_back(fcepstral);
		fcepstral.resize(0);
	}
	/*for(int i=0;i<fi.size();i++)
	{
		cout<<"Average : Cepstral Coefficient of "<<i+1<<"th Frame"<<endl;
		for(int j=0;j<12;j++)
		{
			cout<<fi[i][j]<<" ";
		}
		cout<<endl;
	}*/

	ofstream outputfile("r.txt",ios::out | ios::app | ios::ate);
	for(int i=0;i<fi.size();i++)
	{
		for(int j=0;j<12;j++)
		{
			outputfile<<fi[i][j]<<"\n";
		}
	}
	cout<<"Update the code replace a , e , i , o , u  in '2014101051_a_' string every time for reference file creation in this order only for testing correctly"<<endl;
}



