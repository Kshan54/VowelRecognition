// testingvowel.cpp : Defines the entry point for the console application.
#include "stdafx.h"
#include<iostream>
#include<fstream>
#include <math.h> 
#include<sstream>
#include<vector>
#include <string>
#include<climits>
using namespace std;
long double weights[12]={1.0,3.0,7.0,13.0,19.0,22.0,25.0,33.0,42.0,50.0,56.0,61.0};
// Tokhura Distance calculation
void tokhura(vector<vector<long double>>ref,vector<vector<long double>>test ){
	long double min_dist=0;
	long int min_index;
	long double distance=0.0;
	for(int i=0;i<ref.size();i++)
	{
		int index=(i+1)%5;
		if(index==0)
		{
			for(int j=0;j<12;j++)
			{
				distance+=weights[j]*((test[index][j+1]-ref[i][j])*(test[index][j+1]-ref[i][j]));
			}
			distance=distance/5.0;
			if(min_dist==0)
			{
				min_dist=distance;
				min_index=i;
			}
			else if(min_dist!=0)
			{
				if(distance<min_dist)
				{
					min_dist=distance;
					min_index=i;
				}
			}
			distance=0;
		}
		else
		{
			for(int j=0;j<12;j++)
			{
				distance+=weights[j]*((test[index][j+1]-ref[i][j])*(test[index][j+1]-ref[i][j]));
			}
		}
	}
	// using index for detecting a / e / i / o / u that is reference file is also created in this order to make it work correctly
	if(min_index<5)
	{
		cout<<"A Detected"<<" DISTANCE : "<<min_dist<<endl;
		cout<<endl;
	}
	else if(min_index>=5 && min_index<10)
	{
		cout<<"E Detected"<<" DISTANCE : "<<min_dist<<endl;
		cout<<endl;
	}
	else if(min_index>=10 && min_index<15)
	{
		cout<<"I Detected"<<" DISTANCE : "<<min_dist<<endl;
		cout<<endl;
	}
	else if(min_index>=15 && min_index<20)
	{
		cout<<"O detected"<<" DISTANCE : "<<min_dist<<endl;
		cout<<endl;
	}
	else
	{
		cout<<"U detected"<<" DISTANCE : "<<min_dist<<endl;
		cout<<endl;
	}
}
// Steady Frames around Maximum short term energy (5 frames)
vector<long double> steadyf(vector<long double> &v){
	long double max=0;
	long double sum=0;
	long long int max_ind=0;
	long long int skip=320*10; // skipping 10 frames start to avoid error during live testing (*due to mic spike)
	long long int i=skip;
	for(;i<v.size();i++)
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
	//cout<<"Max index "<<max_ind<<" Max amp "<<max<<endl;
	vector<long double> vst;
	long long int start_ind=max_ind-(320)*3;
	for(int i=start_ind;i<start_ind+(320)*5;i++)
	{
		vst.push_back(v[i]);
	}
	return vst;
}
// Cepstral Ceoefficients Calculation without inverting ai's
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
// This Functions caculates lpc coefficients pushes it in vector of vector lpc_a
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
// Normalize whole signal with max amplitude 10000 as reference
void normalization(vector<long double> &v1 ){
	long double max_amp=0.0; // for storing +ve max amplitude
	long double nmax_amp=0.0;// for storing -ve max amplitude
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
// Appyling Hamming Window and Calculating Ri's Threshold here is of silence 
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
	// Loading and storing data from reference File r2.txt
	ifstream Filer;
	Filer.open("r2.txt");
	vector<long double> vref;
	long double k;
	while(Filer >> k){
		vref.push_back(k);
	}
	Filer.close();
	vector<long double> init;
	vector<vector<long double>> ref;
	int j=0;
	for(int i=0;i<vref.size();i++)
	{
		int z=i%12;
		if(z==0 && i!=0)
		{
			ref.push_back(init);
			init.resize(0);
		}
		init.push_back(vref[i]);
	}
	ref.push_back(init);
	vector< vector<long double>> cepstral;
	cout<<"MENU ::"<<endl;
	cout<<"1).Enter the Vowel for testing"<<endl;
	cout<<"2).Enter r to start Live testing "<<endl;
	cout<<"3).Enter z to stop Testing"<<endl;
	while(1)
	{
		string stv;
		cout<<"Enter the Vowel for testing a , e , i , o , u in  small caps / Enter 'z' to stop Testing / Enter 'r' for live testing : "<<endl;
		cin>>stv;
		if(stv=="z")
		{
			break;
		}
		else if(stv=="a" || stv=="r" || stv=="e" || stv=="i" || stv=="o" || stv=="u")
		{
			int flag=0; // To know whether Live testing is being done or not
			if(stv=="r")
			{
				cout<<"Recording module started"<<endl;
				flag=1;
			}
			int total=10;
			int count=11;
			while(total--)
			{
				ifstream File;
				vector<long double> v1;
				// Convert Number to string and generating File name
				if(flag==0)
				{
					string s1 ="204101051_";
					s1+=stv;
					s1+="_";
					stringstream ss;
					ss << count;
					string str = ss.str();
					s1+=str;
					string s3=".txt";
					s1+=s3;
					File.open(s1);
					cout<<s1<<endl;
					cout<<endl;
				}
				// else Part For Live testing
				else
				{
					system("Recording_Module.exe 5 input_file.wav input_file.txt");
					File.open("input_file.txt");
					cout<<"Record File name: input_file.txt"<<endl;
					flag=0;
					total=0;
				}
				
				long double k;
				while(File >> k){
					v1.push_back(k);
				}
				File.close();
				for(long long int i=0;i<v1.size();i++)
				{
					v1[i]=v1[i]+0.0220778;
				}
				// normalization call
				normalization(v1);
				// steady Frames Calculation
				v1=steadyf(v1);
				vector<vector<long double>> v_r;
				for(int i=0;i<v1.size();i+=320)
				{
					if(i+320<=v1.size())
					{
						// threshold 61312 for silence use for autocorellation so Ri's calculation are Skipped for silence
						vector<long double> result=hamming(v1,i,61312);
						if(result.size()!=0)
						{
							v_r.push_back(result);
						}
						result.resize(0,0);
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
					// cep function call for calculating cepstarl coeffiecients
					for(int j=0;j<lpc_a[i].size();j++)
					{
						cep(ceps,lpc_a[i]);
					}
					cepstral.push_back(ceps);
					ceps.resize(0);
				}
				// Applying Raised sine window
				for(int i=0;i<cepstral.size();i++)
				{
					for(int j=1;j<=12;j++)
					{
						long int k1=j;
						cepstral[i][j]=(1+6*sin((3.141*k1)/12))*cepstral[i][j];
					}
				}
				// tokhura function call for distance calculation and detecting vowel
				tokhura(ref,cepstral);
				lpc_a.clear();
				v_r.clear();
				cepstral.clear();
				v1.clear();
				count++;
		}
	 }
		else
		{
			cout<<"*********************************************"<<endl;
			cout<<"Enter The Correct Choice"<<endl;
			cout<<"*********************************************"<<endl;
		}
	}
}

