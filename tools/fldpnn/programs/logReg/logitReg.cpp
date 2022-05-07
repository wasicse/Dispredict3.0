#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <cmath>
//#include <cstdlib>

const int maxline=10000;
const double meanpl=0.180,minpl=0.005,maxpl=0.625;
const double meanpr=0.151,minpr=0.000,maxpr=0.991;
const double meanpd=0.245,minpd=0.003,maxpd=0.884;
const double meanpp=0.799,minpp=0.041,maxpp=0.980;
const double meanpm=0.5,minpm=0.0407,maxpm=0.9163;
const double pfun_threshold=0.535912;//0.374011;
const double pfun2_threshold=0.464220;//0.334392;
const int localwindowsize=5;
const int longwindowsize=15;
const double THRESHOLD=0.3533715368902425;

using namespace std;

class FEATURE{
public:
	int idx;
	char res;
	double pln,prn,pdn,ppn,pmn,pdisln,pdissn;
	int plb,prb,pdb,ppb,pmb,pdislb,pdissb;
	double pssm[20];
	double pcon;
	double pss[3];
	int pssb[3];
	double pfunn,pfun2n;
	int pfunb, pfun2b;
friend ostream& operator<<(ostream& os,const FEATURE& f);	
};

ostream& operator<<(ostream& os,const FEATURE& f)
{
	os<<f.pln<<"\t"<<f.plb<<"\t";
        os<<f.prn<<"\t"<<f.prb<<"\t";
        os<<f.pdn<<"\t"<<f.pdb<<"\t";
        os<<f.ppn<<"\t"<<f.ppb<<"\t";
        os<<f.pmn<<"\t"<<f.pmb<<"\t";
        os<<f.pcon<<"\t";
        os<<f.pdisln<<"\t"<<f.pdislb<<"\t";
        os<<f.pdissn<<"\t"<<f.pdissb<<"\t";
        os<<f.pfunn<<"\t"<<f.pfunb<<"\t";
        os<<f.pfun2n<<"\t"<<f.pfun2b<<"\t";
        for(int j=0;j<3;j++)
                os<<f.pss[j]<<"\t";
        for(int j=0;j<3;j++)
                os<<f.pssb[j]<<"\t";
        for(int j=0;j<20;j++)
                os<<f.pssm[j]<<"\t";
	return os;
}

class FEATURE_SET{
public:
        int idx;
        char res;
        double pln,prn,pdn,ppn,pmn,pdisln,pdissn;
        double plb,prb,pdb,ppb,pmb,pdislb,pdissb;
        double pssm[20];
        double pcon;
        double pss[3];
        double pssb[3];
        double pfunn,pfun2n;
        double pfunb, pfun2b;
	void add(FEATURE&);
	void divid(double);
	FEATURE_SET(); 
friend ostream& operator<<(ostream& os,const FEATURE_SET& f);
};

FEATURE_SET::FEATURE_SET()
{
	pln=0;
	prn=0;
	pdn=0;
	ppn=0;
	pmn=0;
	pdisln=0;
	pdissn=0;
	pfunn=0;
	pfun2n=0;

	plb=0;
	prb=0;
	pdb=0;
	ppb=0;
	pmb=0;
	pdislb=0;
	pdissb=0;
	pfunb=0;
	pfun2b=0;
	
	pcon=0;

	for(int i=0;i<20;i++)
		pssm[i]=0;
	for(int i=0;i<3;i++)
	{
		pss[i]=0;
		pssb[i]=0;
	}
}

void FEATURE_SET::add(FEATURE& f)
{
	pln+=f.pln;
	prn+=f.prn;
	pdn+=f.pdn;
	ppn+=f.ppn;
	pmn+=f.pmn;
	pdisln+=f.pdisln;
	pdissn+=f.pdissn;
	pfunn+=f.pfunn;
	pfun2n+=f.pfun2n;

	plb+=f.plb;
	prb+=f.prb;
	pdb+=f.pdb;
	ppb+=f.ppb;
	pmb+=f.pmb;
	pdislb+=f.pdislb;
	pdissb+=f.pdissb;
	pfunb+=f.pfunb;
	pfun2b+=f.pfun2b;
	
	pcon+=f.pcon;

	for(int i=0;i<20;i++)
		pssm[i]+=f.pssm[i];
	for(int i=0;i<3;i++)
	{
		pss[i]+=f.pss[i];
		pssb[i]+=f.pssb[i];
	}
} 

void FEATURE_SET::divid(double d)
{
	pln/=d;
	prn/=d;
	pdn/=d;
	ppn/=d;
	pmn/=d;
	pdisln/=d;
	pdissn/=d;
	pfunn/=d;
	pfun2n/=d;

	plb/=d;
	prb/=d;
	pdb/=d;
	ppb/=d;
	pmb/=d;
	pdislb/=d;
	pdissb/=d;
	pfunb/=d;
	pfun2b/=d;
	
	pcon/=d;

	for(int i=0;i<20;i++)
		pssm[i]/=d;
	for(int i=0;i<3;i++)
	{
		pss[i]/=d;
		pssb[i]/=d;
	}
}

ostream& operator<<(ostream& os,const FEATURE_SET& f)
{
	os<<f.pln<<"\t"<<f.plb<<"\t";
        os<<f.prn<<"\t"<<f.prb<<"\t";
        os<<f.pdn<<"\t"<<f.pdb<<"\t";
        os<<f.ppn<<"\t"<<f.ppb<<"\t";
        os<<f.pmn<<"\t"<<f.pmb<<"\t";
        os<<f.pcon<<"\t";
        os<<f.pdisln<<"\t"<<f.pdislb<<"\t";
        os<<f.pdissn<<"\t"<<f.pdissb<<"\t";
        os<<f.pfunn<<"\t"<<f.pfunb<<"\t";
        os<<f.pfun2n<<"\t"<<f.pfun2b<<"\t";
        for(int j=0;j<3;j++)
                os<<f.pss[j]<<"\t";
        for(int j=0;j<3;j++)
                os<<f.pssb[j]<<"\t";
        for(int j=0;j<20;j++)
                os<<f.pssm[j]<<"\t";
	return os;
}

class SAMPLE{
public:
	FEATURE single[5];
	FEATURE_SET local;
	FEATURE_SET global;
	int label;
	double pred;
friend ostream &operator<<(ostream &os,const SAMPLE &sample);
};

ostream &operator<<(ostream &os,const SAMPLE &sample)
{
	for(int i=0;i<5;i++)
		os<<sample.single[i];
	os<<sample.local<<sample.global;
	return os;
}


bool read_dfl(vector<FEATURE>&,string,bool clear=true);
bool read_rdpbind(vector<FEATURE>&,string,bool clear=false);
bool read_fmorf(vector<FEATURE>&,string,bool clear=false);
bool read_pssm(vector<FEATURE>&,string,bool clear=false);
double conServation(double*);
bool read_ss(vector<FEATURE>&,string,bool clear=false);
bool read_iupred(vector<FEATURE>&,string,int,bool clear=false);
bool getpfun(vector<FEATURE>&);
bool generateAllFeatures(vector<SAMPLE>&,vector<FEATURE>&);

void logReg(vector<SAMPLE>&);
void OutputPred(vector<SAMPLE>&,string);
void OutputFeatures(vector<FEATURE>&);
void OutputSamples(vector<SAMPLE>&,string,bool label_flag=false);

int main(int argc,char* argv[])
{
	string abspath=argv[1];
	string id=argv[2];
	vector<FEATURE> vfeature;
	vector<SAMPLE> vsample;
	
	read_dfl(vfeature,abspath+"/"+id+".dfl");
	read_rdpbind(vfeature,abspath+"/"+id+".rdp");
	read_fmorf(vfeature,abspath+"/"+id+".fmorf");
	read_pssm(vfeature,abspath+"/"+id+".pssm");
	read_ss(vfeature,abspath+"/"+id+".ss2");
	read_iupred(vfeature,abspath+"/"+id+".long",0);
	read_iupred(vfeature,abspath+"/"+id+".short",1);
	
	getpfun(vfeature);
	
	generateAllFeatures(vsample,vfeature);

//	OutputFeatures(vfeature);
	OutputSamples(vsample,abspath+"/"+id+".score");
	logReg(vsample);
	OutputPred(vsample,abspath+"/"+id+".log.pred");
	return 1;
}

double logitModel(double* x)
{
	double result=0;
	result=1+exp(7.681907791611347 + 0.0006812454718442421*x[1] + 0.0005060513300529592*x[2] + 0.0008023570887506658*x[3] + 0.24411500713666315*x[6] + 0.3763641202957227*x[22] + 0.346118811336949*x[112] + 0.8229551115306304*x[233] - 3.815283489641216*x[236] - 0.7455368304741901*x[240] - 0.38105945379014156*x[241] + 0.8119750693990343*x[243] - 7.279572075805828*x[247] - 0.17919425885684406*x[249] + 0.5900489573772846*x[254] + 0.17637759372428125*x[260] - 0.3456384769371399*x[262] - 0.24563612161717147*x[266] - 0.10415645976012089*x[269] - 0.14294363106209418*x[270] - 0.2253670755996932*x[271] + 0.13843900321751187*x[274] + 3.5005469734289067*x[279] - 0.1932111305613408*x[281] - 13.652167878006654*x[283] - 0.3893428527431058*x[284] + 7.342381721490256*x[286] + 1.6296885761145483*x[287] - 11.353653394696403*x[288] + 4.168089616881257*x[292] - 1.5766839694759467*x[293] + 0.8309957086930082*x[301] + 0.21305413049206307*x[305] - 1.1686052587461102*x[311] + 0.6496420259831437*x[312] - 0.8803314290419961*x[313] + 0.07673940600630128*x[315]);
	result=1/result;
	return result;
}

void convert(double* x,SAMPLE sam,int len)
{
	int id=sam.single[2].idx;
	x[1]=id<(len+1-id)?id:(len+1-id);
	x[2]=len+1-id;
	x[3]=id;
	int i=0;
	for(i=0;i<5;i++)
	{
		x[5+i*45]=sam.single[i].pln;
		x[6+i*45]=sam.single[i].plb;
		x[7+i*45]=sam.single[i].prn;
		x[8+i*45]=sam.single[i].prb;
		x[9+i*45]=sam.single[i].pdn;
		x[10+i*45]=sam.single[i].pdb;
		x[11+i*45]=sam.single[i].ppn;
		x[12+i*45]=sam.single[i].ppb;
		x[13+i*45]=sam.single[i].pmn;
		x[14+i*45]=sam.single[i].pmb;
		x[15+i*45]=sam.single[i].pcon;
		x[16+i*45]=sam.single[i].pdisln;
		x[17+i*45]=sam.single[i].pdislb;
		x[18+i*45]=sam.single[i].pdissn;
		x[19+i*45]=sam.single[i].pdissb;
		x[20+i*45]=sam.single[i].pfunn;
		x[21+i*45]=sam.single[i].pfunb;
		x[22+i*45]=sam.single[i].pfun2n;
		x[23+i*45]=sam.single[i].pfun2b;
		for(int j=0;j<3;j++)
			x[24+j+i*45]=sam.single[i].pss[j];
		for(int j=0;j<3;j++)
			x[27+j+i*45]=sam.single[i].pssb[j];
		for(int j=0;j<20;j++)
			x[30+j+i*45]=sam.single[i].pssm[j];
	}
	i=5;
	{
		x[5+i*45]=sam.local.pln;
		x[6+i*45]=sam.local.plb;
		x[7+i*45]=sam.local.prn;
		x[8+i*45]=sam.local.prb;
		x[9+i*45]=sam.local.pdn;
		x[10+i*45]=sam.local.pdb;
		x[11+i*45]=sam.local.ppn;
		x[12+i*45]=sam.local.ppb;
		x[13+i*45]=sam.local.pmn;
		x[14+i*45]=sam.local.pmb;
		x[15+i*45]=sam.local.pcon;
		x[16+i*45]=sam.local.pdisln;
		x[17+i*45]=sam.local.pdislb;
		x[18+i*45]=sam.local.pdissn;
		x[19+i*45]=sam.local.pdissb;
		x[20+i*45]=sam.local.pfunn;
		x[21+i*45]=sam.local.pfunb;
		x[22+i*45]=sam.local.pfun2n;
		x[23+i*45]=sam.local.pfun2b;
		for(int j=0;j<3;j++)
			x[24+j+i*45]=sam.local.pss[j];
		for(int j=0;j<3;j++)
			x[27+j+i*45]=sam.local.pssb[j];
		for(int j=0;j<20;j++)
			x[30+j+i*45]=sam.local.pssm[j];
	}
	i=6;
	{
		x[5+i*45]=sam.global.pln;
		x[6+i*45]=sam.global.plb;
		x[7+i*45]=sam.global.prn;
		x[8+i*45]=sam.global.prb;
		x[9+i*45]=sam.global.pdn;
		x[10+i*45]=sam.global.pdb;
		x[11+i*45]=sam.global.ppn;
		x[12+i*45]=sam.global.ppb;
		x[13+i*45]=sam.global.pmn;
		x[14+i*45]=sam.global.pmb;
		x[15+i*45]=sam.global.pcon;
		x[16+i*45]=sam.global.pdisln;
		x[17+i*45]=sam.global.pdislb;
		x[18+i*45]=sam.global.pdissn;
		x[19+i*45]=sam.global.pdissb;
		x[20+i*45]=sam.global.pfunn;
		x[21+i*45]=sam.global.pfunb;
		x[22+i*45]=sam.global.pfun2n;
		x[23+i*45]=sam.global.pfun2b;
		for(int j=0;j<3;j++)
			x[24+j+i*45]=sam.global.pss[j];
		for(int j=0;j<3;j++)
			x[27+j+i*45]=sam.global.pssb[j];
		for(int j=0;j<20;j++)
			x[30+j+i*45]=sam.global.pssm[j];
	}	
}
void logReg(vector<SAMPLE>& vsample)
{
	int n=vsample.size();
	double data[400]={0};
	
	for(int i=0;i<n;i++)
	{
		convert(data,vsample[i],n);
		//cout<<data[1]<<"\t"<<data[2]<<"\t"<<data[3]<<"\t"<<data[319]<<endl;
		vsample[i].pred=logitModel(data);	
		////cout<<vsample[i].pred<<endl;	
	}
}

void OutputPred(vector<SAMPLE>& vsample,string filename)
{
	ofstream out(filename.c_str());
	int n=vsample.size();
	out.setf(ios::fixed);
	for(int i=0;i<n;i++)
	{
		out<<vsample[i].single[2].idx<<"\t"<<vsample[i].single[2].res<<"\t"<<vsample[i].pred<<"\t";
		if(vsample[i].pred>=THRESHOLD)
			out<<"1"<<endl;
		else
			out<<"0"<<endl;
	}
}

bool generateAllFeatures(vector<SAMPLE>& vsample,vector<FEATURE>& vf)
{
	int n=vf.size();
	int halfwin=floor(localwindowsize/2);
	int long_halfwin=floor(longwindowsize/2);
	SAMPLE msample;
	FEATURE_SET glo;
	FEATURE_SET loc;

	vsample.clear();
	
	for(int i=0;i<n;i++)
		glo.add(vf[i]);
	glo.divid((double)n);
//	cout<<glo<<endl;
	for(int i=0;i<n;i++)
	{
		for(int j=-halfwin;j<halfwin+1;j++)
		{
			if(i+j<0)
				msample.single[j+halfwin]=vf[0];
			else if(i+j>=n)
				msample.single[j+halfwin]=vf[n-1];
			else
				msample.single[j+halfwin]=vf[i+j];
		}
		int num=0;
		FEATURE_SET loc;
		for(int j=-long_halfwin;j<long_halfwin+1;j++)
		{
			if(i+j<0||i+j>=n)
				continue;
			loc.add(vf[i+j]);
			num++;	
		}
		loc.divid((double)num);
	//	cout<<loc<<endl;
		msample.local=loc;
		msample.global=glo;
		vsample.push_back(msample);
	}
	return 1;
}

void OutputFeatures(vector<FEATURE>& vf)
{
	int n=vf.size();
	for(int i=0;i<n;i++)
	{
		cout<<vf[i].idx<<"\t"<<vf[i].res<<"\t";
		cout<<vf[i];
		cout<<endl;
	}
}

void OutputSamples(vector<SAMPLE>& vsample,string filename,bool label_flag)
{
	ofstream out(filename.c_str());
	int n=vsample.size();
	out.setf(ios::fixed);
//	out.precision(8);
//	out.width(12);
	for(int i=0;i<n;i++)
	{
		out<<vsample[i].single[2].idx<<"\t"<<vsample[i].single[2].res<<"\t";
		out<<vsample[i];
		if(label_flag==true)
			out<<vsample[i].label;
		out<<endl;
	}
}

bool read_dfl(vector<FEATURE>& vf,string filename,bool clear)
{
	ifstream in(filename.c_str());
	char tmp[maxline];
	char seq[maxline];
	double score[maxline]={0};
	int pred[maxline]={0};
	in.getline(tmp,maxline);
	in.getline(seq,maxline);
	int seqn=strlen(seq);
	for(int i=0;i<seqn;i++)
		in>>score[i];
	
	for(int i=0;i<seqn;i++)
	{
		if(seq[i]>='A'&&seq[i]<='Z')
	       		pred[i]=1;
	       	else
	       		pred[i]=0;
	}
	if(clear==true)
	{
		vf.clear();
		FEATURE node;
		for(int i=0;i<seqn;i++)
		{
			node.idx=i+1;
			node.res=toupper(seq[i]);
			node.pln=score[i];
			node.plb=pred[i];
			vf.push_back(node);
		}
	}
	else
	{
		if(vf.size()!=seqn)
		{
			cout<<filename<<" length is different"<<endl;
			return false;
		}
		for(int i=0;i<seqn;i++)
		{
			vf[i].pln=score[i];
			vf[i].plb=pred[i];
		}
	}
	return 1;
}

bool read_rdpbind(vector<FEATURE>& vf,string filename,bool clear)
{
	ifstream in(filename.c_str());
	char tmp[maxline];
	char seq[maxline];
	double scorer[maxline]={0};
	double scored[maxline]={0};
	double scorep[maxline]={0};
	char predc[maxline]={0};
	int predr[maxline]={0};
	int predd[maxline]={0};
	int predp[maxline]={0};

	in.getline(tmp,maxline);
	in.getline(seq,maxline);
	in.getline(predc,maxline);
	int seqn=strlen(seq);
	for(int i=0;i<seqn;i++)
	{
		if(predc[i]=='1')
			predr[i]=1;
		else
			predr[i]=0;
	}
	for(int i=0;i<seqn;i++)
		in>>scorer[i];

	in.getline(tmp,100);//garantee change line;
	in.getline(predc,maxline);
	for(int i=0;i<seqn;i++)
	{
		if(predc[i]=='1')
			predd[i]=1;
		else
			predd[i]=0;
	}
	for(int i=0;i<seqn;i++)
		in>>scored[i];
	
	in.getline(tmp,100);//garantee change line;
	in.getline(predc,maxline);
	for(int i=0;i<seqn;i++)
	{
		if(predc[i]=='1')
			predp[i]=1;
		else
			predp[i]=0;
	}
	for(int i=0;i<seqn;i++)
		in>>scorep[i];

	if(clear==true)
	{
		vf.clear();
		FEATURE node;
		for(int i=0;i<seqn;i++)
		{
			node.idx=i+1;
			node.res=toupper(seq[i]);
			node.prn=scorer[i];
			node.prb=predr[i];
			node.pdn=scored[i];
			node.pdb=predd[i];
			node.ppn=scorep[i];
			node.ppb=predp[i];
			vf.push_back(node);
		}
	}
	else
	{
		if(vf.size()!=seqn)
		{
			cout<<filename<<" length is different"<<endl;
			return false;
		}
		for(int i=0;i<seqn;i++)
		{
			vf[i].prn=scorer[i];
			vf[i].prb=predr[i];
			vf[i].pdn=scored[i];
			vf[i].pdb=predd[i];
			vf[i].ppn=scorep[i];
			vf[i].ppb=predp[i];
		}
	}
	return 1;
}

bool read_fmorf(vector<FEATURE>& vf,string filename,bool clear)
{
	ifstream in(filename.c_str());
	char tmp[maxline];
	char seq[maxline];
	double score[maxline]={0};
	char predc[maxline];
	int pred[maxline]={0};
	in.getline(tmp,maxline);
	in.getline(seq,maxline);
	int seqn=strlen(seq);
	for(int i=0;i<seqn;i++)
		in>>score[i];
	in.getline(tmp,100);
	in.getline(predc,maxline);
	in.getline(tmp,maxline);
		
	for(int i=0;i<seqn;i++)
	{
		if(predc[i]=='1')
	       		pred[i]=1;
	       	else
	       		pred[i]=0;
	}

	if(clear==true)
	{
		vf.clear();
		FEATURE node;
		for(int i=0;i<seqn;i++)
		{
			node.idx=i+1;
			node.res=toupper(seq[i]);
			node.pmn=score[i];
			node.pmb=pred[i];
			vf.push_back(node);
		}
	}
	else
	{
		if(vf.size()!=seqn)
		{
			cout<<filename<<" length is different"<<endl;
			return false;
		}
		for(int i=0;i<seqn;i++)
		{
			vf[i].pmn=score[i];
			vf[i].pmb=pred[i];
		}
	}
	return 1;
}

bool read_pssm(vector<FEATURE>& vf,string filename,bool clear)
{
        ifstream in(filename.c_str());

	double pssm[20]={0};
	int id;
	char res;
	double freq[20]={0};
	double a,b;

	if(clear==true)
		vf.clear();

	while(!in.eof())
        {
                in>>id;
                if(in.eof())
                        break;
                in>>res;
                for(int i=0;i<20;i++)
                        in>>pssm[i];
                for(int i=0;i<20;i++)
                        in>>freq[i];
                in>>a>>b;
		if(clear==true)
		{
			FEATURE node;
			node.idx=id;
			node.res=res;
			for(int i=0;i<20;i++)
			{
				node.pssm[i]=pssm[i];
			}
			node.pcon=conServation(freq);
			vf.push_back(node);
		}
		else
		{
			vf[id-1].idx=id;
			vf[id-1].res=res;
			for(int i=0;i<20;i++)
			{
				vf[id-1].pssm[i]=pssm[i];
			}
			vf[id-1].pcon=conServation(freq);
		}
        }
	return 1;
}
double conServation(double* freq)
{
        double bkf[20]={0.08259245,0.05537017,0.04060056,0.05462622,0.01380998,0.03932678,0.06737165,0.07079676,0.02275314,0.05927524,0.09657261,0.05819463,0.02415950,0.03865212,0.04730270,0.06622815,0.05353031,0.01097168,0.02919945,0.06866589};
        double conS=0;
        double sumf=0;
        for(int i=0;i<20;i++)
        {
                sumf+=freq[i];
        }
        for(int i=0;i<20;i++)
        {
                if(freq[i]!=0&sumf!=0)
                {
                        double p_i=freq[i]/sumf;
                        conS+=p_i*log(p_i/bkf[i])/log(20);
                }
        }
        return conS;
}

bool read_ss(vector<FEATURE>& vf,string filename,bool clear)
{
	ifstream in(filename.c_str());
	//skip the first two lines
	char tmp[maxline]={0};
	in.getline(tmp,maxline);
	in.getline(tmp,maxline);
	char pred;
	FEATURE node;	
	
	if(clear==true)
		vf.clear();

	int id=0;
	while(!in.eof())
	{
		in>>node.idx>>node.res>>pred>>node.pss[0]>>node.pss[1]>>node.pss[2];
		if(in.eof())
			break;
		for(int i=0;i<3;i++)
			node.pssb[i]=0;
		if(pred=='H')
			node.pssb[1]=1;
		else if(pred=='E')
			node.pssb[2]=1;
		else//'C'
			node.pssb[0]=1;
		if(clear==true)
			vf.push_back(node);
		else
		{
			for(int i=0;i<3;i++)
			{
				vf[id].pss[i]=node.pss[i];
				vf[id].pssb[i]=node.pssb[i];
			}
			id++;
		}
	}
	return 1;
}

//flag=0:long mode;flag=1:short mode
bool read_iupred(vector<FEATURE>& vf,string filename,int flag,bool clear)
{
	ifstream in(filename.c_str());
	double score;
	int id=0;
	if(clear==true)
		vf.clear();
	FEATURE node;

	while(!in.eof())
	{
		in>>node.idx>>node.res>>score;
		if(in.eof())
			break;
		if(flag==0)
		{
			if(clear==true)
			{
				node.pdisln=score;
				if(score>=0.5)
					node.pdislb=1;
				else
					node.pdislb=0;
				vf.push_back(node);
			}
			else
			{	
				vf[id].pdisln=score;
				if(score>=0.5)
					vf[id].pdislb=1;
				else
					vf[id].pdislb=0;
			}	
		}
		else if(flag==1)
		{
			if(clear==true)
			{
				node.pdissn=score;
				if(score>=0.5)
					node.pdissb=1;
				else
					node.pdissb=0;
				vf.push_back(node);
			}
			else
			{	
				vf[id].pdissn=score;
				if(score>=0.5)
					vf[id].pdissb=1;
				else
					vf[id].pdissb=0;
			}	
		}
		id++;
	}	
	return 1;
}

bool getpfun(vector<FEATURE>& vf)
{
	int nn=vf.size();
	double scores[5];
	double mean[5];
	double min[5];
	double max[5];
	mean[0]=meanpl,mean[1]=meanpr,mean[2]=meanpd,mean[3]=meanpp,mean[4]=meanpm;
	min[0]=minpl,min[1]=minpr,min[2]=minpd,min[3]=minpp,min[4]=minpm;
	max[0]=maxpl,max[1]=maxpr,max[2]=maxpd,max[3]=maxpp,max[4]=maxpm;
	
	for(int i=0;i<nn;i++)
        {
               	scores[0]=vf[i].pln;
               	scores[1]=vf[i].prn;
               	scores[2]=vf[i].pdn;
               	scores[3]=vf[i].ppn;
               	scores[4]=vf[i].pmn;
		
		for(int j=0;j<5;j++)
	       	{
	       		 scores[j]-=mean[j];
               		 if(scores[j]<0)
               		         scores[j]=0.5*(scores[j]+mean[j]-min[j])/(mean[j]-min[j]);
               		 else
               		         scores[j]=0.5*(scores[j]+max[j]-mean[j])/(max[j]-mean[j]);
		}
		double maxval=-1;
		double secondmax=-1;
		int maxid=0;
		for(int j=0;j<5;j++)
		{
			if(scores[j]>maxval)
			{
				maxval=scores[j];
				maxid=j;
			}
		}
		for(int j=0;j<5;j++)
		{
			if(j!=maxid&&scores[j]>secondmax)
				secondmax=scores[j];
		}
		vf[i].pfunn=maxval;
		if(vf[i].pfunn>=pfun_threshold)
			vf[i].pfunb=1;
		else
			vf[i].pfunb=0;
		
		vf[i].pfun2n=(maxval+secondmax)/2;
		if(vf[i].pfun2n>=pfun2_threshold)
			vf[i].pfun2b=1;
		else
			vf[i].pfun2b=0;

	}
}
