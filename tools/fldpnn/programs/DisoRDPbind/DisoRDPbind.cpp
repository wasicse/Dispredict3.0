#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

#include <sys/time.h>
#include <unistd.h>

using namespace std;

#define MAXSEQLEN 10000
#define PN_LEN 64
#define indB -1.0E10
#define infoNum 3
#define eValTh 0.1

void printHelp();
long checkFileFormat(string inFN, string outFN);

void allBindTypePred(string fastaFN, string predFN);
void RNA_pred(string protSeq, string ssPred, int ws, string blastPred, vector<string>&rnaPredRes);
void DNA_pred(string protSeq, int ws, string blastPred, vector<string>&dnaPredRes);
void Prot_pred(string protSeq, string IUpred, int ws, string blastPred, vector<string>&protPredRes);

void formatDB(string blastDB_FD, string blastDB, string outDBflag);
void blast(string curDir, string protFN, string alignFN, string rstDir);
void transferAnnotation(string alignFN, string trFunFN, string protSeq, string *balstPred);
string trasferAnnotationPerAlign(string queryID, string sbjctID, string realAnnFN, fstream &alignFile, string line, string *predAnn);
void getRealAnn(string realAnnFN, string sbjct, int start, string *realAnn);
void transferRealAnnPerAlign(string *alignInfo, int start_S, int end_S, string *realAnn, string *predAnn);
string insersion(string ref, string str, char inChar);
string transferAnnotation_AlignPart_ID(string que, string sbj, string realAnn_alignPart);
string transferAnnotation_AlignPart_SI(string que, string alignPart, string realAnn_alignPart);
void updatePredAnn(int start_Q, int end_Q, string *predAnn_alignPart, string *predAnn);
string str_logicOR(string bin_1, string bin_2);
int is_right(string *alignInfo);
int getAlignStartpoint(string str);
void initilization(string *str, int len, char inilab);

string runpsipred_single(string curDir, string protName);
string readOutSSpred(string ssPredFN);
string IUPredGlob(string curDir, string protName);
string readOutGlobDoamin(string predFN);

void getAAindVal(string aaIND_Flag, double val[20]);
int replaceAAseqwithAAindex(string aaSeq, double aaindex[], double aaToInd[]);
int AAtoNum(char ch);

double avg_aaind(double *real, int start, int end);
double getContent(string claSeq, int start, int end, char lab);
double getAvgDiff(double *real, int preStart, int start, int end, int postEnd);
double my_mean(double x[], int len);
double getContentDiff(string claSeq, int preStart, int start, int end, int postEnd, char lab);

double logisticRegression(double *feat, double *coef, int featNum);

string getRidOfCarryChar(string str);
string getCurDir();
void filePerProt(string protName, string protSeq, string seqFN);
string replace_str(string str, string chFlag);
void removeTempFiles(string curDir, string protName);

string updateProtSeqBasedOnPred(string rawSeq, string allBinary);

int main(int argc, char *argv[])
{
	if(argc!=3)
	{
		cout << "missing input file or output file" << endl;
		printHelp();
		exit(0);
	}
	string inFN = argv[1];
	string outFN = argv[2];

	string stdFN = inFN + ".fa"; // example.fasta.fa
	cout << "check the format of the input file ......" << endl;
	long pNum = checkFileFormat(inFN, stdFN);

	cout << "perform prediction for the input file ......" << endl;
	
	struct timeval start, end;
	long mtime, seconds, useconds;
	gettimeofday(&start, NULL); // record the strat time

	allBindTypePred(stdFN, outFN);
	
	gettimeofday(&end, NULL);
	seconds  = end.tv_sec  - start.tv_sec;
	useconds = end.tv_usec - start.tv_usec; // record the ending time
	mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5; // used time in milliseconds

	system(("rm \""+stdFN+"\"").c_str());

	double sec = mtime*1.0/1000;
	double mi = sec/60;
	double hour = mi/60;

	cout << std::fixed;
	cout << "\non " << argv[1] << "\n" << pNum  << " proteins" << endl;
	if(sec < 1)
		cout << mtime << " milliseconds for SLIDER prediction" << endl;
	else if(mi < 1)
		cout << std::setprecision(1) << sec << " seconds for DisoRDPbind prediction" << endl;
	else if(hour < 1)
		cout << std::setprecision(1) << mi << " minutes for DisoRDPbind prediction" << endl;
	else
		cout << std::setprecision(1) << hour << " hours for DisoRDPbind prediction" << endl;

	return 0;
}

void allBindTypePred(string fastaFN, string predFN)
{
	fstream dbFile, out;

	out.open(predFN.c_str(), ios::out);
	if(out.is_open() == 0)
	{
		cout << "cannot open the file: " << predFN << endl;
		exit(0);
	}
	dbFile.open(fastaFN.c_str(), ios::in);
	if(dbFile.is_open() == 0)
	{
		cout << "cannot open the file: " << fastaFN << endl;
		exit(0);
	}

	string line;
	string protName, protSeq;
	string *blastPred = new string[3];
	string realAnnFD = "data/db/realAnn";
	vector<string>rnaPredRes;
	vector<string>dnaPredRes;
	vector<string>protPredRes;
	string allBinary;

	int pNum = 0;

	string curDir = getCurDir();
	string rstDir = predFN.substr(0, predFN.find_last_of ("/"));
	//formatDB(curDir+"/data/db", "train.fasta", "blastDB");
	//exit(0);

	getline(dbFile, line);
	while(!dbFile.eof())
	{
		if(line.at(0) == '>')
		{
			pNum++;

			protName = line.substr(1);
			// if(protName.length() > PN_LEN)
			// 	protName = protName.substr(0, PN_LEN);
			cout << pNum << ":" << protName << endl;

			out << ">" << protName << endl;
			getline(dbFile, line);
			protSeq = line;
			//out << protSeq << endl;
			if(protSeq.length() >= MAXSEQLEN)
			{
				cout << pNum << "\t" << protName << "\t" << protSeq.length() << "\twarning: >= " << MAXSEQLEN << endl;
				out << "WARNING: cannot predict the proteins with size >=" << MAXSEQLEN << endl;
				getline(dbFile, line);
				continue;
			}

			filePerProt(protName, protSeq, rstDir+"/"+protName+".seq");
			

			allBinary = "";
			for(int i=0; i<protSeq.length(); i++)
				allBinary.push_back('0');
			for(int ti=0; ti<4; ti++)
			{				
				if(ti==0)	{
					blast(curDir, protName+".seq", protName+".align", rstDir);
					transferAnnotation(rstDir+"/"+protName+".align", realAnnFD, protSeq, blastPred);
				}
				if(ti==1)	{
					string ssPred = runpsipred_single(rstDir, protName);//curDir-->rstDir
					RNA_pred(protSeq, ssPred, 55, blastPred[0], rnaPredRes);

					allBinary = str_logicOR(allBinary, rnaPredRes.at(1));
				}
				if(ti==2)	{
					DNA_pred(protSeq, 21, blastPred[1], dnaPredRes); 

					allBinary = str_logicOR(allBinary, dnaPredRes.at(1));
				}
				if(ti==3)	{
					//out << "protein-binding:";
					string IUpred = IUPredGlob(rstDir, protName);//curDir-->rstDir
					Prot_pred(protSeq, IUpred, 33, blastPred[2], protPredRes); 

					allBinary = str_logicOR(allBinary, protPredRes.at(1));
				}
			}
			blastPred[0] = blastPred[1] = blastPred[2] = "";
			
			//removeTempFiles(rstDir, protName);//curDir-->rstDir
		}
				
		//update protein sequence
		if(allBinary.find('1') != -1)
			protSeq = updateProtSeqBasedOnPred(protSeq, allBinary);

		//output prediction results
		out << protSeq << endl;
		out << "RNA-binding residues:" << rnaPredRes.at(1) << endl;
		out << "RNA-binding propensity:" << rnaPredRes.at(0) << endl;
		out << "DNA-binding residues:" << dnaPredRes.at(1) << endl;
		out << "DNA-binding propensity:" << dnaPredRes.at(0) << endl;
		out << "protein-binding residues:" << protPredRes.at(1) << endl;
		out << "protein-binding propensity:" << protPredRes.at(0) << endl;


		getline(dbFile, line);
	}
	delete []blastPred;

	dbFile.close();
	out.close();
}

void formatDB(string blastDB_FD, string blastDB, string outDBflag)
{
	string tempFD = blastDB_FD + "/" + outDBflag;
	string command = "mkdir \"" + tempFD + "\"";
	system(command.c_str());

	command  = "cp \"" + blastDB_FD+"/"+blastDB + "\" \"" + tempFD+"/"+outDBflag+"\"";
	system(command.c_str());

	string formatdbDir = "blast-2.2.24/bin";
	command = "cd \"" + formatdbDir + "\"\n./formatdb -i \"" + tempFD+"/"+outDBflag + "\"  -o T  -t " + outDBflag;	
	system(command.c_str());

	command = "rm \"" + tempFD+"/"+outDBflag+"\"";
	system(command.c_str());
}
void blast(string curDir, string protFN, string alignFN, string rstDir)
{
	string dbFN = curDir + "/data/db/blastDB/blastDB";
	string blastDir = "blast-2.2.24/bin";
	string blastRstDir = "";
	string command = "cd \"" + blastDir + "\"\n./blastpgp -e 1 -b 1000 -j 1 -i \"" + rstDir+"/"+protFN + "\" -d \"" + dbFN + "\" -o \"" + rstDir+"/"+alignFN+"\"";
	system(command.c_str());
}
void transferAnnotation(string alignFN, string realAnnFD, string protSeq, string *balstPred)
{
	fstream alignFile;
	alignFile.open(alignFN.c_str(), ios::in);
	if(alignFile.is_open() == 0)
	{
		cout <<  "cannot open the file " << alignFN << endl;
		exit(0);
	}
	
	string line;
	string sbjctID,	queryID;
	string funAnnFN;

	getline(alignFile, line);	
	while(!alignFile.eof())
	{	
		line = getRidOfCarryChar(line);

		if(line.find("Query= ") == 0)	{
			queryID = line.substr(line.find_last_of(" ")+1);

			getline(alignFile, line);
			line = getRidOfCarryChar(line);
			line = line.substr(line.find("(")+1);
			line = line.substr(0, line.find(" "));
			if(line.find(",") != -1)
				line = replace_str(line, ",");
			int len = atoi(line.c_str());
			if(len != protSeq.length())	{
				cout << "error in protSeq or blast file for " << queryID << "(v.s " << sbjctID << ")" << endl;
				exit(0);
			}
			
			initilization(balstPred, len, '0');
		}

		if(line.find(">") == 0)
		{
			sbjctID = line.substr(1);
			if(sbjctID.find(" ") != -1)
				sbjctID = sbjctID.substr(0, sbjctID.find(" "));

			//cout << queryID << " vs " << sbjctID << endl;

			while(!alignFile.eof())
			{
				funAnnFN = realAnnFD + "/" + sbjctID + ".real.ann";
				line = trasferAnnotationPerAlign(queryID, sbjctID, funAnnFN, alignFile, line, balstPred);
				if(line.find(">")==0 || line.find("Matrix: BLOSUM62")==0)
					break;
				if(line.find("Score =") != -1)
					continue;
			}

			continue;
		}
		
		getline(alignFile, line);
	}
	alignFile.close();
}
void initilization(string *str, int len, char inilab)
{
	string iniVec;
	for(int i=0; i<len; i++)	
	{
		iniVec.push_back(inilab);	
	}

	for(int i=0; i<3; i++)
		str[i] = iniVec;
}

string trasferAnnotationPerAlign(string queryID, string sbjctID, string realAnnFN, fstream &alignFile, string line, string *predAnn)
{
	string ass, beginAss, endAss;	
	long start_Q = 1000000, start_S = 1000000;
	long end_Q = 0, end_S = 0;
	double evalue = 0;
	int numAlgin = 0;
	string *alignInfo = new string[infoNum];//0->query; 1->pairAlign; 2->sbjct;

	if(line.find(">") == 0)
		getline(alignFile, line);

	while(!alignFile.eof())
	{		
		if(line.find(">")==0 || line.find("Matrix: BLOSUM62")==0)			
			break;
		if(line.find("Score = ") != -1)	{
			numAlgin++;
			if(numAlgin > 1)	break;
		}			
		if(line.find("Expect = ") != -1)	{
			line = line.substr(line.find("Expect = ")+9);
			line = line.substr(0, line.find(","));
			evalue = atof(line.c_str());
			if(evalue > eValTh)
			{
				for(;;)
				{
					if(line.find(">")==0 || line.find("Matrix: BLOSUM62")==0 || line.find("Score = ")!=-1)
						break;
					getline(alignFile, line);
				}
				break;
			}
		}
		if(line.find("Query:") == 0)	{	
			long a = 6 + getAlignStartpoint(line.substr(6));
			long b = a + line.substr(a).find(" ");

			//0->query; 1->pairAlign; 2->sbjct;
			for(int i=0; i<infoNum; i++)
			{	
				beginAss = replace_str(line.substr(6, a-6), " ");
				endAss = line.substr(line.find_last_of(" ")+1);
				ass = line.substr(a, b-a);
				if(i==0) 
				{
					if(start_Q > atoi(beginAss.c_str()))
						start_Q = atoi(beginAss.c_str());
					if(end_Q < atoi(endAss.c_str()))
						end_Q = atoi(endAss.c_str());					
				}
				if(i==2) 
				{
					if(start_S > atoi(beginAss.c_str()))
						start_S = atoi(beginAss.c_str());
					if(end_S < atoi(endAss.c_str()))
						end_S = atoi(endAss.c_str());
				}
				alignInfo[i] += ass; // information for alignments;

				getline(alignFile, line);
			}

			continue;
		}

		getline(alignFile, line);
	}	
	
	if(is_right(alignInfo) == 0)	
	{
		cout << "error: reading out alignment for pair of proteins: Que:" <<  queryID << " v.s sbj:" << sbjctID << endl;		
		exit(0);
	}

	if(alignInfo[0].length() > 0)
	{
		string *realAnn = new string[3];
		getRealAnn(realAnnFN, alignInfo[2], start_S, realAnn);

		string *predAnn_alignPart = new string[3];
		transferRealAnnPerAlign(alignInfo, start_S, end_S, realAnn, predAnn_alignPart);
		updatePredAnn(start_Q, end_Q, predAnn_alignPart, predAnn);

		delete []predAnn_alignPart;
		delete[]realAnn;
	}	
	delete []alignInfo;

	return line;
}
void getRealAnn(string realAnnFN, string sbjct, int start, string *realAnn)
{
	fstream annFile;
	annFile.open(realAnnFN.c_str(), ios::in);
	if(annFile.is_open() == 0)
	{
		cout <<  "cannot open the file " << realAnnFN << endl;
		exit(0);
	}

	string line;
	string protName, ass;
	int diff;
	int countLine = 0;

	getline(annFile, line);
	while(!annFile.eof())
	{
		line = getRidOfCarryChar(line);
		if(countLine == 1) //sbjSeq
		{
			ass = replace_str(sbjct, "-");
			diff = line.substr(start-1).find(ass);
			if(diff != 0)
			{
				cout << "error" << endl;
				exit(0);
			}	
		}
		if(countLine > 2)
		{
			if(countLine-3 >= 3)
			{
				cout << "error in read out real functional annotations" << endl;
				exit(0);
			}
			realAnn[countLine-3] = line;
		}

		countLine++;
		getline(annFile, line);
	}
	annFile.close();
}
void transferRealAnnPerAlign(string *alignInfo, int start_S, int end_S, string *realAnn, string *predAnn_alignPart)
{
	string *realAnn_alignPart = new string[3];
	string query = alignInfo[0];
	string alignPart = alignInfo[1];
	string sbjct = alignInfo[2];
	string ass;

	for(int i=0; i<3; i++)
	{
		ass = realAnn[i].substr(realAnn[i].find(":")+1);
		realAnn_alignPart[i] = ass.substr(start_S-1, end_S-start_S+1);

		//sbjct = "----" + sbjct.substr(0, 10)+"-----"+sbjct.substr(10, sbjct.length()-10)+"-----------";
		if(sbjct.find("-") != -1) // means the alignment has some insersion for sbjct, for these transfered annotation should be negative '0';
			realAnn_alignPart[i] = insersion(sbjct, realAnn_alignPart[i], '-');

		//predAnn_alignPart[i] = transferAnnotation_AlignPart_ID(query, sbjct, realAnn_alignPart[i]);
		predAnn_alignPart[i] = transferAnnotation_AlignPart_SI(query, alignPart, realAnn_alignPart[i]);
	}

	delete []realAnn_alignPart;
}
string insersion(string ref, string str, char inChar)
{
	string ass;
	for(;;)
	{
		if(ref.find(inChar) == -1)
			break;
		ass = str.substr(ref.find(inChar));
		str = str.substr(0, ref.find(inChar));
		str.push_back('0');
		str = str + ass;

		ref.replace(ref.find(inChar), 1, "*");
	}

	return str;
}
string transferAnnotation_AlignPart_ID(string que, string sbj, string realAnn_alignPart)
{
	if(que.length()!=sbj.length() || realAnn_alignPart.length()!=sbj.length())
	{
		cout << "error inputs for  transferAnnotation_AlignPart_ID" << endl;
		exit(0);
	}
	int len = que.length();
	string predAnn_alignPart;

	for(int i=0; i<len; i++)
	{
		if(que.at(i) == '-')
			continue;

		if(que.at(i) == sbj.at(i)) // only transfer annotation for identical residues
			predAnn_alignPart.push_back(realAnn_alignPart.at(i));
		else
			predAnn_alignPart.push_back('0');
	}

	return predAnn_alignPart;
}
string transferAnnotation_AlignPart_SI(string que, string alignPart, string realAnn_alignPart)
{
	if(que.length()!=alignPart.length() || realAnn_alignPart.length()!=alignPart.length())
	{
		cout << "error inputs for  transferAnnotation_AlignPart_ID" << endl;
		exit(0);
	}
	int len = que.length();
	string predAnn_alignPart;

	for(int i=0; i<len; i++)
	{
		if(que.at(i) == '-')
			continue;

		if(alignPart.at(i) != ' ') // transfer annotation for all similar residues
			predAnn_alignPart.push_back(realAnn_alignPart.at(i));
		else
			predAnn_alignPart.push_back('0');
	}

	return predAnn_alignPart;
}
void updatePredAnn(int start_Q, int end_Q, string *predAnn_alignPart, string *predAnn)
{
	for(int i=0; i<3; i++)
	{
		if(predAnn_alignPart[i].length() != end_Q-start_Q+1)
		{
			cout << "error inputs for updatePredAnn" << endl;
			exit(0);
		}
	}
	
	for(int i=0; i<3; i++)
	{
		if(predAnn_alignPart[i].find("1") != -1) // only update the annotation for positives by logical_or
		{
			string preStr = predAnn[i].substr(0, start_Q-1);
			string postStr = predAnn[i].substr(end_Q);
			string logicOr = str_logicOR(predAnn[i].substr(start_Q-1, end_Q-start_Q+1), predAnn_alignPart[i]);
			predAnn[i] = preStr + logicOr + postStr;
		}
	}
}
string str_logicOR(string bin_1, string bin_2)
{
	if(bin_1.length()==0 || bin_2.length()==0)
	{
		cout << "error: (str_logicOR) invalid inputs" << endl;
		exit(0);
	}		
	if(bin_1.length() != bin_2.length())
	{
		cout << "error in binary anntation" << endl;
		exit(0);
	}

	string logicSum;
	if(bin_1.find('1')==-1 || bin_2.find('1')==-1)
	{
		if(bin_1.find('1') == -1)
			logicSum = bin_2;
		else 
			logicSum = bin_1;
	}
	else
	{
		int len = bin_1.length();
		char ch;
		for(int i=0; i<len; i++)
		{
			if(bin_1.at(i)=='1' || bin_2.at(i)=='1')
				ch = '1';
			else
				ch = '0';

			logicSum.push_back(ch);
		}
	}

	return logicSum;
}
int getAlignStartpoint(string str)
{
	int i = 0;
	for(i=0; i<str.length(); i++)
	{
		if(str.at(i)>='A' && str.at(i)<='Z' || str.at(i) == '-')
			break;
	}

	return i;
}
int is_right(string *alignInfo)
{
	int flag = 1;
	int alingLen = alignInfo[0].length();
	for(int i=0; i<infoNum; i++)
	{
		if(alignInfo[i].length() != alingLen)
		{
			flag = 0;
			break;
		}
	}
	return flag;
}

string runpsipred_single(string curDir, string protName)
{
	string sspredDir = "psipred";
	string command = "cd \"" + sspredDir + "\"\n./runpsipred_single \"" +  curDir+"/"+protName+".seq\" >\""+curDir+"/"+protName+".screenOUT\"";
	system(command.c_str());

	string ssPred = readOutSSpred(sspredDir + "/" + protName + ".horiz");
        command = "mv ./" + sspredDir+ "/"+ protName + ".* \"" + curDir + "/\"";
	system(command.c_str()); 	
	return ssPred;
}
string readOutSSpred(string ssPredFN)
{
	fstream predFile;
	predFile.open(ssPredFN.c_str(), ios::in);
	if(predFile.is_open() == 0)
	{
		cout << "cannot open the file: " << ssPredFN << endl;
		exit(0);
	}

	string line, pred;
	
	getline(predFile, line);
	while(!predFile.eof())
	{
		if(line.find("Pred: ") == 0)
			pred = pred + line.substr(line.find("Pred: ")+6);

		getline(predFile, line);
	}
	predFile.close();
	return pred;
}

void RNA_pred(string protSeq, string ssPred, int ws, string blastPred, vector<string>&rnaPredRes)
{
	if(protSeq.length()!=ssPred.length() || protSeq.length()!=blastPred.length())
	{
		cout << "error inputs for RNA_pred" << endl;
		cout << protSeq << endl;
		cout << ssPred << endl;
		cout << blastPred << endl;
		exit(0);
	}

	int protLen = protSeq.length();

	string aaindTag[9] = {"AURR980103", "AURR980120", "BUNA790102", "CHOP780206", "WILM950103", "CHOP780215", "KLEP840101", "NAKH900113", "QIAN880113"};
	double coef[12] = {-17.4107, 4.8655, 2.2867, -12.0037, -0.7179, 35.253, -3.0753, 0.9816, 5.7314, 1.2239, -1.3558, 11.6733};
	double *aaIND = new double[20];
	double *aaToInd = new double[protLen];
	double *feat = new double[11];	
	int start, preStart, end, postEnd;
	int featNum = 0;
	char *ss = new char[10];

	rnaPredRes.clear();
	rnaPredRes.push_back(""); //for propensity scores
	rnaPredRes.push_back(""); //for rna-binding residues
	for(int i=0; i<protLen; i++)
	{
		start = 0 > (i-ws/2) ? 0 : (i-ws/2);
		end = (i+ws/2) < protLen-1 ? (i+ws/2) : protLen-1;
		preStart = 0 > (start-ws/2) ? 0 : (start-ws/2);
		postEnd = (end+ws/2) < protLen-1 ? (end+ws/2) : protLen-1;

		featNum = 0;
		for(int ind=0; ind<9; ind++)
		{
			getAAindVal(aaindTag[ind], aaIND);
			int n = replaceAAseqwithAAindex(protSeq, aaIND, aaToInd);
			if(protLen != n)
			{
				cout << "error in function: replaceAAseqwithAAindex " << endl;
				exit(0);
			}

			if(ind<5) // avg
			{
				feat[featNum] = avg_aaind(aaToInd, start, end);
				featNum++;
			}
			else //diff
			{
				feat[featNum] = getAvgDiff(aaToInd, preStart, start, end, postEnd);
				featNum++;
			}
		}// feature based on AA indices

		feat[featNum] = getContentDiff(protSeq, preStart, start, end, postEnd, 'G'); // diff_comp_G;
		featNum++;

		//diff_con_H
		feat[featNum] = getContentDiff(ssPred, preStart, start, end, postEnd, 'H');
		featNum++;

		if(featNum!= 11)
		{
			cout << "error in RNA_pred" << endl;
			exit(0);
		}	

		double logProb = logisticRegression(feat, coef, featNum); //logProb[i]
		double finalProb = 0;
		if(blastPred.at(i)=='1')
			finalProb = (1+logProb)/2; //finalProb[i]
		else
			finalProb = logProb; //finalProb[i]
		//out << finalProb << ",";
		sprintf(ss, "%.3f,", finalProb);
		rnaPredRes.at(0) = rnaPredRes.at(0) + ss;

		//cutoff 0.104 updated on 2015/10/13
		if(finalProb>=0.151)
			rnaPredRes.at(1).push_back('1');
		else
			rnaPredRes.at(1).push_back('0');
	}
	//out << endl;

	delete []aaIND;
	delete []aaToInd;
	delete []feat;
	delete []ss;
}
void DNA_pred(string protSeq, int ws, string blastPred, vector<string>&dnaPredRes)
{
	if(protSeq.length()!=blastPred.length())
	{
		cout << "error inputs for DNA_pred" << endl;
		cout << protSeq << endl;
		cout << blastPred << endl;
		exit(0);
	}

	int protLen = protSeq.length();

	string aaindTag[5] = {"KLEP840101", "QIAN880139", "RACS820103", "ROSM880103", "WERD780103"};
	double coef[8] = {-0.2713, -3.0389, -0.7549, -3.4309, -3.4163, 4.6694, -10.1431, 1.2784};
	double *aaIND = new double[20];
	double *aaToInd = new double[protLen];
	double *feat = new double[11];	
	int start, end;
	int featNum = 0;
	char *ss = new char[10];

	//string binary = "";
	//out << "DNA-binding propensity:";
	dnaPredRes.clear();
	dnaPredRes.push_back(""); //for propensity scores
	dnaPredRes.push_back(""); //for dna-binding residues
	for(int i=0; i<protLen; i++)
	{
		start = 0 > (i-ws/2) ? 0 : (i-ws/2);
		end = (i+ws/2) < protLen-1 ? (i+ws/2) : protLen-1;

		featNum = 0;
		for(int ind=0; ind<5; ind++)
		{
			getAAindVal(aaindTag[ind], aaIND);
			int n = replaceAAseqwithAAindex(protSeq, aaIND, aaToInd);
			if(protLen != n)
			{
				cout << "error in function: replaceAAseqwithAAindex " << endl;
				exit(0);
			}
			
			feat[featNum] = avg_aaind(aaToInd, start, end);
			featNum++;
		}// feature based on AA indices

		feat[featNum] = getContent(protSeq, start, end, 'K'); // compostition of K
		featNum++;

		feat[featNum] = getContent(protSeq, start, end, 'W'); // compostition of W
		featNum++;

		if(featNum!= 7)
		{
			cout << "error in DNA_pred" << endl;
			exit(0);
		}	

		double logProb = logisticRegression(feat, coef, featNum); //logProb[i]
		double finalProb = 0;
		if(blastPred.at(i)=='1')
			finalProb = (1+logProb)/2; //finalProb[i]
		else
			finalProb = logProb; //finalProb[i]
		//out << finalProb << ",";
		sprintf(ss, "%.3f,", finalProb);
		dnaPredRes.at(0) = dnaPredRes.at(0) + ss;

		if(finalProb>=0.245)
			dnaPredRes.at(1).push_back('1');
		else
			dnaPredRes.at(1).push_back('0');
	}
	//out << endl;
	//out << "DNA-binding residues:" << binary << endl;

	delete []aaIND;
	delete []aaToInd;
	delete []feat;
	delete []ss;
}
void Prot_pred(string protSeq, string IUpred, int ws, string blastPred, vector<string>&protPredRes)
{
	if(protSeq.length()!=IUpred.length() || protSeq.length()!=blastPred.length())
	{
		cout << "error inputs for RNA_pred" << endl;
		cout << protSeq << endl;
		cout << IUpred << endl;
		cout << blastPred << endl;
		exit(0);
	}

	int protLen = protSeq.length();
	
	string	aaindTag[5] = {"KOEP990102", "PALJ810108", "QIAN880113", "SNEP660103", "DAWD720101"};
	double coef[8] = {1.3314, 1.3696, -5.1621, -5.42, 6.0624, 0.0875, -1.0194, -1.4813};
	double *aaIND = new double[20];
	double *aaToInd = new double[protLen];
	double *feat = new double[11];	
	int start, preStart, end, postEnd;
	int featNum = 0;
	char *ss = new char[10];

	//string binary = "";
	//out << "protein-binding propensity:";
	protPredRes.clear();
	protPredRes.push_back(""); //for propensity scores
	protPredRes.push_back(""); //for rna-binding residues	
	for(int i=0; i<protLen; i++)
	{
		start = 0 > (i-ws/2) ? 0 : (i-ws/2);
		end = (i+ws/2) < protLen-1 ? (i+ws/2) : protLen-1;
		preStart = 0 > (start-ws/2) ? 0 : (start-ws/2);
		postEnd = (end+ws/2) < protLen-1 ? (end+ws/2) : protLen-1;

		featNum = 0;
		for(int ind=0; ind<5; ind++)
		{
			getAAindVal(aaindTag[ind], aaIND);
			int n = replaceAAseqwithAAindex(protSeq, aaIND, aaToInd);
			if(protLen != n)
			{
				cout << "error in function: replaceAAseqwithAAindex " << endl;
				exit(0);
			}

			if(ind<3) // avg
			{
				feat[featNum] = avg_aaind(aaToInd, start, end);
				featNum++;
			}
			else if(ind == 3) //avg and avg_diff "SNEP660103"
			{
				feat[featNum] = avg_aaind(aaToInd, start, end);
				featNum++;

				feat[featNum] = getAvgDiff(aaToInd, preStart, start, end, postEnd);
				featNum++;
			}
			else //diff
			{
				feat[featNum] = getAvgDiff(aaToInd, preStart, start, end, postEnd);
				featNum++;
			}
		}// feature based on AA indices

		//IUPredG_content
		feat[featNum] = getContent(IUpred, start, end, '1');
		featNum++;

		if(featNum!= 7)
		{
			cout << "error in Prot_pred" << endl;
			exit(0);
		}	

		double logProb = logisticRegression(feat, coef, featNum); //logProb[i]
		double finalProb = 0;
		if(blastPred.at(i)=='1')
			finalProb = (1+logProb)/2; //finalProb[i]
		else
			finalProb = logProb; //finalProb[i]
		//out << finalProb << ",";
		sprintf(ss, "%.3f,", finalProb);
		protPredRes.at(0) = protPredRes.at(0) + ss;

		if(finalProb>=0.799)
			protPredRes.at(1).push_back('1');
		else
			protPredRes.at(1).push_back('0');
	}
	//out << endl;
	//out << "protein-binding residues:" << binary << endl;

	delete []aaIND;
	delete []aaToInd;
	delete []feat;
	delete []ss;
}
string IUPredGlob(string curDir, string protName)
{
	string IUPredDir= "IUPRED";
	string predFN = protName + ".IUpred";
	string command = "cd \"" + IUPredDir+ "\"\n" + "./iupred \"" + curDir+"/"+protName+".seq\" glob >\"" + curDir+"/"+predFN+"\"";
	system(command.c_str());

	string IUpred = readOutGlobDoamin(curDir+"/"+predFN);
	
	return IUpred;
}
string readOutGlobDoamin(string predFN)	
{
	fstream predFile;

	predFile.open(predFN.c_str(), ios::in);
	if(predFile.is_open() == 0)
	{
		cout << "Cannot open the file: " << predFN << endl;
		exit(0);
	}

	string line, seq, pred;
	int diff = 'a'-'A';
	
	pred = "";
	getline(predFile, line);
	while(!predFile.eof())
	{
		if(line.at(0) == '>')
		{			
			getline(predFile, line);
			while(!predFile.eof())
			{
				for(int i=0; i<line.length(); i++)
				{
					if(line.at(i)>='a' && line.at(i)<='z')
					{
						pred = pred + "0";
						seq.push_back(line.at(i) - diff);
					}
					if(line.at(i)>='A' && line.at(i)<='Z')
					{
						pred = pred + "1";
						seq.push_back(line.at(i));
					}
				}			
				getline(predFile, line);
			}
			break;
		}

		getline(predFile, line);
	}
	predFile.close();	

	return pred;
}

double avg_aaind(double *real, int start, int end)
{
	double *x = new double[end-start+1];
	for(int i=start; i<=end; i++)
	{
		if(fabs(real[i]+1.0E5)<=1.0-40)
		{
			cout << "error in passing the parameter real" << endl;
			exit(0);
		}
		x[i-start] = real[i];
	}

	double avg = my_mean(x, end-start+1);
	delete []x;

	return avg;
}
double getContent(string claSeq, int start, int end, char lab)
{
	double content = 0;
	for(int i=start; i<=end; i++)
	{
		if(claSeq.at(i) == lab)
		{
			content++;
		}
	}

	content = content/(end-start+1);
	return content;
}
double getAvgDiff(double *real, int preStart, int start, int end, int postEnd)
{
	double *cenVal = new double[end-start+1];
	double *sidVal = new double[start-preStart+postEnd-end];
	int numSid = 0;
	for(int i=preStart; i<=postEnd; i++)
	{
		if(fabs(real[i]+1.0E5)<=1.0-40)
		{
			cout << "error in passing the parameter real" << endl;
			exit(0);
		}

		if(i>=start && i<=end)
			cenVal[i-start] = real[i];
		else
		{
			sidVal[numSid] = real[i];
			numSid++;
		}
	}
	if(numSid != start-preStart+postEnd-end)
	{
		cout << "error to get the preStart and postEnd" << endl;
		exit(0);
	}

	double avgCen = my_mean(cenVal, end-start+1);

	double avgSid = 0;
	if(numSid > 0)
		avgSid = my_mean(sidVal, numSid);
	double diff = avgCen - avgSid;
	delete []cenVal;
	delete []sidVal;

	return diff;
}
double getContentDiff(string claSeq, int preStart, int start, int end, int postEnd, char lab)
{
	double sidCon = 0;
	double cenCon = 0;

	int numSid = 0;
	for(int i=preStart; i<=postEnd; i++)
	{
		if(claSeq.at(i) == lab)
		{
			if(i>=start && i<=end)
				cenCon++;
			else
				sidCon++;
		}

		if(i<start || i>end)
			numSid++;
	}
	if(numSid != (start-preStart+postEnd-end))
	{
		cout << "error to get the preStart and postEnd" << endl;
		exit(0);
	}

	double diff;
	cenCon = cenCon/(end-start+1);

	if(numSid > 0)
		sidCon = sidCon/numSid;

	diff = cenCon - sidCon;

	return diff;
}

void getAAindVal(string aaIND_Flag, double val[20])
{
	if(aaIND_Flag.compare("AURR980103")==0) { val[0]=1.05; val[1]=0.81; val[2]=0.91; val[3]=1.39; val[4]=0.6; val[5]=0.87; val[6]=1.11; val[7]=1.26; val[8]=1.43; val[9]=0.95; val[10]=0.96; val[11]=0.97; val[12]=0.99; val[13]=0.95; val[14]=1.05; val[15]=0.96; val[16]=1.03; val[17]=1.06; val[18]=0.94; val[19]=0.62; }
	else if(aaIND_Flag.compare("AURR980120")==0) { val[0]=0.71; val[1]=1.09; val[2]=0.95; val[3]=1.43; val[4]=0.65; val[5]=0.87; val[6]=1.19; val[7]=1.07; val[8]=1.13; val[9]=1.05; val[10]=0.84; val[11]=1.1; val[12]=0.8; val[13]=0.95; val[14]=1.7; val[15]=0.65; val[16]=0.086; val[17]=1.25; val[18]=0.85; val[19]=1.12; }
	else if(aaIND_Flag.compare("BUNA790102")==0) { val[0]=4.349; val[1]=4.396; val[2]=4.755; val[3]=4.765; val[4]=4.686; val[5]=4.373; val[6]=4.295; val[7]=3.972; val[8]=4.63; val[9]=4.224; val[10]=4.385; val[11]=4.358; val[12]=4.513; val[13]=4.663; val[14]=4.471; val[15]=4.498; val[16]=4.346; val[17]=4.702; val[18]=4.604; val[19]=4.184; }
	else if(aaIND_Flag.compare("CHOP780206")==0) { val[0]=0.7; val[1]=0.34; val[2]=1.42; val[3]=0.98; val[4]=0.65; val[5]=0.75; val[6]=1.04; val[7]=1.41; val[8]=1.22; val[9]=0.78; val[10]=0.85; val[11]=1.01; val[12]=0.83; val[13]=0.93; val[14]=1.1; val[15]=1.55; val[16]=1.09; val[17]=0.62; val[18]=0.99; val[19]=0.75; }
	else if(aaIND_Flag.compare("WILM950103")==0) { val[0]=-1.64; val[1]=-3.28; val[2]=0.83; val[3]=0.7; val[4]=9.3; val[5]=-0.04; val[6]=1.18; val[7]=-1.85; val[8]=7.17; val[9]=3.02; val[10]=0.83; val[11]=-2.36; val[12]=4.26; val[13]=-1.36; val[14]=3.12; val[15]=1.59; val[16]=2.31; val[17]=2.61; val[18]=2.37; val[19]=0.52; }
	else if(aaIND_Flag.compare("CHOP780215")==0) { val[0]=0.058; val[1]=0.085; val[2]=0.091; val[3]=0.081; val[4]=0.128; val[5]=0.098; val[6]=0.064; val[7]=0.152; val[8]=0.054; val[9]=0.056; val[10]=0.07; val[11]=0.095; val[12]=0.055; val[13]=0.065; val[14]=0.068; val[15]=0.106; val[16]=0.079; val[17]=0.167; val[18]=0.125; val[19]=0.053; }
	else if(aaIND_Flag.compare("NAKH900113")==0) { val[0]=1.61; val[1]=0.4; val[2]=0.73; val[3]=0.75; val[4]=0.37; val[5]=0.61; val[6]=1.5; val[7]=3.12; val[8]=0.46; val[9]=1.61; val[10]=1.37; val[11]=0.62; val[12]=1.59; val[13]=1.24; val[14]=0.67; val[15]=0.68; val[16]=0.92; val[17]=1.63; val[18]=0.67; val[19]=1.3; }
	else if(aaIND_Flag.compare("QIAN880113")==0) { val[0]=-0.08; val[1]=0.05; val[2]=-0.08; val[3]=-0.24; val[4]=-0.25; val[5]=-0.28; val[6]=-0.19; val[7]=-0.1; val[8]=0.29; val[9]=-0.01; val[10]=0.28; val[11]=0.45; val[12]=0.11; val[13]=0; val[14]=-0.42; val[15]=0.07; val[16]=-0.33; val[17]=0.36; val[18]=0; val[19]=-0.13; }
	else if(aaIND_Flag.compare("KLEP840101")==0) { val[0]=0; val[1]=1; val[2]=0; val[3]=-1; val[4]=0; val[5]=0; val[6]=-1; val[7]=0; val[8]=0; val[9]=0; val[10]=0; val[11]=1; val[12]=0; val[13]=0; val[14]=0; val[15]=0; val[16]=0; val[17]=0; val[18]=0; val[19]=0; }
	else if(aaIND_Flag.compare("QIAN880139")==0) { val[0]=0.08; val[1]=-0.01; val[2]=-0.06; val[3]=0.04; val[4]=0.37; val[5]=0.48; val[6]=0.36; val[7]=-0.02; val[8]=-0.45; val[9]=0.09; val[10]=0.24; val[11]=-0.27; val[12]=0.16; val[13]=0.34; val[14]=0.16; val[15]=-0.35; val[16]=-0.04; val[17]=-0.06; val[18]=-0.2; val[19]=0.18; }
	else if(aaIND_Flag.compare("RACS820103")==0) { val[0]=0.82; val[1]=2.6; val[2]=2.07; val[3]=2.64; val[4]=0; val[5]=0; val[6]=2.62; val[7]=1.63; val[8]=0; val[9]=2.32; val[10]=0; val[11]=2.86; val[12]=0; val[13]=0; val[14]=0; val[15]=1.23; val[16]=2.48; val[17]=0; val[18]=1.9; val[19]=1.62; }
	else if(aaIND_Flag.compare("ROSM880103")==0) { val[0]=0.4; val[1]=0.3; val[2]=0.9; val[3]=0.8; val[4]=0.5; val[5]=0.7; val[6]=1.3; val[7]=0; val[8]=1; val[9]=0.4; val[10]=0.6; val[11]=0.4; val[12]=0.3; val[13]=0.7; val[14]=0.9; val[15]=0.4; val[16]=0.4; val[17]=0.6; val[18]=1.2; val[19]=0.4; }
	else if(aaIND_Flag.compare("WERD780103")==0) { val[0]=0.15; val[1]=-0.37; val[2]=0.69; val[3]=-0.22; val[4]=-0.19; val[5]=-0.06; val[6]=0.14; val[7]=0.36; val[8]=-0.25; val[9]=0.02; val[10]=0.06; val[11]=-0.16; val[12]=0.11; val[13]=1.18; val[14]=0.11; val[15]=0.13; val[16]=0.28; val[17]=-0.12; val[18]=0.19; val[19]=-0.08; }
	else if(aaIND_Flag.compare("DAWD720101")==0) { val[0]=2.5; val[1]=7.5; val[2]=5; val[3]=2.5; val[4]=3; val[5]=6; val[6]=5; val[7]=0.5; val[8]=6; val[9]=5.5; val[10]=5.5; val[11]=7; val[12]=6; val[13]=6.5; val[14]=5.5; val[15]=3; val[16]=5; val[17]=7; val[18]=7; val[19]=5; }
	else if(aaIND_Flag.compare("KOEP990102")==0) { val[0]=-0.12; val[1]=0.34; val[2]=1.05; val[3]=1.12; val[4]=-0.63; val[5]=1.67; val[6]=0.91; val[7]=0.76; val[8]=1.34; val[9]=-0.77; val[10]=0.15; val[11]=0.29; val[12]=-0.71; val[13]=-0.67; val[14]=0; val[15]=1.45; val[16]=-0.7; val[17]=-0.14; val[18]=-0.49; val[19]=-0.7; }
	else if(aaIND_Flag.compare("PALJ810108")==0) { val[0]=1.34; val[1]=0.91; val[2]=0.83; val[3]=1.06; val[4]=1.27; val[5]=1.13; val[6]=1.69; val[7]=0.47; val[8]=1.11; val[9]=0.84; val[10]=1.39; val[11]=1.08; val[12]=0.9; val[13]=1.02; val[14]=0.48; val[15]=1.05; val[16]=0.74; val[17]=0.64; val[18]=0.73; val[19]=1.18; }
	else if(aaIND_Flag.compare("SNEP660103")==0) { val[0]=-0.11; val[1]=0.079; val[2]=-0.136; val[3]=-0.285; val[4]=-0.184; val[5]=-0.067; val[6]=-0.246; val[7]=-0.073; val[8]=0.32; val[9]=0.001; val[10]=-0.008; val[11]=0.049; val[12]=-0.041; val[13]=0.438; val[14]=-0.016; val[15]=-0.153; val[16]=-0.208; val[17]=0.493; val[18]=0.381; val[19]=-0.155; }
	else {for(int i=0; i<20; i++) val[i]=indB;}
}
int replaceAAseqwithAAindex(string aaSeq, double aaindex[], double aaToInd[])
{
	int i = 0;
	for(i = 0; i<aaSeq.length(); i++)
	{
		int flag = AAtoNum(aaSeq[i]);
		if(flag == -1)
			aaToInd[i] = 0;
		else
			aaToInd[i] = aaindex[flag];
	}

	return i;
}
int AAtoNum(char ch)
{
	//I    A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V 
	string legalAAs = "ARNDCQEGHILKMFPSTWYV"; //ARNDCQEGHILKMFPSTWYV
	int d = -1;

	for(int i=0; i<20; i++)
	{
		if(legalAAs[i] == ch)
		{
			d = i;
			break;
		}
	}

	return d;
}
double my_mean(double x[], int len)
{
	if(len==0)
		return 0;

	double avg = 0;
	int i;
	for(i=0; i<len; i++)
	{
		avg = avg + x[i];
	}
	return avg/len;
}

double logisticRegression(double *feat, double *coef, int featNum)
{
	double sum = 0;
	for(int i=0; i<featNum; i++)
	{
		sum = sum + feat[i]*coef[i];
	}
	sum = sum + coef[featNum];

	double prob = exp(sum)/(1+exp(sum));

	return prob;
}
string getRidOfCarryChar(string str)
{
	if(str.length() > 0)
	{
		if(str.at(str.length()-1) == 13)
			str.erase(str.length()-1);
	}
		
	return str;
}
string getCurDir()
{
  char buff[FILENAME_MAX];
  getcwd( buff, FILENAME_MAX );
  std::string current_working_dir(buff);
//cout << "curDir( from: DisoRDPbind.cpp) " <<  current_working_dir << endl;
  return current_working_dir;
/*
	system("pwd >curDir");

	fstream temp;
	temp.open("curDir", ios::in);
	if(temp.is_open() == 0)
	{
		cout << "Cannot open output file " << "curDir" << endl;
		exit(0);
	}
	string curDir;

	getline(temp, curDir);
	temp.close();

	system("rm curDir");

	return curDir;
*/
}
void filePerProt(string protName, string protSeq, string seqFN)
{
	fstream out;
	out.open(seqFN.c_str(), ios::out);
	if(out.is_open() == 0)
	{
		cout << "cannot open the file: " << seqFN << endl;
		exit(0);
	}
	out << ">" << protName << endl;
	out << protSeq << endl;
	out.close();
}
string replace_str(string str, string chFlag)
{
	while(str.find(chFlag) != -1)
	{
		str.replace(str.find(chFlag), 1, "");
	}

	return str;
}
void removeTempFiles(string curDir, string protName)
{
	string command = "cd \"" + curDir + "\"\nrm " + protName + ".*";
	system(command.c_str());
	
	string sspredDir = "psipred";
	command = "cd \"" + sspredDir + "\"\nrm " + protName + ".*";
	system(command.c_str());
}

void printHelp()
{
	cout << "\n============= HELP =============" << endl;
	cout << "command line: ./DisoRDPbind inputFile outputFile" << endl;
	cout << "where inputFile has to follow fasta format\n" << endl;
}
long checkFileFormat(string inFN, string outFN)
{
	fstream dbFile, out;

	out.open(outFN.c_str(), ios::out);
	if(out.is_open() == 0)
	{
		cout << "cannot open the file: " << outFN << endl;
		exit(0);
	}

	dbFile.open(inFN.c_str(), ios::in);
	if(dbFile.is_open() == 0)
	{
		cout << "cannot open the file: " << inFN << endl;
		exit(0);
	}	

	string line, protName, protSeq;
	string std_AAs = "ARNDCQEGHILKMFPSTWYV";
	char ch[2];
	long pNum = 0;

	getline(dbFile, line);
	while(!dbFile.eof())
	{
		if(line.at(0) == '>')
		{
			pNum++;

			line = getRidOfCarryChar(line);
			protName = line.substr(1);
			// if(protName.length() > PN_LEN)
			// 	protName = protName.substr(0, PN_LEN) + "...";
			
			protSeq = "";
			getline(dbFile, line);
			while(!dbFile.eof() && line.at(0) != '>')
			{
				line = getRidOfCarryChar(line);
				if(line.length() == 0)
				{
					getline(dbFile, line);
					continue;
				}

				for(int i=0; i<line.length(); i++)
				{
					if(line.at(i)>='A' && line.at(i)<='Z')
						ch[0] = line.at(i);
					else if(line.at(i)>='a' && line.at(i)<='z')
						ch[0] = line.at(i) - ('a'-'A');
					else
					{
						cout << "the sequence of the protein " << protName << " comprises the unexpected character " << line.at(i) << endl;						
						exit(0);
					}
					ch[1] = '\0';

					if(std_AAs.find(ch) == -1)
					{
						cout << "warning: the amino acid " << ch[0] << " is replaced by U" << endl;
						ch[0] = 'U';
					}
					
					protSeq = protSeq + ch;
				}

				getline(dbFile, line);
			}			

			//out << ">" << protName << "|" << protSeq.length() << " AAs" << endl;
			out << ">" << protName << endl;
			out << protSeq << endl;
			out.flush();
			continue;
		}
		else
		{
			cout << "the input data has to be fasta format" << endl;
			exit(0);
		}

		getline(dbFile, line);
	}	

	dbFile.close();
	out.close();

	return pNum;
}
string updateProtSeqBasedOnPred(string rawSeq, string allBinary)
{
	if(rawSeq.length()!=allBinary.length())
	{
		cout << "error inputs for function: updateProtSeqBasedOnPred" << endl;
		exit(0);
	}
	
	string newSeq;
	char AA;
	for(int i=0; i<rawSeq.length(); i++)
	{
		AA = rawSeq.at(i);
		if(allBinary.at(i) == '1')
			AA = AA + ('a' - 'A');
		else if(allBinary.at(i) == '0')
			;
		else {
			cout << "error annotations (updateProtSeqBasedOnPred)" << endl;
			exit(0);
		}

		newSeq.push_back(AA);
	}

	return newSeq;
}
