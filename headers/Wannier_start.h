//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <complex>
//#include <strings.h>
//#include <armadillo>
#include <fftw3.h>   //for fourier transform
//#include <iomanip>   //essential for setw
//#include <stdio.h>      /* printf, scanf, puts, NULL */
//#include <stdlib.h>     /* srand, rand */
//#include <time.h>       /* time */
//
//using namespace std;
//using namespace arma;
//
//#include <Constants.h>
//#include <BoostArrays.h>
//#include <typedef.h>
//#include <Coordinate.h>
//-------------------------Variables





//v-------------------variables
#ifndef Variables_HPP
#define Variables_HPP
int NumberOfWannierFunctions;
vector<int> Nk(3);

vec2d a;
vector<int> DegeneracyOfR;
vector<vector<double>> k;
vector<vector<int>> RCoordinates;
vector<vector<int>> RSpace_U;
vector<vector<vector<complex<double>>>> H;
vector<vector<vector<complex<double>>>> U;
vector<vector<vector<complex<double>>>> U_interpolated;
vector<vector<vector<complex<double>>>> UR;
vector<vector<vector<complex<double>>>> XiX;
vector<vector<vector<complex<double>>>> XiY;
vector<vector<vector<complex<double>>>> XiZ;
double threshold=1.e-02;
double RadiusOfSingularities = 0.1;
vector<int> Rindex_H;
vector<int> Rindex_XiX;
vector<int> Rindex_XiY;
vector<int> Rindex_XiZ;
vector<int> Rindex_U;
vector<vector<vector<int>>> Cluster;
vector<vector<vector<complex<double>>>> Difference;


#endif





//--------------------------Functions
vector<string> ReadFile(string& filename);
vector<string> SplitStringInWords(string& String);
vector<vector<double>> SetKSpace(vector<int>& Nk);
void ReadTB(string& filename);
void HermitianDFT(vector<vector<vector<complex<double>>>>& Matrix);
void ResumeParameters();
void ResumeParametersRSpace();
vector<vector<complex<double>>> DFT(vector<vector<vector<complex<double>>>>& MatrixR, vector<double>& k, vector<int>& Rindex, vector<vector<int>>& R);
vector<vector<complex<double>>> UMatrix(vector<double>& k);
void PrintDatas(vector<int>& Rindex, vector<vector<int>>& R, vector<vector<vector<complex<double>>>>& Matrix,  const string& string);

void FixPhase(vector<vector<vector<complex<double>>>>& U);
vector<int> CleanRSpace(vector<vector<vector<complex<double>>>>& Matrix);
vector<vector<int>> SetRSpace(vector<int>& NR);
vector<vector<vector<complex<double>>>> matrixFFT(vector<vector<vector<complex<double>>>>& Uk,vector<int>& Nk);
void PrintDatask(vector<vector<double>>& k, vector<vector<vector<complex<double>>>>& Matrix,  const string& string);
vector<vector<vector<complex<double>>>>	CalculateDifference(vector<vector<vector<complex<double>>>>& U, vector<vector<vector<complex<double>>>>& U_interpolated);

vector<vector<bool>> IsItHigh(vector<vector<vector<complex<double>>>>& Difference);
vector<int> NumberOfTrues(vector<vector<bool>>& Singularity);
bool IsInside(int ikneighbor1D, vector<int>& iksToInspect);
void Clusterize(vector<vector<bool>>& Singularity);
vector<int> FindMaximumChord(vector<int>& Cluster);			


void ConvertRToStaticArray(vector<int>& Rindex_H, vector<vector<int>>& RCoordinates, vec2i& R);
void ConvertToStaticArray(vector<int>&  Rindex_H, vector<vector<vector<complex<double>>>>& H, vec3x& MatrixR);

//main function
void Wannier::Wannier_start(string& name)
{	
	system("mkdir WannierInformations");
	Nk[0] = 200;      Nk[1] = 200;      Nk[2] = 1;
	//reading file, making sure all the observables are hermitian
	ReadTB(name);
	Coordinate::set_crys_to_cart(a);
	HermitianDFT(H);
	HermitianDFT(XiX);
	HermitianDFT(XiY);
	HermitianDFT(XiZ);
	
	ResumeParameters();
	
	k = SetKSpace(Nk);

	for(vector<double> Kpoint : k)
		U.push_back(UMatrix(Kpoint));
	PrintDatask(k, U,            "WannierInformations/Uk_randomphase");

	FixPhase(U);		

	UR = matrixFFT(U,Nk);
	RSpace_U   = SetRSpace(Nk);

	//print files with acquired data 
	vector<int> Rindex_trivial;
	for(int iR=0; iR<RCoordinates.size(); iR++) Rindex_trivial.push_back(iR);
	PrintDatas(Rindex_trivial, RCoordinates, H,   "WannierInformations/H_tb.txt");
	PrintDatas(Rindex_trivial, RCoordinates, XiX, "WannierInformations/XiX_tb.txt");
	PrintDatas(Rindex_trivial, RCoordinates, XiY, "WannierInformations/XiY_tb.txt");
	PrintDatas(Rindex_trivial, RCoordinates, XiZ, "WannierInformations/XiZ_tb.txt");
	PrintDatas(Rindex_trivial, RCoordinates, XiZ, "WannierInformations/XiZ_tb.txt");
	Rindex_trivial.erase(Rindex_trivial.begin(), Rindex_trivial.end());
	
	for(int iR=0; iR<RSpace_U.size(); iR++) Rindex_trivial.push_back(iR);
	PrintDatas(Rindex_trivial, RSpace_U, UR, "WannierInformations/U_tb.txt");	
	Rindex_trivial.erase(Rindex_trivial.begin(), Rindex_trivial.end());

	//Clean R space
	Rindex_H   = CleanRSpace(H);
	Rindex_XiX = CleanRSpace(XiX);
	Rindex_XiY = CleanRSpace(XiY);
	Rindex_XiZ = CleanRSpace(XiZ);
	Rindex_U   = CleanRSpace(UR);
	//print files with cleaned datas
	PrintDatas(Rindex_H,   RCoordinates, H,   "WannierInformations/H_cleaned.txt");
	PrintDatas(Rindex_XiX, RCoordinates, XiX, "WannierInformations/XiX_cleaned.txt");
	PrintDatas(Rindex_XiY, RCoordinates, XiY, "WannierInformations/XiY_cleaned.txt");
	PrintDatas(Rindex_XiZ, RCoordinates, XiZ, "WannierInformations/XiZ_cleaned.txt");
	PrintDatas(Rindex_U,   RSpace_U,     UR,  "WannierInformations/U_cleaned.txt");
	
	ResumeParametersRSpace();

	//Find singularities
	vector<vector<vector<complex<double>>>> Hk;
	vector<vector<vector<complex<double>>>> Hdiag;
	for(vector<double> kpt : k)
		U_interpolated.push_back(DFT(UR, kpt, Rindex_U, RSpace_U)); //Rindex_U, RSpace_U));

	Hdiag.resize(U.size());
	for(int ik=0; ik<k.size(); ik++)
	{
		Hk.push_back(DFT(H,k[ik],Rindex_H,RCoordinates));
		Hdiag[ik].resize(U[0].size());
		for(int row=0; row<U[0].size(); row++)
		{
			Hdiag[ik][row].assign(U[0].size(), complex<double>(0,0));
			for(int col=0; col<U[0].size(); col++)
				for(int col1=0; col1<U[0].size(); col1++)
					for(int col2=0; col2<U[0].size(); col2++)
						Hdiag[ik][row][col] += U[ik][row][col1]*Hk[ik][col1][col2]*conj(U[ik][col][col2]);
		}
	}

	PrintDatask(k, Hk,            "WannierInformations/Hk");
	PrintDatask(k, Hdiag,         "WannierInformations/Energy");
	PrintDatask(k, U,             "WannierInformations/Uk");
	PrintDatask(k, U_interpolated,"WannierInformations/UkCleaned");
       

	Difference = CalculateDifference(U, U_interpolated);
	PrintDatask(k, Difference, "WannierInformations/Difference");


	vector<vector<bool>> SingularityRow;
	SingularityRow = IsItHigh(Difference);
	Singularity.resize(SingularityRow.size());
	Singularity.fill(false);

	for(int irow=0; irow<SingularityRow[0].size(); irow++)
		for(int ik=0; ik<SingularityRow.size(); ik++)
			Singularity[ik] = Singularity[ik] || SingularityRow[ik][irow];

	ConvertRToStaticArray(Rindex_H,   RCoordinates, RvvH);
	ConvertRToStaticArray(Rindex_XiX, RCoordinates, Rvvx);
	ConvertRToStaticArray(Rindex_XiY, RCoordinates, Rvvy);
	ConvertRToStaticArray(Rindex_XiZ, RCoordinates, Rvvz);
	ConvertRToStaticArray(Rindex_U,   RSpace_U, RvvU);

	ConvertToStaticArray(Rindex_H, H, Hvv);
	ConvertToStaticArray(Rindex_XiX, XiX, xvv);
	ConvertToStaticArray(Rindex_XiY, XiY, yvv);
	ConvertToStaticArray(Rindex_XiZ, XiZ, zvv);
	ConvertToStaticArray(Rindex_U, UR,Uvv);
}





vector<string> ReadFile(string& filename)
{
	/*********************************************************************************
	*      INPUT  -> name of the file to read                                        *
	*      WannierInformations -> A vector of string, each element represent a line of the file   *
	*                                                                                *
	*      This function creates a vector of strings, each component is a line       *
	*      inside the file.                                                          *
	**********************************************************************************/
	vector<string> lines;
	ifstream fp_input;
    fp_input.open(filename.c_str());

    if(!fp_input) exit(1);    
    int i=0;
    lines.push_back("");
    while(getline(fp_input,lines[i]))   {lines.push_back(""); i++;}  
    fp_input.close();
    return lines;
}


vector<string> SplitStringInWords(string& String)
{
	/*********************************************************************************
	*      INPUT  -> A string we want to separate in words                           *
	*      WannierInformations -> A vector of string, each element is a different word of the inpu*
	*                                                                                *
	**********************************************************************************/

	vector<string> Words;
	Words.push_back("");

	for (char const &c : String) 
	{
		if(c != ' ')
			Words[Words.size()-1] += c;
		else if(Words.back() != "")
			Words.push_back("");
	}
	
	if(Words.back()=="")
		Words.pop_back();
	return Words;
}




void ReadTB(string& filename)
{
	/*********************************************************************************
	*      INPUT  -> The seedname of the _tb.dat file                                *
    *                                                                                *
	*      This function stores in global variables the values stored in TB file     *                                                                                *
	**********************************************************************************/
	filename += "_tb.dat";
	vector<string> LinesInFile;
	vector<vector<string>> SplittedLines;

	LinesInFile = ReadFile(filename);
	SplittedLines.resize(LinesInFile.size());

	for(size_t index=0; index < SplittedLines.size(); index++)
		SplittedLines[index] = SplitStringInWords(LinesInFile[index]);

	//From splitted lines to local variables
	a.resize(3,3);
	for(size_t aIndex=0; aIndex<3; aIndex++)
		for(size_t coor=0; coor<3; coor++)
			a[aIndex][coor] = atof(SplittedLines[aIndex+1][coor].c_str())*space_A_au;

	NumberOfWannierFunctions = atoi(SplittedLines[4][0].c_str());
	RCoordinates.resize(atoi(SplittedLines[5][0].c_str()));
	for(size_t i=0; i<RCoordinates.size(); i++)
		RCoordinates[i].resize(3);

	int iEndOfDegeneracySection = 6;
	do
		iEndOfDegeneracySection++;
	while(SplittedLines[iEndOfDegeneracySection].size() != 0);

	for(size_t line = 6; line <= iEndOfDegeneracySection; line++)
		for(string deg : SplittedLines[line])
			DegeneracyOfR.push_back(atoi(deg.c_str()));

	H.resize(RCoordinates.size());
	XiX.resize(RCoordinates.size());
	XiY.resize(RCoordinates.size());
	XiZ.resize(RCoordinates.size());
	for(int iR=0; iR<H.size(); iR++)
	{
		H[iR].resize(NumberOfWannierFunctions);
		XiX[iR].resize(NumberOfWannierFunctions);
		XiY[iR].resize(NumberOfWannierFunctions);
		XiZ[iR].resize(NumberOfWannierFunctions);
		for(int iWannier=0; iWannier<NumberOfWannierFunctions; iWannier++)
		{
			H[iR][iWannier].resize(NumberOfWannierFunctions);
			XiX[iR][iWannier].resize(NumberOfWannierFunctions);
			XiY[iR][iWannier].resize(NumberOfWannierFunctions);
			XiZ[iR][iWannier].resize(NumberOfWannierFunctions);
		}
	}

	size_t iLine      = iEndOfDegeneracySection+1;
	size_t iLineEnd   = iLine + RCoordinates.size()*(2 + NumberOfWannierFunctions*NumberOfWannierFunctions);
	size_t iR = 0;
	while(iLine < iLineEnd)
	{
		for(int i=0; i<3; i++)
			RCoordinates[iR][i] = atoi(SplittedLines[iLine][i].c_str());

		iLine++;
		while(SplittedLines[iLine].size() != 0)
		{
			size_t row = atoi(SplittedLines[iLine][0].c_str())-1;
			size_t column = atoi(SplittedLines[iLine][1].c_str())-1;
			double RealPart = atof(SplittedLines[iLine][2].c_str())*energy_eV_au;
			double ImagPart = atof(SplittedLines[iLine][3].c_str())*energy_eV_au;			
			H[iR][row][column] = complex<double>(RealPart, ImagPart);
			iLine++;
		}
		iLine++;
		iR++;
	}

	iR = 0;
	while(iLine < SplittedLines.size()-1)
	{
		iLine++;
		while(SplittedLines[iLine].size() != 0)
		{
			size_t row = atoi(SplittedLines[iLine][0].c_str())-1;
			size_t column = atoi(SplittedLines[iLine][1].c_str())-1;
			double XRealPart = atof(SplittedLines[iLine][2].c_str())*space_A_au;
			double XImagPart = atof(SplittedLines[iLine][3].c_str())*space_A_au;			
			double YRealPart = atof(SplittedLines[iLine][4].c_str())*space_A_au;
			double YImagPart = atof(SplittedLines[iLine][5].c_str())*space_A_au;			
			double ZRealPart = atof(SplittedLines[iLine][6].c_str())*space_A_au;
			double ZImagPart = atof(SplittedLines[iLine][7].c_str())*space_A_au;			
			XiX[iR][row][column] = complex<double>(XRealPart, XImagPart);
			XiY[iR][row][column] = complex<double>(YRealPart, YImagPart);
			XiZ[iR][row][column] = complex<double>(ZRealPart, ZImagPart);
			iLine++;
		}
		iLine++;
		iR++;
	}


}




void ResumeParameters()
{
	printf("\n*********************************************************************************************\n");
	printf(  "*                              Summary of parameters                                        *\n");
	printf(  "*********************************************************************************************\n");  
	printf(  "*  a1:%30s(%8.4f,%8.4f,%8.4f) angstrom%19s*\n"," ",a[0][0],a[0][1],a[0][2]," ");
	printf(  "*  a2:%30s(%8.4f,%8.4f,%8.4f) angstrom%19s*\n"," ",a[1][0],a[1][1],a[1][2]," ");
	printf(  "*  a3:%30s(%8.4f,%8.4f,%8.4f) angstrom%19s*\n"," ",a[2][0],a[2][1],a[2][2]," ");
	printf(  "*  Number of R points:%25s%4i%41s*\n"," ",int(RCoordinates.size())," ");
	printf(  "*  First and last R:  %15s(%2i,%2i,%2i)%4s(%2i,%2i,%2i)%31s*\n"," ",RCoordinates[0][0], RCoordinates[0][1],RCoordinates[0][2]," ",RCoordinates[RCoordinates.size()-1][0],RCoordinates[RCoordinates.size()-1][1],RCoordinates[RCoordinates.size()-1][2]," ");
	printf(  "*  First and last deg:%24s%2i %2i%41s*\n"," ",DegeneracyOfR[0],DegeneracyOfR[DegeneracyOfR.size()-1]," ");
	printf(  "*********************************************************************************************\n");  
}	

vector<vector<double>> SetKSpace(vector<int>& Nk)
{
	/*********************************************************************************
	*      INPUT  -> The number of k points in each direction                        *
    *      WannierInformations -> The k space in fractional coordinates                           *                                                   
    *                                                                                *
	*      This function stores in a 2d array the k space                            *                                                    
	*      This function stores in a 2d array the k space                            *    
	*      WARNING: in order for the program to work the best is to put the k space  *    
	*      in the region [0,1] instead of [-.5,.5]                                   *    
	**********************************************************************************/

	vector<vector<double>> k;
	k.resize(Nk[0]*Nk[1]*Nk[2]);
	vector<int> ik3D(3);
	for(size_t ik=0; ik<k.size(); ik++)
	{
		ik3D[0] = (ik/(Nk[2]*Nk[1]))%Nk[0];
	    ik3D[1] = (ik/Nk[2])%Nk[1];
    	ik3D[2] = ik%Nk[2];

		k[ik].resize(3);

		for(int icoor=0; icoor<3; icoor++)
			//k[ik][icoor] = -.5 + (ik3D[icoor]+.5)/double(Nk[icoor]);
			k[ik][icoor] = ik3D[icoor]/double(Nk[icoor]);
	}
	return k;
}


vector<vector<complex<double>>> DFT(vector<vector<vector<complex<double>>>>& MatrixR, vector<double>& k, vector<int>& Rindex, vector<vector<int>>& R)
{
	/*********************************************************************************
	*      INPUT  -> The matrix function in space R, the point where we calculate    *
    *      WannierInformations -> The matrix function at k                                        *                                                    
    *                                                                                *
	*      This function stores in a 3d array the Matrix function at k               *                                                                                *
	**********************************************************************************/
	vector<vector<complex<double>>> MatrixK;
	MatrixK.assign(MatrixR[0].size(), vector<complex<double>>(MatrixR[0].size(), 0.));


	for(int iR=0; iR<Rindex.size(); iR++)
	{	
		double KdotR = 0.;
		for(int coor=0; coor<3; coor++)
			KdotR += R[Rindex[iR]][coor]*k[coor];
		for(size_t row=0; row<MatrixR[Rindex[iR]].size(); row++)
		{
			for(size_t column=0; column<MatrixR[Rindex[iR]][row].size(); column++)
				MatrixK[row][column] += exp(2.*pi*c1*KdotR)*MatrixR[Rindex[iR]][row][column];
		}
	}
	return MatrixK;	
}





vector<vector<complex<double>>> UMatrix(vector<double>& k)
{
	/*********************************************************************************
	*      This function calculates the U matrix from the H                          *                                                                                *
    *      with a diagonalization.                                                   *
	**********************************************************************************/

	vector<int> Rindex_trivial;
	for(int iR=0; iR<RCoordinates.size(); iR++) Rindex_trivial.push_back(iR);
	
	vector<vector<complex<double>>> HK;
	HK = DFT(H,k,Rindex_trivial,RCoordinates);

    arma::vec EigenValues;
    arma::cx_mat Uarma;
    arma::cx_mat HKarma; 
    EigenValues.zeros(HK.size());
    Uarma.zeros(HK.size(), HK.size());
    HKarma.zeros(HK.size(), HK.size());
    //copy the Ek matrix in an armadillo matrix
    for (size_t row=0; row<HK.size(); row++)
        for (size_t column=0; column<HK.size(); column++)
            HKarma(row, column) = HK[row][column];

    arma::eig_sym(EigenValues, Uarma, HKarma);

    vector<vector<complex<double>>> UK;
    UK.resize(HK.size());
    for (size_t row=0; row<HK.size(); row++)
    {   
    	UK[row].resize(HK.size());
        for (size_t column=0; column<HK.size(); column++)
 			UK[row][column] = conj(Uarma(column,row));
    }
    return UK;
}







void FixPhase(vector<vector<vector<complex<double>>>>& U)
{
	/*********************************************************************************
	*      INPUT  -> a 3D array                                                      *                                                                                *
    *      WannierInformations -> the same array with the phase fixed                             *
    *                                                                                *
    *      ALGORITHM:                                                                *
    *      - Calculate: - the distance between two contiguous eigenvectors           *     
    *                   - the distance between these with one changed by sign        *
    *      - Choose the sign as the one that minimizes the distance                  *
    *                                                                                *
	*      Note that the eigenvectors are the (conj) of the row of U, so we fix the  *
	*      phase of the row.                                                         *
	**********************************************************************************/

	for(size_t ik=1; ik<U.size(); ik++)
	{	
		for(size_t row=0; row<U[ik].size(); row++)
		{
			double Distance1  = 0.;
			double Distance2  = 0.; 
	
			for(size_t column=0; column<U[ik].size(); column++)
			{
					Distance1  += real((U[ik][row][column]-U[ik-1][row][column])*conj(U[ik][row][column]-U[ik-1][row][column]));
					Distance2  += real((U[ik][row][column]+U[ik-1][row][column])*conj(U[ik][row][column]+U[ik-1][row][column]));
			}
			
			if(Distance1 > Distance2)
				for(size_t column=0; column<U[ik].size(); column++)
					U[ik][row][column] *= -1.;
		}
	}
}




vector<int> CleanRSpace(vector<vector<vector<complex<double>>>>& Matrix)
{
	/*********************************************************************************
	*      INPUT  -> a 3D array with noise                                           *                                                                                *
    *      WannierInformations -> the same array with less noise                                  *
    *                                                                                *
    *      ALGORITHM:                                                                *
    *      - Calculate the maximum of each matrix element                            *     
    *      - We keep only the values that are high enough (>threshold)               *             
	**********************************************************************************/

	vector<int> Rindex;
//	vector<vector<double>> MaximumMatrix;
//	MaximumMatrix.assign(Matrix[0].size(), vector<double>(Matrix[0].size(), 0.));
	vector<double> MaximumMatrix;
	MaximumMatrix.assign(Matrix[0].size(), 0.);

	//Find maximum of each matrix element
	for(size_t iR=0; iR<Matrix.size(); iR++)
		for(size_t row=0; row<Matrix[0].size(); row++)
			for(size_t column=0; column<Matrix[0].size(); column++)
//				if(abs(Matrix[iR][row][column]) > MaximumMatrix[row][column])
//					MaximumMatrix[row][column] = abs(Matrix[iR][row][column]);
				if(abs(Matrix[iR][row][column]) > MaximumMatrix[row])
					MaximumMatrix[row] = abs(Matrix[iR][row][column]);

	//Fill Rindex vector
	for(size_t iR=0; iR<Matrix.size(); iR++)
	{
		bool Filled=false;
		for(size_t row=0; row<Matrix[0].size(); row++)
			for(size_t column=0; column<Matrix[0].size(); column++)
//				if(abs(abs(Matrix[iR][row][column])-MaximumMatrix[row][column])/MaximumMatrix[row][column] > threshold && !Filled)
//				{
//					Rindex.push_back(int(iR));
//					Filled = true;
//				}
				if(abs(Matrix[iR][row][column])/MaximumMatrix[row] > threshold && !Filled)
				{
					Rindex.push_back(int(iR));
					Filled = true;
				}
	}
	return Rindex;
}


void ResumeParametersRSpace()
{
	printf(  "*********************************************************************************************\n");  
	printf(  "*  Number of R points H:                       %3i                                          *\n", int(Rindex_H.size()));
	printf(  "*  Number of R points XiX:                     %3i                                          *\n", int(Rindex_XiX.size()));
	printf(  "*  Number of R points XiY:                     %3i                                          *\n", int(Rindex_XiY.size()));
	printf(  "*  Number of R points XiZ:                     %3i                                          *\n", int(Rindex_XiZ.size()));
	printf(  "*  Number of R points U:                       %3i                                          *\n", int(Rindex_U.size()));
	printf(  "*********************************************************************************************\n");  
}





vector<vector<int>> SetRSpace(vector<int>& NR)
{
	/*********************************************************************************
	*      INPUT  -> The number of points in each direction                          *                                                                                *
    *      WannierInformations -> The R space                                                     *
    *                                                                                *
    *      This function can calculate the R space with the indeces required by      *
    *      the fft routine.                                                          *      
    **********************************************************************************/

	vector<vector<int>> R;
	R.assign(NR[0]*NR[1]*NR[2], vector<int>(3,0));
	vector<int> maxR(3);
	for(size_t i=0; i<3; i++)
		maxR[i] = NR[i]/2;

	//defining vectors in R space: 
    for(size_t iR1=0; iR1<NR[0];iR1++)
    	for(size_t iR2=0; iR2<NR[1]; iR2++)
    		for(size_t iR3=0; iR3<NR[2]; iR3++)
    		{
    			int index1d = iR3+NR[2]*iR2+NR[1]*NR[2]*iR1;
    			R[index1d][0] = -maxR[0]+iR1;
        		R[index1d][1] = -maxR[1]+iR2;
        		R[index1d][2] = -maxR[2]+iR3;
    		}
    return R;
}





vector<vector<vector<complex<double>>>> matrixFFT(vector<vector<vector<complex<double>>>>& Uk,vector<int>& Nk)
{
	/*********************************************************************************
	*      INPUT  -> Matrix in k space                                               *                                                                                *
    *      WannierInformations -> Matrix in R space                                               *
    *                                                                                *
    *      This function calculate the FFT of a matrix in k space. Then it saves     *
    *      it with the right indeces.                                                *      
	*      cfr. index1d_t that we use here and index1d of the Rspace          .      *
	*      Remember that fftw put in out vector first the positive R values and then *
	*      the negative R value                                                      *
	*      Note that we need several loops:                                          *
	*      -We go through the matrix and we fourier transform each element           *
	*      -We go through R space and we save in the matrix UR the values we found   *
    **********************************************************************************/
	fftw_complex* in  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Uk.size());
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Uk.size());

	vector<vector<vector<complex<double>>>> UR;
	UR.assign(Uk.size(), vector<vector<complex<double>>>(Uk[0].size(), vector<complex<double>>(Uk[0].size())));

	for(size_t row=0; row<Uk[0].size(); row++)
	{
		for(size_t column=0; column<Uk[0].size(); column++)
		{
		    fftw_plan p;
			for(int ik=0; ik<U.size(); ik++)
			{
				in[ik][0]=Uk[ik][row][column].real();
				in[ik][1]=Uk[ik][row][column].imag();
			}    	
	
    		//executing 3d fft
    		p = fftw_plan_dft_3d(Nk[0], Nk[1], Nk[2], in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    		fftw_execute(p);
    		fftw_destroy_plan(p);

    		//save the element of the matrix
			vector<int> maxR(3);
			for(size_t i=0; i<3; i++)
				maxR[i] = Nk[i]/2;

		    for(int iR1=0; iR1<Nk[0];iR1++)
    		{
    			int sign1;
    			if(iR1<maxR[0]) sign1=+1;
    			else sign1 = -1;
		
    			for(int iR2=0; iR2<Nk[1]; iR2++)
    			{
    				int sign2;
    				if(iR2<maxR[1]) sign2=+1;
    				else sign2 = -1;
		
	    				for(int iR3=0; iR3<Nk[2]; iR3++)
    				{
		    			int sign3;
    					if(iR3<maxR[2]) sign3=+1;
    					else sign3 = -1;
		    			int index1d = iR3+Nk[2]*iR2+Nk[1]*Nk[2]*iR1;
	   					int index1d_t = Nk[1]*Nk[2]*(iR1+sign1*maxR[0]) + Nk[2]*(iR2+sign2*maxR[1]) + iR3+sign3*maxR[2];
						UR[index1d][row][column] = complex<double> (out[index1d_t][0]/Uk.size(), out[index1d_t][1]/Uk.size());
					}//end iR3
				}//end iR2
			}//end iR1
		}//end column
	}//end row

	fftw_free(in);
	fftw_free(out);
	return UR;
}



void HermitianDFT(vector<vector<vector<complex<double>>>>& Matrix)
{
	/*********************************************************************************
	*      INPUT  -> Matrix in R space                                               *                                                                                *
    *      WannierInformations -> The same hermitian matrix                                       *
    *                                                                                *
    *      This function makes sure that the matrices in WannierInformations of wannier90 are     *
    *      Hermitian. If not then:                                                   *      
    *      Matrix_nm(R) = .5(Matrix_nm(R)+Matrix_mn(-R)*)                            *      
    *      We dont write .5 because we eject all the R points that are above the middle      
    *      i(-R) = #R-1-iR                                                           *      
    **********************************************************************************/
	for(size_t iR=0; iR<Matrix.size(); iR++)
	{	
		int iR2=RCoordinates.size()-1-iR;
		for(size_t row=0; row<Matrix[0].size(); row++)
			for(size_t column=0; column<Matrix[0].size(); column++)
				Matrix[iR][row][column] = .5*(Matrix[iR][row][column] + conj(Matrix[iR2][column][row]));
	}



}



void PrintDatas(vector<int>& Rindex, vector<vector<int>>& R, vector<vector<vector<complex<double>>>>& Matrix,  const string& string)
{
	ofstream WannierInformations;
	WannierInformations.open(string.c_str());
	for(int iR=0; iR<Rindex.size(); iR++)
	{
		Coordinate Rpoint;
		Rpoint.setcrys(R[Rindex[iR]][0], R[Rindex[iR]][1], R[Rindex[iR]][2]);
		for(int row =0; row<Matrix[0].size(); row++)
			for(int column =0; column<Matrix[0].size(); column++)
				WannierInformations << setprecision(6) << Rpoint.cart[0] << " " << Rpoint.cart[1] << " " << Rpoint.cart[2] << " " << sqrt(Rpoint.dot(Rpoint)) << " " << abs(Matrix[Rindex[iR]][row][column]) << endl;
	}
	WannierInformations.close();

}


void PrintDatask(vector<vector<double>>& k, vector<vector<vector<complex<double>>>>& Matrix,  const string& seedname)
{
    multivec2D<ofstream> WannierInformations; 
    WannierInformations.resize(U[0].size(),U[0].size());
    
    for(int row=0; row<U[0].size(); row++)	
 	{	
 		for(int column=0; column<U[0].size(); column++)
		{
			string NameWannierInformationsFile = seedname + to_string(row) + to_string(column) + ".txt";
			WannierInformations[row][column].open(NameWannierInformationsFile.c_str());
		}
	}

	for(int ik=0; ik<k.size(); ik++)
		for(int row =0; row<Matrix[0].size(); row++)
			for(int column =0; column<Matrix[0].size(); column++)
				WannierInformations[row][column] << setprecision(6) << k[ik][0] << " " << k[ik][1] << " " << k[ik][2] << " " <<  real(Matrix[ik][row][column]) << " " << imag(Matrix[ik][row][column]) << endl;
	
    for(int row=0; row<U[0].size(); row++)	
 		for(int column=0; column<U[0].size(); column++)
			WannierInformations[row][column].close();

}




vector<vector<vector<complex<double>>>>	CalculateDifference(vector<vector<vector<complex<double>>>>& U, vector<vector<vector<complex<double>>>>& U_interpolated)
{
	vector<vector<vector<complex<double>>>> Difference;
	Difference.assign(U.size(), vector<vector<complex<double>>>(U[0].size(), vector<complex<double>>(U[0][0].size())));

	for(size_t ik=0; ik<U.size(); ik++)
		for(size_t row =0; row<U[0].size(); row++)
			for(size_t column=0; column<U[0][0].size(); column++)
				Difference[ik][row][column] = (U[ik][row][column] - U_interpolated[ik][row][column])/abs(U[ik][row][column]);
	return Difference;
}



vector<vector<bool>> IsItHigh(vector<vector<vector<complex<double>>>>& Difference)
{
	vector<vector<bool>> Singularity;
	Singularity.assign(Difference.size(), vector<bool>(Difference[0].size(), false));
	for(size_t ik=0; ik<Singularity.size(); ik++)
		for(size_t row =0; row<Difference[0].size(); row++)
			for(size_t column=0; column<Difference[0][0].size(); column++)
				Singularity[ik][row] = Singularity[ik][row] || ( abs(Difference[ik][row][column]) > 10*threshold ); 	
	return Singularity;			
}





vector<int> NumberOfTrues(vector<vector<bool>>& Singularity)
{
	vector<int> Counter(Singularity[0].size(),0);
	for(size_t ik=0; ik<Singularity.size(); ik++)
		for(size_t row=0; row<Singularity[0].size(); row++)
			Counter[row] += Singularity[ik][row];
	//printing step
	printf("\n*********************************************************************************************\n");
	for(size_t row=0; row<Singularity[0].size(); row++)
		printf(  "* Eigenvector %3i     ---->           %10i singularities                              *\n",int(row),Counter[row]);
	printf("\n*********************************************************************************************\n");

	return Counter;
}



bool IsInside(int ikneighbor1D, vector<int>& iksToInspect)
{
	bool IsInside = false;
	for(int ik=0; ik<iksToInspect.size(); ik++)
		if(iksToInspect[ik] == ikneighbor1D) 
		{
			IsInside = true;
			break;
		}
	return IsInside;
}





void Clusterize(vector<vector<bool>>& Singularity)
{
	/*********************************************************************************
	*      ALGORITHM:                                                                *      
    *      for each eigenvector (row) we count the number of point in k space in whic*
    *      we expect the singularity.                                                *
    *      In Singularity we have the k space with the singularities indicated with true
    *      state.                                                                    *
    *      Therefore for each row we make a search randomly in k space. when we find *
    *      the first singularity we inspect all the neighbors points and the neighbors
    *      of the neighbors, until we don't close the loop                           *
    *      We stop when we expected all the points with singularities.               * 
    **********************************************************************************/

	vector<int> CounterSingularityPoints;
	CounterSingularityPoints = NumberOfTrues(Singularity);
	Cluster.resize(CounterSingularityPoints.size());
  	srand (time(NULL));

  	/* generate secret number between 1 and 10: */

	for(int row=0; row < CounterSingularityPoints.size(); row++)
	{
		while(CounterSingularityPoints[row] > 0)
		{

			Cluster[row].resize(Cluster[row].size()+1);
			//initialize k point
			int ik;
			do	ik = rand() % (Nk[0]*Nk[1]*Nk[2]);
			while (!Singularity[ik][row]);
			vector<int> iksToInspect;
			vector<int> ik3D(3);
			//grow cluster
			iksToInspect.push_back(ik);

			while(iksToInspect.size() > 0)
			{
				ik = iksToInspect[iksToInspect.size()-1];

				iksToInspect.pop_back();
				if(!Singularity[ik][row]) continue;

				Cluster[row][Cluster[row].size()-1].push_back(ik);
				Singularity[ik][row] = false;

			    ik3D[0] = (ik/(Nk[2]*Nk[1]))%Nk[0];
	    		ik3D[1] = (ik/Nk[2])%Nk[1];
    			ik3D[2] = ik%Nk[2];

				vector<int> ikneighbor(3);

				for(int sign1=-1; sign1<=1; sign1++)
				{
					ikneighbor[0] = (ik3D[0] + sign1)%Nk[0];

					if (ikneighbor[0] < 0) ikneighbor[0] += Nk[0];
					for(int sign2=-1; sign2<=1; sign2++)
					{
						ikneighbor[1] = (ik3D[1] + sign2)%Nk[1];
						if (ikneighbor[1] < 0) ikneighbor[1] += Nk[1];
						for(int sign3=-1; sign3<=1; sign3++)
						{
							ikneighbor[2] = (ik3D[2] + sign3)%Nk[2];
							if (ikneighbor[2] < 0) ikneighbor[2] += Nk[2];

							int ikneighbor1D = ikneighbor[2] + Nk[2]*ikneighbor[1] + Nk[2]*Nk[1]*ikneighbor[0];
							if(!IsInside(ikneighbor1D, iksToInspect) && ikneighbor1D != ik) iksToInspect.push_back(ikneighbor1D);
						}
					}
				}
			}
			CounterSingularityPoints[row] -= Cluster[row][Cluster[row].size()-1].size();
		}
	}
}






//void PrintClusters()
//{
//	ofstream Singul;
//	Singul.open("Singularity.txt");
//
//	for(size_t ik=0; ik<Singularity.size(); ik++)
//		for(size_t row =0; row<Difference[0].size(); row++)
//			Singul << k[ik][0] << " " << k[ik][1] << " " << k[ik][2] << " " << Singularity[ik][row] << endl;
//	Singul.close();
//	vector<vector<ofstream>> clusterWannierInformations(Cluster.size());
//	
//	for(int i=0; i<Cluster.size(); i++)
//	{	
//		clusterWannierInformations[i].resize(Cluster[i].size());
//		for(int j=0; j<Cluster[i].size(); j++)
//		{
//			string NameWannierInformationsFile = "Cluster"+ to_string(i) + to_string(j) + ".txt";
//			clusterWannierInformations[i][j].open(NameWannierInformationsFile.c_str());
//		}
//	}
//
//	for(int i=0; i<Cluster.size(); i++)
//		for(int j=0; j<Cluster[i].size(); j++)				
//			for(int ik=0; ik<Cluster[i][j].size(); ik++)
//				clusterWannierInformations[i][j] << k[Cluster[i][j][ik]][0] << " " << k[Cluster[i][j][ik]][1] << " " << k[Cluster[i][j][ik]][2] << endl; 
//}






/*PRINT U in cartesian coordinates
	ofstream U_cartesian_K;
	U_cartesian_K.open("U_fast.txt");
	for(int ikx=-1; ikx<2; ikx++)
	{
		for(int iky=-1; iky<2; iky++)
		{
			double aAngstrom=2.46;
			vector<double> Kpoint(3);
			double kx = 0. + 0.2/0.5*double(ikx)/3.;
			double ky = 4.*pi/(3.*aAngstrom) + 0.2/0.5*double(iky)/3.;
			Kpoint[0]=0;
			Kpoint[1] = kx*aAngstrom/(2.*pi)*cos(pi/6.) + ky*aAngstrom/(2.*pi)*sin(pi/6.);
			Kpoint[2] = kx*aAngstrom/(2.*pi)*cos(-pi/6.) + ky*aAngstrom/(2.*pi)*sin(-pi/6.);
			vector<vector<complex<double>>> UK;
			UK = UMatrix(Kpoint);
			U_cartesian_K << kx << " " << ky << " " << UK[0][0].real() << " " << UK[0][0].imag() << " " <<UK[0][1].real() << " " << UK[0][1].imag()<<endl;
		}
	}
	U_cartesian_K.close();
*/



void ConvertRToStaticArray(vector<int>& Rindex, vector<vector<int>>& RCoordinates, vec2i& R)
{
	R.resize(Rindex.size(), 3);
	for(size_t iR=0; iR<Rindex.size(); iR++)
		for(size_t icoor=0; icoor<3; icoor++)
			R[iR][icoor] = RCoordinates[Rindex[iR]][icoor];
}



void ConvertToStaticArray(vector<int>&  Rindex, vector<vector<vector<complex<double>>>>& MatrixVector, vec3x& MatrixR)
{
	MatrixR.resize(MatrixVector[0].size(), MatrixVector[0].size(), Rindex.size());

	for(int irow=0; irow<MatrixVector[0].size(); irow++)
		for(int icolumn=0; icolumn<MatrixVector[0].size(); icolumn++)
			for(int iR=0; iR<Rindex.size(); iR++)
				MatrixR[irow][icolumn][iR] = MatrixVector[Rindex[iR]][irow][icolumn];
}
