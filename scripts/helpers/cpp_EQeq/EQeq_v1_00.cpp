// The EQeq Method v. 1.0 //////////////////////////////////////////////////////////////////////////////////////////////
// Author: Christopher E. Wilmer & Randall Q. Snurr (advisor)                             //////////////////////////////
// Date: Mar. 19, 2012                                                                    //////////////////////////////
// First published in paper "An Extended Charge Equilibration Method"                     //////////////////////////////
// If you have questions, please e-mail c.wilmer@gmail.com OR snurr@northwestern.edu      //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// This source code is released under the FreeBSD license. 
// What does that mean? 
// Short answer: THIS SOFTWARE IS FREE TO USE AS YOU LIKE (commercial, academic, or fun)
// Slightly longer answer: See below.
//
// Copyright (c) <2012>, < Christopher E. Wilmer & Randall Q. Snurr >
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met: 
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer. 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution. 
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies, 
// either expressed or implied, of the FreeBSD Project.

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DISCLAIMER:                                                                            //////////////////////////////
// This is a development snapshot of the source code that was used in conjuction with     //////////////////////////////
// the published paper. It is not guaranteed to be free of errors nor is it guaranteed    //////////////////////////////
// to run on any operating system or hardware. For updated source code, with new features //////////////////////////////
// and bug fixes (as they are found) please contact the Snurr group or Chris Wilmer.      //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Running the program:                                                                   //////////////////////////////
// Program expects two input files "ionization.dat" and "chargecenters.dat" as well as    //////////////////////////////
// a CIF file for an input structure. To run the program, pass the input CIF file path    //////////////////////////////
// as the first argument to the executable (i.e: \EQeq_v1_00.exe MyDirectory/myfile.cif ) //////////////////////////////
// Additional input parameters are optional. Please look at source code below to see      //////////////////////////////
// what the other optional inputs are for (should be mostly self-explanatory).            //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The source code in this program demonstrates the charge equilibration method described //////////////////////////////
// in the accompanying paper. The purpose of the source code provided is to be            //////////////////////////////
// minimalistic and do "just the job" described. In practice, you may wish to add various //////////////////////////////
// features to the source code to fit the particular needs of your project.               //////////////////////////////
//                                                                                        //////////////////////////////
// Major highlights of program:                                                           //////////////////////////////
//      - Obtains charges for atoms in periodic systems without iteration                 //////////////////////////////
//      - Can use non-neutral charge centers for more accurate point charges              //////////////////////////////
//      - Designed for speed (but without significant code optimizations)                 //////////////////////////////
//                                                                                        //////////////////////////////
// Features not implemented but that you may want to consider adding:                     //////////////////////////////
//      - Spherical cut-offs (for both real-space and reciprocal-space sums)              //////////////////////////////
//      - An iterative loop that guesses the appropriate charge center (so the user does  //////////////////////////////
//                                                                      not have to guess)//////////////////////////////
//      - Ewald parameter auto-optimization                                               //////////////////////////////
// 		- Various code optimizations                                                      //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>		// To read files
#include <fstream>		// To output files
#include <string>
#include <vector>
#include <map>			// For string enumeration (C++ specific)
#include <cmath>		// For basic math functions
#include <cstdlib>
using namespace std;

#define TABLE_OF_ELEMENTS_SIZE 84
#define PI 3.1415926535897932384626433832795	// 32 digits of PI

// This is a clumsy way to enable switch statements with atom labels in C++
enum StringAtomLabels {
	ev_H,
	ev_He,
	ev_Li,
	ev_Be,
	ev_B,
	ev_C,
	ev_N,
	ev_O,
	ev_F,
	ev_Ne,
	ev_Na,
	ev_Mg,
	ev_Al,
	ev_Si,
	ev_P,
	ev_S,
	ev_Cl,
	ev_Ar,
	ev_K ,
	ev_Ca,
	ev_Sc,
	ev_Ti,
	ev_V ,
	ev_Cr,
	ev_Mn,
	ev_Fe,
	ev_Co,
	ev_Ni,
	ev_Cu,
	ev_Zn,
	ev_Ga,
	ev_Ge,
	ev_As,
	ev_Se,
	ev_Br,
	ev_Kr,
	ev_Rb,
	ev_Sr,
	ev_Y ,
	ev_Zr,
	ev_Nb,
	ev_Mo,
	ev_Tc,
	ev_Ru,
	ev_Rh,
	ev_Pd,
	ev_Ag,
	ev_Cd,
	ev_In,
	ev_Sn,
	ev_Sb,
	ev_Te,
	ev_I ,
	ev_Xe,
	ev_Cs,
	ev_Ba,
	ev_La,
	ev_Ce,
	ev_Pr,
	ev_Nd,
	ev_Pm,
	ev_Sm,
	ev_Eu,
	ev_Gd,
	ev_Tb,
	ev_Dy,
	ev_Ho,
	ev_Er,
	ev_Tm,
	ev_Yb,
	ev_Lu,
	ev_Hf,
	ev_Ta,
	ev_W ,
	ev_Re,
	ev_Os,
	ev_Ir,
	ev_Pt,
	ev_Au,
	ev_Hg,
	ev_Tl,
	ev_Pb,
	ev_Bi,
	ev_Po
};

// Map to associate the strings with the enum values
std::map<std::string, StringAtomLabels> s_mapStringAtomLabels;

class Coordinates {
	public:
		// Constructor
		Coordinates();

		double x;
		double y;
		double z;
};

class IonizationDatum {
	public:
		IonizationDatum();

		// TODO: Mass, radii and other properties can be added here if that would help for some reason
		string Label;
		vector<bool> isDataAvailable; // True or false
		vector<double> ionizationPotential; // The first 8 ionization potentials and the electron affinity
		int chargeCenter;
};

// EQeq function headers (alphabetical order)
void DetermineReciprocalLatticeVectors();
double GetJ(int i, int j);
void InitializeStringAtomLabelsEnumeration();
void LoadIonizationData(string filename);
void LoadChargeCenters(string filename);
void LoadCIFFile(string filename); // Reads in CIF files, periodicity can be switched off
void OutputCIFFormatFile(string filename); 
void OutputPDBFormatFile(string filename);
void OutputMOLFormatFile(string filename); // Outputs 'RASPA' MOL file
void Qeq();
void RoundCharges(int digits); // Make *slight* adjustments to the charges for nice round numbers

// Algebra helper functions (alphaAnglebetical order)
vector<double> Cross(vector<double> a, vector<double> b);
double Dot(vector<double> a, vector<double> b);
double Mag(vector<double> a);
double Round(double num);
vector<double> Scalar(double a, vector<double> b);
vector<double> SolveMatrix(vector<vector<double> > A, vector<double> b);

// Global variables
bool isPeriodic = true;
bool useEwardSums = true; // will use direct sums if false
double aLength; double bLength; double cLength;
double alphaAngle; double betaAngle; double gammaAngle;
double unitCellVolume;
vector<double> aV(3); vector<double> bV(3); vector<double> cV(3); // Real-space vectors
vector<double> hV(3); vector<double> jV(3); vector<double> kV(3); // Reciprocal-lattice vectors
int numAtoms; // To be read from input file
double Qtot; // To be read in from file
vector<Coordinates> Pos; // Array of atom positions
vector<double> J; // Atom "hardness"
vector<double> X; // Atom electronegativity
vector<double> Q; // Partial atomic charge
vector<string> Label; // Atom labels (e.g., "C1" "C2" "ZnCation" "dummyAtom")
vector<string> Symbol; // Atom symbols (e.g., "C" "O" "Zn")
vector<IonizationDatum> IonizationData(TABLE_OF_ELEMENTS_SIZE);

// Parameters and constants
double k = 14.4; // Physical constant: the vacuum permittivity 1/(4pi*epsi) [units of Angstroms * electron volts]
double eta = 50; // Ewald splitting parameter
double lambda = 1.2; // Coulomb scaling parameter
float hI0 = -2.0; // Default value used in paper
float hI1 = 13.598; // This is the empirically mesaured 1st ionization energy of hydrogen
int chargePrecision = 3; // Number of digits to use for point charges
int mR = 2;  int mK = 2;
int aVnum = mR; int bVnum = mR; int cVnum = mR; // Number of unit cells to consider in per. calc. ("real space")
int hVnum = mK; int jVnum = mK; int kVnum = mK; // Number of unit cells to consider in per. calc. ("frequency space")
/*****************************************************************************/
/*****************************************************************************/
int main (int argc, char *argv[]) {
	string inputFilename,method;

	// The only mandatory parameter is the input file parameter
	if (argc <= 1) { cout << "Error, invalid input!" << endl; exit(1); }
	if (argc > 1) inputFilename = argv[1];
	if (argc > 2) lambda = atof(argv[2]); // The dielectric screening parameter (optional, default value above)
	if (argc > 3) hI0 = atof(argv[3]); // The electron affinity of hydrogen (optional, default value above)
	if (argc > 4) chargePrecision = atoi(argv[4]); // Num of digits to use for charges (optional, default value above)
	if (argc > 5) {
		method = argv[5];
		if ((method == "NonPeriodic")||(method == "Nonperiodic")||(method == "nonperiodic")) isPeriodic = false;
		else
		if ((method == "Ewald")||(method == "ewald")) useEwardSums = true;
		else
		useEwardSums = false; // Direct sums are used if the 5th argument is misspelled
	}
	if (argc > 6) mR = atoi(argv[6]);
	if (argc > 7) mK = atoi(argv[7]);
	if (argc > 8) eta = atof(argv[8]);
	
	InitializeStringAtomLabelsEnumeration();  // Part of the clumsy way to enable string-based switch statements
	LoadIonizationData("ionizationdata.dat");
	LoadChargeCenters("chargecenters.dat");
	LoadCIFFile(inputFilename);
	
	cout << "==================================================" << endl;
	cout << "===== Calculating charges... please wait. ========" << endl;
	Qtot = 0; // Can be non-zero for non-periodic structures
	Qeq();
	RoundCharges(chargePrecision);
	cout << "===== ... done!                           ========" << endl;
	cout << "==================================================" << endl;
	
	char buffer[50];
	if (useEwardSums) method = "Ewald"; else method = "Direct";
	if (!isPeriodic) method = "NonPeriodic";
	sprintf(buffer,"_EQeq_%s_%4.2f_%4.2f",method.c_str(),lambda,hI0);
	
	OutputCIFFormatFile(inputFilename+buffer+".cif"); 
	OutputMOLFormatFile(inputFilename+buffer+".mol"); 
	OutputPDBFormatFile(inputFilename+buffer+".pdb"); 

	return 0;
}
/*****************************************************************************/
/*****************************************************************************/
Coordinates::Coordinates() {
  x = 0; y = 0; z = 0; // default coordinates
}
/*****************************************************************************/
IonizationDatum::IonizationDatum() {
	isDataAvailable.resize(9,false);
	ionizationPotential.resize(9,0);
}
/*****************************************************************************/
void DetermineReciprocalLatticeVectors() {
	vector<double> crs;
	double pf; // pf => PreFactor
	
	crs = Cross(bV, cV);
	pf = 2*PI / Dot(aV, crs);
	hV[0] = pf * crs[0];
	hV[1] = pf * crs[1];
	hV[2] = pf * crs[2];

	crs = Cross(cV, aV);
	pf = 2*PI / Dot(bV, crs);
	jV[0] = pf * crs[0];
	jV[1] = pf * crs[1];
	jV[2] = pf * crs[2];

	crs = Cross(aV, bV);
	pf = 2*PI / Dot(cV, crs);
	kV[0] = pf * crs[0];
	kV[1] = pf * crs[1];
	kV[2] = pf * crs[2];
}
/*****************************************************************************/
void InitializeStringAtomLabelsEnumeration() {
	s_mapStringAtomLabels["H "] = ev_H;	// 1
	s_mapStringAtomLabels["He"] = ev_He;// 2
	s_mapStringAtomLabels["Li"] = ev_Li;// 3
	s_mapStringAtomLabels["Be"] = ev_Be;// 4
	s_mapStringAtomLabels["B "] = ev_B;	// 5
	s_mapStringAtomLabels["C "] = ev_C;	// 6
	s_mapStringAtomLabels["N "] = ev_N;	// 7
	s_mapStringAtomLabels["O "] = ev_O;	// 8
	s_mapStringAtomLabels["F "] = ev_F;	// 9
	s_mapStringAtomLabels["Ne"] = ev_Ne;//10
	s_mapStringAtomLabels["Na"] = ev_Na;//11
	s_mapStringAtomLabels["Mg"] = ev_Mg;//12
	s_mapStringAtomLabels["Al"] = ev_Al;//13
	s_mapStringAtomLabels["Si"] = ev_Si;//14
	s_mapStringAtomLabels["P "] = ev_P;	//15
	s_mapStringAtomLabels["S "] = ev_S;	//16
	s_mapStringAtomLabels["Cl"] = ev_Cl;//17
	s_mapStringAtomLabels["Ar"] = ev_Ar;//18
	s_mapStringAtomLabels["K "] = ev_K ;//19
	s_mapStringAtomLabels["Ca"] = ev_Ca;//20
	s_mapStringAtomLabels["Sc"] = ev_Sc;//21
	s_mapStringAtomLabels["Ti"] = ev_Ti;//22
	s_mapStringAtomLabels["V "] = ev_V ;//23
	s_mapStringAtomLabels["Cr"] = ev_Cr;//24
	s_mapStringAtomLabels["Mn"] = ev_Mn;//25
	s_mapStringAtomLabels["Fe"] = ev_Fe;//26
	s_mapStringAtomLabels["Co"] = ev_Co;//27
	s_mapStringAtomLabels["Ni"] = ev_Ni;//28
	s_mapStringAtomLabels["Cu"] = ev_Cu;//29
	s_mapStringAtomLabels["Zn"] = ev_Zn;//30
	s_mapStringAtomLabels["Ga"] = ev_Ga;//31
	s_mapStringAtomLabels["Ge"] = ev_Ge;//32
	s_mapStringAtomLabels["As"] = ev_As;//33
	s_mapStringAtomLabels["Se"] = ev_Se;//34
	s_mapStringAtomLabels["Br"] = ev_Br;//35
	s_mapStringAtomLabels["Kr"] = ev_Kr;//36
	s_mapStringAtomLabels["Rb"] = ev_Rb;//37
	s_mapStringAtomLabels["Sr"] = ev_Sr;//38
	s_mapStringAtomLabels["Y "] = ev_Y ;//39
	s_mapStringAtomLabels["Zr"] = ev_Zr;//40
	s_mapStringAtomLabels["Nb"] = ev_Nb;//41
	s_mapStringAtomLabels["Mo"] = ev_Mo;//42
	s_mapStringAtomLabels["Tc"] = ev_Tc;//43
	s_mapStringAtomLabels["Ru"] = ev_Ru;//44
	s_mapStringAtomLabels["Rh"] = ev_Rh;//45
	s_mapStringAtomLabels["Pd"] = ev_Pd;//46
	s_mapStringAtomLabels["Ag"] = ev_Ag;//47
	s_mapStringAtomLabels["Cd"] = ev_Cd;//48
	s_mapStringAtomLabels["In"] = ev_In;//49
	s_mapStringAtomLabels["Sn"] = ev_Sn;//50
	s_mapStringAtomLabels["Sb"] = ev_Sb;//51
	s_mapStringAtomLabels["Te"] = ev_Te;//52
	s_mapStringAtomLabels["I "] = ev_I ;//53
	s_mapStringAtomLabels["Xe"] = ev_Xe;//54
	s_mapStringAtomLabels["Cs"] = ev_Cs;//55
	s_mapStringAtomLabels["Ba"] = ev_Ba;//56
	s_mapStringAtomLabels["La"] = ev_La;//57
	s_mapStringAtomLabels["Ce"] = ev_Ce;//58
	s_mapStringAtomLabels["Pr"] = ev_Pr;//59
	s_mapStringAtomLabels["Nd"] = ev_Nd;//60
	s_mapStringAtomLabels["Pm"] = ev_Pm;//61
	s_mapStringAtomLabels["Sm"] = ev_Sm;//62
	s_mapStringAtomLabels["Eu"] = ev_Eu;//63
	s_mapStringAtomLabels["Gd"] = ev_Gd;//64
	s_mapStringAtomLabels["Tb"] = ev_Tb;//65
	s_mapStringAtomLabels["Dy"] = ev_Dy;//66
	s_mapStringAtomLabels["Ho"] = ev_Ho;//67
	s_mapStringAtomLabels["Er"] = ev_Er;//68
	s_mapStringAtomLabels["Tm"] = ev_Tm;//69
	s_mapStringAtomLabels["Yb"] = ev_Yb;//70
	s_mapStringAtomLabels["Lu"] = ev_Lu;//71
	s_mapStringAtomLabels["Hf"] = ev_Hf;//72
	s_mapStringAtomLabels["Ta"] = ev_Ta;//73
	s_mapStringAtomLabels["W "] = ev_W ;//74
	s_mapStringAtomLabels["Re"] = ev_Re;//75
	s_mapStringAtomLabels["Os"] = ev_Os;//76
	s_mapStringAtomLabels["Ir"] = ev_Ir;//77
	s_mapStringAtomLabels["Pt"] = ev_Pt;//78
	s_mapStringAtomLabels["Au"] = ev_Au;//79
	s_mapStringAtomLabels["Hg"] = ev_Hg;//80
	s_mapStringAtomLabels["Tl"] = ev_Tl;//81
	s_mapStringAtomLabels["Pb"] = ev_Pb;//82
	s_mapStringAtomLabels["Bi"] = ev_Bi;//83
	s_mapStringAtomLabels["Po"] = ev_Po;//84
}
/*****************************************************************************/
double GetJ(int i, int j) {
	// Note to reader - significant consolidation of code may be possible in this function
	if (isPeriodic == false) {
		//////////////////////////////////////////////////////////////////////
		//  NonPeriodic                                                     //
		//////////////////////////////////////////////////////////////////////
		if (i == j) {
			return J[i]; // Return the hardness/idempotential
		} else {
			double dx = Pos[i].x - Pos[j].x;
			double dy = Pos[i].y - Pos[j].y;
			double dz = Pos[i].z - Pos[j].z;
			double RabSq = dx*dx + dy*dy + dz*dz;
			double Rab = sqrt(RabSq); 
			
			double Jij = sqrt(J[i] * J[j]);
			double a = Jij / k;
			double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab); // Other functional forms are OK too
			
			double Jab = lambda * (k/2) * ((1/Rab) + orbitalOverlapTerm);
			
			return Jab;
		}
		//////////////////////////////////////////////////////////////////////
		//////////////////////////////////////////////////////////////////////
	} else 
	if (isPeriodic == true) {
		aVnum = mR; bVnum = mR; cVnum = mR; // Number of unit cells to consider in per. calc. (in "real space")
		hVnum = mK; jVnum = mK; kVnum = mK; // Number of unit cells to consider in per. calc. (in "frequency space")
		if (useEwardSums == false) {
			//////////////////////////////////////////////////////////////////////
			// Direct sums                                                      //
			//////////////////////////////////////////////////////////////////////
			if (i == j) {
			
				double sigmaStar = 0;
				for (int u = -aVnum; u <= aVnum; u++) {
					for (int v = -bVnum; v <= bVnum; v++) {
						for (int w = -cVnum; w <= cVnum; w++) {
							if (!((u==0)&&(v==0)&&(w==0))) {
								double dx = u*aV[0] + v*bV[0] + w*cV[0];
								double dy = u*aV[1] + v*bV[1] + w*cV[1];
								double dz = u*aV[2] + v*bV[2] + w*cV[2];
								double Rab = sqrt(dx*dx + dy*dy + dz*dz);
								double RabSq = dx*dx + dy*dy + dz*dz;
								
								double Jij = sqrt(J[i] * J[j]);
								double a = Jij / k;
								double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab); 
								// Other functional forms (for orbital overlap) are OK too
								
								sigmaStar += (1/Rab) + orbitalOverlapTerm;
							}
						}
					}
				}
				return J[i] + lambda * (k/2)*sigmaStar;
				
			} else {
			
				double sigma = 0;
				for (int u = -aVnum; u <= aVnum; u++) {
					for (int v = -bVnum; v <= bVnum; v++) {
						for (int w = -cVnum; w <= cVnum; w++) {
							double dx = Pos[i].x - Pos[j].x + u*aV[0] + v*bV[0] + w*cV[0];
							double dy = Pos[i].y - Pos[j].y + u*aV[1] + v*bV[1] + w*cV[1];
							double dz = Pos[i].z - Pos[j].z + u*aV[2] + v*bV[2] + w*cV[2];
							double Rab = sqrt(dx*dx + dy*dy + dz*dz);
							double RabSq = dx*dx + dy*dy + dz*dz;
				
							double Jij = sqrt(J[i] * J[j]);
							double a = Jij / k;
							double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab);
							// Other functional forms (for orbital overlap) are OK too							
							
							sigma += (1/Rab) + orbitalOverlapTerm;
						}
					}
				}
				
				return lambda * (k/2) * sigma;
			}
		} else {
			//////////////////////////////////////////////////////////////////////
			// Ewald sums                                                       //
			//////////////////////////////////////////////////////////////////////
			if (i == j) {
				// Orbital energy term
				double orbital = 0;
				for (int u = -aVnum; u <= aVnum; u++) {
					for (int v = -bVnum; v <= bVnum; v++) {
						for (int w = -cVnum; w <= cVnum; w++) {
							if ((u==0) && (v==0) && (w==0)) {
								// do nothing
							} else {
								double dx = u*aV[0] + v*bV[0] + w*cV[0];
								double dy = u*aV[1] + v*bV[1] + w*cV[1];
								double dz = u*aV[2] + v*bV[2] + w*cV[2];
								double Rab = sqrt(dx*dx + dy*dy + dz*dz);
								double RabSq = dx*dx + dy*dy + dz*dz;
					
								double Jij = sqrt(J[i] * J[j]);
								double a = Jij / k;
								double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab);
								// Other functional forms (for orbital overlap) are OK too							
								
								orbital += orbitalOverlapTerm;
							}
						}
					}
				}
			
				// Real-space Coulomb component
				double alphaStar = 0;
				for (int u = -aVnum; u <= aVnum; u++) {
					for (int v = -bVnum; v <= bVnum; v++) {
						for (int w = -cVnum; w <= cVnum; w++) {
							if ((u==0) && (v==0) && (w==0)) {
								// do nothing
							} else {
								double dx = u*aV[0] + v*bV[0] + w*cV[0];
								double dy = u*aV[1] + v*bV[1] + w*cV[1];
								double dz = u*aV[2] + v*bV[2] + w*cV[2];
								double Rab = sqrt(dx*dx + dy*dy + dz*dz);
					
								alphaStar += erfc( Rab / eta ) / Rab;
							}
						}
					}
				}
				
				// K-space component
				double betaStar = 0;
				double h = 0; double b = 0; 
				vector<double> RLV(3); // reciprocal lattice vector
				for (int u = -hVnum; u <= hVnum; u++) {
					for (int v = -jVnum; v <= jVnum; v++) {
						for (int w = -kVnum; w <= kVnum; w++) {
							if ((u==0) && (v==0) && (w==0)) {
								// do nothing
							} else {
								RLV[0] = u*hV[0] + v*jV[0] + w*kV[0];
								RLV[1] = u*hV[1] + v*jV[1] + w*kV[1];
								RLV[2] = u*hV[2] + v*jV[2] + w*kV[2];
								
								h = Mag(RLV);
								b = 0.5 * h * eta;
								
								//beta += cos( RLV[0]*dx + RLV[1]*dy + RLV[2]*dz ) / (h*h) * exp(-b*b);
								betaStar += 1 / (h*h) * exp(-b*b);
							}
						}
					}
				}
				betaStar *= 4*PI / unitCellVolume; 
				
				return J[i] + lambda * (k/2) * (alphaStar + betaStar + orbital - 2/(eta*sqrt(PI)));
				
			} else {
				// Orbital energy term
				double orbital = 0;
				for (int u = -aVnum; u <= aVnum; u++) {
					for (int v = -bVnum; v <= bVnum; v++) {
						for (int w = -cVnum; w <= cVnum; w++) {
							double dx = Pos[i].x - Pos[j].x + u*aV[0] + v*bV[0] + w*cV[0];
							double dy = Pos[i].y - Pos[j].y + u*aV[1] + v*bV[1] + w*cV[1];
							double dz = Pos[i].z - Pos[j].z + u*aV[2] + v*bV[2] + w*cV[2];
							double Rab = sqrt(dx*dx + dy*dy + dz*dz);
							double RabSq = dx*dx + dy*dy + dz*dz;
				
							double Jij = sqrt(J[i] * J[j]);
							double a = Jij / k;
							double orbitalOverlapTerm = exp(-(a*a*RabSq))*(2*a - a*a*Rab - 1/Rab); 
							// Other functional forms (for orbital overlap) are OK too							
							
							orbital += orbitalOverlapTerm;
						}
					}
				}
			
				// Real-space Coulomb component
				double alpha = 0;
				for (int u = -aVnum; u <= aVnum; u++) {
					for (int v = -bVnum; v <= bVnum; v++) {
						for (int w = -cVnum; w <= cVnum; w++) {
							double dx = Pos[i].x - Pos[j].x + u*aV[0] + v*bV[0] + w*cV[0];
							double dy = Pos[i].y - Pos[j].y + u*aV[1] + v*bV[1] + w*cV[1];
							double dz = Pos[i].z - Pos[j].z + u*aV[2] + v*bV[2] + w*cV[2];
							double Rab = sqrt(dx*dx + dy*dy + dz*dz);
				
							alpha += erfc( Rab / eta ) / Rab;
						}
					}
				}
				
				// K-space component
				double beta = 0;
				double h = 0; double b = 0; 
				vector<double> RLV(3); // reciprocal lattice vector
				for (int u = -hVnum; u <= hVnum; u++) {
					for (int v = -jVnum; v <= jVnum; v++) {
						for (int w = -kVnum; w <= kVnum; w++) {
							if ((u==0) && (v==0) && (w==0)) {
								// do nothing
							} else {
								RLV[0] = u*hV[0] + v*jV[0] + w*kV[0];
								RLV[1] = u*hV[1] + v*jV[1] + w*kV[1];
								RLV[2] = u*hV[2] + v*jV[2] + w*kV[2];
								
								h = Mag(RLV);
								b = 0.5 * h * eta;
								
								double dx = Pos[i].x - Pos[j].x;
								double dy = Pos[i].y - Pos[j].y;
								double dz = Pos[i].z - Pos[j].z;
								
								beta += cos( RLV[0]*dx + RLV[1]*dy + RLV[2]*dz ) / (h*h) * exp(-b*b);
							}
						}
					}
				}
				beta *= 4*PI / unitCellVolume; 
				
				return lambda * (k/2) * (alpha + beta + orbital);
			}
		}
	} else {
		cout << "Serious error specifying periodic boundary conditions. Exiting" << endl;
		exit(1);
	}
}
/*****************************************************************************/
void LoadChargeCenters(string filename) {
	// Loads charge centers to be used, atoms are assumed to be 

	ifstream fileInput(filename.c_str(),ios::in);
	string tmp, tStr; 
	int sInd, Z;

	if(!fileInput) { // Error checking
		printf("%s is not a valid filename\n\n", filename.c_str());
		exit(1);
	}

	while(!fileInput.eof()) { 
		getline(fileInput, tmp); // Read line-by-line
		
		// Read atom symbol
		sInd = tmp.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", 0);
		tStr = tmp.substr(sInd, 2);
		if (tStr[1] == '\t') tStr[1] = ' '; // Convert tabs to spaces
		Z = s_mapStringAtomLabels[tStr]; // Get Z number from label
		 
		// Read charge center (must be a positive integer)
		sInd = tmp.find_first_of("0123456789", sInd);
		tStr = tmp.substr(sInd, 1);
		IonizationData[Z].chargeCenter = atoi(tStr.c_str());
	}

}
/*****************************************************************************/
void LoadIonizationData(string filename) {
	// Loads ionization data into a global vector called IonizationData

	ifstream fileInput(filename.c_str(),ios::in);
	string data, tmp, cStr, tStr;
	int sInd = 0;
	int eInd = 0;

	if(!fileInput) { // Error checking
		printf("%s is not a valid filename\n\n", filename.c_str());
		exit(1);
	}

	while(!fileInput.eof()) { // Read file into a gigantic string
		getline(fileInput, tmp);
		data += tmp;
		data += "\n";
	}

	for (int i = 0; i < TABLE_OF_ELEMENTS_SIZE; i++) {
		// The Atom Label
		sInd = data.find("\t",sInd) + 1;
		eInd = data.find("\t",sInd);
		cStr = data.substr(sInd, eInd - sInd);
		if (cStr.length() == 1) cStr += " ";
		IonizationData[i].Label = cStr;
		
		// The Data Status
		sInd = data.find("\t",sInd) + 1;
		eInd = data.find("\t",sInd);
		cStr = data.substr(sInd, eInd - sInd);

		// The Electron Affinity
		sInd = data.find("\t",sInd) + 1;
		eInd = data.find("\t",sInd);
		cStr = data.substr(sInd, eInd - sInd);
		if ((cStr.find("<0.5") >= 0) && (cStr.find("<0.5") < 100)) {
			IonizationData[i].isDataAvailable[0] = true;
			IonizationData[i].ionizationPotential[0] = 0.5;
		} else
		if ((cStr.find("na") >= 0) && (cStr.find("na") < 100)){
			IonizationData[i].isDataAvailable[0] = false;
			IonizationData[i].ionizationPotential[0] = 0.0;
		} else {
			IonizationData[i].isDataAvailable[0] = true;
			IonizationData[i].ionizationPotential[0] = atof( cStr.c_str() );
		}
		
		int trueCount = 0;
		for (int j  = 1; j < 9; j++) {
			// The J'th Ionization Potential
			sInd = data.find("\t",sInd) + 1;
			eInd = data.find("\t",sInd);
			cStr = data.substr(sInd, eInd - sInd);
			if ((cStr.find("na") >= 0) && (cStr.find("na") < 100)){
				IonizationData[i].isDataAvailable[j] = false;
				IonizationData[i].ionizationPotential[j] = 0.0;
			} else
			if ((cStr.find("np") >= 0) && (cStr.find("np") < 100)){
				IonizationData[i].isDataAvailable[j] = false;
				IonizationData[i].ionizationPotential[j] = 0.0;
			} else {
				IonizationData[i].isDataAvailable[j] = true;
				IonizationData[i].ionizationPotential[j] = atof( cStr.c_str() );
				trueCount++;
			}
		}
		
		sInd = data.find("\n",sInd) + 1; // Go to next line
	}
}
/*****************************************************************************/
void LoadCIFFile(string filename) {
	// Two string index variables used for generating substrings from larger strings

	ifstream fileInput(filename.c_str(),ios::in);
	string data, tmp;

	if(!fileInput) { // Error checking
		printf("%s is not a valid filename\n\n", filename.c_str());
		exit(1);
	}

	while(!fileInput.eof()) { // Read file into a gigantic string
		getline(fileInput, tmp);
		data += tmp;
		data += "\n";
	}
	
	string cStr; // current string
	string tStr; // temp string
	int sInd = 0;
	int eInd = 0;

	// Read in unit cell dimensions
	sInd = data.find("_cell_length_a") + 15;
	eInd = data.find("\n", sInd);
	cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
	aLength = atof( cStr.c_str() );
	
	sInd = data.find("_cell_length_b") + 15;
	eInd = data.find("\n", sInd);
	cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
	bLength = atof( cStr.c_str() );

	sInd = data.find("_cell_length_c") + 15;
	eInd = data.find("\n", sInd);
	cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
	cLength = atof( cStr.c_str() );

	// Read in unit cell angles
	sInd = data.find("_cell_angle_alpha") + 18;
	eInd = data.find("\n", sInd);
	cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
	alphaAngle = atof( cStr.c_str() );

	sInd = data.find("_cell_angle_beta") + 17;
	eInd = data.find("\n", sInd);
	cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
	betaAngle = atof( cStr.c_str() );

	sInd = data.find("_cell_angle_gamma") + 18;
	eInd = data.find("\n", sInd);
	cStr = data.substr(sInd, eInd - sInd); // Read in the number of atoms in the file
	gammaAngle = atof( cStr.c_str() );

	// Convert to radians
	alphaAngle *= (PI / 180.0);
	betaAngle *= (PI / 180.0);
	gammaAngle *= (PI / 180.0);

	// Initialize unit cell vectors from |a|,|b|,|c| and alphaAngle, betaAngle, gammaAngle information
	// Here we are applying the A along x-axis, B in xy plane convention
	aV[0] = aLength; aV[1] = 0; aV[2] = 0;
	bV[0] = bLength*cos(gammaAngle); bV[1] = bLength*sin(gammaAngle); bV[2] = 0;
	cV[0] = cLength*cos(betaAngle);
	cV[1] = (cLength*bLength*cos(alphaAngle) - bV[0]*cV[0])/bV[1]; 
	cV[2] = sqrt(cLength*cLength - cV[0]*cV[0] - cV[1]*cV[1]);
	
	if (useEwardSums == true) DetermineReciprocalLatticeVectors();
	
	// Unitcell Volume
	vector<double> crs;
	crs = Cross(bV,cV);
	unitCellVolume = fabs( aV[0]*crs[0] + aV[1]*crs[1] + aV[2]*crs[2] ); // Volume of a parallelipiped

	// Find first line that does not contain underscore
	bool underscoreFound = true;
	int eInd2 = eInd; // we need another index
	while (underscoreFound == true) {
		sInd = eInd2; // End of the previous line
		eInd2 = data.find("\n", eInd2 + 1); // End of the next line
		cStr = data.substr(sInd, eInd2 - sInd); // The line
		if ((cStr.find("_",0) >=0) && (cStr.find("_",0) < cStr.size())) {
			underscoreFound = true; // Under score found, skip to the next line
		} else {
			// Underscore not found, we are on a legitimate line of data
			underscoreFound = false;
		}
	}
	
	cout << "==================================================" << endl;
	cout << "========= Atom types - X & J values used =========" << endl;
	cout << "==================================================" << endl;

	// Read in atom positions, symbols, and names
	Coordinates tempAtom;
	while (underscoreFound == false) {
		if ((cStr.find("_",0) >=0) && (cStr.find("_",0) < cStr.size())) {
			underscoreFound = true; // Under score found, skip to the next line
		} else {
			// Underscore not found, we are on a legitimate line of data
			underscoreFound = false;
			
			//Read atom label
			sInd = cStr.find_first_of(" \t", 1);
			tStr = cStr.substr(1, sInd-1);
			//tempAtom.label = tStr;
			Label.push_back(tStr);

			// Read atom symbol
			sInd = cStr.find_first_of("ABCDEFGHIJKLMNOPQRSTUVWXYZ", sInd);
			eInd = sInd + 1; 
			tStr = cStr.substr(sInd, eInd - sInd + 1);
			Symbol.push_back(tStr);
			
			// Find first "x" coordinate
			sInd = cStr.find(".",sInd) - 2;
			eInd = cStr.find_first_of(" \t",sInd + 2);
			tStr = cStr.substr(sInd, eInd - sInd);
			tempAtom.x = atof( tStr.c_str() );	// X Position
			
			// Find first "y" coordinate
			sInd = cStr.find(".",eInd) - 2;
			eInd = cStr.find_first_of(" \t",sInd + 2);
			tStr = cStr.substr(sInd, eInd - sInd);
			tempAtom.y = atof( tStr.c_str() );	// Y Position

			// Find first "z" coordinate
			sInd = cStr.find(".",eInd) - 2;
			eInd = cStr.find_first_of(" \t\n",sInd + 2);
			tStr = cStr.substr(sInd, eInd - sInd);
			tempAtom.z = atof( tStr.c_str() );	// Z Position
			
			// Change from fractional to cartesian:
			tempAtom.x = tempAtom.x * aV[0] + tempAtom.y * bV[0] + tempAtom.z * cV[0];
			tempAtom.y = tempAtom.x * aV[1] + tempAtom.y * bV[1] + tempAtom.z * cV[1];
			tempAtom.z = tempAtom.x * aV[2] + tempAtom.y * bV[2] + tempAtom.z * cV[2];

			Pos.push_back(tempAtom);
			
			int i = Symbol.size() - 1;
			int Z = s_mapStringAtomLabels[Symbol[i]]; // Get Z number from label
			
			if (Symbol[i] == "H ") {
				X.push_back(0.5*(hI1 + hI0));
				J.push_back(hI1 - hI0);
			} else {
				int cC = IonizationData[Z].chargeCenter;
				X.push_back(0.5*(IonizationData[Z].ionizationPotential[cC+1] + 
					IonizationData[Z].ionizationPotential[cC]));
				J.push_back(IonizationData[Z].ionizationPotential[cC+1] - 
					IonizationData[Z].ionizationPotential[cC]);
				X[i] -= cC*(J[i]);
			}
			
			bool beenDone = false;
			for (int j = 0; j < i; j++) {
				if (Symbol[i] == Symbol[j]) beenDone = true;
			}
			if (beenDone == false) {
				cout << Symbol[i] << "\t";
				cout << "Z: " << Z+1 << "\t";
				cout << "Ch. Cent: " << IonizationData[Z].chargeCenter << "\t";
				cout << "X: " << X[i] << "\t";
				cout << "J: " << J[i] << "\t" << endl;
			}
		}
		sInd = eInd2; // End of the previous line
		eInd2 = data.find("\n", eInd2 + 1); // End of the next line
		cStr = data.substr(sInd, eInd2 - sInd); // The line
	}

	numAtoms = Pos.size();
	
	Q.resize(numAtoms, 0); // initialize charges to zero
}

// #include <sstream>
// #include <cctype>

// static inline std::string trim(const std::string& s) {
//     size_t b = s.find_first_not_of(" \t\r\n");
//     if (b == std::string::npos) return "";
//     size_t e = s.find_last_not_of(" \t\r\n");
//     return s.substr(b, e - b + 1);
// }

// void LoadCIFFile(string filename) {
//     ifstream fileInput(filename.c_str(), ios::in);
//     if(!fileInput) {
//         printf("%s is not a valid filename\n\n", filename.c_str());
//         exit(1);
//     }

//     // Slurp whole file (we still do this for cell params; later we stream lines)
//     string data, tmp;
//     while (getline(fileInput, tmp)) { data += tmp; data += "\n"; }

//     // -------------------------
//     // Unit cell lengths/angles
//     // -------------------------
//     auto fetch_after = [&](const string& key)->double{
//         size_t s = data.find(key);
//         if (s == string::npos) { cerr << "CIF key not found: " << key << endl; exit(1); }
//         s += key.size();
//         size_t e = data.find("\n", s);
//         string v = trim(data.substr(s, e - s));
//         return atof(v.c_str());
//     };

//     aLength = fetch_after("_cell_length_a");
//     bLength = fetch_after("_cell_length_b");
//     cLength = fetch_after("_cell_length_c");
//     alphaAngle = fetch_after("_cell_angle_alpha");
//     betaAngle  = fetch_after("_cell_angle_beta");
//     gammaAngle = fetch_after("_cell_angle_gamma");

//     // radians
//     alphaAngle *= (PI/180.0);
//     betaAngle  *= (PI/180.0);
//     gammaAngle *= (PI/180.0);

//     // Lattice vectors (A along x; B in xy plane)
//     aV[0] = aLength; aV[1] = 0.0; aV[2] = 0.0;
//     bV[0] = bLength * cos(gammaAngle);
//     bV[1] = bLength * sin(gammaAngle);
//     bV[2] = 0.0;
//     cV[0] = cLength * cos(betaAngle);
//     cV[1] = (cLength*bLength*cos(alphaAngle) - bV[0]*cV[0]) / bV[1];
//     cV[2] = sqrt(max(0.0, cLength*cLength - cV[0]*cV[0] - cV[1]*cV[1]));

//     if (useEwardSums == true) DetermineReciprocalLatticeVectors();

//     // Unit cell volume
//     vector<double> crs = Cross(bV, cV);
//     unitCellVolume = fabs(aV[0]*crs[0] + aV[1]*crs[1] + aV[2]*crs[2]);

//     // -------------------------
//     // Re-scan line-by-line for atom loop
//     // -------------------------
//     istringstream iss(data);
//     string line;

//     // We will start reading atoms after the header loop of _atom_site_* labels ends.
//     bool in_atom_section = false;
//     vector<string> atom_headers;

//     cout << "==================================================\n";
//     cout << "========= Atom types - X & J values used =========\n";
//     cout << "==================================================\n";

//     while (std::getline(iss, line)) {
//         string L = trim(line);
//         if (L.empty()) continue;           // skip blanks
//         if (L[0] == '#') continue;         // skip comments

//         // Header lines:
//         if (L.rfind("loop_", 0) == 0) {    // starts with "loop_"
//             atom_headers.clear();
//             in_atom_section = false;
//             continue;
//         }
//         if (L.size() > 0 && L[0] == '_') {
//             // Could be _atom_site_* header
//             if (L.rfind("_atom_site_", 0) == 0) {
//                 atom_headers.push_back(L);
//                 in_atom_section = true;    // we are collecting headers
//             }
//             continue;
//         }

//         // If we have collected _atom_site_ headers and now see a non-underscore line,
//         // treat it as a data row.
//         if (in_atom_section) {
//             // Tokenize row (CIF allows quoted strings; most simple MOF CIFs are whitespace-separated)
//             std::istringstream ls(L);
//             std::vector<std::string> f; std::string tok;
//             while (ls >> tok) f.push_back(tok);

//             // We need at least label, type_symbol, fract_x, fract_y, fract_z
//             if (f.size() < 5) continue; // skip malformed lines gracefully

//             string label = f[0];
//             string sym   = f[1];
//             // normalize symbol to map keys used by s_mapStringAtomLabels
//             if (sym.size() == 1) sym += " ";
//             else sym = sym.substr(0,2); // e.g., "Cl", "Zn"

//             // Accept numbers that might have uncertainties, e.g. 0.123(4)
//             auto parse_frac = [](const std::string& s)->double{
//                 // strip any parenthetical uncertainty
//                 size_t p = s.find('(');
//                 string core = (p == string::npos) ? s : s.substr(0, p);
//                 return atof(core.c_str());
//             };

//             Coordinates tempAtom;
//             tempAtom.x = parse_frac(f[2]);
//             tempAtom.y = parse_frac(f[3]);
//             tempAtom.z = parse_frac(f[4]);

//             // fractional -> cartesian with temporaries
//             double fx = tempAtom.x, fy = tempAtom.y, fz = tempAtom.z;
//             double cx = fx * aV[0] + fy * bV[0] + fz * cV[0];
//             double cy = fx * aV[1] + fy * bV[1] + fz * cV[1];
//             double cz = fx * aV[2] + fy * bV[2] + fz * cV[2];
//             tempAtom.x = cx; tempAtom.y = cy; tempAtom.z = cz;

//             Label.push_back(label);
//             Symbol.push_back(sym);
//             Pos.push_back(tempAtom);

//             int i = (int)Symbol.size() - 1;
//             // lookup Z; if unknown, fail early with message
//             if (s_mapStringAtomLabels.find(Symbol[i]) == s_mapStringAtomLabels.end()) {
//                 cerr << "Unknown element symbol in CIF: '" << Symbol[i]
//                      << "' (raw: '" << f[1] << "'). Fix symbol mapping." << endl;
//                 exit(1);
//             }
//             int Z = s_mapStringAtomLabels[Symbol[i]];

//             if (Symbol[i] == "H ") {
//                 X.push_back(0.5*(hI1 + hI0));
//                 J.push_back(hI1 - hI0);
//             } else {
//                 int cC = IonizationData[Z].chargeCenter;
//                 X.push_back(0.5*(IonizationData[Z].ionizationPotential[cC+1] +
//                                  IonizationData[Z].ionizationPotential[cC]));
//                 J.push_back(IonizationData[Z].ionizationPotential[cC+1] -
//                             IonizationData[Z].ionizationPotential[cC]);
//                 X[i] -= cC * (J[i]);
//             }

//             bool seen = false;
//             for (int j = 0; j < i; j++) if (Symbol[i] == Symbol[j]) { seen = true; break; }
//             if (!seen) {
//                 cout << Symbol[i] << "\t";
//                 cout << "Z: " << (int)Z + 1 << "\t";
//                 cout << "Ch. Cent: " << IonizationData[Z].chargeCenter << "\t";
//                 cout << "X: " << X[i] << "\t";
//                 cout << "J: " << J[i] << "\t" << endl;
//             }
//         }
//     }

//     numAtoms = (int)Pos.size();
//     if (numAtoms == 0) {
//         cerr << "No atoms parsed from CIF. Check that your CIF has an _atom_site loop." << endl;
//         exit(1);
//     }
//     Q.resize(numAtoms, 0.0);
// }

/*****************************************************************************/
void OutputCIFFormatFile(string filename) {
	FILE *out; string str = "";
	int fileFormat = 2;
	out = fopen(filename.c_str(),"wt");

	str += "data_functionalizedCrystal"; str	+= "\n";
	str += "_audit_creation_method\t"; str += "'EQeq! by Chris Wilmer'"; str += "\n";
	str += "_symmetry_space_group_name_H-M\t"; str += "'P1'"; str += "\n"; 
	str += "_symmetry_Int_Tables_number\t"; str += "1"; str += "\n"; 
	str += "_symmetry_cell_setting\t"; str += "triclinic"; str += "\n"; 
	str += "loop_"; str += "\n";
	str += "_symmetry_equiv_pos_as_xyz"; str += "\n";
	str += "  x,y,z"; str += "\n";
	fprintf(out,str.c_str());
	fprintf(out,"_cell_length_a\t%f\n",aLength);
	fprintf(out,"_cell_length_b\t%f\n",bLength);
	fprintf(out,"_cell_length_c\t%f\n",cLength);
	fprintf(out,"_cell_angle_alpha\t%f\n",alphaAngle*(180 / PI));
	fprintf(out,"_cell_angle_beta\t%f\n",betaAngle*(180 / PI));
	fprintf(out,"_cell_angle_gamma\t%f\n",gammaAngle*(180 / PI));
	str = "";
	str += "loop_"; str += "\n";
	str += "_atom_site_label"; str += "\n";
	str += "_atom_site_type_symbol"; str += "\n";
	str += "_atom_site_fract_x"; str += "\n";
	str += "_atom_site_fract_y"; str += "\n";
	str += "_atom_site_fract_z"; str += "\n";
	str += "_atom_site_charge"; str += "\n";
	fprintf(out,str.c_str());
    
	int k = 0;
    
	// For all atoms
	for (int i = 0; i < numAtoms ; i++) {
		k++;
		// Determine the fractional coordinates
		double dx = Pos[i].x;
		double dy = Pos[i].y;
		double dz = Pos[i].z;
        
		// Convert to fractional coordinates (below is the "inverse transform matrix") 
		double a = (bV[2]*cV[1]*dx - bV[1]*cV[2]*dx - bV[2]*cV[0]*dy + bV[0]*cV[2]*dy + bV[1]*cV[0]*dz - bV[0]*cV[1]*dz)/
		           (aV[2]*bV[1]*cV[0] - aV[1]*bV[2]*cV[0] - aV[2]*bV[0]*cV[1] + 
				    aV[0]*bV[2]*cV[1] + aV[1]*bV[0]*cV[2] - aV[0]*bV[1]*cV[2]);
		double b = (aV[2]*cV[1]*dx - aV[1]*cV[2]*dx - aV[2]*cV[0]*dy + aV[0]*cV[2]*dy + aV[1]*cV[0]*dz - aV[0]*cV[1]*dz)/
		           (-(aV[2]*bV[1]*cV[0]) + aV[1]*bV[2]*cV[0] + aV[2]*bV[0]*cV[1] - 
				   aV[0]*bV[2]*cV[1] - aV[1]*bV[0]*cV[2] + aV[0]*bV[1]*cV[2]);
		double c = (aV[2]*bV[1]*dx - aV[1]*bV[2]*dx - aV[2]*bV[0]*dy + aV[0]*bV[2]*dy + aV[1]*bV[0]*dz - aV[0]*bV[1]*dz)/
		           (aV[2]*bV[1]*cV[0] - aV[1]*bV[2]*cV[0] - aV[2]*bV[0]*cV[1] + 
				   aV[0]*bV[2]*cV[1] + aV[1]*bV[0]*cV[2] - aV[0]*bV[1]*cV[2]);
        
		fprintf(out,"cg");
		fprintf(out,Symbol[i].c_str());
		fprintf(out,"\t");
		fprintf(out,Symbol[i].c_str());
		fprintf(out,"\t%f\t%f\t%f\t%f\n",a,b,c,Q[i]);
	}
	
	fprintf(out,"_end\n");
	fclose(out);
}
/*****************************************************************************/
void OutputPDBFormatFile(string filename) {
	FILE *out;
	out = fopen(filename.c_str(),"wt");
	
	fprintf(out,"TITLE       YourMoleculeNameHere            \n");
	fprintf(out,"REMARK   4\n");
	fprintf(out,"REMARK   4      COMPLIES WITH FORMAT V. 2.2, 16-DEC-1996\n");
	if (isPeriodic == true) {
		fprintf(out,"CRYST1    %5.2f    %5.2f    %5.2f  %3.2f  %3.2f  %3.2f P1\n",
			aLength,bLength,cLength,alphaAngle*180/PI,betaAngle*180/PI,gammaAngle*180/PI);
	} else {
		// Do nothing
	}
	for (int i = 0; i < numAtoms; i++) {
		fprintf(out,"ATOM    %3d %s   MOL A   0     % 7.3f % 7.3f % 7.3f % 5.2f                %s\n",
			i+1,Symbol[i].c_str(),Pos[i].x,Pos[i].y,Pos[i].z,Q[i],Symbol[i].c_str());
	}
	
	fclose(out);
}
/*****************************************************************************/
void OutputMOLFormatFile(string filename) {
	FILE *out;
	out = fopen(filename.c_str(),"wt");
	
	fprintf(out," Molecule_name: hypotheticalMOF\n"); // This should be updated
	fprintf(out,"\n");
	fprintf(out,"  Coord_Info: Listed Cartesian None\n");
	fprintf(out,"        %d\n",numAtoms);
	
	for (int i = 0; i < numAtoms; i++) {
		fprintf(out,"  %4d  % 8.4f % 8.4f % 8.4f  Mof_%s   % 6.3f  0  0\n",
			i+1,Pos[i].x,Pos[i].y,Pos[i].z,Symbol[i].c_str(),Q[i]);
	}
	
	fprintf(out,"\n");
	fprintf(out,"\n");
	fprintf(out,"\n");
	fprintf(out,"  Fundcell_Info: Listed\n");
	fprintf(out,"        %8.4f      %8.4f      %8.4f\n",aLength,bLength,cLength);
	fprintf(out,"        %8.4f      %8.4f      %8.4f\n",alphaAngle*180/PI,betaAngle*180/PI,gammaAngle*180/PI);
	fprintf(out,"        0.00000        0.00000       0.00000\n");
	fprintf(out,"        %8.4f      %8.4f      %8.4f\n",aLength,bLength,cLength);
	
	fclose(out);
}
/*****************************************************************************/
void Qeq() {
	int i, j; // generic counter;

	// Formulate problem in the form of A x = b
	vector<double> dummyRow(numAtoms, 0); // Is doing this necessary?
	vector<vector<double> > A(numAtoms, dummyRow);
	vector<double> b(numAtoms,0);

	// First row of A is all ones
	for (int i = 0; i < numAtoms; i++) {
		A[0][i] = 1;
	}

	// First element in b is the total charge
	b[0] = Qtot;

	// Rest of elements in b are the differences in electronegativity
	for (int i = 1; i < numAtoms; i++) {
		b[i] = X[i] - X[i-1];
	}

	// Fill in 2nd to Nth rows of A
	for (int i = 1; i < numAtoms; i++) {
		cout << ".";
		for (int j = 0; j < numAtoms; j++) {
			A[i][j] = GetJ(i-1, j) - GetJ(i, j);
		}
	}

	Q = SolveMatrix(A,b);
}
/*****************************************************************************/
void RoundCharges(int digits) {
	
	double qsum = 0;
	double factor = pow((double)10,digits);
	
	for(int i=0; i < numAtoms; i++) {
		Q[i] = Round(Q[i]*factor)/factor;
		qsum += Q[i];
	}
	
	if (qsum == 0) { // Great, rounding worked on the first try!
		// do nothing
	} else { // There is a small excess charge from rounding, adjust it
		int numAtomsToAdjust = (int)(abs(qsum * factor) + 0.5); // Weird double-to-int conversion tricks
		cout << " adjusting the charge of " << numAtomsToAdjust << " atoms!" << endl;
		
		int sign; if (qsum > 0) sign = -1; else sign = 1;
		for (int i=0; i < numAtomsToAdjust; i++) { // Adjust
			Q[i] += sign*(1/factor);
		}
	}
	
}
/*****************************************************************************/
vector<double> Cross(vector<double> a, vector<double> b) {

	vector<double> c(3);

	c[0] = a[1]*b[2] - a[2]*b[1];
	c[1] = a[2]*b[0] - a[0]*b[2];
	c[2] = a[0]*b[1] - a[1]*b[0];

	return c;
}
/*****************************************************************************/
double Dot(vector<double> a, vector<double> b) {
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
/*****************************************************************************/
double Mag(vector<double> a) {
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}
/*****************************************************************************/
double Round(double num) {
	return (num > 0.0) ? floor(num + 0.5) : ceil(num - 0.5);
}
/*****************************************************************************/
vector<double> Scalar(double a, vector<double> b) {
	vector<double> c(3);
	c[0] = a*b[0]; c[1] = a*b[1]; c[2] = a*b[2];
	return c;
}
/*****************************************************************************/
vector<double> SolveMatrix(vector<vector<double> > A, vector<double> b) {
	// Assumptions: A x = b, A is a MxN matrix, M = rows, N = cols, x is a vector, b is vector
	// matrix has more rows than columns
	// number of rows of matrix is equal to size of the vector x

    // Initialize x = b
    vector<double> x;
    x = b;

    int i, j, k;
    int N = A.size();
    int M = A[0].size();

    vector<double> d (N);

    /* Perform Householder transformation */
    for (i = 0; i < N; i++) {
        const double aii = A[i][i];
        double alef, f, ak;
        double max_norm = 0.0;
        double r = 0.0;

        for (k = i; k < M; k++) {
          r += A[k][i] * A[k][i];
        }

        if (r == 0) {
          cout << "Error! Matrix is rank deficient." << endl;
          // return -1;
        }

        if (A[i][i] < 0)
          alef = (-1)*sqrt(r);
        else
          alef = sqrt(r);

        ak = 1.0 / (r + alef * A[i][i]);

        A[i][i] +=  alef;

        d[i] = -alef;

        for (k = i + 1; k < N; k++) {
          double norm = 0.0;
          f = 0.0;

          for (j = i; j < M; j++) {
            norm += A[j][k] * A[j][k];
            f += A[j][k] * A[j][i];
          }

          max_norm = max(max_norm, norm);

          f *= ak;

          for (j = i; j < M; j++) {
            A[j][k] -= f * A[j][i];
          }
        }

        if (fabs(alef) < 0.00001) {
          cout << "Apparent singularity in matrix." << endl;
          // return -1;
        }

        f = 0.0;

        for (j = i; j < M; j++) {
          f += x[j] * A[j][i];
        }

        f *= ak;

        for (j = i; j < M; j++) {
          x[j] -= f * A[j][i];
        }
    }

    /* Perform back-substitution */

    for (i = N-1; i >= 0; i--) {
      double sum = 0.0;

      for (k = i + 1; k < N; k++) {
        sum += A[i][k] * x[k];
      }

      x[i] = (x[i] - sum) / d[i] ;
    }

    return x;
}
/*****************************************************************************/