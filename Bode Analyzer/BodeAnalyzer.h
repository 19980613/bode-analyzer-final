#ifndef _BODEANALYZER_H // Einmalige Einbindung sicherstellen
#define _BODEANALYZER_H

#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <math.h>
using namespace std;

vector<double> readCoefficients(int grad, const string& polynomName);
vector<double> calculateFrequencyPoints(double fmin, double fmax, int fpointsVar);
double adjustPhaseContinuity(double currentPhase, double previousPhase);
complex<double> hornerScheme(const vector<double>& currentCoeff, complex<double> s);
bool criterionHurwitzCondition1(const vector<double>& currentCoeff);
vector<vector<double>> criterionHurwitzCondition2(const vector<double>& coeffPolynom);
vector<double> calculateDeterminante(const vector<vector<double>>& hurwitzMatrix, const vector<double>& koeffNenner);
void generateGnuplotScript(const string& csvFile, const string& zaehlerFormel, const string& nennerFormel,
    double phaseMargin, double magMargin);
string vectorToString(const vector<double>& coeffs);

#endif
