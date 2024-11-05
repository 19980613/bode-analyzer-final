#ifndef _BODEANALYZER_H // Einmalige Einbindung sicherstellen
#define _BODEANALYZER_H

#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <math.h>
using namespace std;

double adjustPhaseContinuity(double currentPhase, double previousPhase);
complex<double> hornerScheme(const vector<double>& koeffZaehler, complex<double> s);
bool criterionHurwitzCondition1(const vector<double>& koeffNenner);
vector<vector<double>> criterionHurwitzCondition2(const vector<double>& koeffNenner);
vector<double> calculateDeterminante(const vector<vector<double>>& hurwitzMatrix, const vector<double>& koeffNenner);


#endif
