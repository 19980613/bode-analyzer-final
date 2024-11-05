#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <math.h>
#include "BodeAnalyzer.h"
using namespace std;


/* Phasensprünge berechnen */

double adjustPhaseContinuity(double currentPhase, double previousPhase) {
    
	// Berechne den Phasensprung
    double phaseDiff = currentPhase - previousPhase;

    // Prüfe auf Sprünge und passe an
    if (phaseDiff > 180.0) {
        currentPhase -= 360.0;
    }
    else if (phaseDiff < -180.0) {
        currentPhase += 360.0;
    }
    return currentPhase;
}


/* Horner Schema zur Auswertung des Polynoms */

complex<double> hornerScheme(const vector<double>& koeffZaehler, complex<double> s) {
    complex<double> result = (0.0, 0.0);
    // Multipliziere das aktuelle Ergebnis mit x und 
    // addiere den nächsten Koeffizienten
    for (int i = koeffZaehler.size() - 1; i >= 0; i--) {
        result = result * s + koeffZaehler[i];
    }
    return result;
}


/* Prüfung auf gleiche Vorzeichen */

bool criterionHurwitzCondition1(const vector<double>& koeffNenner) {

    /* Bedingung 1: Alle Koeff haben das gleiche Vorzeichen
    */

    // Variable ist true, wenn Vorzeichen negativ ist
    bool erstesVorzeichenNegativ = signbit(koeffNenner[0]);

    for (int i = 0; i < koeffNenner.size() - 1; i++) {
        if (signbit(koeffNenner[i]) != erstesVorzeichenNegativ) return false;
    }
    return true;
}


/* Hurwitz Matrix aufstellen */

vector<vector<double>> criterionHurwitzCondition2(const vector<double>& koeffNenner) {

	int anzahlCoeff = koeffNenner.size();
	int dimension = 3;

	vector<vector<double>> hurwitzMatrix(dimension, vector<double>(dimension, 0.0));

	// Matrix füllen 

	// Zeile 1
	if (anzahlCoeff >= 2)	hurwitzMatrix[0][0] = koeffNenner[1];
	else hurwitzMatrix[0][0] = 0;

	if (anzahlCoeff >= 4)	hurwitzMatrix[0][1] = koeffNenner[3];
	else hurwitzMatrix[0][1] = 0;

	if (anzahlCoeff >= 6)    hurwitzMatrix[0][2] = koeffNenner[5];
	else hurwitzMatrix[0][2] = 0;

	// Zeile 2
	if (anzahlCoeff >= 0)	hurwitzMatrix[1][0] = koeffNenner[0];
	else hurwitzMatrix[1][0] = 0;

	if (anzahlCoeff >= 3)	hurwitzMatrix[1][1] = koeffNenner[2];
	else hurwitzMatrix[1][1] = 0;

	if (anzahlCoeff >= 5)	hurwitzMatrix[1][2] = koeffNenner[4];
	else hurwitzMatrix[1][2] = 0;

	// Zeile 3
	if (anzahlCoeff >= 0)	hurwitzMatrix[2][0] = 0;

	if (anzahlCoeff >= 2)	hurwitzMatrix[2][1] = koeffNenner[1];
	else hurwitzMatrix[2][1] = 0;

	if (anzahlCoeff >= 4)	hurwitzMatrix[2][2] = koeffNenner[3];
	else hurwitzMatrix[2][2] = 0;


	return hurwitzMatrix;
}


/* Determinanten berechnen */

vector<double> calculateDeterminante(const vector<vector<double>>& hurwitzMatrix, const vector<double>& koeffNenner) {

	// Berechnet die Hauptminoren

	int dimension = 2;
	if (koeffNenner.size() >= 4) dimension++;
	vector<double> determinanten(dimension, 0.0);

	determinanten[0] = hurwitzMatrix[0][0];

	determinanten[1] = hurwitzMatrix[0][0] * hurwitzMatrix[1][1] -
		hurwitzMatrix[0][1] * hurwitzMatrix[1][0];

	if (koeffNenner.size() >= 4) {
		determinanten[2] = hurwitzMatrix[0][0] * (hurwitzMatrix[1][1] * hurwitzMatrix[2][2] - hurwitzMatrix[1][2] * hurwitzMatrix[2][1]) -
			hurwitzMatrix[0][1] * (hurwitzMatrix[1][0] * hurwitzMatrix[2][2] - hurwitzMatrix[1][2] * hurwitzMatrix[2][0]) +
			hurwitzMatrix[0][2] * (hurwitzMatrix[1][0] * hurwitzMatrix[2][1] - hurwitzMatrix[1][1] * hurwitzMatrix[2][0]);
	}

	return determinanten;
}