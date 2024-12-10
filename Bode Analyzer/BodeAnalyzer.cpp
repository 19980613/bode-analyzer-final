#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <math.h>
#include "BodeAnalyzer.h"
using namespace std;


/*
 * Einlesen der Zähler- und Nennerpolynome
 * grad: Der Grad des jeweiligen Polynoms
 * polynomName: "Zähler" oder "Nenner"
 * Rückgabewert: Ein Vektor mit den Koeffizienten des Polynoms
 */
vector<double> readCoefficients(int grad, const string& polynomName) {
	vector<double> coefficients(grad + 1);
	cout << "Geben Sie die Koeffizienten des " << polynomName << "-Polynoms an:" << endl;

	for (int i = grad; i >= 0; i--) {
		cout << "Koeffizient für s^" << i << ": ";
		cin >> coefficients[i];
	}
	return coefficients;
}

/*
 * Berechnet Frequenzpunkte, die gleichmäßig im Bereich [fmin, fmax] verteilt sind.
 * fmin: Die minimale Frequenz
 * fmax: Die maximale Frequenz
 * fpointsVar: Die Anzahl der Frequenzpunkte
 * Rückgabewert: Ein Vektor mit den berechneten Frequenzpunkten.
 */
vector<double> calculateFrequencyPoints(double fmin, double fmax, int fpointsVar) {
	vector<double> fpoints(fpointsVar);

	for (int i = 0; i < fpointsVar; ++i) {
		fpoints[i] = fmin + static_cast<double>(i) * (fmax - fmin) / (fpointsVar - 1);
	}

	return fpoints;
}

/* 
 * Phasensprünge berechnen 
 * currentPhase: Ergebnis der aktuellen Phasenberechnung 
 * previousPhase: Ergebnis der vorherigen Phasenberechnung
 * Rückgabewert: Berechnete Phase
 */
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


/* 
 * Horner Schema zur Auswertung des Polynoms 
 * currentCoeff: Koeffizienten des jeweiligen Polynoms
 * s: Komplexe Variable
 * Rückgabewert: Auswertung des Polynoms an einer Stelle
 */
complex<double> hornerScheme(const vector<double>& currentCoeff, complex<double> s) {
    complex<double> result = (0.0, 0.0);
    // Multipliziere das aktuelle Ergebnis mit x und 
    // addiere den nächsten Koeffizienten
    for (int i = currentCoeff.size() - 1; i >= 0; i--) {
        result = result * s + currentCoeff[i];
    }
    return result;
}


/* 
 * Prüfung auf gleiche Vorzeichen 
 * currentCoeff: Koeffizienten des jeweiligen Polynoms
 * Rückgabewert: True/False 
 */
bool criterionHurwitzCondition1(const vector<double>& currentCoeff) {

    /* Bedingung 1: Alle Koeff haben das gleiche Vorzeichen */
    
    // Variable ist true, wenn Vorzeichen negativ ist
    bool erstesVorzeichenNegativ = signbit(currentCoeff[0]);

    for (int i = 0; i < currentCoeff.size() - 1; i++) {
        if (signbit(currentCoeff[i]) != erstesVorzeichenNegativ) return false;
    }
    return true;
}


/* 
 * Hurwitz Matrix aufstellen 
 * coeffPolynom: Koeffizienten des Nennerpolynoms
 * Rückgabewert: aufgestellte Hurwitz-Matrix
 */
vector<vector<double>> criterionHurwitzCondition2(const vector<double>& coeffPolynom) {

	int anzahlCoeff = coeffPolynom.size();
	int dimension = 3;

	vector<vector<double>> hurwitzMatrix(dimension, vector<double>(dimension, 0.0));

	// Zeile 1
	if (anzahlCoeff >= 2)	hurwitzMatrix[0][0] = coeffPolynom[1];
	else hurwitzMatrix[0][0] = 0;

	if (anzahlCoeff >= 4)	hurwitzMatrix[0][1] = coeffPolynom[3];
	else hurwitzMatrix[0][1] = 0;

	if (anzahlCoeff >= 6)    hurwitzMatrix[0][2] = coeffPolynom[5];
	else hurwitzMatrix[0][2] = 0;

	// Zeile 2
	if (anzahlCoeff >= 0)	hurwitzMatrix[1][0] = coeffPolynom[0];
	else hurwitzMatrix[1][0] = 0;

	if (anzahlCoeff >= 3)	hurwitzMatrix[1][1] = coeffPolynom[2];
	else hurwitzMatrix[1][1] = 0;

	if (anzahlCoeff >= 5)	hurwitzMatrix[1][2] = coeffPolynom[4];
	else hurwitzMatrix[1][2] = 0;

	// Zeile 3
	if (anzahlCoeff >= 0)	hurwitzMatrix[2][0] = 0;

	if (anzahlCoeff >= 2)	hurwitzMatrix[2][1] = coeffPolynom[1];
	else hurwitzMatrix[2][1] = 0;

	if (anzahlCoeff >= 4)	hurwitzMatrix[2][2] = coeffPolynom[3];
	else hurwitzMatrix[2][2] = 0;


	return hurwitzMatrix;
}


/* 
 * Determinanten berechnen
 * hurwitzMatrix: Übergabe der Parameter der HM
 * koeffNenner: Übergabe der Anzahl der Koeffizienten des Nennerpolynoms
 * Rückgabewert: Determinanten
 */
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


/* 
 * Erstellt ein GNUplot-Skript zur Visualisierung des Bode-Diagramms
 * csvFile: Name der CSV-Datei mit den Daten
 * zaehlerFormel: Zähler der Transferfunktion 
 * nennerFormel: Nenner der Transferfunktion 
 * phaseMargin: Phasenreserve
 * magMargin: Amplitudenreserve
 */
void generateGnuplotScript(const string& csvFile, const string& zaehlerFormel, const string& nennerFormel,
    double phaseMargin, double magMargin) {
    ofstream script("bode_plot.gp");

    // Gnuplot-Terminal und Ausgabe-Datei definieren
    script << "set terminal png size 900,600\n";
    script << "set output 'bode_plot.png'\n";

    // Logarithmische x-Achse und Gitter
    script << "set logscale x\n";
    script << "set grid\n";

    // Achsenbeschriftungen
    script << "set xlabel 'Frequency (Hz)'\n";
    script << "set ylabel 'Amplitude (dB)'\n";
    script << "set y2label 'Phase (Degrees)'\n";
    script << "set ytics nomirror\n";
    script << "set y2tics\n";

    // Diagrammränder anpassen
    script << "set lmargin at screen 0.1\n";      // Linker Rand
    script << "set rmargin at screen 0.9\n";      // Rechter Rand
    script << "set tmargin at screen 0.8\n";      // Oberer Rand
    script << "set bmargin at screen 0.1\n";      // Unterer Rand

    // Transferfunktion beschreiben
    script << "set label 1 'G(s) = " << zaehlerFormel << " / " << nennerFormel << "' at screen 0.25, 0.9 center\n";

    // Phasen- und Amplitudenreserve anzeigen
    script << "set label 2 'Phasenreserve (Grad) = " << phaseMargin << "' at screen 0.6, 0.93 left\n";
    script << "set label 3 'Amplitudenreserve (dB) = " << magMargin << "' at screen 0.6, 0.89 left\n";

    // Daten plotten
    script << "plot '" << csvFile << "' using 1:2 with lines title 'Amplitude (dB)', \\\n";
    script << "     '" << csvFile << "' using 1:3 axes x1y2 with lines title 'Phase (Degrees)'\n";
    script.close();
}


/*
 * Vektor von Koeffizienten -> String-Darstellung eines Polynoms
 * coeffs: Die Koeffizienten des Polynoms
 * Rückgabewert:  Polynom als String
 */
string vectorToString(const vector<double>& coeffs) {
    ostringstream stream;
    for (int i = coeffs.size() - 1; i >= 0; i--) {
        stream << coeffs[i];
        if (i > 0) stream << "s^" << i << " + ";
    }
    return stream.str();
}