
/* Bibliotheken einbinden */
#include <iostream>
#include <vector>
#include <complex>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "BodeAnalyzer.h"
using namespace std;

/* Konstanten definieren */
#define PI 3.14159265359



int main() {

    /* Deklaration: Fmin, Fmax, Fpoints */
    double fmin;
    double fmax;
    double fpointsVar;

    /* Deklaration der Variablen für Phasen- und Amplitudenreserve */
    bool magCrossed = false;
    bool phaseCrossed = false;
    double magMargin = 0.0;
    double phaseMargin = 0.0;
    double magCrossFrequency = 0.0;
    double phaseCrossFrequency = 0.0;
    double previousPhase = 0.0;

    /* Variablen für den Grad der Zähler- und Nennerpolynome */
    int gradZaehler;
    int gradNenner;

    /* Einlesen der Koeffizienten */
    cout << "Geben Sie den Grad des Zaehlerpolynoms an: ";
    cin >> gradZaehler;
    vector<double> koeffZaehler(gradZaehler + 1);
    for (int i = gradZaehler; i >= 0; i--) {
        cout << "Koeffizient s^" << i << ": ";
        cin >> koeffZaehler[i];
     }

    cout << "\nGeben Sie den Grad des Nennerpolynoms an (max. 5. Grades): ";
    cin >> gradNenner;

    /* Fehlerbehandlung bei grad>5 */
    if (gradNenner > 5) return 0;

    vector<double> koeffNenner(gradNenner + 1);
    for (int i = gradNenner; i >= 0; i--) {
        cout << "Koeffizient s^" << i << ": ";
        cin >> koeffNenner[i];
    }

    /* Einlesen der Variablen: Fmin/Fmax/Fpoints */
    cout << "\n" << "In welchem Frequenzbereich soll die Berechnung erfolgen?" << "\n";
    cout << "fmin = ";
    cin >> fmin;
    cout << "fmax = ";
    cin >> fmax;
    cout << "fpoints = ";
    cin >> fpointsVar;

    /* Speicherung der Frequenzpunkte im Vector durch */ 
    vector<double> fpoints{};
    for (int i = 0; i < fpointsVar; i++) {
        double hilf = fmin + (double)i * (fmax - fmin) / (fpointsVar - 1);
        fpoints.push_back(hilf);
    }

    /* CSV-Datei für die Ausgabe der Amplituden- und Phasenwerte */ 
    ofstream outputFile("bode_plot_data.csv");


    /* Amplitude und Phase berechnen und in csv.File ablegen */ 
    for (int i = 0; i < fpointsVar; i++) {

        complex<double> s(0.0, fpoints[i]);   // s = jw

        complex<double> valueZaehler = hornerScheme(koeffZaehler, s);
        //cout << valueZaehler << "\n";
        complex<double> valueNenner = hornerScheme(koeffNenner, s);
        //cout << valueNenner << "\n";

        complex<double> H = valueZaehler / valueNenner;
        //cout << H << "\n";

        double mag = 20.0 * log10(abs(H));
        double phase = arg(H) * 180.0 / PI;

        phase = adjustPhaseContinuity(phase, previousPhase);
        previousPhase = phase;


        /* Suche Amplitudengrenzfrequenz */ 
        if (!magCrossed && mag <= 0.0) {
            magCrossed = true;
            magCrossFrequency = fpoints[i];
            phaseMargin = 180.0 + phase;    // Phasenreserve
        }

        /* Suche Phasengrenzfrequenz */
        if (!phaseCrossed && phase <= -180.0) {
            phaseCrossed = true;
            phaseCrossFrequency = fpoints[i];
            magMargin = -mag;               // Amplitudenreserve
        }

        /* Schreiben der Amplituden und Phasenwerte in die csv.Datei */
        outputFile << fpoints[i] << "\t" << mag << "\t" << phase << "\n";

    }

    /* Ausgabe der Amplituden- und Phasenreserve */
    cout << "\n" << "---------------------------------------------" << "\n";
    cout << "Amplitudenreserve (db): " << magMargin << endl;
    cout << "Phasengrenzfrequenz (Hz): " << phaseCrossFrequency << "\n" << endl;

    cout << "Phasenreserve (Grad): " << phaseMargin << endl;
    cout << "Amplitudengrenzfrequenz (Hz): " << magCrossFrequency << endl;
    cout << "---------------------------------------------" << "\n";

    /* Schließen der csv.Datei */
    outputFile.close();

    /* Stabilitätsprüfung nach Hurwitz */
    bool condition1 = criterionHurwitzCondition1(koeffNenner);
    vector<vector<double>> hurwitzMatrix = criterionHurwitzCondition2(koeffNenner);
    vector<double> det = calculateDeterminante(hurwitzMatrix, koeffNenner);

    /* Ausgabe der Hurwitz-Matrix */
    cout << "\n" << "Hurwitz-Matrix" << endl;
    cout << "------------------" << endl;
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; ++j) {
            cout << hurwitzMatrix[i][j] << "\t";
        }
        cout << endl;
    }

    /* Ausgabe der Determinanten */
    cout << "\n" << "Hauptdeterminaten\n";
    cout << "------------------" << endl;
    bool allDetPos = false;
    for (int i = 0; i < det.size(); i++) {
        if (!allDetPos && det[i] > 0) allDetPos = false;
        else allDetPos = true;
        cout << "D" << i + 1 << ": " << det[i] <<  endl;
    }

    /* Ausgabe System stabil/instabil */
    if ((condition1 == true) && (!allDetPos)) {
        cout << "\n";
        cout << "Alle Koeffizienten haben das gleiche Vorzeichen!" << endl;
        cout << "Alle Hauptdeterminanten sind positiv!" << endl;
        cout << "\nDas System ist STABIL!" << "\n" << endl;
    }
    else {
        cout << "\n";
        cout << "Alle Koeffizienten haben NICHT das gleiche Vorzeichen oder mind. eine Hauptdeterminante ist NICHT >0!\n";
        cout << "\nDas System ist NICHT STABIL!" << "\n" << endl;
    }

    /* Ausgabe des Bode-Diagramms mittels GNU-Plot */  
    string csvFile = "bode_plot_data.csv";
    ofstream script("bode_plot.gp");

    script << "set terminal png size 900,600\n";
    script << "set output 'bode_plot.png'\n";
    script << "set logscale x\n";               
    script << "set grid\n";
    script << "set xlabel 'Frequency (Hz)'\n";
    script << "set ylabel 'Amplitude (dB)'\n";
    script << "set y2label 'Phase (Degrees)'\n";
    script << "set ytics nomirror\n";
    script << "set y2tics\n";
    
    /* Ränder des Plots anpassen */
    script << "set lmargin at screen 0.1\n";      // Linker Rand des Diagramms
    script << "set rmargin at screen 0.9\n";      // Rechter Rand des Diagramms
    script << "set tmargin at screen 0.8\n";       // Oberer Rand des Diagramms
    script << "set bmargin at screen 0.1\n";      // Unterer Rand des Diagramms
    
    /* Zählerstring generieren */
    ostringstream zaehlerStream;
    for (int i = gradZaehler; i >= 0; i--) {
        zaehlerStream << koeffZaehler[i] ;
        if (i > 0) zaehlerStream << "s^" << i << " + ";
    }
    string zaehlerFormel = zaehlerStream.str();
    
    /* Nennerstring generieren */
    ostringstream nennerStream;
    for (int i = gradNenner; i >= 0; i--) {
        nennerStream << koeffNenner[i];
        if (i > 0) nennerStream << "s^" << i << " + ";
    }
    string nennerFormel = nennerStream.str();

    /* Transferfunktion ausgeben */
    script << "set label 1 'G(s) = " << zaehlerFormel << " / " << nennerFormel << "'at screen 0.25, 0.9 center\n";

    /* Phasen- und Amplitudenreserve ausgeben */
    script << "set label 2 'Phasenreserve (Grad) = " << phaseMargin << "'at screen 0.6, 0.93 left\n";
    script << "set label 3 'Amplitudenreserve (dB) = " << magMargin << "'at screen 0.6, 0.89 left\n";

    script << "plot '" << csvFile << "' using 1:2 with lines title 'Amplitude (dB)', \\\n";
    script << "     '" << csvFile << "' using 1:3 axes x1y2 with lines title 'Phase (Degrees)'\n";
    script.close();

    /* Gnuplot-Skript ausführen */ 
    system("bode_plot.gp");
    cout << "Bode-Diagramm erstellt und als 'bode_plot.png' gespeichert." << endl;

    /* Bode-Plot öffnen */
    system("start bode_plot.png");
    cout << "\n\n";

 
    system("pause");
    return 0;
}