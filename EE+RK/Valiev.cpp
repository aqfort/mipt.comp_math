#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>

#include <ctime>

#define _L 10.0
#define _LMBD 0.3
#define _G 10.0
#define _U_0 1.0
#define _V_0 0.5

#define _PI 3.0

#define _N 150
#define _PERIODS 3

using namespace std;

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void LaTeX(const int &N, double *AN, double *EXP, double *EE, double *RK) {
    ofstream tex("Report.tex");
    tex << "\\documentclass[12pt, a4paper, final]{article}" << endl << endl;

    tex << "\\usepackage{cmap}" << endl;
    tex << "\\usepackage[T2A]{fontenc}" << endl;
    tex << "\\usepackage{mathtext}" << endl;
    tex << "\\usepackage[utf8]{inputenc}" << endl;
    tex << "\\usepackage[english,russian]{babel}" << endl;
    tex << "\\usepackage[left=2cm,right=2cm,top=2cm,bottom=2cm]{geometry}" << endl;
    tex << "\\usepackage{graphicx}" << endl;
    tex << "\\usepackage{color}" << endl;
    tex << "\\usepackage{float}" << endl;
    tex << "\\usepackage{wrapfig}" << endl;
    tex << "\\usepackage{pgfplots, pgfplotstable}" << endl;
    tex << "\\pgfplotsset{compat=1.9}" << endl << endl;
    tex << "\\usepackage{amsfonts}" << endl;
    tex << "\\usepackage{amsmath}" << endl;
    tex << "\\usepackage{amssymb}" << endl;
    tex << "\\usepackage{amsthm}" << endl;
    tex << "\\usepackage{mathtools}" << endl;
    tex << "\\usepackage{cancel}" << endl;
    tex << "\\usepackage{mathrsfs}" << endl;
    tex << "\\usepackage{icomma}" << endl;

    tex << "\\mathtoolsset{showonlyrefs=true}" << endl;
    tex << "\\pagestyle{empty}" << endl << endl;

    tex << "\\begin{document}" << endl << endl;

    //----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\\begin{center}" << endl;
    tex << "    \\textsc{\\LARGE Московский Физико-Технический Институт}\\\\[1,5cm]" << endl;
    tex << "    \\textsc{\\Large Вычислительная математика}\\\\[0,5cm]" << endl;
    tex << "    \\textsc{\\large Компьютерное задание \\textnumero 1}\\\\[1cm]" << endl << endl;

    tex << "    \\noindent\\rule{\\textwidth}{1pt}" << endl;
    tex << "    \\\\[0.5cm]" << endl;
    tex << "    { \\huge \\bfseries Практическое применение явной схемы Эйлера и схемы Рунге-Кутты 4-го порядка}" << endl;
    tex << "    \\\\[0.1cm]" << endl;
    tex << "    \\noindent\\rule{\\textwidth}{1pt}" << endl;
    tex << "\\end{center}" << endl << endl;

    tex << "Модель свободных колебаний физического маятника с затуханием, описываемая уравнением:" << endl;
    tex << "\\[" << endl;
    tex << "l u''(t) + 2 \\lambda u'(t) + g \\sin u(t) = 0" << endl;
    tex << "\\]" << endl;
    tex << "\\[" << endl;
    tex << "u(0) = u_0 ~~~~~~~ u'(0) = v_0" << endl;
    tex << "\\]" << endl;
    tex << "Период колебаний считается равным:" << endl;
    tex << "\\[" << endl;
    tex << "T = 2 \\pi \\sqrt{\\dfrac{l}{g}}" << endl;
    tex << "\\]" << endl;
    tex << "Параметры $ l, \\lambda, g, u_0, v_0 $ и период разбиения задаются внутри программы. Для расчетов используются явная схема Эйлера (EE) и схема Рунге-Кутты 4-го порядка (RK).\\\\" << endl << endl;

    tex << "Для малых колебаний используется формула точного решения (AN):" << endl;
    tex << "\\[" << endl;
    tex << "u(t) = Ae^{-\\gamma t}\\sin(\\omega t + \\varphi)" << endl;
    tex << "\\]" << endl;
    tex << "\\[" << endl;
    tex << "\\varphi = \\arctan \\dfrac{\\omega u_0}{v_0 + \\gamma u_0} ~~~~~ A = \\dfrac{u_0}{\\sin \\varphi} ~~~~~ \\gamma = \\dfrac{\\lambda}{l} ~~~~~ \\omega = \\sqrt{\\omega_0^2 - \\gamma^2}" << endl;
    tex << "\\]\\\\" << endl << endl;

    tex << "В работе рассматривается случай: $ \\omega_0^2 > \\gamma^2 $." << endl << endl;

    tex << "\\vfill" << endl << endl;

    tex << "\\begin{minipage}[b]{0.33\\textwidth}" << endl;
    tex << "    \\textit{Выполнил:}\\\\" << endl;
    tex << "    Р.Р. Валиев, 715 гр.\\\\\\\\" << endl;
    tex << "    \\textit{Проверил:}\\\\" << endl;
    tex << "    Н.Б. Явич" << endl;
    tex << "\\end{minipage}" << endl << endl;

    tex << "\\newpage" << endl << endl;

    //----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    double *DATA[] = {AN, EE, RK};
    const char *NAME[] = {"Exact solution (AN)",
                    "Explicit Euler method (EE)",
                    "Runge-Kutta 4th order method (RK)"};

    for(int k = 0; k < 3; k++) {
        tex << "\\begin{figure}[H]" << endl;
        tex << "    \\centering" << endl;
        tex << "    \\begin{tikzpicture}" << endl << endl;

        tex << "    \\pgfplotstableread{" << endl;
        tex << "    X Y" << endl;
        for(int i = 0; i < N; i++) {
            tex << "    " << i << " " << DATA[k][i] << endl;
        }
        tex << "}{\\DATA}" << endl << endl;

        tex << endl;

        tex << "    \\pgfplotstableread{" << endl;
        tex << "    X Y" << endl;
        for(int i = 0; i < N; i++) {
            tex << "    " << i << " " << EXP[i] << endl;
        }
        tex << "}{\\EXP}" << endl << endl;

        tex << "\\begin{axis}[" << endl;
        tex << "    width = \\textwidth," << endl;
        tex << "    xlabel = $ n $," << endl;
        tex << "    ylabel = $ U(t) \\equiv \\widetilde{U}(n) $," << endl;
        //tex << "    ymin = -5," << endl;
        //tex << "    ymax = 5," << endl;
        tex << "    xmin = 0," << endl;
        tex << "    xmax = " << _N * _PERIODS << "," << endl;
        tex << "    grid = major" << endl;
        tex << "    ]" << endl << endl;

        tex << "    \\addplot[" << endl;
        tex << "    only marks," << endl;
        tex << "    color = red," << endl;
        tex << "    mark = *" << endl;
        tex << "    ]" << endl;
        tex << "    table[" << endl;
        tex << "    x," << endl;
        tex << "    y" << endl;
        tex << "    ] {\\DATA};" << endl << endl;

        tex << "    \\addplot[" << endl;
        tex << "    only marks," << endl;
        tex << "    color = blue," << endl;
        tex << "    mark = *" << endl;
        tex << "    ]" << endl;
        tex << "    table[" << endl;
        tex << "    x," << endl;
        tex << "    y" << endl;
        tex << "    ] {\\EXP};" << endl << endl;

        tex << "    \\end{axis}" << endl;
        tex << "    \\end{tikzpicture}" << endl;
        tex << "    \\caption{" << NAME[k] << "}" << endl;
        tex << "    \\label{graph_" << k + 1 << "}" << endl;
        tex << "\\end{figure}" << endl << endl;

        tex << "\\newpage" << endl << endl;
    }

    int i = 0;
    int j = 0;
    while(i < N) {
        tex << "\\begin{table}[H]" << endl;
        tex << "    \\centering" << endl;
        tex << "    \\begin{tabular}{|c|ccc|}" << endl;
        tex << "        \\hline" << endl;
        tex << "        $ n $ & $ AN $ & $ EE $ & $ RK $ \\\\ \\hline" << endl;
        for(j += 50; (i < j) && (i < N); i++) {
            tex << "        " << i << " & " << AN[i] << " & " << EE[i] << " & " << RK[i] << " \\\\" << endl;
        }
        tex << "        \\hline" << endl;
        tex << "    \\end{tabular}" << endl;
        tex << "    \\caption{Point comparison}" << endl;
        //tex << "    \\label{TABULAR}" << endl;
        tex << "\\end{table}" << endl << endl;
    }

    tex << "\\end{document}" << endl;

    system("pdflatex Report.tex");
    system("rm -f Report.aux Report.log Report.tex");
}

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----

int main(int argc, char **argv, char **env) {
    try {
        clock_t _T_0, _T_1;
        _T_0 = clock();

        const double _OMEGA_0 = sqrt(_G / _L),
                     _T = 2 * _PI / _OMEGA_0,
                     _TAU = _T / _N;
        const int _MAX_N = _PERIODS * _N;

        if(_G / _L <= _LMBD * _LMBD / _L / _L) {
            throw invalid_argument("invalid arguments: _G, _L, _LMBD");
        }

        //exact solution (AN)

        double _U_AN[_MAX_N];
        //u = A * exp(-lmbd * t / l) * sin(omega * t + phi)
        double _OMEGA = sqrt(_G / _L - _LMBD * _LMBD / _L / _L);
        double _PHI = atan(_OMEGA * _U_0 / (_V_0 + _LMBD / _L * _U_0));
        double _A = _U_0 / sin(_PHI);
        for(int i = 0; i < _MAX_N; i++) {
            _U_AN[i] = _A * exp(-_LMBD / _L * _TAU * i) * sin(_OMEGA * _TAU * i + _PHI);
        }


        //Envelope

        double _EXP[_MAX_N];
        for(int i = 0; i < _MAX_N; i++) {
            _EXP[i] = _A * exp(-_LMBD / _L * _TAU * i);
        }

        //Explicit Euler method (EE)

        //1) u' = p
        //2) p' + 2 * (lmbd / l) * p + omega_0 * omega_0 * sin(u)
        double _U_EE[_MAX_N];
        double _P_EE[_MAX_N];
        _U_EE[0] = _U_0;
        _P_EE[0] = _V_0;
        
        
        for (int i = 0; i < _MAX_N - 1; i++) {
            _U_EE[i + 1] = _U_EE[i] + _TAU * _P_EE[i];
            _P_EE[i + 1] = _P_EE[i] - _TAU * (2 * _LMBD / _L * _P_EE[i] + _G / _L * sin(_U_EE[i]));
        }

        //Runge-Kutta 4th order method (RK)

        double _U_RK[_MAX_N];
        double _P_RK[_MAX_N];
        _U_RK[0] = _U_0;
        _P_RK[0] = _V_0;
        
        double _K_1 = 0, _K_2 = 0, _K_3 = 0, _K_4 = 0;
        
        for(int i = 0; i < _MAX_N - 1; i++) {
            _K_1 = -_G / _L * sin(_U_RK[i]) - 2 * _LMBD / _L * (_P_RK[i]);
            _K_2 = -_G / _L * sin(_U_RK[i]) - 2 * _LMBD / _L * (_P_RK[i] + _TAU * _K_1 / 2);
            _K_3 = -_G / _L * sin(_U_RK[i]) - 2 * _LMBD / _L * (_P_RK[i] + _TAU * _K_2 / 2);
            _K_4 = -_G / _L * sin(_U_RK[i]) - 2 * _LMBD / _L * (_P_RK[i] + _TAU * _K_3);
            _P_RK[i + 1] = _P_RK[i] + _TAU * (_K_1 + 2 * _K_2 + 2 * _K_3 + _K_4) / 6;

            _U_RK[i + 1] = _U_RK[i] + _TAU * _P_RK[i];
        }

        //Creating pdf-file

        LaTeX(_MAX_N, _U_AN, _EXP, _U_EE, _U_RK);

        _T_1 = clock();

        cout << endl;
        cout << "---------- ---------- ---------- ---------- ----------" << endl;
        cout << "Spent time: " << (double) (_T_1 - _T_0) / CLOCKS_PER_SEC << " sec" << endl;
        cout << "---------- ---------- ---------- ---------- ----------" << endl;
        cout << endl;
        cout << "---------- ---------- ---------- ---------- ----------" << endl;
        cout << "Done, check \"Report.pdf\" file" << endl;
        cout << "---------- ---------- ---------- ---------- ----------" << endl;
        cout << endl;

    } catch(const exception &ex) {
        cout << "U have got an exception::: " << ex.what() << endl;
        return 0;
    }

    return 0;
}

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
