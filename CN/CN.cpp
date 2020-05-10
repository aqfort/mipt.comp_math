#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <fstream>

#define _USE_MATH_DEFINES
#define H 50
#define T 20

using namespace std;

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void PRINT_TO_FILE(const vector<vector<double>> & DATA,
                   const string & NAME);
void LaTeX();

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

double FUNCTION_EXACT(const double & x, const double & t);
double KAPPA(const double & x);
double PHI(const double & x);

void GET_DATA_EXACT(vector<vector<double>> & DATA_EXACT);

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void GET_F_EXPLICIT(vector<vector<double>> & fooo_EXPLICIT);

void EXPLICIT(vector<vector<double>> & DATA_EXPLICIT,
              vector<vector<double>> & fooo_EXPLICIT);

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void GET_F_CN(vector<vector<double>> & fooo_CN);

void PROGONKA(const int & N, // H + 1
              vector<double> A, // A_CN
              vector<double> B, // B_CN
              vector<double> C, // C_CN
              vector<double> D, // D_CN
              vector<double> & X);

void CN(vector<vector<double>> & DATA_CN,
        vector<vector<double>> & fooo_CN);

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

int main(int argc, char ** argv, char ** env) {
    cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
    cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░STARTING░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
    cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";

    vector<vector<double>> DATA_EXACT(H + 1, vector<double>(T + 1));

    GET_DATA_EXACT(DATA_EXACT);

    PRINT_TO_FILE(DATA_EXACT, "DATA_EXACT.dat");

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    vector<vector<double>> DATA_EXPLICIT(H + 1, vector<double>(T + 1));
    vector<vector<double>> fooo_EXPLICIT(H + 1, vector<double>(T + 1));

    GET_F_EXPLICIT(fooo_EXPLICIT);
    PRINT_TO_FILE(fooo_EXPLICIT, "fooo_EXPLICIT.dat");

    EXPLICIT(DATA_EXPLICIT, fooo_EXPLICIT);

    PRINT_TO_FILE(DATA_EXPLICIT, "DATA_EXPLICIT.dat");

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    vector<vector<double>> DATA_CN(H + 1, vector<double>(T + 1));
    vector<vector<double>> fooo_CN(H + 1, vector<double>(T + 1));

    GET_F_CN(fooo_CN);
    PRINT_TO_FILE(fooo_CN, "fooo_CN.dat");

    CN(DATA_CN, fooo_CN);

    PRINT_TO_FILE(DATA_CN, "DATA_CN.dat");



    LaTeX();



    cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
    cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░DONE░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";
    cout << "░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░\n";

    return 0;
}

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void PRINT_TO_FILE(const vector<vector<double>> & DATA,
                   const string & NAME) {
    ofstream out(NAME);
    long flag = out.precision();
    out << fixed << setprecision(7);

    for(int t = 0; t <= T; t++) {
        for(int x = 0; x <= H; x++) {
            out << x << ' ' << t  << ' ' << DATA.at(x).at(t) << endl;
        }
        out << endl;
    }

    out.precision(flag);
    out << resetiosflags(ios::fixed);

    out.close();
}

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void LaTeX() {
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

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\\begin{center}" << endl << endl;
    tex << "    \\textsc{\\LARGE Московский Физико-Технический Институт}\\\\[1,5cm]" << endl;
    tex << "    \\textsc{\\Large Вычислительная математика}\\\\[0,5cm]" << endl;
    tex << "    \\textsc{\\large Компьютерное задание \\textnumero 3}\\\\[1cm]" << endl << endl;

    tex << "    \\noindent\\rule{\\textwidth}{1pt}" << endl;
    tex << "    \\\\[0.5cm]" << endl;
    tex << "    { \\huge \\bfseries Уравнение теплопроводности}" << endl;
    tex << "    \\\\[0.1cm]" << endl;
    tex << "    \\noindent\\rule{\\textwidth}{1pt}" << endl;
    tex << "\\end{center}" << endl << endl;

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\\vspace{1cm}\n" << endl;

    tex << "В данной работе предлагается ознакомиться с реализацией явной схемы (ЯС) и схемы Кранка-Никольсона (КН) для нахождения приближенного решения уравнения теплопроводности:" << endl;
    tex << "\\begin{equation*}" << endl;
    tex << "\\begin{cases}" << endl;
    tex << "\\dfrac{\\partial u}{\\partial t} - \\dfrac{\\partial}{\\partial x} \\left[ \\varkappa(x) \\dfrac{\\partial u}{\\partial x} \\right] = f(x, t)\\\\" << endl;
    tex << "u(x, 0) = \\varphi(x)\\\\" << endl;
    tex << "u(0, t) = u(1, t) = 0" << endl;
    tex << "\\end{cases}" << endl;
    tex << "~~~~~" << endl;
    tex << "\\begin{cases}" << endl;
    tex << "0 < x < 1\\\\" << endl;
    tex << "0 < t < 1\\\\" << endl;
    tex << "\\varkappa(x) = \\left| x - 0.5 \\right|\\\\" << endl;
    tex << "\\varphi(x) = \\sin(\\pi x)" << endl;
    tex << "\\end{cases}" << endl;
    tex << "\\end{equation*}" << endl;
    tex << "" << endl;
    tex << "Для решения системы на верхнем слое в схеме Кранка-Никольсона использовался метод прогонки.\\\\" << endl;
    tex << "" << endl;
    tex << "Работа алгоритма приведена для нахождения решения $ u(x, t) = \\sin(\\pi (t + 0.5)) \\sin(\\pi x) $, сетка выбрана с шагами: $ h = 1 / " << H << " $ по пространству, $ \\tau = 1 / " << T << " $ по времени.\\\\" << endl;
    tex << "" << endl;
    tex << "Программа выводит графики: искомой функции $ u(x, t) = \\sin(\\pi (t + 0.5)) \\sin(\\pi x) $, результат работы ЯС и схемы КН.\\\\" << endl;

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\\vfill" << endl << endl;

    tex << "\\begin{minipage}[b]{0.33\\textwidth}" << endl;
    tex << "    \\textit{Выполнил:}\\\\" << endl;
    tex << "    Р.Р. Валиев, 715 гр.\\\\\\\\" << endl;
    tex << "    \\textit{Проверил:}\\\\" << endl;
    tex << "    Н.Б. Явич" << endl;
    tex << "\\end{minipage}" << endl << endl;

    tex << "\\newpage" << endl << endl;

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    // fooo_EXPLICIT and fooo_CN

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\n";

    tex << "\\subsection*{Сравнение функции справа $ f(x, t) $ для ЯС и схемы КН.}\n";

    tex << "\n";

    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[view={120}{10}, grid=both, xlabel={$h$\\text{-отсчеты}}, ylabel={$\\tau$\\text{-отсчеты}}, title={\\text{Функция справа (ЯС)}}]\n";
    tex << "\\addplot3[surf, mesh/rows=" << T + 1 << "] file {fooo_EXPLICIT.dat};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";

    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[view={120}{10}, grid=both, xlabel={$h$\\text{-отсчеты}}, ylabel={$\\tau$\\text{-отсчеты}}, title={\\text{Функция справа (КН)}}]\n";
    tex << "\\addplot3[surf, mesh/rows=" << T + 1 << "] file {fooo_CN.dat};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";

    tex << "\n";

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    // real points

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\n";

    tex << "\\subsection*{Искомая функция для сравнения с результатами работы.}\n";

    tex << "\n";

    tex << "\\begin{center}\n";
    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[\n";
    tex << "view={120}{10},\n";
    tex << "xmin=0,\n";
    tex << "xmax=1,\n";
    tex << "ymin=0,\n";
    tex << "ymax=1,\n";
    tex << "grid=both, xlabel={x}, ylabel={t}, zlabel={$ u(x, t) = \\sin(\\pi (t + 0.5)) \\sin(\\pi x) $}, title={\\text{Искомая функция}}]\n";
    tex << "\\addplot3[surf,samples=50,domain=0:1]\n";
    tex << "{sin(180 * (y + 0.5)) * sin(180 * x)};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";
    tex << "\\end{center}\n";

    tex << "\n";

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    // DATA_EXACT and DATA_EXPLICIT

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\n";

    tex << "\\subsection*{Сравнение результата работы явной схемы с точными значениями искомой функции в узлах сетки.}\n";

    tex << "\n";

    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[view={120}{10}, grid=both, xlabel={$h$\\text{-отсчеты}}, ylabel={$\\tau$\\text{-отсчеты}}, title={\\text{Точные значения}}]\n";
    tex << "\\addplot3[surf, mesh/rows=" << T + 1 << "] file {DATA_EXACT.dat};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";

    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[view={120}{10}, grid=both, xlabel={$h$\\text{-отсчеты}}, ylabel={$\\tau$\\text{-отсчеты}}, title={\\text{Результат работы ЯС}}]\n";
    tex << "\\addplot3[surf, mesh/rows=" << T + 1 << "] file {DATA_EXPLICIT.dat};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";

    tex << "\n";

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    // DATA_EXACT and DATA_CN

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\n";

    tex << "\\subsection*{Сравнение результата работы схемы Кранка-Никольсона с точными значениями искомой функции в узлах сетки.}\n";

    tex << "\n";

    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[view={120}{10}, grid=both, xlabel={$h$\\text{-отсчеты}}, ylabel={$\\tau$\\text{-отсчеты}}, title={\\text{Точные значения}}]\n";
    tex << "\\addplot3[surf, mesh/rows=" << T + 1 << "] file {DATA_EXACT.dat};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";

    tex << "\\begin{tikzpicture}\n";
    tex << "\\begin{axis}[view={120}{10}, grid=both, xlabel={$h$\\text{-отсчеты}}, ylabel={$\\tau$\\text{-отсчеты}}, title={\\text{Результат работы схемы КН}}]\n";
    tex << "\\addplot3[surf, mesh/rows=" << T + 1 << "] file {DATA_CN.dat};\n";
    tex << "\\end{axis}\n";
    tex << "\\end{tikzpicture}\n";

    tex << "\n";

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex << "\\end{document}" << endl;

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    tex.close();

    system("pdflatex Report.tex");
    system("rm -f Report.aux Report.log Report.tex");
    system("rm -f DATA_CN.dat DATA_EXACT.dat DATA_EXPLICIT.dat fooo_CN.dat fooo_EXPLICIT.dat");
}

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

double FUNCTION_EXACT(const double & x, const double & t) {
    return sin(M_PI * t / ((double) T) + M_PI * 0.5) * sin(M_PI * x / ((double) H));
}

double KAPPA(const double & x) {
    return abs(x / ((double) H) - 0.5);
}

double PHI(const double & x) {
    return sin(M_PI * x / ((double) H));
}

void GET_DATA_EXACT(vector<vector<double>> & DATA_EXACT) {
    for(int x = 0; x <= H; x++) {
        for(int t = 0; t <= T; t++) {
            DATA_EXACT.at(x).at(t) = FUNCTION_EXACT(x, t);
        }
    }
}

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void GET_F_EXPLICIT(vector<vector<double>> & fooo_EXPLICIT) {
    for(int t = 0; t <= T; t++) {
        for(int x = 0; x <= H; x++) {
            fooo_EXPLICIT.at(x).at(t) = T * (FUNCTION_EXACT(x, t + 1) - FUNCTION_EXACT(x, t)) -
                                        H * (KAPPA(x + 0.5) * H * (FUNCTION_EXACT(x + 1, t) - FUNCTION_EXACT(x, t)) -
                                             KAPPA(x - 0.5) * H * (FUNCTION_EXACT(x, t) - FUNCTION_EXACT(x - 1, t)));
        }
    }
}

void EXPLICIT(vector<vector<double>> & DATA_EXPLICIT,
              vector<vector<double>> & fooo_EXPLICIT) {
    // initial conditions
    for(int x = 0; x <= H; x++) {
        DATA_EXPLICIT.at(x).at(0) = PHI(x);
    }

    // boundary conditions
    for(int t = 0; t <= T; t++) {
        DATA_EXPLICIT.at(0).at(t) = 0;
        DATA_EXPLICIT.at(H).at(t) = 0;
    }

    // internal points
    for(int t = 0; t < T; t++) {
        for(int x = 1; x < H; x++) {
            DATA_EXPLICIT.at(x).at(t + 1) = DATA_EXPLICIT.at(x).at(t) + 1.0 / T *
                                            (fooo_EXPLICIT.at(x).at(t) + H *
                                                (KAPPA(x + 0.5) * H * (FUNCTION_EXACT(x + 1, t) - FUNCTION_EXACT(x, t)) -
                                                 KAPPA(x - 0.5) * H * (FUNCTION_EXACT(x, t) - FUNCTION_EXACT(x - 1, t))));
        }
    }
}

// ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

void GET_F_CN(vector<vector<double>> & fooo_CN) {
    for(int t = 0; t <= T; t++) {
        for(int x = 0; x <= H; x++) {
            fooo_CN.at(x).at(t) = T * (FUNCTION_EXACT(x, t + 1) - FUNCTION_EXACT(x, t)) -
                                  KAPPA(x) * H * H * 0.5 *
                                  (FUNCTION_EXACT(x + 1, t + 1) - 2 * FUNCTION_EXACT(x, t + 1) + FUNCTION_EXACT(x - 1, t + 1) +
                                   FUNCTION_EXACT(x + 1, t) - 2 * FUNCTION_EXACT(x, t) + FUNCTION_EXACT(x - 1, t));
        }
    }
}

    //https://ru.wikibooks.org/wiki/Реализации_алгоритмов/Метод_прогонки

    /**
     * N - число уравнений (строк матрицы)
     * C - диагональ, лежащая над главной (нумеруется: [0; N - 2])
     * B - главная диагональ матрицы A (нумеруется: [0; N - 1])
     * A - диагональ, лежащая под главной (нумеруется: [1; N - 1])
     * D - правая часть (столбец)
     * X - решение, массив X будет содержать ответ
     */
void PROGONKA(const int & N, // H + 1
              vector<double> A, // A_CN
              vector<double> B, // B_CN
              vector<double> C, // C_CN
              vector<double> D, // D_CN([0, H], t)
              vector<double> & X) {
    double TEMP;
    for (int i = 1; i < N; i++) {
        TEMP = A.at(i) / B.at(i - 1);
        B.at(i) = B.at(i) - TEMP * C.at(i - 1);
        D.at(i) = D.at(i) - TEMP * D.at(i - 1);
    }

    X.at(N - 1) = D.at(N - 1) / B.at(N - 1);

    for (int i = N - 2; i >= 0; i--) {
        X.at(i) = (D.at(i) - C.at(i) * X.at(i + 1)) / B.at(i);
    }
}

void CN(vector<vector<double>> & DATA_CN,
        vector<vector<double>> & fooo_CN) {
    // initial conditions
    for(int x = 0; x <= H; x++) {
        DATA_CN.at(x).at(0) = PHI(x);
    }

    // boundary conditions
    for(int t = 0; t <= T; t++) {
        DATA_CN.at(0).at(t) = 0;
        DATA_CN.at(H).at(t) = 0;
    }

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    vector<double> A_CN(H + 1);
    vector<double> B_CN(H + 1);
    vector<double> C_CN(H + 1);
    vector<double> D_CN(H + 1);

    for(int x = 0; x <= H; x++) {
        A_CN.at(x) = -KAPPA(x) * H * H * 0.5;
        B_CN.at(x) = T + KAPPA(x) * H * H;
        C_CN.at(x) = -KAPPA(x) * H * H * 0.5;
    }

    // ----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    for(int t = 1; t <= T; t++) {
        for(int x = 1; x < H; x++) {
            D_CN.at(x) = fooo_CN.at(x).at(t - 1) +
                         KAPPA(x) * H * H * 0.5 * DATA_CN.at(x - 1).at(t - 1) +
                         (T - KAPPA(x) * H * H) * DATA_CN.at(x).at(t - 1) +
                         KAPPA(x) * H * H * 0.5 * DATA_CN.at(x + 1).at(t - 1);
        }

        vector<double> X(H + 1);

        PROGONKA(H + 1, A_CN, B_CN, C_CN, D_CN, X);

        for(int x = 0; x <= H; x++) {
            DATA_CN.at(x).at(t) = X.at(x);
        }
    }

    // initial conditions
    for(int x = 0; x <= H; x++) {
        DATA_CN.at(x).at(0) = PHI(x);
    }

    // boundary conditions
    for(int t = 0; t <= T; t++) {
        DATA_CN.at(0).at(t) = 0;
        DATA_CN.at(H).at(t) = 0;
    }
}
