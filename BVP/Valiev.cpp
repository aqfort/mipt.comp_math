#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>

#define N 100

using namespace std;

void LaTeX(double *Y_EX, double *Y_NO, const double &NORM) {
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

tex << "\\begin{center}" << endl << endl;
tex << "    \\textsc{\\LARGE Московский Физико-Технический Институт}\\\\[1,5cm]" << endl;
tex << "    \\textsc{\\Large Вычислительная математика}\\\\[0,5cm]" << endl;
tex << "    \\textsc{\\large Компьютерное задание \\textnumero 2}\\\\[1cm]" << endl << endl;

tex << "    \\noindent\\rule{\\textwidth}{1pt}" << endl;
tex << "    \\\\[0.5cm]" << endl;
tex << "    { \\huge \\bfseries Решение краевой задачи методом конечных разностей}" << endl;
tex << "    \\\\[0.1cm]" << endl;
tex << "    \\noindent\\rule{\\textwidth}{1pt}" << endl;
tex << "\\end{center}" << endl << endl;

tex << "В данной работе предлагается ознакомиться с решением краевой задачи:" << endl;
tex << "\\begin{equation*}" << endl;
tex << "\\begin{cases}" << endl;
tex << "\\dfrac{\\text{d}}{\\text{d}x} \\left( k(x) \\dfrac{\\text{d}y}{\\text{d}x} \\right) - p(x) = f(x)\\\\" << endl;
tex << "\\left.\\dfrac{\\text{d}y}{\\text{d}x}\\right|_{x = 0} = a\\\\" << endl;
tex << "\\left.\\dfrac{\\text{d}y}{\\text{d}x}\\right|_{x = 1} = b" << endl;
tex << "\\end{cases}" << endl;
tex << "~~~~~ x \\in [0; 1]" << endl;
tex << "\\end{equation*}" << endl << endl;

tex << "В качестве примера использованы функции:" << endl;
tex << "\\begin{equation*}" << endl;
tex << "\\begin{cases}" << endl;
tex << "k(x) = -x\\\\" << endl;
tex << "p(x) = -x^2 - 1" << endl;
tex << "\\end{cases}" << endl;
tex << "\\end{equation*}" << endl << endl;

tex << "Также в работе используется равномерная сетка из~$ N = " << N << " $~интервалов. Окончательный результат вычисляется методом прогонки.\\\\" << endl << endl;

tex << "Функция~$ f(x) $ ищется в узлах сетки для известной нам функции~$ y_0(x) = \\sqrt{x + 1} $. Так же определяются коэффициенты $ a $~и~$ b $ начальных условий. Таким образом, в конце работы результат $ y(x) $ будет сравниваться с истинной функцией~$ y_0(x) $.\\\\" << endl << endl;

tex << "В качестве результата приводятся графики функций: истинная~$ y_0(x) $ и вычисленная~$ y(x) $, а также таблица значений функций в узлах сетки для поточечного сравнения.\\\\" << endl << endl;

tex << "Спойлер: кубическая норма разницы истинного и вычисленного функций равна:\\\\" << endl;

long flag = tex.precision();

tex << "$$ \\left\\| y_0(x_n) - y(x_n) \\right\\|_3 = \\max_{0 \\le n \\le 100,~n \\in \\mathbb{N}} \\left| y_0(x_n) - y(x_n) \\right| = " << fixed << setprecision(7) << NORM << "$$";

tex.precision(flag);

tex << resetiosflags(ios::fixed);

    tex << "\\vfill" << endl << endl;

    tex << "\\begin{minipage}[b]{0.33\\textwidth}" << endl;
    tex << "    \\textit{Выполнил:}\\\\" << endl;
    tex << "    Р.Р. Валиев, 715 гр.\\\\\\\\" << endl;
    tex << "    \\textit{Проверил:}\\\\" << endl;
    tex << "    Н.Б. Явич" << endl;
    tex << "\\end{minipage}" << endl << endl;

    tex << "\\newpage" << endl << endl;

    //----- ----- ----- ----- ----- ----- ----- ----- ----- -----

    double *DATA[] = {Y_EX, Y_NO};
    const char *NAME[] = {"Истинное решение",
                          "Результат вычислений"};

    for(int k = 0; k < 2; k++) {
        tex << "\\begin{figure}[H]" << endl;
        tex << "    \\centering" << endl;
        tex << "    \\begin{tikzpicture}" << endl << endl;

        tex << "    \\pgfplotstableread{" << endl;
        tex << "    X Y" << endl;
        for(int i = 0; i <= N; i++) {
            tex << "    " << ((double) i) / ((double) N) << " " << DATA[k][i] << endl;
        }
        tex << "}{\\DATA}" << endl << endl;

        tex << endl;

        // tex << "    \\pgfplotstableread{" << endl;
        // tex << "    X Y" << endl;
        // for(int i = 0; i < N; i++) {
        //     tex << "    " << i << " " << EXP[i] << endl;
        // }
        // tex << "}{\\EXP}" << endl << endl;

        tex << "\\begin{axis}[" << endl;
        tex << "    width = \\textwidth," << endl;
        tex << "    xlabel = $ x $," << endl;
        tex << "    ylabel = $ y(x) $," << endl;
        tex << "    ymin = 1," << endl;
        tex << "    ymax = 1.5," << endl;
        tex << "    xmin = 0," << endl;
        tex << "    xmax = 1," << endl;
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

        // tex << "    \\addplot[" << endl;
        // tex << "    only marks," << endl;
        // tex << "    color = blue," << endl;
        // tex << "    mark = *" << endl;
        // tex << "    ]" << endl;
        // tex << "    table[" << endl;
        // tex << "    x," << endl;
        // tex << "    y" << endl;
        // tex << "    ] {\\EXP};" << endl << endl;

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
        tex << "    \\begin{tabular}{|c|cc|}" << endl;
        tex << "        \\hline" << endl;
        tex << "        x & Истинное значение & Результат вычислений \\\\ \\hline" << endl;

        long FFF = tex.precision();

        tex << fixed;

        for(j += 50; (i < j) && (i < N); i++) {
            tex << "        " << setprecision(2) << ((double) i) / ((double) N) << setprecision(7) << " & " << Y_EX[i] << " & " << Y_NO[i] << " \\\\" << endl;
        }
        
        tex.precision(FFF);

		tex << resetiosflags(ios::fixed);
        
        tex << "        \\hline" << endl;
        tex << "    \\end{tabular}" << endl;
        tex << "    \\caption{Поточечное сравнение}" << endl;
        //tex << "    \\label{TABULAR}" << endl;
        tex << "\\end{table}" << endl << endl;
    }

    tex << "\\end{document}" << endl;

    system("pdflatex Report.tex");
    system("rm -f Report.aux Report.log Report.tex");
}

double Y_N(const double &n) { //calculating y(x_n)
    return sqrt(1 + ( n / ((double) N)));
}

double K_N(const double &n) { //calcilating k(x_n)
    return (-(n / ((double) N)));
}

double P_N(const double &n) { //calculating p(x_n)
    return (-(n * n / ((double) (N * N)) + 1));
}

void Init(double *F, double &a, double &b) { //initialize f(x_n) and coef-s a & b from exact solution
    for(int i = 0; i <= N; i++) {
        F[i] = -1 / (2 * Y_N(i)) - K_N(i) / (4 * Y_N(i) * Y_N(i) * Y_N(i)) - P_N(i) * Y_N(i);
    }
    a = 0.5;
    b = 0.5 / sqrt(2);
}

int main(int argc, char **argv, char **env) {
    double F[N + 1];
    double a, b;

    double Y_EX[N + 1];

    for(int i = 0; i <= N; i++) { //exact solution
        Y_EX[i] = Y_N(i);
    }

    Init(F, a, b); //initialize f(x_n)

    double A[N + 1];
    double B[N + 1];
    double C[N + 1];

    for(int i = 0; i <= N; i++) { //calculating tridiogonal coef-s a_n, b_n, c_n (in our case d_n == f_n)
        A[i] = N * N * K_N(((double) i) - 0.5);
        B[i] = N * N * (K_N(((double) i) + 0.5) + K_N(((double) i) - 0.5)) + P_N(i);
        C[i] = N * N * K_N(((double) i) + 0.5);
    }

    double Y_NO[N + 1];

    double P[N + 1];
    double Q[N + 1];

    P[1] = 1; //init cond-s
    Q[1] = -a / ((double) N); //init cond-s

    Y_NO[0] = F[0] + a; //init cond-s

    for(int i = 2; i <= N; i++) { //conversion coef-s
        P[i] = C[i - 1] / (B[i - 1] - A[i - 1] * P[i - 1]);
        Q[i] = (A[i - 1] * Q[i - 1] - F[i - 1]) / (B[i - 1] - A[i - 1] * P[i - 1]);
    }

    for(int i = 1; i <= N; i++) { //calculating y(x_n)
        Y_NO[i] = (Y_NO[i - 1] - Q[i]) / P[i];
    }

    double NORM = 0;
    for(int i = 0; i <= N; i++) { //calculating max deviation
        if(abs(Y_EX[i] - Y_NO[i]) > NORM)
            NORM = abs(Y_EX[i] - Y_NO[i]);
    }

    LaTeX(Y_EX, Y_NO, NORM); //compile report

    return 0;
}


