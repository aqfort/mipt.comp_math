#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>

#include <ctime>

#define T_SZ long
#define T_EL double

#define MY_FILE "input.txt"
#define MY_FILE_OUT "output.txt"
#define EPS 1e-6

#define OMEGA 1.5

using namespace std;

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----

class Matrix {
public:
    Matrix() {
        ROWS = 0;
        COLS = 0;
        BuildMatrix();
    }

    Matrix(const T_SZ &new_SIZE) {
        if(new_SIZE < 0) {
            throw out_of_range("arguments are negative");
        }
        ROWS = new_SIZE;
        COLS = new_SIZE;
        BuildMatrix();
    }

    Matrix(const T_SZ &new_ROWS, const T_SZ &new_COLS) {
        AreSizeArgsNegative(new_ROWS, new_COLS);
        ROWS = new_ROWS;
        COLS = new_COLS;
        BuildMatrix();
    }

    T_SZ GetNumRows() const {
        return ROWS;
    }

    T_SZ GetNumColumns() const {
        return COLS;
    }

    void Reset(const T_SZ &new_ROWS, const T_SZ &new_COLS) {
        AreSizeArgsNegative(new_ROWS, new_COLS);
        ROWS = new_ROWS;
        COLS = new_COLS;
        BuildMatrix();
    }

    T_EL At(const T_SZ &row, const T_SZ &col) const {
        AreSizeArgsNegative(row, col);
        AreSizeArgsOutOfRange(row, col);
        return DATA.at(row).at(col);
    }

    T_EL &At(const T_SZ &row, const T_SZ &col) {
        AreSizeArgsNegative(row, col);
        AreSizeArgsOutOfRange(row, col);
        return DATA.at(row).at(col);
    };

    friend bool operator==(const Matrix &lhs, const Matrix &rhs);
    friend Matrix operator+(const Matrix &lhs, const Matrix &rhs);
    friend Matrix operator-(const Matrix &lhs, const Matrix &rhs);

    Matrix &operator=(const Matrix &matrix) {
        if(this == &matrix) {
            return *this;
        }
        ROWS = matrix.GetNumRows();
        COLS = matrix.GetNumColumns();
        for(T_SZ row = 0; row < ROWS; row++) {
            for(T_SZ element = 0; element < COLS; element++) {
                this->At(row, element) = matrix.At(row, element);
            }
        }
        return *this;
    }

    T_EL CubeNorm() const {
        T_EL result = this->At(0,0), temp = 0;
        for(T_SZ row = 0; row < ROWS; row++) {
            for(T_SZ col = 0; col < COLS; col++) {
                temp += abs(this->At(row, col));
            }
            if(temp > result) {
                result = temp;
            }
            temp = 0;
        }
        return result;
    }

private:
    T_SZ ROWS;   //number of rows of matrix
    T_SZ COLS;   //number of columns of matrix
    vector<vector<T_EL>> DATA;   //data of matrix

    void AreSizeArgsNegative(const T_SZ &row, const T_SZ &col) const {
        if(row < 0 || col < 0) {
            throw out_of_range("arguments are negative");
        }
    }

    void AreSizeArgsOutOfRange(const T_SZ &row, const T_SZ &col) const {
        if(row > ROWS || col > COLS) {
            throw out_of_range("arguments are out of range");
        }
    }

    void BuildMatrix() {
        DATA.clear();
        DATA.resize(ROWS);
        for(int row = 0; row < ROWS; row++) {
            DATA[row].resize(COLS);
        }
    }

    bool AreMatrixSizesEqual(const Matrix &matrix) const {
        return (ROWS == matrix.GetNumRows() && 
                COLS == matrix.GetNumColumns());
    }
};

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----

istream &operator>>(istream &in, Matrix &matrix) {
    T_SZ new_ROWS(0), new_COLS(0);
    string line;
    getline(in, line);
    stringstream line_stream(line);
    line_stream >> new_ROWS >> ws; // eat up any leading white spaces
    if(line_stream.peek() == EOF) {
        new_COLS = new_ROWS;
    } else {
        line_stream >> new_COLS;
    }
    //in >> new_ROWS >> new_COLS;
    matrix.Reset(new_ROWS, new_COLS);
    for(T_SZ row = 0; row < new_ROWS; row++) {
        for(T_SZ element = 0; element < new_COLS; element++) {
            in >> matrix.At(row, element);
        }
    }
    return in;
}

ostream &operator<<(ostream &out, const Matrix &matrix) {
    out << "ROWS: " << matrix.GetNumRows() << '\n' << "COLS: " << matrix.GetNumColumns() << endl;

    long flag = out.precision();

    out << setprecision(15);

    for(T_SZ row = 0; row < matrix.GetNumRows(); row++) {
        for(T_SZ element = 0; element < matrix.GetNumColumns(); element++) {
            out << setw(20) << matrix.At(row, element) << ' ';
        }
        if((row + 1) < matrix.GetNumRows()) {
            out << endl;
        }
    }

    out.precision(flag);

    return out;
}

bool operator==(const Matrix &lhs, const Matrix &rhs) {
    if(!lhs.AreMatrixSizesEqual(rhs)) {
        return false;
    }
    for(T_SZ row = 0; row < lhs.GetNumRows(); row++) {
        for(T_SZ element = 0; element < lhs.GetNumColumns(); element++) {
            if(lhs.At(row, element) != rhs.At(row, element)) {
                return false;
            }
        }
    }
    return true;
}

Matrix operator+(const Matrix &lhs, const Matrix &rhs) {
    if(!lhs.AreMatrixSizesEqual(rhs)) {
        throw invalid_argument("different sizes");
    }
    Matrix result(lhs.GetNumRows(), lhs.GetNumColumns());
    for(T_SZ row = 0; row < lhs.GetNumRows(); row++) {
        for(T_SZ element = 0; element < lhs.GetNumColumns(); element++) {
            result.At(row, element) = lhs.At(row, element) + rhs.At(row, element);
        }
    }
    return result;
}

Matrix operator-(const Matrix &lhs, const Matrix &rhs) {
    if(!lhs.AreMatrixSizesEqual(rhs)) {
        throw invalid_argument("different sizes");
    }
    Matrix result(lhs.GetNumRows(), lhs.GetNumColumns());
    for(T_SZ row = 0; row < lhs.GetNumRows(); row++) {
        for(T_SZ element = 0; element < lhs.GetNumColumns(); element++) {
            result.At(row, element) = lhs.At(row, element) - rhs.At(row, element);
        }
    }
    return result;
}

Matrix operator*(const Matrix &lhs, const Matrix &rhs) {
    if(lhs.GetNumColumns() != rhs.GetNumRows()) {
        throw invalid_argument("not matched for multiplication");
    }
    T_SZ N = lhs.GetNumColumns();
    Matrix result(lhs.GetNumRows(), rhs.GetNumColumns());
    for(T_SZ row = 0; row < lhs.GetNumRows(); row++) {
        for(T_SZ element = 0; element < rhs.GetNumColumns(); element++) {
            for(T_SZ item = 0; item < N; item++) {
                result.At(row, element) += lhs.At(row, item) * rhs.At(item, element);
            }
        }
    }
    return result;
}

Matrix operator*(const T_EL &lhs, const Matrix &rhs) {
    Matrix result(rhs.GetNumRows(), rhs.GetNumColumns());
    for(T_SZ row = 0; row < rhs.GetNumRows(); row++) {
        for(T_SZ element = 0; element < rhs.GetNumColumns(); element++) {
            result.At(row, element) = lhs * rhs.At(row, element);
        }
    }
    return result;
}

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----

Matrix ZeydelStep(const Matrix &A, const Matrix &b, const Matrix &x) {
    T_SZ N = A.GetNumRows(), M = A.GetNumColumns();
    Matrix result(N, 1);
    T_EL temp = 0;
    for(T_SZ row = 0; row < N; row++) {
        for(T_SZ col = 0; col < M; col++) {
            if(row != col) {
                temp += A.At(row, col) * x.At(col, 0);
            }
        }
        if(A.At(row, row) == 0) {
            throw invalid_argument("diagonal element equals zero");
        }
        result.At(row, 0) = (b.At(row, 0) - temp) / A.At(row, row);
        temp = 0;
    }
    return result;
}

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----

int main() {
    Matrix A;
    ifstream input(MY_FILE);
    ofstream out(MY_FILE_OUT);
    input >> A;

    clock_t t0, t1, t2, t3;

    T_SZ N = A.GetNumRows(), M = A.GetNumColumns();

    Matrix b(N, 1), x(N, 1), y(N, 1), NEV(N, 1);
    for(T_SZ element = 0; element < N; element++) {
        b.At(element, 0) = 1;
        x.At(element, 0) = 0;
        y.At(element, 0) = 0;
    }

    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";
    out << "MATRIX A:\n" << A << endl;
    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";
    out << "VECTOR b:\n" << b << endl;
    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";
    out << "STARTER VECTOR x:\n" << x << endl;
    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";

    out << "\n\nGauss–Seidel\n\n\n";
    t0 = clock();
    NEV = b - A * x;
    out << "diff = " << NEV.CubeNorm() << endl;
    while(NEV.CubeNorm() > EPS) {
        x = ZeydelStep(A, b, x);
        NEV = b - A * x;
        out << "diff = " << NEV.CubeNorm() << endl;
    }
    t1 = clock();

    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";
    out << "VECTOR x:\n" << x << endl;
    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";



    for(T_SZ element = 0; element < N; element++) {
        x.At(element, 0) = 0;
        y.At(element, 0) = 0;
    }

    out << "\n\nSuccessive over-relaxation\n\n\n";
    t2 = clock();
    NEV = b - A * x;
    out << "diff = " << NEV.CubeNorm() << endl;
    while(NEV.CubeNorm() > EPS) {
        y = ZeydelStep(A, b, x);
        x = (1 - OMEGA) * x + OMEGA * y;
        NEV = b - A * x;
        out << "diff = " << NEV.CubeNorm() << endl;
    }
    t3 = clock();

    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";
    out << "VECTOR x:\n" << x << endl;
    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";

    if((double)(t0 - t1) < (double)(t3 - t2)) {
        out << "\n\nGauss–Seidel is faster and difference = " << ((double)(t3 - t2) - (double)(t1 - t0)) / CLOCKS_PER_SEC << " sec!!!\n\n\n";
    } else {
        out << "\n\nSuccessive over-relaxation is faster and difference = " << ((double)(t1 - t0) - (double)(t3 - t2)) / CLOCKS_PER_SEC << " sec!!!\n\n\n";
    }

    out << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";

    cout << "\n----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n";
    cout << "\nDONE! CHECK \"" << MY_FILE_OUT << "\"!\n\n";
    cout << "----- ----- ----- ----- ----- ----- ----- ----- ----- -----\n\n";

    return 0;
}

//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
//----- ----- ----- ----- ----- ----- ----- ----- ----- -----