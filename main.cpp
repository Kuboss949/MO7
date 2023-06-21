#include <iostream>
#include <fstream>

using namespace std;

void wypiszWektor(double* v, int size){
    for (int i = 0; i < size; ++i) {
        cout << v[i] << " ";
    }
}

double * residuum(double* v, double** A, double* b, int n){
    auto *r = new double[n];
    for (int i = 0; i < n; ++i) {
        r[i] = b[i];
        for (int j = 0; j < n; ++j) {
            r[i] -= A[i][j] * v[j];
        }
    }
    return r;
}

double* errorEstimator(const double* prev, const double* next, int n){
    auto *e = new double[n];
    for (int i = 0; i < n; ++i) {
        e[i] = fabs(next[i] - prev[i]);
    }
    return e;
}

bool checkEstimator(double* e){
    if(e[0]<0.00001 && e[1]<0.00001 && e[2]<0.00001 && e[3]<0.00001)
        return true;
    return false;
}

bool checkResiduum(double* e){
    if(e[0]<0.00001 && e[1]<0.00001 && e[2]<0.00001 && e[3]<0.00001)
        return true;
    return false;
}

void jacobi(int n, double** A, double *b, double *xo, int iterations){
    auto* u=new double[n];
    auto* x=new double[n];
    for(int i=0; i<n; ++i){
        x[i]=xo[i];
    }
    double* estimator, *res;
    double sum;
    for (int k = 0; k < iterations; ++k) {
        for (int i = 0; i < n; ++i) {
            sum=0;
            for (int j = 0; j < n; ++j) {
                if(i!=j)
                    sum+=A[i][j]*x[j];
            }
            u[i]=(b[i]-sum)/A[i][i];
        }

        estimator = errorEstimator(x, u, n);
        res = residuum(u, A, b, n);

        cout << "Iteracja nr: " << k+1 << "\nPrzyblizenie wyniku: ";
        wypiszWektor(u, n);
        cout << "\nEstymator bledu rozwiazania: ";
        wypiszWektor(estimator, n);
        cout << "\nResiduum: ";
        wypiszWektor(res, n);
        cout << "\n\n";

        for (int i = 0; i < n; ++i) {
            x[i]=u[i];
        }
        if(k+1==iterations){
            cout << "Przerwano ze wzgledu na wyczerpanie limitu iteracji.";
        }
        if(checkEstimator(estimator) && checkResiduum(res)){
            cout << "Przerwano ze wzgledu na kryterium dokladnosci wyznaczenia xn\noraz kryterium wiarygodnosci xn jako przyblizenia pierwiastka\nRozwiazanie: ";
            wypiszWektor(x, n);
            break;
        }
    }
}

void gaussSeidel(int n, double** A, double *b, double *xo, int iterations){
    auto* x1=new double[n];
    auto* x=new double[n];
    for(int i=0; i<n; ++i){
        x1[i]=xo[i];
        x[i]=xo[i];
    }
    double* estimator, *res;
    double sum;
    for (int k = 0; k < iterations; ++k) {
        for (int i = 0; i < n; ++i) {
            sum=0;
            for (int j = 0; j < n; ++j) {
                if(i!=j)
                    sum+=A[i][j]*x1[j];
            }
            x1[i]=(b[i]-sum)/A[i][i];
        }

        estimator = errorEstimator(x, x1, n);
        res = residuum(x1, A, b, n);

        cout << "Iteracja nr: " << k+1 << "\nPrzyblizenie wyniku: ";
        wypiszWektor(x1, n);
        cout << "\nEstymator bledu rozwiazania: ";
        wypiszWektor(estimator, n);
        cout << "\nResiduum: ";
        wypiszWektor(res, n);
        cout << "\n\n";

        for (int i = 0; i < n; ++i) {
            x[i]=x1[i];
        }
        if(k+1==iterations){
            cout << "Przerwano ze wzgledu na wyczerpanie limitu iteracji.";
        }
        if(checkEstimator(estimator) && checkResiduum(res)){
            cout << "Przerwano ze wzgledu na kryterium dokladnosci wyznaczenia xn\noraz kryterium wiarygodnosci xn jako przyblizenia pierwiastka\nRozwiazanie: ";
            wypiszWektor(x, n);
            break;
        }
    }
}

void sor(int n, double** A, double *b, double *xo, double omega, int iterations){
    auto* x1=new double[n];
    auto* x=new double[n];
    for(int i=0; i<n; ++i){
        x1[i]=xo[i];
        x[i]=xo[i];
    }
    double* estimator, *res;
    double sum;
    for (int k = 0; k < iterations; ++k) {
        for (int i = 0; i < n; ++i) {
            sum = 0;
            for (int j = 0; j < n; ++j) {
                if (i != j)
                    sum += A[i][j] * x1[j];
            }
            x1[i] = (1 - omega) * x[i] + omega * (b[i] - sum) / A[i][i];
        }

        estimator = errorEstimator(x, x1, n);
        res = residuum(x1, A, b, n);

        cout << "Iteracja nr: " << k+1 << "\nPrzyblizenie wyniku: ";
        wypiszWektor(x1, n);
        cout << "\nEstymator bledu rozwiazania: ";
        wypiszWektor(estimator, n);
        cout << "\nResiduum: ";
        wypiszWektor(res, n);
        cout << "\n\n";

        for (int i = 0; i < n; ++i) {
            x[i]=x1[i];
        }
        if(k+1==iterations){
            cout << "Przerwano ze wzgledu na wyczerpanie limitu iteracji.";
        }
        if(checkEstimator(estimator) && checkResiduum(res)){
            cout << "Przerwano ze wzgledu na kryterium dokladnosci wyznaczenia xn\noraz kryterium wiarygodnosci xn jako przyblizenia pierwiastka\nRozwiazanie: ";
            wypiszWektor(x, n);
            break;
        }
    }
}

int main() {
    ifstream plik("matrix.txt");

    int n;
    plik >> n;
    auto **matrix = new double*[n];
    auto *b = new double[n];
    auto *x = new double[n];

    for (int i = 0; i < n; ++i) {
        matrix[i]=new double[n];
        for (int j = 0; j < n; ++j) {
            plik >> matrix[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        plik >> b[i];
    }
    for (int i = 0; i < n; ++i) {
        plik >> x[i];
    }
    plik.close();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << matrix[i][j] << " ";
        }
        cout << '\n';
    }
    cout << '\n';
    for (int i = 0; i < n; ++i) {
        cout << b[i] << " ";
    }

    cout << "Metoda Jacobiego:\n";
    jacobi(n, matrix, b, x, 1000);
    cout << "\n-----------------------------------------------------------\n";
    cout << "Metoda Gaussa-Seidela:\n";
    gaussSeidel(n, matrix, b, x, 1000);
    cout << "\n-----------------------------------------------------------\n";
    cout << "Metoda SOR z parametrem omega=1/2:\n";
    sor(n, matrix, b, x, 0.5, 1000);


    return 0;
}
