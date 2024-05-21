#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

int main(){
    //random complex number

    complex<double> comp(1.0, 2.0);

    complex<double> comp_star(conj(comp));

    comp_star = real(comp) - imag(comp) ;

    cout << "comp = " << comp << endl;
    cout << "comp_star = " << comp_star << endl;

    return 0;
}
