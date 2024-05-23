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
typedef vector<complex<double>> vec_cmplx;
const double pi = 3.14159265358979323846; //LUIZA : Définition de pi


// Nil : Integrales pas OK
// NIL : Pour toutes les fonctions intégrales, il faut return la partie réelle de la valeur de l'intégrale (car la partie imaginaire c'est juste 0 mais pas vraiment avec la précision machine)

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(vector<T> const& diag,  vector<T> const& lower, vector<T> const& upper,
                 vector<T> const& rhs, vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (int i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

// LUIZA : Définition du potentiel V(x) : 
// NIL : ok
double V(double V0, double n_v, double x, double xL, double xR)
{
    return (0.5*V0*(1+cos(2*pi*n_v*(x-xL)/(xR-xL))));
}

// @TODO compute the folliwing quantities
// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - E:    calcule son energie moyenne,

// NIL : j'ai pas fait les calculs pour la méthode des trapèzes mais je te fais confiance

// LUIZA : Calcul de la probabilité que la particule se trouve à gauche de la barrière, avec xc le maximum local du potentiel :
//NIL : ok
double prob_left(double xL, double xR, double n_v, double dx, vec_cmplx psi, vector<double> x)
{//Définition des variables dont on a besoin
	double xc;
	int c;
	complex<double> integrale(0.00); 
	// Calcul des variables 
	xc = xL + (xR - xL)/(n_v);
	// Calcul de l'intégrale en utilisant la règle des trapèzes
	int i = 0;
	while (x[i] < xc)
	{
		
		// double psi_i_module = pow(real(psi[i]),2)+pow(imag(psi[i]),2);
		// double psi_iplus_module = pow(real(psi[i+1]),2)+pow(imag(psi[i+1]),2);

		complex<double> psi_i_module = conj(psi[i])*psi[i];
		complex<double> psi_iplus_module = conj(psi[i+1])*psi[i+1];

		integrale += dx*(psi_i_module + psi_iplus_module)/2.0;
		i++;
	}
    return real(integrale);
}
//LUIZA : Calcul de la probabilité qu'on trouve la particule à droite de la barrière, avec xc le maximum local du potentiel
//NIL : Ok
double prob_right(double xL, double xR, double n_v, int Npoints, double dx, vec_cmplx psi, vector<double> x)
{
	//Définition des variables dont on a besoin
	double xc;
	int c;
	complex<double> integrale(0.00) ; 
	// Calcul des variables
	xc = xL + (xR - xL)/(n_v);
	int i = 0;
	while (x[i] < xc)
	{
		i++;
	}
	for (int j = i; j < (Npoints-1); ++j)
	{
		complex<double> psi_i_module = conj(psi[j])*psi[j];
		complex<double> psi_iplus_module = conj(psi[j+1])*psi[j+1];

		integrale += dx*(psi_i_module + psi_iplus_module)/2.0;
	}
	return real(integrale);
}	

double Proba_totale(double xL, double xR, double n_v, int Npoints, double dx, vec_cmplx psi)
{
	//Définition des variables dont on a besoin
	double integrale(0.00) ; 
	// Calcul de l'intégrale en utilisant la règle des trapèzes
	for(int i = 0 ; i < (Npoints-1) ; ++i){
		integrale += (pow(abs(psi[i]),2) + pow(abs(psi[i+1]),2))/2;
	}
	integrale *= dx ; 
	return integrale;
}

// LUIZA : Définition de l'énergie de la particule, moyenne de l'Hamiltonien
// NIL : pas ok, voir commentaire dans la fonction
//NIL : re pas ok, le psi_star était mal défini 

//UTILISER LA MATRICE
double E(vec_cmplx psi, vector<double> x, double V0, double m, double n_v, double xL, double xR, double dx, double hbar, int Npoints, vec_cmplx dH, vec_cmplx aH, vec_cmplx cH)
{
	// // Définition des variables dont on a besoin :
	// vec_cmplx psi_star(Npoints) ;
	// complex<double> integrale(0.00); 
	// vec_cmplx derivee2_psi(Npoints);
	// //Calcul du conjugué complexe de psi
	// for(int i=0; i < Npoints ; ++i){ 
	// 	psi_star[i] = conj(psi[i]) ; 
	// }
	// //Calcul de la dérivée spatiale de psi :
	// //Pour les points d'extremité : 
	// derivee2_psi[0] = 0 ; 
	// derivee2_psi[Npoints-1] = 0;
	// //Pour les points intérieurs :
	// for(int i=1 ; i < (Npoints - 1) ; ++i){
	// 	derivee2_psi[i] = (psi[i+1]-2.0*psi[i]+psi[i-1])/pow(dx,2) ; 
	// }
	// // Calcul de l'intégrale en utilisant la règle du trapèze : 
	// for(int i = 0; i < (Npoints-1) ; ++i){
	// 	//paranthèse mal fermée
	// 	// ton code
	// 	// integrale += (psi_star[i]*(-(pow(hbar,2.0)/(2.0*m))*derivee2_psi[i]+V(V0,n_v,x[i],xL,xR)*psi[i]) + psi_star[i+1]*(-(pow(hbar,2.0)/(2.0*m))*derivee2_psi[i+1] + V(V0, n_v, x[i+1], xL, xR)*psi[i+1]))/2.0 ; 
	// 	//mon code :
	// 	integrale += (psi_star[i]*(-(pow(hbar,2.0)/(2.0*m))*derivee2_psi[i]+V(V0,n_v,x[i],xL,xR))*psi[i] + psi_star[i+1]*(-(pow(hbar,2.0)/(2.0*m))*derivee2_psi[i+1] + V(V0, n_v, x[i+1], xL, xR))*psi[i+1])/2.0 ;
	// }
    // integrale *= dx ; 

	// we have to compute psi^* H psi

	// compute H psi
	vec_cmplx H_psi(Npoints, 0.);
	for (int i(0); i < Npoints; ++i) {
		H_psi[i] = dH[i] * psi[i];
		if (i > 0)
			H_psi[i] += aH[i-1] * psi[i - 1];
		if (i < Npoints - 1)
			H_psi[i] += cH[i] * psi[i + 1];
	}

	// compute psi^* H psi
	complex<double> integrale(0.);
	for (int i(0); i < Npoints; ++i) {
		integrale += conj(psi[i]) * H_psi[i];
	}

	integrale *= dx;

    return real(integrale);
}

// LUIZA : Calcul de la position moyenne de la particule 
//NIL : pas ok, psi star
double xmoy(vec_cmplx psi, vector<double> x, double dx, int Npoints)
{
	// Définition des variables dont on a besoin : 
	vec_cmplx psi_star(Npoints) ;
	complex<double> integrale(0.00); 
	//Calcul du conjugué complexe de psi
	for(int i=0; i < Npoints ; ++i){ 
		psi_star[i] = conj(psi[i]) ; 
	}
	// Calcul de l'intégrale en utilisant la règle des trapèzes
	for(int i=0 ; i < (Npoints-1) ; ++i){
		integrale += ((psi_star[i]*x[i]*psi[i] + psi_star[i+1]*x[i+1]*psi[i+1])/2.0) ; 
	}
	integrale *= dx ;
		
    return real(integrale);
}
//LUIZA : Calcul de la position au carré moyenne de la particule : 
//NIL : pas ok, psi_star
double x2moy(vec_cmplx psi, vector<double> x, double dx, int Npoints)
{
	// Définition des variables dont on a besoin : 
	vec_cmplx psi_star(Npoints) ;
	complex<double> integrale(0.00); 
	//Calcul du conjugué complexe de psi
	for(int i=0; i < Npoints ; ++i){ 
		psi_star[i] = conj(psi[i]) ; 
	}
	// Calcul de l'intégrale en utilisant la règle des trapèzes
	for(int i = 0 ; i < (Npoints-1) ; ++i){
		integrale += ((psi_star[i]*pow(x[i],2)*psi[i] + psi_star[i+1]*pow(x[i+1],2)*psi[i+1])/2.0) ; 
	}
	integrale *= dx ; 
	
    return real(integrale);
}

// LUIZA : Définition de la quantité de mouvement moyenne de la particule 
//NIL : pas ok, psi_star
double pmoy(vec_cmplx psi, complex<double> complex_i, double dx, double hbar, int Npoints)
{
	// Définition des variables dont on a besoin :
	vec_cmplx psi_star(Npoints) ;
	complex<double> integrale(0.00); 
	vec_cmplx derivee_psi(Npoints);
	//Calcul du conjugué complexe de psi
	for(int i=0; i < Npoints ; ++i){ 
		psi_star[i] = conj(psi[i]) ; 
	}
	//Calcul de la dérivée spatiale de psi :
	//Pour le point d'extrémité gauche, en utilisant les différences finies "forward" :
	derivee_psi[0] = ((psi[1]-psi[0])/dx) ; 
	//Pour le point d'extrémité droite, en utilisant les différences finies "backward" : 
	derivee_psi[Npoints-1] = ((psi[Npoints-1]-psi[Npoints-2])/dx);
	//Calcul de la dérivée spatiale de psi pour les points intérieurs:
	for(int i=1 ; i < (Npoints-1) ; ++i){
		derivee_psi[i] = ((psi[i+1]-psi[i-1])/(2*dx)) ; 
	}
	//Calcul de l'intégrale en utilisant la règle du trapèze : 
	for(int i = 0; i < (Npoints-1) ; ++i){
		integrale += ((psi_star[i]*(-complex_i*hbar*derivee_psi[i]) + psi_star[i+1]*(-complex_i*hbar*derivee_psi[i+1]))/2.0) ;
	}	
	integrale *= dx ; 
	
    return real(integrale);
}

//LUIZA : Définition de la quantité de mouvement de la particule au carré moyenne :
//NIL : ok
double p2moy(vec_cmplx psi, double dx, double hbar, int Npoints)
{
	// Définition des variables dont on a besoin :
	vec_cmplx psi_star(Npoints) ;
	complex<double> integrale(0.00); 
	vec_cmplx derivee2_psi(Npoints);
	//Calcul du conjugué complexe de psi
	for(int i=0; i < Npoints ; ++i){ 
		psi_star[i] = conj(psi[i]) ; ; 
	}
	//Calcul de la dérivée spatiale de psi :
	//Pour les points d'extremité : 
	derivee2_psi[0] = 0 ; 
	derivee2_psi[Npoints-1] = 0;
	//Pour les points intérieurs :
	for(int i=1 ; i < (Npoints - 1) ; ++i){
		derivee2_psi[i] = ((psi[i+1]-2.0*psi[i]+psi[i-1])/(pow(dx,2))) ; 
	}
	// Calcul de l'intégrale en utilisant la règle du trapèze : 
	for(int i = 0 ; i < (Npoints - 1) ; ++i){
		integrale += ((psi_star[i]*(-pow(hbar,2.0)*derivee2_psi[i]) + psi_star[i+1]*(-pow(hbar,2.0)*derivee2_psi[i+1]))/2.0);
	}
	integrale *= dx ; 
    return real(integrale);
}
//LUIZA : Définition d'une fonction qui normalise psi : 
//NIL : pas ok, voir commentaire dans la fonction
vec_cmplx normalize(vec_cmplx psi, double dx, int Npoints)
{
	//Définition des variables dont on a besoin
    vec_cmplx psi_norm(Npoints);
    double norme(0.00) ;
    //Calcul de la "norme" de psi en utilisant la règle des trapèzes
    for(int i=0; i < (Npoints -1); ++i){
		norme += ((pow(abs(psi[i]),2)+pow(abs(psi[i+1]),2))/2);
	}
	norme *= dx ; 

	//ton code :
	// for(int i=0; i < Npoints ; ++i){
	// 	psi_norm[i] = (psi[i]/norme) ;
	// }

	//mon code :
	for (int i = 0; i < Npoints; ++i) {
		psi_norm[i] = psi[i] / sqrt(norme);
	}

    return psi_norm;
}
//Définition d'une autre fonction mieux 
//void normalize_l(vec_cmplx& psi, double dx, int Npoints)
//{
//	double norme(0.00);
//	for(int i=0; i<(Npoints-1) ; ++i){
//		norme += ((pow(abs(psi[i]),2) + pow(abs(psi[i+1]),2))/2) ; 
//	}
//	norme *= dx ; 
//	for(int i = 0; i< Npoints ; ++i){
//		psi[i] /= norme ;
//	}
//}

//LUIZA : définition de l'incertitude de la position : 
//NIL : ok
double xincertitude(vec_cmplx psi, vector<double> x, double dx, int Npoints)
{
	return (sqrt(x2moy(psi,x,dx,Npoints)-pow(xmoy(psi,x,dx,Npoints),2.0)));
}
//LUIZA : Définition de l'incertitude sur la quantité de mouvement 
//NIL : ok
double pincertitude(vec_cmplx psi, complex<double> complex_i, double dx, double hbar, int Npoints)
{
	return (sqrt(p2moy(psi,dx,hbar,Npoints)- pow(pmoy(psi, complex_i, dx, hbar, Npoints),2.0)));
}	

int
main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double V0 = configFile.get<double>("V0");
    double n_v = configFile.get<double>("n_v");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double

    // Parametres numeriques :

    double dt = configFile.get<double>("dt");
    int Nintervals = configFile.get<int>("Nintervals");
    int Npoints = Nintervals + 1;
    double dx = (xR - xL) / Nintervals;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    // LUIZA : Définition de la grille x : 
    x[0] = xL ; //Initialisationn de x_0
    x[Npoints-1] = xR ; //Initialisation de x_N ; 
    for(int i=1 ; i < (Npoints-1) ; ++i){
		x[i] = x[i-1] + dx ; 
	 }

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((xR*0.5 - xL)/(xR-xL)*Npoints); //chosen xR*0.5 since top of potential is at half x domain
  
    double x0 = configFile.get<double>("x0");
    double k0 = 2 * M_PI * n / (xR - xL);
    double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);
   
    //LUIZA : Initialisation du paquet d'onde
    // Caclul de C comme l'inverse de la "norme" de psi, i.e. le coefficient de normalisation
	//NIL : pas ok, tu normalises déjà psi à la ligne 361, donc on a pas besoin de le faire avant, d'autant plus qu'on ne peut pas normaliser psi avant de l'avoir défini
    // double norme(0.00) ; 
    // for(int i=0; i < (Npoints -1); ++i){
	// 	//ton code
	// 	// norme += ((pow(abs(psi[i]),2) + pow(abs(psi[i+1]),2))/2);
	// 	//
	// 	norme += (pow(real(psi[i]),2)+pow(imag(psi[i]),2)+pow(real(psi[i+1]),2)+pow(imag(psi[i+1]),2))/2;
	// }
	// norme *= dx ; 
	// double C(1.0/norme);
	//Nil : t'as aussi oublié de mettre au carré (x-x0)
	//Calcul du paquet d'onde initial
	for(int i=0; i < Npoints ; ++i){
		psi[i] = exp(complex_i*k0*x[i])*exp(-pow(x[i]-x0,2)/(2*pow(sigma0,2)));
	}

    //LUIZA : Modification des valeurs aux bords : 
    psi[0] = 0 ; 
    psi[Npoints-1]=0;
    
    // LUIZA : Normalisation  
    psi = normalize(psi,dx, Npoints);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a =
      complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a de l'equation (4.100)

	//NIL : j'ai rajt ça 
	complex<double> b = complex_i * dt / (2.*hbar); // Coefficient complexe b de l'equation (4.100)

	//LUIZA : Calcul des éléments de la matrice A : 
	// Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales supérieures et inférieures
	//NIL : pas ok, j'ai changé avec a et b pour simplifier le code
	//Aussi, t'as oublié V dans les éléments de la matrice A et B, mais c pcq je l avais oublié dans la matrice sur overleaf
	dA[0] = 1 ; 
	dA[Npoints-1] = 1 ; 
	cA[0]=0 ;
	cA[Nintervals-1]= 0 ; 
	aA[0] = 0;
	aA[Nintervals-1]=0;
	for(int i = 1; i < (Nintervals-1);++i){
		aA[i] = -a;
		cA[i] = -a;
	}
	for(int i = 1; i < (Npoints -1) ; ++i){
		dA[i] = (1.0 + 2.0*a + b*V(V0, n_v, x[i], xL, xR));
	}
	//LUIZA : Calcul des éléments de la matrice B :
	dB[0] = 1;
	dB[Npoints-1]=1;
	cB[0] = 0;
	cB[Nintervals-1] = 0 ; 
	aB[0] = 0;
	aB[Nintervals-1]=0;
	for(int i = 1; i < (Nintervals-1);++i){
		aB[i] = a;
		cB[i] = a;
	}
	for(int i = 1; i < (Npoints-1) ; ++i){
		dB[i] = (1.0 - 2.0*a - b*V(V0, n_v, x[i], xL, xR)) ; 
	}
	// LUIZA : Définition des éléments de la matrice H :
	//NIL : ok
	dH[0] = 1 ; 
	dH[Npoints-1] = 1 ; 
	cH[0] = 0;
	cH[Nintervals-1]=0;
	aH[0]=0;
	aH[Nintervals-1]=0;
	for(int i = 1; i < (Npoints -1) ; ++i){
		dH[i] = ((pow(hbar,2)/(m*pow(dx,2))) + V(V0, n_v, x[i], xL, xR));
	 }
	 for(int i = 1; i < (Nintervals-1);++i){
		 cH[i] = -(pow(hbar,2)/(2*m*pow(dx,2))) ; 
		 aH[i] = -(pow(hbar,2)/(2*m*pow(dx,2)));
	}
	
	


    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel(("./outputs/" + output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(V0, n_v, x[i], xL, xR) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi(("./outputs/" + output + "_psi2.out").c_str());
    fichier_psi.precision(15);

    ofstream fichier_observables(("./outputs/" + output + "_obs.out").c_str());
    fichier_observables.precision(15);
	
    // t0 writing
    for (int i(0); i < Npoints; ++i){
		// cout << "i = " << i;
		// cout << "psi[i] = " << psi[i] << endl;

        fichier_psi << pow(abs(psi[i]), 2)  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;


	// cout << "prob_left = " << prob_left(xL, xR, n_v, dx, psi) << endl;
	// cout << "prob_right = " << prob_right(xL, xR, n_v, Npoints, dx, psi) << endl;
	// cout << "E = " << E(psi, x, V0, m, n_v, xL, xR, dx, hbar, Npoints) << endl;
	// cout << "xmoy = " << xmoy(psi, x, dx, Npoints) << endl;
	// cout << "x2moy = " << x2moy(psi, x, dx, Npoints) << endl;
	// cout << "pmoy = " << pmoy(psi, complex_i, dx, hbar, Npoints) << endl;
	// cout << "p2moy = " << p2moy(psi, dx, hbar, Npoints) << endl;
	// cout << "xincertitude = " << xincertitude(psi, x, dx, Npoints) << endl;
	// cout << "pincertitude = " << pincertitude(psi, complex_i, dx, hbar, Npoints) << endl;

    // Ecriture des observables :
    fichier_observables << t << " " << prob_left(xL, xR, n_v, dx, psi,x) << " " << prob_right(xL, xR, n_v, Npoints, dx, psi,x) << " " << E(psi, x, V0, m, n_v, xL, xR, dx, hbar, Npoints,dH,aH,cH) << " " << xmoy(psi, x, dx, Npoints) << " "  
                << x2moy(psi, x, dx, Npoints) << " " << pmoy(psi, complex_i, dx, hbar, Npoints) << " " << p2moy(psi, dx, hbar, Npoints) << " " << xincertitude(psi, x, dx, Npoints) << " " 
                << pincertitude(psi, complex_i, dx, hbar, Npoints) << endl; 
    // Boucle temporelle :    
    while (t < tfin) {

        // TODO Calcul du membre de droite :
        //LUIZA : Je me casse la tête dessus depuis jsp cb de temps, j'ai strictement aucune idée de ce qu'il veut que je fasse 
        vec_cmplx psi_tmp(Npoints, 0.);
		// NIL :
		//enft, on veut résoudre l'équation (4.100) pour psi_tmp, donc on a besoin de calculer le membre de droite de l'équation
		//qui est donné par  B * psi
		//on a déjà calculé les éléments de la matrice B, donc on peut juste faire un produit matrice vecteur
		
		//ce que j ai ajt
		for (int i(0); i < Npoints; ++i) {
			psi_tmp[i] = dB[i] * psi[i];
			if (i > 0)
				psi_tmp[i] += aB[i-1] * psi[i - 1];
			if (i < Npoints - 1)
				psi_tmp[i] += cB[i] * psi[i + 1];
		}

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

		double psi_smodule;
        // t0 writing
        for (int i(0); i < Npoints; ++i){
			// NIL : même remarque pour le module de psi au carré
			psi_smodule = pow(real(psi[i]),2)+pow(imag(psi[i]),2); 
            fichier_psi << psi_smodule  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
            }
        fichier_psi << endl;

        // Ecriture des observables :
        fichier_observables << t << " " << prob_left(xL,xR,n_v,dx,psi,x) << " " << prob_right(xL, xR, n_v, Npoints, dx, psi,x)
                    << " " << E(psi, x, V0, m, n_v, xL, xR, dx, hbar, Npoints,dH,aH,cH) << " " << xmoy(psi, x, dx, Npoints) << " "  
                    << x2moy(psi, x, dx, Npoints) << " " << pmoy(psi, complex_i, dx, hbar, Npoints) << " " << p2moy(psi, dx, hbar, Npoints) << " " 
                    << xincertitude(psi, x, dx, Npoints) << " " << pincertitude(psi, complex_i, dx, hbar, Npoints) << endl; 
    } // Fin de la boucle temporelle



    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
