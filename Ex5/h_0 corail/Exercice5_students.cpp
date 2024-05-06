#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;

//LUIZA : Définition des conditions au bord : 
//Les conditions étaient sur x_1 et x_2, ainsi que x_{Nx} et x_{Nx -1}, je suis partie du principe que x_1 
// correspondait à l'indice 0, x_2 à l'indice 1, x_{Nx} à l'indice N et x_{Nx-1} à l'indice N-1
// Pour beta, comme c'est un vecteur de beta 2 je suis partie du principe que le beta utilisé dans la théorie c'était beta[i]
void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, \
		double const& t,double const& dt, \
		vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
  // Condition au bord gauche : 
  // NIL : OK
	if(bc_l == "fixe"){
		fnext[0] = fnow[0] ;
	};
  // NIL : OK
	if( bc_l == "libre"){
		fnext[0] = fnext[1] ;  
	}; 
  // NIL : t'as fait les calculs pour obtenir ça ? je te fais confiance
	if(bc_l == "sortie"){
		fnext[0] = fnow[0] - sqrt(beta2[0])*(fnow[0] - fnow[1]) ; 
	};
	//Condition au bord droit
  // NIL : OK
	if(bc_r == "fixe"){
		fnext[N] = fnow[N] ; 
	};
  // NIL : OK
	if(bc_r == "libre"){
		fnext[N] = fnext[N-1]; 
	};
  //NIL : OK
	if(bc_r == "sortie"){
		fnext[N] = fnow[N] - sqrt(beta2[N]) *(fnow[N] - fnow[N-1]) ;
	}
}

// LUIZA : Définition des modes propres : 
// NIL : OK
double finit_mode(double x, double xL, double n_init, double xR)
{
  double finit_(0.0);
  const double PI = 3.1415926535897932384626433832795028841971e0;
  finit_ = sin(((n_init*PI)/(xR - xL))*x) ; 
  return finit_;
}

//NIL : je rajoute pour le cas ou on initialise pas avec mode propre
double finit_eq4(double x, double x1, double x2, double A)
{
  double finit_(0.0);
  const double PI = 3.1415926535897932384626433832795028841971e0;
  if (x <= x1)
    finit_ = 0;
  else if (x < x2)
    finit_ = (A / 2)*(1 - cos(2*PI*((x-x1)/(x2-x1))));
  else
    finit_ = 0;

  return finit_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());
  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if(len > 0)
    o << v[len-1];
  return o;
}

//
// Main
//
int main(int argc, char* argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g  = 9.81;
  double dx;
  double dt;
  double t;
  double Nsteps;
  int stride(0);

  string inputPath("input_example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres de simulation :
  double tfin    = configFile.get<double>("tfin");
  int nx          = configFile.get<int>("nx"); // nb d'intervalles
  int N = nx+1;                                // nb de pts de maillage
  double CFL     = configFile.get<double>("CFL");
  double nsteps  = configFile.get<double>("nsteps");
  double A       = configFile.get<double>("A");
  double n_init   = configFile.get<double>("n_init");
  double hL      = configFile.get<double>("hL");
  double hR      = configFile.get<double>("hR");
  double hC      = configFile.get<double>("hC");
  double h00     = configFile.get<double>("h00"); // profondeur, cas uniforme
  double x1      = configFile.get<double>("x1");
  double x2      = configFile.get<double>("x2");
  double xa      = configFile.get<double>("xa");
  double xb      = configFile.get<double>("xb");
  double xc      = configFile.get<double>("xc");
  double xd      = configFile.get<double>("xd");
  double xL      = configFile.get<double>("xL");
  double xR      = configFile.get<double>("xR");
  int n_stride(configFile.get<int>("n_stride"));
 
// Conditions aux bords:
  string bc_l           = configFile.get<string>("cb_gauche");
  string bc_r           = configFile.get<string>("cb_droite");

// Type de forme initiale de la vague: selon donnée Eq.(4) ou mode propre
// (par exemple 'mode' pour mode propre, autrement Eq.(4))
  string initialization = configFile.get<string>("initialization"); 

// Onde partant initialement vers la gauche ou vers la droite ou statique
// (par exemple 'left', 'right', 'static')
  string initial_state = configFile.get<string>("initial_state");

// Selecteur pour le cas h0 uniforme:
  bool v_uniform        = configFile.get<bool>("v_uniform");

// Selecteur pour choisir le pas de temps:
// true --> dt=tfin/nsteps; t final est exactement tfin
// false --> dt tel que beta_CFL=1; attention, t final n'est pas exactement tfin
  bool impose_nsteps    = configFile.get<bool>("impose_nsteps");

  vector<double> h0(N) ;   // profondeur aux points de maillage
  vector<double> vel2(N) ; // u^2 = g*h_0 aux points de maillage
  vector<double> x(N) ;    // positions des points de maillage
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = (xR - xL) / (N-1);
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non

 // Eq.(1) ou Eq.(2) [ou Eq.(6) (facultatif)]: Eq1, Eq2 ou Eq6
  string equation_type = configFile.get<string>("equation_type");  
  
// LUIZA : définition des positions des points de maillage
// Ca n'a pas été fait dans le code, or on en a vraiment besoin donc je l'ai fait 
// NIL : OK
  x[0] = xL ; 
  for(int i(1); i<=N; ++i){
	  x[i] = x[i-1] + dx ; 
  };

// LUIZA : Initialisation des indices : 
// Donc ici le principe est qu'on a beaucoup de conditions pour les initialisations sur les x : 
// entre x_a et x_b, entre x_b et x_c etc... Donc sous le conseil plus que douteux d'un assistant
// j'ai traduit ces conditions sur x par des conditions sur les indices, un peu d'algèbre et on trouve 
// que x_a correspond à l'indice a comme défini ici (dans un maillage uniforme) : 

// NIL : Pas besoin, voir autres commentaires
  // int a((xa - xL)/dx); // indice correspondant à xa
  // int b((xb - xL)/dx); // indice correspondant à xb
  // int c((xc - xL)/dx); // indice correspondant à xc
  // int d((xd - xL)/dx); // indice correspondant à xd
  // int i_1((x1 - xL)/dx) ; //indice correspondant à x1 
  // int i_2 ((x2 - xL)/dx); //indice correspondant à x2 
 
// LUIZA : Initialisation de la profondeur h0 et de la vitesse au carré vel2
// ici j'ai utilisé la théorie du cours et de l'énoncé
// NIL : OK
  for(int i(0); i<N; ++i){ 
    if(v_uniform){
		h0[i] = h00 ;  
	}else{
    //NIL : il suffit d evaluer en x[i] enft, donc pas besoin d indices
    //Alexis: J'ai changé les and et des conditions
		if((x[i] >= xL) && (x[i] <= xa)){
			h0[i] = hL ;
		};
		if((x[i] > xa) && (x[i] < xb)){
			h0[i] = 0.5*(hL + hC) + 0.5*(hL - hC)*cos(PI*((x[i]-xa)/(xb -xa))) ; 
		};
		if((x[i] >= xb) && (x[i] <= xc)){
			h0[i] = hC ; 
		};
		if((x[i] > xc) && (x[i] < xd)){
			h0[i] = 0.5*(hR + hC) - 0.5*(hR - hC)*cos(PI*((x[i]-xc)/(xd - xc))) ; 
		};
		if((x[i] >= xd) && (x[i] <= xR)){
			h0[i] = hR ; 
		};
	};
  //NIL : OK
	vel2[i] = g*h0[i] ; 
  };

  // cout << vel2 << endl;

  auto max_vel2_iter = std::max_element(vel2.begin(), vel2.end());
  double max_vel2 = *max_vel2_iter;

  // LUIZA : Définition de dt avec le CFL 
  // Ici ce que j'ai fait c'est définir la norme de la vitesse u, par rapport à la vitesse au carré vel2
  // pour définir dt avec le CFL, comme demandé dans le code initial. 
  // Je vois pas trop d'autres moyens
  // NIL : pas OK, j'ai retema l'énoncé et en fait faut définir dt avec le CFL et max_vel2

  // ton code
  // double somme(0); 
  // double norme_u(0) ; 
  // for(int i(0); i<N; ++i){
	// somme += vel2[i] ;
  // }
  // norme_u = sqrt(somme) ;

  // mon code
  //ici on veut que dt soit tel que max(beta_CFL) = 1, donc on a dt = dx/sqrt(max_vel2) 
  dt = dx/sqrt(max_vel2) ; 
  
  // LUIZA : définition de dt et CLF quand on veut fixer le nombre de nsteps 
  // J'ai fait la même chose que juste au dessus mais dans le sens inverse. 
  // NIL : PAS OK, voir commentaire ci-dessous
  if(impose_nsteps){ 
    dt = tfin/nsteps; // je propse ça plutôt
    // CFL = max_vel2 * dt / dx; // pas nécessaire, jsp a quoi sert le paramètre CFL enft
  }

  

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_x((output + "_x").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output + "_v").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output + "_f").c_str());
  fichier_f.precision(15);


  // Initialisation des tableaux du schema numerique :

  //LUIZA : Initialization de fnow et fpast : 
  // ATTENTION : J'ai fait les conditions avec les noms en input en anglais
  // PROBLEME : Initialization de beta2 qui fait pas beaucoup de sens, a revoir
  for(int i(0); i<N; ++i) {
  // NIL : OK
  //Forme initiale de la vague
	if(initialization == "mode"){ //Initialisation avec les modes propres
		fnow[i] = finit_mode(x[i], xL, n_init, xR) ; 
	}else{						// Initialisation avec l'equation (4) 
		// if(i <= i_1){
		// 	fnow[i] = 0;
		// };
		// if((i > i_1) and (i < i_2)){
		// 	fnow[i] = (A / 2)*(1 - cos(2*PI*(x[i]-x1)/(x2-x1)));
		// };
		// if((i>=i_2) and (i<N)){
		// 	fnow[i] = 0 ; 

    //NIL : pas besoin de faire avec des indices !! (voir fct que g déf. finit_eq4)
    fnow[i] = finit_eq4(x[i], x1, x2, A);
		
	}
  // NIL : OK
  //mode de propagation initial de la vague
	if(initial_state == "static"){ //Onde initialement au repos
		fpast[i] = fnow[i] ; 
	};
	if(initial_state == "right"){ // Onde se propageant vers la droite
    // NIL : le truc avec l'indice marche pas, car la diff est trop petite, donc g def une fct eq4 qui le fait de maniere continue
		fpast[i] = finit_eq4(x[i] + sqrt(vel2[i])*dt,x1,x2,A); 
	};
	if(initial_state == "left"){ // Onde se propageant vers la gauche
    // NIL : same
		fpast[i] = finit_eq4(x[i] - sqrt(vel2[i])*dt,x1,x2,A); 
	};
  // NIL : OK
  beta2[i] = vel2[i]*pow(dt/dx, 2);
  
  };
	

  // cout <<"beta2[0] is "<<beta2[0]<<endl;
  // cout <<"dt is "<< dt <<endl;


  // Boucle temporelle :
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
    } else if (t==tfin-dt) {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
    }

    ++stride;

    // Evolution :
    for(int i(1); i<N-1; ++i)
    {
      // LUIZA : Ecriture de l'expression pour fnext 
      // Fait avec l'expression des notes du cours
      // NIL : OK
      if (equation_type == "Eq1")
      {
        
        fnext[i] = 2*(1 - beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1] + fnow[i-1]);
      }
      else if (equation_type == "Eq2")
      {
        fnext[i] = 2*(1 - beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1] + fnow[i-1]) + 0.25*(beta2[i+1] - beta2[i-1])*(fnow[i+1] - fnow[i-1]);

        // cout << "Eq2" << endl;
        // fnext[i] = 2*(1 - beta2[i])*fnow[i] - fpast[i] + (beta2[i+1]/4 + beta2[i] - beta2[i-1]/4)*(fnow[i+1] + fnow[i-1]);
        // cout << "fnext[i] = " << fnext[i] << endl;
        // if (t==0.)
        // {
        //   cout << "beta2[i+1] = " << beta2[i+1] << endl;
        // cout << "beta2[i] = " << beta2[i] << endl;
        // }
      }
    }

    //Application des conditions au bord : 
    boundary_condition(fnext, fnow, A, t, dt, beta2, bc_l, bc_r, N);

    //LUIZA : Mise à jour
    // NIL : OK
    fpast = fnow;
    fnow  = fnext;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl;
  fichier_x << x << endl;
  fichier_v << vel2 << endl;

  fichier_f.close();
  fichier_x.close();
  fichier_v.close();

  return 0;
}
