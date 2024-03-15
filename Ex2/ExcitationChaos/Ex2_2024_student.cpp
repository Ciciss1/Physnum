#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> // Se usi fun
using namespace std;

class Exercice2
{

private:
  double t, dt, tFin;
  double m, g, L, alpha;
  double d, Omega, kappa;
  double theta, thetadot;
  int N_excit, nsteps_per;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      // Définition de l'énergie mécanique et de la puissance à l'aide des fonctions définies ultérieurement
      double emec = Emec(theta, thetadot, t);
      double pnc = Pnonc(theta, thetadot, t);

      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

    // Définition des fonctions d'accélération a_1 et a_2
  valarray<double> acceleration(double theta_, double thetadot_, double t_)
  {
    valarray<double> acc = valarray<double>(2);

    acc[0] = - g*sin(theta_)/length(t_); // accélération dépendant uniquement de theta
    acc[1] = (-2/length(t_))*lendot(t_)*thetadot_ ; // accélération dépendant uniquement de la dérivée temporelle de theta

    return acc;
  }

  // Définition de l'énergie mécanique à l'aide de fonctions définies ultérieurement
  double Emec(double theta_, double thetadot_, double t_)
  {
    double E;
    E = 0.5*m*(pow(lendot(t_),2) + pow(length(t_)*thetadot_,2)) - m*g*length(t_)*cos(theta_);
    return  E;
  }

    // Définition de la puissance des forces non conservatives à l'aide de fonctions définies ultérieurement
  double Pnonc(double theta_, double thetadot_, double t_)
  {
    double thetadotdot_ = (-g*sin(theta_) - 2 * lendot(t_) * thetadot_)/length(t_);
    double P = m*(lendot(t_) * lendotdot(t_) + length(t_) * lendot(t_) * pow(thetadot_,2) + thetadot_ * thetadotdot_ * pow(length(t_),2) - g * lendot(t_) * cos(theta_) + g * length(t_) * thetadot_ * sin(theta_));
    return  P;
  }

  // Définition de la longueur du pendule qui dépend du temps
  double length(double t_)
  {
    double l;
    l = L + alpha*t_ + d * sin(Omega*t_);
    return  l;
  }

  // Définition de la première dérivée temporelle de la longueur du pendule
  double lendot(double t_)
  {
    double ldot;
    ldot = alpha + Omega*d*cos(Omega*t_);
    return  ldot;
  }
  // Définition de la deuxième dérivée temporelle de la longueur du pendule
  double lendotdot(double t_)
  {
    double ldotdot;
    ldotdot = -Omega*Omega*d*sin(Omega*t_);
    return ldotdot;
  }


    void step()
  {
    
    // Définition des variables dont on a besoin pour faire le schema numerique
    double theta_old        = theta;
    double thetadot_old     = thetadot;

    valarray<double> a_ = acceleration(theta_old, thetadot_old, t);

    double acc = a_[0] + a_[1];

    theta = theta_old + thetadot_old*dt + 0.5*acc*dt*dt;

    double v_demi = thetadot_old + 0.5*acc*dt;

    valarray<double> acc_demi = acceleration(theta, v_demi, t+0.5*dt);
    valarray<double> acc_plus1 = acceleration(theta, thetadot, t+dt);

    thetadot = thetadot_old + 0.5*(a_[0] + acc_plus1[0])*dt + acc_demi[1]*dt;
  }


public:
  Exercice2(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin     = configFile.get<double>("tFin");      // t final (réécrit si N >0)
    d        = configFile.get<double>("d");         // amplitude du changement periodic de longueur
    Omega    = configFile.get<double>("Omega");     // fréquence du changement periodic de longueur 
    kappa    = configFile.get<double>("kappa");     // coefficient de frottement
    m        = configFile.get<double>("m");         // mass
    g        = configFile.get<double>("g");         // accélération de la gravité
    L        = configFile.get<double>("L");         // longueur
    alpha    = configFile.get<double>("alpha");     // Coefficient linéaire du changement de longueur
    theta    = configFile.get<double>("theta0");    // condition initiale en theta
    thetadot = configFile.get<double>("thetadot0"); // condition initiale en thetadot
    sampling = configFile.get<int>("sampling");     // fréquence d'écriture des données
    N_excit  = configFile.get<int>("N");            // Nombre de périodes d'excitation simulées
    nsteps_per= configFile.get<int>("nsteps");       // nombre de pas de temps par période (si N>0), nombre total de pas de temps si N=0

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);

    if(N_excit>0){
      // TODO Définir le pas temporel et le temp final si N_excit est défini.
      //periode. d'exc. = T = 2pi/omega = pi/omega_0
      double omega_0 = sqrt(g/L);
      double T_exc = pi/omega_0;
      tFin = N_excit*T_exc;
      dt   = T_exc/nsteps_per;
    }
    else{
      // TODO Définir le pas temporel si N_excit n'est pas défini.
      dt = tFin / nsteps_per;
    }
  }

  ~Exercice2()
  {
    outputFile->close();
    delete outputFile;
  };

    void run()
  {
    t = 0.;
    last = 0;
    printOut(true);

    while( t < tFin-0.5*dt )
    {
      step();

      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};

int main(int argc, char* argv[])
{
  Exercice2 engine(argc, argv);
  engine.run();

  return 0;
}
