#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <valarray>
#include <cmath> 
#include <numeric>
using namespace std;

class Exercice3
{

private:
  double t, dt, tFin;
  double r0, r1, v0, a;
  double GM=6.674e-11;
  double mtot=0e0;
  double L2x,L2y;
  valarray<double> m = std::valarray<double>(0.e0, 2);
  int N_excit, nsteps;
  int sampling;
  int last;
  int  nsel_physics;
  bool adapt;
  double Omega;
  double alpha = 0e0;
  double beta = 0e0;
  double tol= 0e0;
  double L2 = 1.51098995045e+11;
  double d = 149.598023e9;
  valarray<double> x0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> x  = std::valarray<double>(0.e0, 4); // Correctly initialized
  ofstream *outputFile;

  void printOut(bool write)
  {
    if((!write && last>=sampling) || (write && last!=1))
    {
      double Energy = compute_energy(x[0],x[1],x[2],x[3]);
      *outputFile << t << " "<< x[0] << " " << x[1] << " "<< x[2] << " " << x[3] << " " \
      << Energy<< " "<< nsteps<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  double norm(double x, double y){
    return sqrt(pow(x,2)+pow(y,2));
  }

  double norm2(double x, double y){
    return pow(x,2)+pow(y,2);
  }

  std::valarray<double> get_f(double t, const std::valarray<double>& x) {
    std::valarray<double> xdot(0.0, 4);
    double r;
    double M_T = m[1];
    double M_S = m[0];
    //TO DO
    if (nsel_physics == 1) {
      r = sqrt(pow(x[0],2) + pow(x[1],2));
      xdot[2] = -GM*m[1]*x[0]/pow(r,3);
      xdot[3] = -GM*m[1]*x[1]/pow(r,3);
    }
    else if (nsel_physics == 2) {
      double x_S = -d*alpha;
      double y_S = 0;
      double x_T = d*beta;
      double y_T = 0;

      double r_S = norm(x[0]-x_S, x[1]-y_S);
      double r_T = norm(x[0]-x_T, x[1]-y_T);

      xdot[2] = -GM*M_S*(x[0]-x_S)/pow(r_S,3) - GM*M_T*(x[0]-x_T)/pow(r_T,3) + pow(Omega, 2)*x[0] + 2*Omega*x[3];
      xdot[3] = -GM*M_S*(x[1]-y_S)/pow(r_S,3) - GM*M_T*(x[1]-y_T)/pow(r_T,3) + pow(Omega, 2)*x[1] - 2*Omega*x[2];
    }
    else{
        cerr << "No dynamics corresponds to this index" << endl;
        return xdot;
    }

    // Velocities
    xdot[0] = x[2];
    xdot[1] = x[3];

    return xdot;
  }  


// Function to compute potential energy per mass in R (nsel_physics=1) or in R'(nsel_physics=2)
double get_Epot(double xx, double yy) {
    //TO DO
    double Epot;
    double M_T = m[1];
    double M_S = m[0];
    if(nsel_physics==1){
      Epot = -GM*M_T/sqrt(pow(xx,2)+pow(yy,2));
  }
  else{
      double x_S = -d*alpha;
      double y_S = 0;
      double x_T = d*beta;
      double y_T = 0;

      double r_S = norm(xx-x_S, yy-y_S);
      double r_T = norm(xx-x_T, yy-y_T);

      Epot = -GM*M_T/r_T - GM*M_S/r_S;

  }
 
    return Epot;
}

// Function to compute mechanical energy per mass in R'
double compute_energy(double xx, double yy, double vx, double vy) {
    double Emec;

    double Epot = get_Epot(xx, yy);
    double Ekin;

    if (nsel_physics==1){
      Ekin = 0.5*norm2(vx, vy);
    }else if (nsel_physics==2)
    {
      Ekin = 0.5*(norm2(vx,vy) + pow(Omega,2)*norm2(xx,yy) + 2*Omega*(xx*vy - yy*vx));
    }
    
    return Epot + Ekin;
}
void initial_condition(void){
  if(nsel_physics==1){
    //TO DO initialize x0
    x0[0] = -r0; //x
    x0[1] = 0; //y
    x0[2] = 0; //vx
    x0[3] = sqrt(2*GM*m[1]*r1/(r0*(r0+r1))); //vy
  }
  else{
    //TO DO initialize x0
    x0[0] = L2; //x
    x0[1] = 0; //y              
    x0[2] = 0; //vx
    x0[3] = -0.1; //vy
    
  }
}


std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double ti, double dt) {
    std::valarray<double> k1, k2, k3, k4, ynew;
    k1 = dt*get_f(ti, yold);
    k2 = dt*get_f(ti+0.5*dt, yold+0.5*k1);
    k3 = dt*get_f(ti+0.5*dt, yold+0.5*k2);
    k4 = dt*get_f(ti+dt, yold+k3);

    ynew = yold + (k1 + 2.0*k2 + 2.0*k3 + k4)/6.0;

    return ynew;
}


public:
  Exercice3(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin         = configFile.get<double>("tFin");            // t final (overwritten if N_excit >0)
  //  dt       = configFile.get<double>("dt");            // time step (overwritten if nsteps_per >0)
    m[0]         = configFile.get<double>("m1");              // mass of the sun
    m[1]         = configFile.get<double>("m2");              // mass of the earth
    r0           = configFile.get<double>("r0");              // r0
    r1           = configFile.get<double>("r1");              // r1
    L2x          = configFile.get<double>("L2x");              // L2x
    L2y          = configFile.get<double>("L2y");              // L2y
    a            = configFile.get<double>("a");               // demi grand-axe (solei-terre hyp MCU)
    nsel_physics = configFile.get<int>("nsel_physics");       //1) one body problem around mass#2 or 2) one body in rotating reference frame of {1,2}
    adapt        = configFile.get<bool>("adapt");             //if 1=true -> adaptive dt, if 0=false -> fixed dt
    tol          = configFile.get<double>("tol");             //tolerance of the adaptive scheme
    sampling     = configFile.get<int>("sampling");     // number of time steps between two writings on file
    nsteps       = configFile.get<int>("nsteps");        // number of time step per period
    Omega        = configFile.get<double>("Omega");       // angular velocity of the rotating frame
    mtot=m[0]+m[1];
    alpha = m[1] / mtot;
    beta = m[0] / mtot;
    //TO DO
    
    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
    //TO DO initialize tFin for nsel_physics=1 and initialize dt for both nsel_physics
    if (nsel_physics==1){
      tFin=pi*sqrt(pow(((r1+r0)/2),3)/(GM*m[1])); //demi periode
      tFin=tFin; //periode
      cout << tFin;
    }else{
      //2years in second
      
    } 
    
    dt=tFin/nsteps;
  }

  ~Exercice3()
  {
    outputFile->close();
    delete outputFile;
  };
  

  void run()
  {
    t = 0.;
    initial_condition();
    x=x0;
    last = 0;
    printOut(true);
    std::valarray<double> y1;
    std::valarray<double> y2;

    if (adapt==false){
      //To DO fixed dt scheme
      while(t<tFin){
        y1 = RK4_do_onestep(x, t, dt);

        x = y1;

        t += dt;
        printOut(false);
      }
    }
    else{
      //TO DO adaptive case
      double d;
      double f(0.9);
      double n(4.0); //ordre de la méthode
      nsteps = 0;
      
      while (t<tFin){

        dt = min(dt, tFin - t);

        y1 = RK4_do_onestep(x, t, dt);
        y2 = RK4_do_onestep(RK4_do_onestep(x, t, dt/2), t+dt/2, dt/2);

        d = sqrt( pow(y1[0]-y2[0],2) + pow(y1[1]-y2[1],2) + pow(y1[2]-y2[2],2) + pow(y1[3]-y2[3],2) );


        if (d>tol) {
          while(d>tol){
            dt = f*dt*pow(tol/d,1/(n+1));

            y1 = RK4_do_onestep(x, t, dt);
            y2 = RK4_do_onestep(RK4_do_onestep(x, t, dt/2), t+dt/2, dt/2);

            d = sqrt( pow(y1[0]-y2[0],2) + pow(y1[1]-y2[1],2) + pow(y1[2]-y2[2],2) + pow(y1[3]-y2[3],2) );
          }
          x = y2;
          t+=dt;
        }else{
          x = y2;
          t+=dt;
          dt = dt*pow(tol/d,1/(n+1));
        }

        nsteps +=1;
        printOut(false);
      }
    };
    
    printOut(true); // ecrire le dernier pas de temps
  };




};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();

  return 0;
}
