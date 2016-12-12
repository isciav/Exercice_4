#include <iostream>  // prova
#include <fstream>
#include <cmath>
#include <iomanip>
#include <valarray>
#include "ConfigFile.tpp" // Fichier .tpp car inclut un template

using namespace std;

//prova iscia

const double G(6.674e-11);

//--------------------DEFINITIONS CLASSES----------------------------
class Exercice4{
	private:
			// Definitions variables
		double t,dt,tFin;
		double m1,m2,m3;
		double epsilon;
		valarray<double> Y;
		int sampling;
		int last;
		ofstream *outputFile;
		
			// PrintOut
		void printOut(bool force){
			if((!force && last>=sampling) || (force && last!=1)){
				energie();
				
				*outputFile << t << " ";
				
				
				for (size_t i(0); i<Y.size(); ++i){
					*outputFile << Y[i] << " ";
				}				
				*outputFile << endl;

				last=1;
			}else{
				last++;
			}
		};
		
			// Fonction qui calcule l'énergie
		void energie(){ 
			
		}
		
			// Fonctions necessaires pour la fonction step()
		valarray<double> acceleration(valarray<double> const& y, double tictac){
			valarray<double> v(0., 6);
			valarray<double> a(0., 6);
			
			v = y[slice(6,6,1)];
			a[0] = G*m2*(y[2]-y[0])/(pow(dist(y[0], y[2], y[1], y[3]),3)) +  G*m3*(y[4]-y[0])/(pow(dist(y[0], y[4], y[1], y[5]),3));	// acceleration sur x de la masse 1
			a[1] = G*m2*(y[3]-y[1])/(pow(dist(y[0], y[2], y[1], y[3]),3)) +  G*m3*(y[5]-y[1])/(pow(dist(y[0], y[4], y[1], y[5]),3));	// acceleration sur y de la masse 1
			a[2] = G*m1*(y[0]-y[2])/(pow(dist(y[0], y[2], y[1], y[3]),3)) +  G*m3*(y[4]-y[2])/(pow(dist(y[2], y[4], y[3], y[5]),3));	// acceleration sur x de la masse 2
			a[3] = G*m1*(y[1]-y[3])/(pow(dist(y[0], y[2], y[1], y[3]),3)) +  G*m3*(y[5]-y[3])/(pow(dist(y[2], y[4], y[3], y[5]),3));	// acceleration sur y de la masse 2
			a[4] = G*m1*(y[0]-y[4])/(pow(dist(y[4], y[0], y[5], y[1]),3)) +  G*m2*(y[2]-y[4])/(pow(dist(y[4], y[2], y[3], y[5]),3));	// acceleration sur x de la masse 3
			a[5] = G*m1*(y[1]-y[5])/(pow(dist(y[0], y[4], y[1], y[5]),3)) +  G*m2*(y[3]-y[5])/(pow(dist(y[2], y[2], y[5], y[3]),3));	// acceleration sur y de la masse 3
			
			valarray<double> f(0., 12);
			f[slice(0,5,1)] = v;
			f[slice(6,11,1)] = a;
			
			return f;
		}
		
		double dist(double x1, double x2, double y1, double y2){
			double d(sqrt(pow(x1-x2,2) + pow(y1-y2,2)));
			if (d == 0) {
				cout << "Division par zero dans la distance!!" << endl;
				return 1;
			}	
			return d;
		}
		
			// step() pour Runge-Kutta 4
		void step(valarray<double> y, double deltat){
			valarray<double> k1(0., 12);
			valarray<double> k2(0., 12);
			valarray<double> k3(0., 12);
			valarray<double> k4(0., 12);
			
			k1 = deltat*acceleration(y, t);
			k2 = deltat*acceleration(y + 0.5*k1, t + deltat*0.5);
			k3 = deltat*acceleration(y + 0.5*k2, t + deltat*0.5);
			k4 = deltat*acceleration(y + k3, t + deltat);
			y = y + (k1 + 2.*k2 + 2.*k3 + k4)/6.;
		}
		
			// restitue le max
		double max(double a, double b, double c){
			if(a>=b && a>=c){
				return a;
			} else { if (b>=a && b>=c){
				return b;
			} else { 
				return c;
			}
			}
		}
		
	
	public:
			// Constructeur Exercice4
		Exercice4(int argc, char* argv[])
		:Y(12)
		{
    
			string inputPath("configuration.in"); // Fichier d'input par defaut
			if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
				inputPath = argv[1];

				ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

				for(int i(2); i<argc; ++i) // Input complementaires ("./Onde config_perso.in input_scan=[valeur]") argc a ete enleve temporairement
				configFile.process(argv[i]);
				
				tFin     = configFile.get<double>("tFin");
				dt       = configFile.get<double>("dt");
				m1       = configFile.get<double>("m1");
				m2       = configFile.get<double>("m2");
				m3       = configFile.get<double>("m3");
				Y[0]     = configFile.get<double>("x01");
				Y[1]     = configFile.get<double>("y01");
				Y[2]     = configFile.get<double>("x02");
				Y[3]     = configFile.get<double>("y02");
				Y[4]     = configFile.get<double>("x03");
				Y[5]     = configFile.get<double>("y03");
				Y[6]     = configFile.get<double>("vx01");
				Y[7]     = configFile.get<double>("vy01");
				Y[8]     = configFile.get<double>("vx02");
				Y[9]     = configFile.get<double>("vy02");
				Y[10]    = configFile.get<double>("vx03");
				Y[11]    = configFile.get<double>("vy03");
				sampling = configFile.get<int>("sampling");
    
			// Ouverture du fichier de sortie
			outputFile = new ofstream(configFile.get<string>("output").c_str());
			outputFile->precision(15);
		};
		
			// Destructeur  
		~Exercice4(){
			outputFile->close();
			delete outputFile;
		};
		
			// run()
		void run(){
			last = 0;
			t = 0;
			printOut(true);
			
			while( t<(tFin-0.5*dt) ) {
				valarray<double> Yprim(Y);
				double d1, d2, d3, d ;
				
				step(Y, dt);
				step(Yprim, dt*0.5);
				step(Yprim, dt*0.5);
					// if
				
				d1 = sqrt(pow(Y[0]-Yprim[0],2) + pow(Y[1]-Yprim[1],2));	
				d2 = sqrt(pow(Y[2]-Yprim[2],2) + pow(Y[3]-Yprim[3],2));	
				d3 = sqrt(pow(Y[4]-Yprim[4],2) + pow(Y[5]-Yprim[5],2));
				
				d = max(d1, d2, d3);
				
				if (d<epsilon){
					dt = dt * pow(epsilon/d,1./5.);
				} else {
					do {	
					dt = 0.99 * dt * pow(epsilon/d,1./5.);
					
					step(Y, dt);
					step(Yprim, dt*0.5);
					step(Yprim, dt*0.5);
					
					d1 = sqrt(pow(Y[0]-Yprim[0],2) + pow(Y[1]-Yprim[1],2));	
					d2 = sqrt(pow(Y[2]-Yprim[2],2) + pow(Y[3]-Yprim[3],2));	
					d3 = sqrt(pow(Y[4]-Yprim[4],2) + pow(Y[5]-Yprim[5],2));
					
					} while (d > epsilon);
				}
				
				t += dt;
				printOut(false);
			}
			
			printOut(true);
		};
};


//-------------------------------MAIN----------------------------------
int main(int argc, char* argv[]) 
{
  Exercice4 planete(argc, argv);
  planete.run();
  return 0;
}
