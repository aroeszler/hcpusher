//
//  main.cpp
//  hcsecondtry.cpp
//
//  Created by Annegret Roeszler on 16.03.20.
//  Copyright © 2020 Annegret Roeszler. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <vector>
#include "vector_add.hpp"

using namespace std;
using Vec3D = vector<double>;
int main ()
{
    double m = 9.1e-31; //Masse in kg
    double q = -1.6022e-19; //Ladung in As
    double c = 299792458;
    //double c = 3e8; // in m/s

    Vec3D u_vector_0 = {0.,0.,0.}; //v/m u bei 0
    Vec3D u_vector_i = u_vector_0; //v/m
    double normu = sqrt(innerprod(u_vector_i,u_vector_i));

   
    Vec3D B = {0.,0.,0.}; //im Test bei 1e-4
    
    
//für B const ungleich 0:
    double normB = sqrt(innerprod(B,B));
    double r = (m*normu)/(q*normB);//((1.13594e-06)/2;)
    //double t = (2.*M_PI*m/(q*normB))/1.; //für B=/=0
    int rounds = 1;
    //int n = 360.*rounds; //Anzahl Zeitschritte ohne Einheit
    //double delta_t = (t/n)*rounds; //Zeit gesamt in s //brauch ich eigentlich gar nicht
    //Vec3D x_0 = {0., 0., -r};  //?

//für E const ungleich 0:
    double delta_t = 0.01; //delta t in s
    int n = 101.;
    double t = n*delta_t;
    Vec3D x_0 = {0,0,0};
    Vec3D x_i = x_0;
    Vec3D E = {(m*c)/q,0.,0.}; //homogenes E-Feld
    Vec3D v_vector = scalmultip(E,(q/(m*c)));

    cout << "m*c/q = " << (m*c)/q << endl << endl;
    
    cout << "u(t=0) = " << vec_print(u_vector_0) << "            ";
    
    cout << "x(t=0) = " << vec_print(x_0) << endl << endl;
    

    for (int i=0; i < n; i++ )
    {

//dann Ortsupdate 1:
   
        Vec3D u_vector_i_by_c = scalmultip(u_vector_i,1./c);
        
        double gamma_i = sqrt(1.+ innerprod(u_vector_i_by_c,u_vector_i_by_c));
    
        Vec3D x_first_half = vec_add( x_i, (scalmultip(u_vector_i,delta_t/(2.*gamma_i))));
        
            //cout << "x(i+1/2) = " << vec_print(x_first_half) << endl;
        
//dann u:
        
        Vec3D u_minus = vec_add(u_vector_i, scalmultip(E,(q*delta_t)/(2.*m))); //ist hier u_vector_i
        
            //cout << "u_minus = " << vec_print(u_minus) << endl;

        double gamma_minus = sqrt(1. + (innerprod(u_minus, u_minus))/(pow(c,2.)));
        
            //cout << "gamma_minus = " << gamma_minus << endl;
    
        Vec3D tau = scalmultip( B, (q * delta_t ) / ( 2. * m ));
        
            //cout << "tau = " << vec_print(tau) << endl;
                
        double sigma = pow(gamma_minus,2.) - innerprod(tau,tau);
        
            //cout << "sigma = " << sigma << endl;
        
        double u_star = innerprod(u_minus,scalmultip(tau,1./c));
        
            //cout << "u_star = " << u_star << endl;
        
            //cout << "sigma^2 = " << pow(sigma,2.) << endl;
        
            //cout << "4tau^2 = " << 4. * (innerprod(tau,tau)) << endl;
                                     
            //cout << "u*^2 = " << pow(u_star,2.) << endl;
        
            //cout << "(sigma^2 + 4tau^2 + u*^2)^(1/2) = " << sqrt( pow(sigma,2.) + 4. * (innerprod(tau,tau)+ pow(u_star,2.))) << endl;

            //cout << "sigma + (sigma^2 + 4tau^2 + u*^2)^(1/2) = " << (sigma + sqrt( pow(sigma,2.) + 4. * (innerprod(tau,tau)+ pow(u_star,2.)))) << endl;
        
        double gamma_plus = sqrt(( sigma + sqrt( pow(sigma,2.) + 4. * (innerprod(tau,tau)+ pow(u_star,2.))))/2.);
        
            //cout << "gamma_plus = " << gamma_plus << endl;
        
        Vec3D t_vector = scalmultip(tau,1./gamma_plus);
        
            //cout << "t_vector = " << vec_print(t_vector) << endl;
        
        double s = 1./(1.+ innerprod(t_vector,t_vector));
        
            //cout << "s = " << s << endl;
        
        Vec3D u_plus = scalmultip(vec_add( vec_add( u_minus , scalmultip( t_vector , innerprod( u_minus , t_vector))), crossprod(u_minus,t_vector)), s);
        
            //cout << "u_plus = " << vec_print(u_plus) << endl;
        
// und jetzt u_vector_i überschreiben
        
        u_vector_i = vec_add(vec_add(u_plus , scalmultip( E, ( q*delta_t) / (2.*m) ) ), crossprod(u_plus,t_vector));
        
//dann ortsupdate 2:

        u_vector_i_by_c = scalmultip(u_vector_i,1./c);
        
        double gamma_final = sqrt(1.+ innerprod(u_vector_i_by_c,u_vector_i_by_c));
        
            //cout << "gamma_final = " << gamma_final << endl;
        
        x_i = vec_add(x_first_half , scalmultip(u_vector_i,delta_t/(2.*gamma_final)));
        //}
        
// jetzt der output:
        if ( i % 10 == 0 )
            
        {
            
        //cout << "Nach der " << (i+360)/360 << ".Runde:" << endl;
            
        cout << "u(t=" << (i+10)*delta_t << ") = " << vec_print(u_vector_i) << "            ";
        
        cout << "x(t=" << (i+10)*delta_t << ") = " << vec_print(x_i) << "            ";
        
        //für E=const ungleich 0:
        //Vec3D x_expected = scalmultip((scalmultip(v_vector, 1./(innerprod(v_vector,v_vector)))),c*(sqrt(1.+(innerprod(v_vector,v_vector))*(pow(((i+1.)*delta_t),2.)))-1.));
            
        //cout << "x_expected = " << vec_print(x_expected) << endl << endl;
        
        //normu = sqrt(innerprod(u_vector_i,u_vector_i));
        
            //cout << "|u| = " << normu << "                ";
        
        //double xnorm = sqrt(innerprod(x_i,x_i));
        
        //cout << "|x| = " << xnorm << endl << endl;
            
        }
    }
    
}
