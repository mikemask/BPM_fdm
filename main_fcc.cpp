#include <cmath>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include "modulus.h"
#include "part.h"
#include "distance.h"
#include "prop_bond.h"
#include "forces_bond_kelvin.h"

int main()
{
    for (int it=1; it<=200; it++)
    {
        double const rho = 7000.;
        double const d = 1.e-3;
        double const E = 2.7e6;
        double const nu = 0.4;
        double const dt = 0.0000001;
        double const eps = 1.e-6;
        double const lambda = 1.5;
        double const gamma = 1.;
        double const pi = 3.14159265;
        int bonds = 0;
        int count = 0;
//        int inc = 2000000;
//        int cnt = 0;
//        int const n = 10000000;
        int n_lays = 3;
        int n_rows = 3;
        int n_cols = 3;
        int n_part = 0;
        double const a = 2.*(d+eps)/sqrt(2.);
        double omega = it*10;
        double const mun = 1.e-1;
        double const mus = 1.;

        double t_tot = 4.*pi*(1./omega);
        double n_double = t_tot/dt;
        int n = (int) n_double;
        int n_out = n/720;

        double const f_ext[3] = {0., 0., 0.};

        std::ofstream out_file;
        out_file.open("freq_sweep_kelvin/fc_cube/mun_1-1_mus_1_ampl_1-6/x-freqsweep_"+ std::to_string(it) +".csv");
//        out_file.open("fcc_out/x-periodic_fcc_kel_test.csv");
//        out_file.open("frequency_sweep_fcc/x-periodic_fcc_"+ std::to_string(it) +".csv");
        out_file << "timestep,ID,x,y,z,u,v,w,fx,fy,fz,omx,omy,omz,tx,ty,tz,omega\n";
//        int const n_out = 10000;
        bool flag = false;

        Part vec_part[100];

        Bond *vec_bond = new Bond[1000];

        /*Initial particle configuration*/

        double const m = rho*(4./3.)*3.14159265*pow((d/2.),3.);

        int id = 0;
        double pos[3]={};
        double iterd = floor(n_cols/2.);
        int iter = int(iterd);

        for (int i=0; i<n_lays; i++)
        {
            for (int j=0; j<n_rows; j++)
            {
                if ((i==0 || i%2==0) && (j==0 || j%2==0))
                {
                    for (int k=0; k<=iter; k++)
                    {
                        pos[0] = k*a;
                        pos[1] = j*0.5*a;
                        pos[2] = i*0.5*a;

                        vec_part[id].setMass(m);
                        vec_part[id].setDiameter(d);
                        vec_part[id].setId(id);

                        vec_part[id].setPos(pos);
                        vec_part[id].setPosIn(pos);
                        id++;
                    }
                }

                else if ((i==0 || i%2==0) && (j%2!=0))
                {
                    for (int k=1; k<=iter; k++)
                    {
                        pos[0] = k*0.5*a;
                        pos[1] = j*0.5*a;
                        pos[2] = i*0.5*a;

                        vec_part[id].setMass(m);
                        vec_part[id].setDiameter(d);
                        vec_part[id].setId(id);

                        vec_part[id].setPos(pos);
                        vec_part[id].setPosIn(pos);
                        id++;
                    }
                }

                else if ((i%2!=0) && (j==0 || j%2==0))
                {
                    for (int k=1; k<=iter; k++)
                    {
                        pos[0] = k*0.5*a;
                        pos[1] = j*0.5*a;
                        pos[2] = i*0.5*a;

                        vec_part[id].setMass(m);
                        vec_part[id].setDiameter(d);
                        vec_part[id].setId(id);

                        vec_part[id].setPos(pos);
                        vec_part[id].setPosIn(pos);
                        id++;
                    }
                }

                else if ((i%2!=0) && (j%2!=0))
                {
                    for (int k=0; k<=iter; k++)
                    {
                        pos[0] = k*a;
                        pos[1] = j*0.5*a;
                        pos[2] = i*0.5*a;

                        vec_part[id].setMass(m);
                        vec_part[id].setDiameter(d);
                        vec_part[id].setId(id);

                        vec_part[id].setPos(pos);
                        vec_part[id].setPosIn(pos);
                        id++;
                    }
                }
            }
        }

        n_part = id;
        for (int i=0; i<n_part; i++)
        {
            std::cout << "particle:" << i << " " << "X:" << vec_part[i].getPos()[0] << ","
                   << vec_part[i].getPos()[1] << "," << vec_part[i].getPos()[2] <<" "
                   << "mass:" << vec_part[i].getMass() << " " << "diameter:" << vec_part[i].getDiameter() << "\n";
        }

        /*Bond initialization*/

        int pp = 1;

        for (int i=0; i<n_part; i++)
        {
            for (int j=pp; j<n_part; j++)
            {
                if (i != j)
                {
                    double *posi, *posj, radi, radj, radi_cr, radj_cr, dist;
                    int IDS[2];

                    posi = vec_part[i].getPos();
                    radi = vec_part[i].getDiameter()/2.;
                    radi_cr = radi*lambda;

                    posj = vec_part[j].getPos();
                    radj = vec_part[j].getDiameter()/2.;
                    radj_cr = radj*lambda;

                    dist = distance(posi, posj);

                    if (dist < (radi_cr+radj_cr))
                    {
                        double radius, inertia[2];
                        double const kn = E/(radi+radj);
                        double const ks = kn/2.5;

                        radius = lambda*std::min(radi,radj);

                        IDS[0] = i;
                        IDS[1] = j;

                        inertia[0] = 0.25*3.14159265*pow(radius,4.);
                        inertia[1] = 2.*inertia[0];

                        vec_bond[bonds].setRad(radius);
                        vec_bond[bonds].setStiffNorm(kn);
                        vec_bond[bonds].setStiffShear(ks);
                        vec_bond[bonds].setIds(IDS);
                        vec_bond[bonds].setInertia(inertia);

                        bonds++;
                    }
                }
            }

            pp++;
        }


        /*Starts the dynamics*/

        for (int i=0; i<n; i++)
        {
            /*first timestep, initialization*/

//            if(i==0)
//            {
//                for (int j=9; j<n_part; j++)
//                {
//                    double *x, *v, *om, *f;
//
//                    x = vec_part[j].getPos();
//                    v = vec_part[j].getVel();
//                    om = vec_part[j].getRot();
////                    vec_part[j].setForce(f_ext);
//
//                    x[0] += 1.e-5;
//
////                    v[0] = 0.1;
////                    v[2] = 0.001;
////                    v[1] = 0.001;
////                    om[0] = 0.001;
//
//                }
//            }
//
//            else
//            {

            /*New positions with integration of EOM due to bond forces*/

                for (int j=5; j<n_part; j++)
                {
                    double *x, *v, *f, *om, *t;
                    double *v_old, *f_old, *t_old;
                    double x_act[3], v_act[3];
                    double m, d, inertia;

                    x = vec_part[j].getPos();
                    v = vec_part[j].getVel();
                    f = vec_part[j].getForce();
                    om = vec_part[j].getRot();
                    t = vec_part[j].getTorque();

                    v_old = vec_part[j].getVelOld();
                    f_old = vec_part[j].getForceOld();
                    t_old = vec_part[j].getTorqueOld();

                    m = vec_part[j].getMass();
                    d = vec_part[j].getDiameter();
                    inertia = (2./5.)*m*pow((d/2.),2.);

                    if (j<9)
                    {
                        if(i==1)
                        {
                            for (int k=0; k<3; k++)
                            {
                                v_act[k] = v[k];

                                x[k] += v[k]*dt;
                                v[k] += (dt/m)*f[k];
                                om[k] += (dt/inertia)*t[k];
                            }
                        }

                        else
                        {
                            for (int k=0; k<3; k++)
                            {
                                v_act[k] = v[k];

                                x[k] += 0.5*dt*(3.*v[k] - v_old[k]);
                                v[k] += 0.5*(dt/m)*(3.*f[k] - f_old[k]);
                                om[k] += 0.5*(dt/inertia)*(3.*t[k] - t_old[k]);
                            }
                        }

                        vec_part[j].setVelOld(v_act);
                        vec_part[j].setForceOld(f);
                        vec_part[j].setTorqueOld(t);
                    }

                  /*Decomment for periodic BC*/
                    else
                    {
                        double *xin;
                        double ampl = 1.e-6;

//                        if (i == inc*cnt)
//                        {
//                            ampl = 1.e-6 + cnt*5.e-7;
//                            cnt += 1;
//                        }

//                        omega = 200.;
                        xin = vec_part[j].getPosIn();

                        x[0] = xin[0] + ampl*sin(omega*dt*i);

                        v[0] = omega*ampl*cos(omega*dt*i);
                        v[1] = 0.;
                        v[2] = 0.;

                        om[0] = 0.;
                        om[1] = 0.;
                        om[2] = 0.;
                    }
                }
//            }

            /*Forces on bonds*/

            for (int j=0; j<bonds; j++)
            {
                int *IDS;
                Bond *bondino;

                IDS = vec_bond[j].getIds();
                bondino = &vec_bond[j];
                forces_bond_kelvin(bondino, vec_part[IDS[0]], vec_part[IDS[1]], dt, mun, mus, i);

            }

            /*Forces on particles*/

            for (int j=0; j<n_part; j++)
            {
                double F_part[3]={}, T_part[3]={};

                /*Bonds contribution*/

                for (int k=0; k<bonds; k++)
                {
                    int *IDS;

                    IDS = vec_bond[k].getIds();

                    if (j == IDS[0])
                    {
                        F_part[0] += -(vec_bond[k].getForce_n()[0] + vec_bond[k].getForce_s()[0]);
                        F_part[1] += -(vec_bond[k].getForce_n()[1] + vec_bond[k].getForce_s()[1]);
                        F_part[2] += -(vec_bond[k].getForce_n()[2] + vec_bond[k].getForce_s()[2]);

                        T_part[0] += -(vec_bond[k].getTorque_s()[0] + vec_bond[k].getTorque_t()[0] + vec_bond[k].getTorque_n()[0]);
                        T_part[1] += -(vec_bond[k].getTorque_s()[1] + vec_bond[k].getTorque_t()[1] + vec_bond[k].getTorque_n()[1]);
                        T_part[2] += -(vec_bond[k].getTorque_s()[2] + vec_bond[k].getTorque_t()[2] + vec_bond[k].getTorque_n()[2]);
                    }

                    else if (j == IDS[1])
                    {
                        F_part[0] += (vec_bond[k].getForce_n()[0] + vec_bond[k].getForce_s()[0]);
                        F_part[1] += (vec_bond[k].getForce_n()[1] + vec_bond[k].getForce_s()[1]);
                        F_part[2] += (vec_bond[k].getForce_n()[2] + vec_bond[k].getForce_s()[2]);

                        T_part[0] += (vec_bond[k].getTorque_s()[0] + vec_bond[k].getTorque_t()[0] + vec_bond[k].getTorque_n()[0]);
                        T_part[1] += (vec_bond[k].getTorque_s()[1] + vec_bond[k].getTorque_t()[1] + vec_bond[k].getTorque_n()[1]);
                        T_part[2] += (vec_bond[k].getTorque_s()[2] + vec_bond[k].getTorque_t()[2] + vec_bond[k].getTorque_n()[2]);
                    }

                }

                /*External forces contribution*/

//                if(j>=9 && i<4000000)
//                {
//                    double const omega = 600.;
//
//                    F_part[0] += f_ext[0];
//                    F_part[1] += f_ext[1];
//                    F_part[2] += f_ext[2];
//                }

                if (std::isnan(T_part[0]) || std::isnan(T_part[1]) || std::isnan(T_part[2]))
                {
                    std::cout << "error at iteration: " << i << "\n";
                    std::cout << "part: " << j << "\n"
                              << "f: " << F_part[0] << "," << F_part[1] << "," << F_part[2] << "\n"
                              << "t: " << T_part[0] << "," << T_part[1] << "," << T_part[2] << "\n";
                }
                vec_part[j].setForce(F_part);
                vec_part[j].setTorque(T_part);
            }

            for(int j=0; j<n_part; j++)
            {
                if (i == 120 || i == 121)
                {
                    std::cout << j << ","
                              << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2]
                              << "," << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2]
                              << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2] << "\n";
                }
            }

            for(int j=0; j<n_part; j++)
            {
                if (i == count*n_out)
                {
                    out_file << i*dt << "," << j << ","
                              << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2]
                              << "," << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2]
                              << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2]
                              << "," << vec_part[j].getRot()[0] << "," << vec_part[j].getRot()[1] << "," << vec_part[j].getRot()[2]
                              << "," << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "," << omega << "\n";
                    flag = true;
                }
            }

            if (flag)
            {
                count++;
                flag = false;
                std::cout << count << "\n";
            }
        }
//        }

        for (int i=0; i<n_part; i++)
        {
            std::cout << i << "," << vec_part[i].getPos()[0] << ","
                                  << vec_part[i].getPos()[1] << ","
                                  << vec_part[i].getPos()[2] << "\n";
        }

        delete [] vec_bond;
//    delete [] vec_part;
    }
    return 0;
}
