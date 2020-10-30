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
#include "forces_bond_undamped.h"
#include "forces_bond_kelvin.h"
#include "forces_bond_maxwell.h"

int main()
{
//    for (int it=2; it<20; it++)
//    {

    double const rho = 7000.;
    double const d = 0.001;
    double const E = 2.7e6;
    double const nu = 0.4;
    double const dt = 0.0000001;
    double const eps = 1.e-6;
    double const lambda = 1.1;
    double const gamma = 1.;
    double const mun = 1.e-4;
    double const mus = 1.e-4;
    int bonds = 0;
    int count = 0;
    int const n = 10000000;
    int const n_part = 2;
    double omega;

//    omega = it;

    std::ofstream out_file;
    out_file.open("test_out/z-displacement_max_test.csv");
    out_file << "timestep,ID,x,y,z,u,v,w,fx,fy,fz,omx,omy,omz,tx,ty,tz\n";
    int const n_out = 10000;
    bool flag = false;

    Part vec_part[n_part];

    Bond vec_bond;

    /*Initial particle configuration*/

    double const m = rho*(4./3.)*3.14159265*pow((d/2.),3.);

    double pos1[3], pos2[3];

    pos1[0] = 0.;
    pos1[1] = 0.;
    pos1[2] = 0.;

    pos2[0] = 0.;
    pos2[1] = 0.;
    pos2[2] = d+eps;

    vec_part[0].setPos(pos1);
    vec_part[0].setPosIn(pos1);
    vec_part[1].setPos(pos2);
    vec_part[1].setPosIn(pos2);

    for (int i=0; i<n_part; i++)
    {
        vec_part[i].setMass(m);
        vec_part[i].setDiameter(d);
    }


//    for (int i=0; i<n_part; i++)
//    {
//        std::cout << "particle:" << i << "  " << "X:" << vec_part[i].getPos()[0] << ","
//               << vec_part[i].getPos()[1] << "," << vec_part[i].getPos()[2] <<"\n" ;
//    }


    /*Bond Initialization*/

    double *posi, *posj, radi, radj, radi_cr, radj_cr, dist;
    int IDS[2];

    posi = vec_part[0].getPos();
    radi = vec_part[0].getDiameter()/2.;
    radi_cr = radi*lambda;

    posj = vec_part[1].getPos();
    radj = vec_part[1].getDiameter()/2.;
    radj_cr = radj*lambda;

    dist = distance(posi, posj);

    if (dist < (radi_cr+radj_cr))
    {
        double radius;
        double const kn = E/(radi+radj);
        double const ks = kn/2.5;
        double inertia[2];
//        double const kn = sqrt(2.)*E*radi/(2.*(1.-2.*nu));
//        double const ks = kn/2.5;

        radius = lambda*std::min(radi,radj);

        IDS[0] = 0;
        IDS[1] = 1;

        inertia[0] = 0.25*3.14159265*pow(radius,4.);
        inertia[1] = 2.*inertia[0];

        vec_bond.setRad(radius);
        vec_bond.setStiffNorm(kn);
        vec_bond.setStiffShear(ks);
        vec_bond.setIds(IDS);
        vec_bond.setInertia(inertia);
    }

    double const f_ext[3] = {0., 0., 0.};

    /*Starts the dynamics*/

    for (int i=0; i<n; i++)
    {
        /*first timestep, initialization*/

        if(i==0)
        {
            double *x, *f, *v, *t, *om;

            x = vec_part[1].getPos();
            v = vec_part[1].getVel();
            f = vec_part[1].getForce();
            om = vec_part[1].getRot();
            t = vec_part[1].getTorque();

//            x[0] += 1.e-5;
            x[2] += 1.e-6;
//            v[0] = 0.001;
//            f[0] = 0.00001;
//            v[2] = 0.001;
//            om[0] = -0.0001;
//            om[1] = 0.001;

        }

      /*New positions with integration of EOM due to interaction and external forces*/

        else
        {
            for (int j=1; j<n_part; j++)
            {
                double *x, *v, *f, *om, *t;
                double *v_old, *f_old, *t_old;
                double v_act[3];
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

//                if (j==0)
//                {
                    if (i==1)
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

              /*Decomment for periodic boundary condition*/
//            else
//            {
             // double const omega = 700.;
//
//                double *x, *xin, *v, *om, v_old[3];
//
//                omega = 200.;
//                x = vec_part[1].getPos();
//                xin = vec_part[1].getPosIn();
//                v = vec_part[1].getVel();
//                om = vec_part[1].getRot();
//
//                v_old[0] = v[0];
//                v_old[1] = v[1];
//                v_old[2] = v[2];
//
//                vec_part[1].setVelOld(v_old);
//
//                x[0] = xin[0] + 1.e-6*sin(omega*dt*i);
//
//                v[0] = 1.e-6*omega*cos(omega*dt*i);
//                v[1] = 0.;
//                v[2] = 0.;
//
//                om[0] = 0.;
//                om[1] = 0.;
//                om[2] = 0.;
//            }
//            }
        }
        /*Forces on bonds*/

        Bond *bondino;

        bondino = &vec_bond;
        forces_bond_maxwell(bondino, vec_part[0], vec_part[1], dt, mun, mus, i);

        if (i==0)
            std::cout << "normal force: " << vec_bond.getForce_n()[0] << "," << vec_bond.getForce_n()[1] << "," << vec_bond.getForce_n()[2] << "\n"
                      << "shear force: " << vec_bond.getForce_s()[0] << "," << vec_bond.getForce_s()[1] << "," << vec_bond.getForce_s()[2] << "\n"
                      << "shear torque: " << vec_bond.getTorque_s()[0] << "," << vec_bond.getTorque_s()[1] << "," << vec_bond.getTorque_s()[2] << "\n";

        /*Forces on particles*/

        for (int j=0; j<n_part; j++)
        {
            double F_part[3], T_part[3];

            /*Bonds contribution*/

            int *IDS;

            IDS = vec_bond.getIds();

            if (j == IDS[0])
            {
                F_part[0] = -(vec_bond.getForce_n()[0] + vec_bond.getForce_s()[0]);
                F_part[1] = -(vec_bond.getForce_n()[1] + vec_bond.getForce_s()[1]);
                F_part[2] = -(vec_bond.getForce_n()[2] + vec_bond.getForce_s()[2]);

                T_part[0] = -(vec_bond.getTorque_s()[0] + vec_bond.getTorque_n()[0] + vec_bond.getTorque_t()[0]);
                T_part[1] = -(vec_bond.getTorque_s()[1] + vec_bond.getTorque_n()[1] + vec_bond.getTorque_t()[1]);
                T_part[2] = -(vec_bond.getTorque_s()[2] + vec_bond.getTorque_n()[2] + vec_bond.getTorque_t()[2]);
            }

            else if (j == IDS[1])
            {
                F_part[0] = (vec_bond.getForce_n()[0] + vec_bond.getForce_s()[0]);
                F_part[1] = (vec_bond.getForce_n()[1] + vec_bond.getForce_s()[1]);
                F_part[2] = (vec_bond.getForce_n()[2] + vec_bond.getForce_s()[2]);

                T_part[0] = vec_bond.getTorque_s()[0] + vec_bond.getTorque_n()[0] + vec_bond.getTorque_t()[0];
                T_part[1] = vec_bond.getTorque_s()[1] + vec_bond.getTorque_n()[1] + vec_bond.getTorque_t()[1];
                T_part[2] = vec_bond.getTorque_s()[2] + vec_bond.getTorque_n()[2] + vec_bond.getTorque_t()[2];
            }

            /*External forces contribution*/

            //if(i==0 && j==1)
            //{
            //    double const omega = 110.;

            //    F_part[0] += f_ext[0];
            //    F_part[1] += f_ext[1];
            //    F_part[2] += f_ext[2];
            //}

            vec_part[j].setForce(F_part);
            vec_part[j].setTorque(T_part);
        }

        for(int j=0; j<n_part; j++)
        {
            if (i == 100 || i == 200)
            {
                std::cout << j << ","
                          << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2]
                          << "," << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2]
                          << "," << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2]
                          << "," << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "\n";
            }
        }

        for(int j=0; j<n_part; j++)
        {
            if (i == count*n_out)
            {
                out_file << i*dt << "," << j << ","
                         << vec_part[j].getPos()[0] << "," << vec_part[j].getPos()[1] << "," << vec_part[j].getPos()[2] << ","
                         << vec_part[j].getVel()[0] << "," << vec_part[j].getVel()[1] << "," << vec_part[j].getVel()[2] << ","
                         << vec_part[j].getForce()[0] << "," << vec_part[j].getForce()[1] << "," << vec_part[j].getForce()[2] << ","
                         << vec_part[j].getRot()[0] << "," << vec_part[j].getRot()[1] << "," << vec_part[j].getRot()[2] << ","
                         << vec_part[j].getTorque()[0] << "," << vec_part[j].getTorque()[1] << "," << vec_part[j].getTorque()[2] << "\n";
                flag = true;
            }
        }

        if (flag)
        {
            count++;
            flag = false;
        }
    }

    for (int i=0; i<n_part; i++)
    {
        std::cout << i << "," << vec_part[i].getPos()[0] << ","
                              << vec_part[i].getPos()[1] << ","
                              << vec_part[i].getPos()[2] << "\n";
    }

    return 0;
}
