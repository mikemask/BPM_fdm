void forces_bond_undamped(Bond *bond, Part parti, Part partj, double dt, int i)
{
    double *omi, *omj, *veli, *velj;
    double *posi, *posj;
    double *posi_in, *posj_in;
    double kn, ks, radius, *inertia, area;
    double *fn, *fs, *ts, *tn, *ttr;

    /*Getting bond values*/

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

    fn = bond -> getForce_n();
    fs = bond -> getForce_s();
    tn = bond -> getTorque_n();
    ts = bond -> getTorque_s();
    ttr = bond -> getTorque_t();

    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;

    /*Getting particles positions, velocities and torques*/

    posi = parti.getPos();
    posj = partj.getPos();

    posi_in = parti.getPosIn();
    posj_in = partj.getPosIn();

    omi = parti.getRot();
    omj = partj.getRot();

    veli = parti.getVel();
    velj = partj.getVel();

    /*Computing relative motion and contact plane versors*/

    double r_f[3], r_0[3], r_c[3], v_rel[3], v_c[3];

    r_f[0] = posi[0] - posj[0];
    r_f[1] = posi[1] - posj[1];
    r_f[2] = posi[2] - posj[2];

    r_0[0] = posi_in[0] - posj_in[0];
    r_0[1] = posi_in[1] - posj_in[1];
    r_0[2] = posi_in[2] - posj_in[2];

    v_rel[0] = veli[0] - velj[0];
    v_rel[1] = veli[1] - velj[1];
    v_rel[2] = veli[2] - velj[2];

//    double x_cont[3];
//
//    x_cont[0] = 0.5*(posi[0] + posj[0]);
//    x_cont[1] = 0.5*(posi[1] + posj[1]);
//    x_cont[2] = 0.5*(posi[2] + posj[2]);
//
//    int e;
//
//    for (int i=0; i<3; i++)
//    {
//        double perm = 0.;
//
//        for (int j=0; j<3; j++)
//        {
//            for (int k=0; k<3; k++)
//            {
//                if (i==j || i==k || j==k || i==j==k)
//                    e = 0;
//                else if ((i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2))
//                    e = 1;
//                else
//                    e = -1;
//
//                perm += e*(omi[j]*(x_cont[k] - posi[k]) - omj[j]*(x_cont[k] - posj[k]));
//            }
//        }
//
////        r_c[i] = r_f[i] + perm*dt;
//        v_c[i] = v_rel[i] + perm;
////        std::cout << "perm: " << perm << "\n";
//        r_c[i] = v_c[i]*dt;
//    }

    double mod_r0, mod_rf, beta;
    double n0[3], n[3], S[3], s[3], T[3], t[3];

    mod_r0 = modulus(r_0);
    mod_rf = modulus(r_f);

    n0[0] = r_0[0]/mod_r0;
    n0[1] = r_0[1]/mod_r0;
    n0[2] = r_0[2]/mod_r0;

    n[0] = r_f[0]/mod_rf;
    n[1] = r_f[1]/mod_rf;
    n[2] = r_f[2]/mod_rf;

    beta = acos((r_0[0]*r_f[0]+r_0[1]*r_f[1]+r_0[2]*r_f[2])/(mod_rf*mod_r0));

    if (std::isnan(beta))
        beta = 0.;

    T[0] = r_0[1]*r_f[2] - r_0[2]*r_f[1];
    T[1] = r_0[2]*r_f[0] - r_0[0]*r_f[2];
    T[2] = r_0[0]*r_f[1] - r_0[1]*r_f[0];

    double mod_T = modulus(T);

    if (mod_T < 1.e-30)
    {
        t[0] = 0.;
        t[1] = 0.;
        t[2] = 0.;
    }

    else
    {
        t[0] = T[0]/mod_T;
        t[1] = T[1]/mod_T;
        t[2] = T[2]/mod_T;
    }

    S[0] = r_f[1]*(r_f[0]*r_0[1]-r_0[0]*r_f[1]) - r_f[2]*(r_0[0]*r_f[2]-r_f[0]*r_0[2]);
    S[1] = r_f[2]*(r_f[1]*r_0[2]-r_0[1]*r_f[2]) - r_f[0]*(r_0[1]*r_f[0]-r_f[1]*r_0[0]);
    S[2] = r_f[0]*(r_f[2]*r_0[0]-r_0[2]*r_f[0]) - r_f[1]*(r_0[2]*r_f[1]-r_f[2]*r_0[1]);

    double mod_S = modulus(S);

    if (mod_S < 1.e-30)
    {
        s[0] = 0.;
        s[1] = 0.;
        s[2] = 0.;
    }

    else
    {
        s[0] = S[0]/mod_S;
        s[1] = S[1]/mod_S;
        s[2] = S[2]/mod_S;
    }

    double delta_theta[3];

    delta_theta[0] = (omi[0] - omj[0])*dt;
    delta_theta[1] = (omi[1] - omj[1])*dt;
    delta_theta[2] = (omi[2] - omj[2])*dt;

    double delta_theta_n, delta_theta_s, delta_theta_t;

    delta_theta_n = delta_theta[0]*n[0] + delta_theta[1]*n[1] + delta_theta[2]*n[2];
    delta_theta_s = delta_theta[0]*s[0] + delta_theta[1]*s[1] + delta_theta[2]*s[2];
    delta_theta_t = delta_theta[0]*t[0] + delta_theta[1]*t[1] + delta_theta[2]*t[2];

    /*Computing forces and torques*/

    double tt[3];

    fn[0] = kn*area*(mod_rf - mod_r0)*n[0];
    fn[1] = kn*area*(mod_rf - mod_r0)*n[1];
    fn[2] = kn*area*(mod_rf - mod_r0)*n[2];

    fs[0] = ks*area*mod_r0*sin(beta)*s[0];
    fs[1] = ks*area*mod_r0*sin(beta)*s[1];
    fs[2] = ks*area*mod_r0*sin(beta)*s[2];

    double mod_fs = modulus(fs);

    tn[0] += ks*inertia[1]*delta_theta_n*n[0];
    tn[1] += ks*inertia[1]*delta_theta_n*n[1];
    tn[2] += ks*inertia[1]*delta_theta_n*n[2];

    ts[0] += kn*inertia[0]*delta_theta_s*s[0];
    ts[1] += kn*inertia[0]*delta_theta_s*s[1];
    ts[2] += kn*inertia[0]*delta_theta_s*s[2];

    ttr[0] += kn*inertia[0]*delta_theta_t*t[0];
    ttr[1] += kn*inertia[0]*delta_theta_t*t[1];
    ttr[2] += kn*inertia[0]*delta_theta_t*t[2];

    tt[0] = ttr[0] + 0.5*mod_rf*mod_fs*t[0];
    tt[1] = ttr[1] + 0.5*mod_rf*mod_fs*t[1];
    tt[2] = ttr[2] + 0.5*mod_rf*mod_fs*t[2];

    bond -> setTorque_t(tt);

    if (std::isnan(fs[0]) || std::isnan(fs[1]) || std::isnan(fs[2]))
    {
        std::cout << "error at iteration: " << i << "\n";
        std::cout << "normal force: " << fn[0] << "," << fn[1] << "," << fn[2] << "\n";
        std::cout << "shear versor: " << s[0] << "," << s[1] << "," << s[2] << "\n";
        std::cout << "beta: " << beta << "\n";
        exit(0);
    }
}
