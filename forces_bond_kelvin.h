void forces_bond_kelvin(Bond *bond, Part parti, Part partj, double dt, double mun, double mus, int i)
{
    double *omi, *omj;
    double *posi, *posj, *veli, *velj;
    double *posi_in, *posj_in;
    double kn, ks, radius, *inertia, area;
    double *fn, *fs, *ts, *tn, *tt, *ttr;

    /*Getting bond values*/

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

    fn = bond -> getForce_n();
    fs = bond -> getForce_s();
    tn = bond -> getTorque_n();
    ts = bond -> getTorque_s();
//    tt = bond -> getTorque_t();
//    ttr = bond -> getTorque_tr();

    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;

    /*Getting particles positions, velocities and torques*/

    posi = parti.getPos();
    posj = partj.getPos();

    posi_in = parti.getPosIn();
    posj_in = partj.getPosIn();

    veli = parti.getVel();
    velj = partj.getVel();

    omi = parti.getRot();
    omj = partj.getRot();

    /*Computing relative motion and contact plane versors*/

    double r_f[3], r_0[3], v_rel[3], v_c[3];

    r_f[0] = posi[0] - posj[0];
    r_f[1] = posi[1] - posj[1];
    r_f[2] = posi[2] - posj[2];

    r_0[0] = posi_in[0] - posj_in[0];
    r_0[1] = posi_in[1] - posj_in[1];
    r_0[2] = posi_in[2] - posj_in[2];

    v_rel[0] = veli[0] - velj[0];
    v_rel[1] = veli[1] - velj[1];
    v_rel[2] = veli[2] - velj[2];

    double x_cont[3];

    x_cont[0] = 0.5*(posi[0] + posj[0]);
    x_cont[1] = 0.5*(posi[1] + posj[1]);
    x_cont[2] = 0.5*(posi[2] + posj[2]);

    int e;

    for (int i=0; i<3; i++)
    {
        double perm = 0.;

        for (int j=0; j<3; j++)
        {
            for (int k=0; k<3; k++)
            {
                if (i==j || i==k || j==k || i==j==k)
                    e = 0;
                else if ((i==1 && j==2 && k==3) || (i==2 && j==3 && k==1) || (i==3 && j==1 && k==2))
                    e = 1;
                else
                    e = -1;

                perm += e*(omi[j]*(x_cont[k] - posi[k]) - omj[j]*(x_cont[k] - posj[k]));
            }
        }
        v_c[i] = v_rel[i] + perm;
    }

    double mod_rf, mod_r0, beta;
    double n[3], S[3], s[3];

    mod_rf = modulus(r_f);
    mod_r0 = modulus(r_0);

    n[0] = r_f[0]/mod_rf;
    n[1] = r_f[1]/mod_rf;
    n[2] = r_f[2]/mod_rf;

    beta = acos((r_0[0]*r_f[0]+r_0[1]*r_f[1]+r_0[2]*r_f[2])/(mod_rf*mod_r0));

    if (std::isnan(beta))
        beta = 0.;

//    T[0] = r_0[1]*r_f[2] - r_0[2]*r_f[1];
//    T[1] = r_0[2]*r_f[0] - r_0[0]*r_f[2];
//    T[2] = r_0[0]*r_f[1] - r_0[1]*r_f[0];
//
//    T[0] = r_0[1]*r_c[2] - r_0[2]*r_c[1];
//    T[1] = r_0[2]*r_c[0] - r_0[0]*r_c[2];
//    T[2] = r_0[0]*r_c[1] - r_0[1]*r_c[0];
//
//    double mod_T = modulus(T);
//
//    if (mod_T < 1.e-30)
//    {
//        t[0] = 0.;
//        t[1] = 0.;
//        t[2] = 0.;
//    }
//
//    else
//    {
//        t[0] = T[0]/mod_T;
//        t[1] = T[1]/mod_T;
//        t[2] = T[2]/mod_T;
//    }

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

    double delta_theta_n, delta_theta_s, delta_theta_t, vc_n, vc_s;

    delta_theta_n = delta_theta[0]*n[0] + delta_theta[1]*n[1] + delta_theta[2]*n[2];
    delta_theta_s = delta_theta[0]*s[0] + delta_theta[1]*s[1] + delta_theta[2]*s[2];
//    delta_theta_t = delta_theta[0]*t[0] + delta_theta[1]*t[1] + delta_theta[2]*t[2];

    vc_n = v_c[0]*n[0] + v_c[1]*n[1] + v_c[2]*n[2];
    vc_s = v_c[0]*s[0] + v_c[1]*s[1] + v_c[2]*s[2];

    /*Computing forces and torques*/

    fn[0] = (area*(mod_rf - mod_r0)*kn + mun*vc_n)*n[0];
    fn[1] = (area*(mod_rf - mod_r0)*kn + mun*vc_n)*n[1];
    fn[2] = (area*(mod_rf - mod_r0)*kn + mun*vc_n)*n[2];

    fs[0] = (ks*area*mod_r0*(beta) + mus*vc_s)*s[0];
    fs[1] = (ks*area*mod_r0*(beta) + mus*vc_s)*s[1];
    fs[2] = (ks*area*mod_r0*(beta) + mus*vc_s)*s[2];

    double mod_fs = modulus(fs), tts[3];

    tn[0] += ks*inertia[1]*delta_theta_n*n[0];
    tn[1] += ks*inertia[1]*delta_theta_n*n[1];
    tn[2] += ks*inertia[1]*delta_theta_n*n[2];

    ts[0] += kn*inertia[0]*delta_theta_s*s[0];
    ts[1] += kn*inertia[0]*delta_theta_s*s[1];
    ts[2] += kn*inertia[0]*delta_theta_s*s[2];

//    tts[0] = 0.5*mod_rf*mod_fs*t[0];
//    tts[1] = 0.5*mod_rf*mod_fs*t[1];
//    tts[2] = 0.5*mod_rf*mod_fs*t[2];
//
//    ttr[0] += kn*inertia[0]*delta_theta_t*t[0];
//    ttr[1] += kn*inertia[0]*delta_theta_t*t[1];
//    ttr[2] += kn*inertia[0]*delta_theta_t*t[2];
//
//    tt[0] = ttr[0] + tts[0];
//    tt[1] = ttr[1] + tts[1];
//    tt[2] = ttr[2] + tts[2];

//    std::cout << "fvn: " << fvn[0] << "," << fvn[1] << "," << fvn[2] << "\n"
//              << "fvs: " << fvs[0] << "," << fvs[1] << "," << fvs[2] << "\n";
//              << "tt: " << ttr[0] << "," << ttr[1] << "," << ttr[2] << "\n";
//
//    std::cout << "dthetan: " << delta_theta_n << "\n"
//              << "dthetas: " << delta_theta_s << "\n"
//              << "dthetat: " << delta_theta_t << "\n";
//    std::cout << "n: " << n[0] << "," << n[1] << "," << n[2] << "\n"
//              << "s: " << s[0] << "," << s[1] << "," << s[2] << "\n"
//              << "t: " << t[0] << "," << t[1] << "," << t[2] << "\n";

    if (std::isnan(fs[0]) || std::isnan(fs[1]) || std::isnan(fs[2]))
    {
        std::cout << "error at iteration: " << i << "\n";
        std::cout << "normal force: " << fn[0] << "," << fn[1] << "," << fn[2] << "\n";
        std::cout << "shear versor: " << s[0] << "," << s[1] << "," << s[2] << "\n";
        std::cout << "beta: " << beta << "\n";
        exit(0);
    }
//    bond -> setForce_n(fn);
//    bond -> setForce_s(fs);
//    bond -> setTorque_st(tst);
//    bond -> setTorque_n(tn);
//    bond -> setTorque_sr(tsr);
//    bond -> setTorque_tr(ttr);
}
