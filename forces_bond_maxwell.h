void forces_bond_maxwell(Bond *bond, Part parti, Part partj, double dt, double mun, double mus, int i)
{
    double *omi, *omj;
    double *posi, *posj, *veli, *velj, *veli_old, *velj_old;
    double *posi_in, *posj_in;
    double kn, ks, radius, *inertia, area;
    double *r_old, *r_old2, *vc_old, *fn, *fs, *ts, *tn, *tt, *ttr;

    /*Getting bond values*/

    double lambdan, lambdas;

    kn = bond -> getStiffNorm();
    ks = bond -> getStiffShear();

    r_old = bond -> getRold();
    r_old2 = bond -> getRold2();
    vc_old = bond -> getVcold();
//    fn_old = bond -> getFn();
//    fs_old = bond -> getFs();

    fn = bond -> getForce_n();
    fs = bond -> getForce_s();
    tn = bond -> getTorque_n();
    ts = bond -> getTorque_s();
//    tt = bond -> getTorque_t();
//    ttr = bond -> getTorque_tr();

    inertia = bond -> getInertia();

    radius = bond -> getRad();
    area = 3.14159265*radius*radius;

    lambdan = mun/(kn*area);
    lambdas = mus/(ks*area);

    /*Getting particles positions, velocities and torques*/

    posi = parti.getPos();
    posj = partj.getPos();

    posi_in = parti.getPosIn();
    posj_in = partj.getPosIn();

    veli = parti.getVel();
    velj = partj.getVel();

    veli_old = parti.getVelOld();
    velj_old = partj.getVelOld();

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

    double mod_rf, mod_r0, mod_rold, mod_rold2, beta, beta_old;
    double n[3], S[3], s[3], S_rot[3], s_rot[3];

    mod_rf = modulus(r_f);
    mod_r0 = modulus(r_0);
    mod_rold = modulus(r_old);
    mod_rold2 = modulus(r_old2);

    n[0] = r_f[0]/mod_rf;
    n[1] = r_f[1]/mod_rf;
    n[2] = r_f[2]/mod_rf;

    S[0] = r_f[1]*(r_f[0]*r_0[1]-r_0[0]*r_f[1]) - r_f[2]*(r_0[0]*r_f[2]-r_f[0]*r_0[2]);
    S[1] = r_f[2]*(r_f[1]*r_0[2]-r_0[1]*r_f[2]) - r_f[0]*(r_0[1]*r_f[0]-r_f[1]*r_0[0]);
    S[2] = r_f[0]*(r_f[2]*r_0[0]-r_0[2]*r_f[0]) - r_f[1]*(r_0[2]*r_f[1]-r_f[2]*r_0[1]);

//    S[0] = r_c[1]*(r_c[0]*r_0[1]-r_0[0]*r_c[1]) - r_c[2]*(r_0[0]*r_c[2]-r_c[0]*r_0[2]);
//    S[1] = r_c[2]*(r_c[1]*r_0[2]-r_0[1]*r_c[2]) - r_c[0]*(r_0[1]*r_c[0]-r_c[1]*r_0[0]);
//    S[2] = r_c[0]*(r_c[2]*r_0[0]-r_0[2]*r_c[0]) - r_c[1]*(r_0[2]*r_c[1]-r_c[2]*r_0[1]);

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

    double delta_theta_n, delta_theta_s, delta_theta_t, vc_n[3], vc_s[3];

    delta_theta_n = delta_theta[0]*n[0] + delta_theta[1]*n[1] + delta_theta[2]*n[2];
    delta_theta_s = delta_theta[0]*s[0] + delta_theta[1]*s[1] + delta_theta[2]*s[2];
//    delta_theta_t = delta_theta[0]*t[0] + delta_theta[1]*t[1] + delta_theta[2]*t[2];

    /*Computing forces and torques*/

//    double fn_rot, fs_rot;
//
//    fn_rot = fn[0]*n[0] + fn[1]*n[1] + fn[2]*n[2];
//    fs_rot = fs[0]*s[0] + fs[1]*s[1] + fs[2]*s[2];

    if (i==0)
    {
//        fn[0] = 0.;
//        fn[1] = 0.;
//        fn[2] = 0.;
//
//        fs[0] = 0.;
//        fs[1] = 0.;
//        fs[2] = 0.;

//        fn[0] = kn*area*(mod_rf - mod_r0)*n[0];
//        fn[1] = kn*area*(mod_rf - mod_r0)*n[1];
//        fn[2] = kn*area*(mod_rf - mod_r0)*n[2];
//
//        fs[0] = ks*area*mod_r0*sin(beta)*s[0];
//        fs[1] = ks*area*mod_r0*sin(beta)*s[1];
//        fs[2] = ks*area*mod_r0*sin(beta)*s[2];

        vc_n[0] = (v_c[0]*n[0] + v_c[1]*n[1] + v_c[2]*n[2])*n[0];
        vc_n[1] = (v_c[0]*n[0] + v_c[1]*n[1] + v_c[2]*n[2])*n[1];
        vc_n[2] = (v_c[0]*n[0] + v_c[1]*n[1] + v_c[2]*n[2])*n[2];

        vc_s[0] = v_c[0] - vc_n[0];
        vc_s[1] = v_c[1] - vc_n[1];
        vc_s[2] = v_c[2] - vc_n[2];

        fn[0] = fn[0]*(1. - dt/lambdan) + kn*area*vc_n[0]*dt;
        fn[1] = fn[1]*(1. - dt/lambdan) + kn*area*vc_n[1]*dt;
        fn[2] = fn[2]*(1. - dt/lambdan) + kn*area*vc_n[2]*dt;

        fs[0] = fs[0]*(1. - dt/lambdas) + ks*area*vc_s[0]*dt;
        fs[1] = fs[1]*(1. - dt/lambdas) + ks*area*vc_s[1]*dt;
        fs[2] = fs[2]*(1. - dt/lambdas) + ks*area*vc_s[2]*dt;

        std::cout << "vc_n: " << vc_n << "\n";
        std::cout << "vc_s: " << vc_s << "\n";
        std::cout << "fn: " << fn[0] << "," << fn[1] << "," << fn[2] << "\n";
        std::cout << "fs: " << fs[0] << "," << fs[1] << "," << fs[2] << "\n";
    }

    else if (i==1)
    {
        beta = acos((r_f[0]*r_0[0] + r_f[1]*r_0[1] + r_f[2]*r_0[2])/(mod_rf*mod_r0));
//        beta = acos((r_c[0]*r_0[0] + r_c[1]*r_0[1] + r_c[2]*r_0[2])/(mod_rc*mod_r0));

        if (std::isnan(beta))
            beta = 0.;

        fn[0] = kn*area*(mod_rf - mod_r0)*(2. - dt/lambdan)*n[0];
        fn[1] = kn*area*(mod_rf - mod_r0)*(2. - dt/lambdan)*n[1];
        fn[2] = kn*area*(mod_rf - mod_r0)*(2. - dt/lambdan)*n[2];

        fs[0] = ks*area*mod_r0*sin(beta)*(2. - dt/lambdas)*s[0];
        fs[1] = ks*area*mod_r0*sin(beta)*(2. - dt/lambdas)*s[1];
        fs[2] = ks*area*mod_r0*sin(beta)*(2. - dt/lambdas)*s[2];
    }

    else
    {
        vc_n[0] = (vc_old[0]*n[0] + vc_old[1]*n[1] + vc_old[2]*n[2])*n[0];
        vc_n[1] = (vc_old[0]*n[0] + vc_old[1]*n[1] + vc_old[2]*n[2])*n[1];
        vc_n[2] = (vc_old[0]*n[0] + vc_old[1]*n[1] + vc_old[2]*n[2])*n[2];

        vc_s[0] = vc_old[0] - vc_n[0];
        vc_s[1] = vc_old[1] - vc_n[1];
        vc_s[2] = vc_old[2] - vc_n[2];

        fn[0] = fn[0]*(1. - dt/lambdan) + kn*area*vc_n[0]*dt;
        fn[1] = fn[1]*(1. - dt/lambdan) + kn*area*vc_n[1]*dt;
        fn[2] = fn[2]*(1. - dt/lambdan) + kn*area*vc_n[2]*dt;

        fs[0] = fs[0]*(1. - dt/lambdas) + ks*area*vc_s[0]*dt;
        fs[1] = fs[1]*(1. - dt/lambdas) + ks*area*vc_s[1]*dt;
        fs[2] = fs[2]*(1. - dt/lambdas) + ks*area*vc_s[2]*dt;
    }

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

//    std::cout << "fn: " << fn[0] << "," << fn[1] << "," << fn[2] << "\n"
//              << "fs: " << fs[0] << "," << fs[1] << "," << fs[2] << "\n";
//              << "tt: " << ttr[0] << "," << ttr[1] << "," << ttr[2] << "\n";
//
//    std::cout << "dthetan: " << delta_theta_n << "\n"
//              << "dthetas: " << delta_theta_s << "\n"
//              << "dthetat: " << delta_theta_t << "\n";
//    std::cout << "n: " << n[0] << "," << n[1] << "," << n[2] << "\n"
//              << "s: " << s[0] << "," << s[1] << "," << s[2] << "\n"
//              << "t: " << t[0] << "," << t[1] << "," << t[2] << "\n";

//    fn_old = fn[0]*n[0] + fn[1]*n[1] + fn[2]*n[2];
//    fs_old = fs[0]*s[0] + fs[1]*s[1] + fs[2]*s[2];

    bond -> setRold(r_f);
    bond -> setRold2(r_old);
    bond -> setVcold(v_c);

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
