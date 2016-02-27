package ctmc.scheme.ionization;

import java.text.DecimalFormat;
import java.util.Random;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 *
 * A = re - rt (R2)
 * B = rp - rt (R1)
 * C = re - rp (R3)
 * C = A - B
 */
public class IonizationScheme {
    private Particle p;
    private Electron e;
    private Particle t;
    private Vector A;
    private Vector dA;
    private Vector B;
    private Vector dB;
    private Vector C;
    private Vector dC;
    private double Eb;
    private double b;
    private double KE;
    private double met;
    private double mep;

    public IonizationScheme() {

    }

    public String toString() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(10);
        df.setMinimumFractionDigits(10);
        StringBuilder sb = new StringBuilder();
        String str = "\nEb = " + df.format(Eb) + " b = " + df.format(b) + " KE = " + df.format(KE);
        sb.append(str);
        sb.append(p.toString());
        sb.append(e.toString());
        sb.append(t.toString());
        str = "A = " + A.toString();
        sb.append(str);
        str = "\nB = " + B.toString();
        sb.append(str);
        str = "\nC = " + C.toString();
        sb.append(str);
        str = "\n|A| = " + A.length();
        sb.append(str);
        str = "\n|B| = " + B.length();
        sb.append(str);
        str = "\n|C| = " + C.length();
        sb.append(str);
        return sb.toString();
    }

    public void setParameters(RecGenerator.Orbital orb, double KE, double b) {
        this.Eb = RecGenerator.getValue(orb);
        this.KE = KE;
        this.b = b;
    }

    public IonizationScheme(RecGenerator.Orbital orb, double KE, double b, Particle.ParticleName carbon_ion) {
        t = new Particle(Particle.ParticleName.WATER);
        p = new Particle(carbon_ion);
        e = new Electron();

        A = new Vector();
        B = new Vector();
        C = new Vector();

        dA = new Vector();
        dB = new Vector();
        dC = new Vector();

        this.setParameters(orb, KE, b);
        this.met = t.m*e.m/(t.m+e.m);
        this.mep = p.m*e.m/(p.m+e.m);
    }

    public double Eet()
    {
        return dA.square()*met*0.5 - t.Vc(A.length());
    }

    public double Eep()
    {
        return dC.square()*mep*0.5 - p.Vc(C.length());
    }

    public double theta() {
        return Math.acos(e.R.z / e.R.length());
    }

    public double calcAzimuthalAngle() {
        return Math.acos(A.z / A.length());
    }

    public void updateABC() {
        A = Vector.subtract(e.R, t.R);
        B = Vector.subtract(p.R, t.R);
        C = Vector.subtract(e.R, p.R);
    }

    public void updateDiffABC() {
        dA = Vector.subtract(e.V, t.V);
        dB = Vector.subtract(p.V, t.V);
        dC = Vector.subtract(e.V, p.V);
    }

    public void initialize(RecGenerator RGen, Random rand1, Random rand2, Random rand3) {
        // Tested OK (20-Jun-15)
        // ---------------------------- CM Frame ---------------------
        //=============Test value==============
                /*double rec = 0.7130;
                double teta = 0.2221;
                double phi = 5.7978;
                double neta = 5.0289;*/
        //=====================================
        //Random a value of rec
        double rec = RGen.nextRec();
        double pec = RGen.Pec(rec);

        //Random spherical coordinates
        double teta = rand1.nextDouble() * Constant.Pi;
        double phi = rand2.nextDouble() * 2.0 * Constant.Pi;
        double neta = rand3.nextDouble() * 2.0 * Constant.Pi;

        //Vector Rec
        //p.33 Thesis Liamsuwan
        Vector Rec = new Vector(Math.sin(teta)*Math.cos(phi), Math.sin(teta)*Math.sin(phi), Math.cos(teta));
        Rec.times(rec);

        //Relative velocity electron vs core
        //Eq.44 Thesis Liamsuwan
        Vector Vec = new Vector();
        Vec.x = -Math.sin(phi) * Math.cos(neta) - Math.cos(phi) * Math.cos(teta) * Math.sin(neta);
        Vec.y = Math.cos(phi) * Math.cos(neta) - Math.sin(phi) * Math.cos(teta) * Math.sin(neta);
        Vec.z = Math.sin(teta) * Math.sin(neta);
        Vec.times(pec/met);

        //------------------------------- LAB Frame ----------------------------
        //============================= INITIAL STATE ==========================
        double Met = (e.m + t.m);
        // >> Electron
        double xe0 = (t.m * Rec.x)/Met;
        double ye0 = (t.m * Rec.y)/Met;
        double ze0 = (t.m * Rec.z)/Met;

        double vex0 = (t.m * Vec.x)/Met;
        double vey0 = (t.m * Vec.y)/Met;
        double vez0 = (t.m * Vec.z)/Met;
        // >> Target
        double xt0 = -(e.m * Rec.x)/Met;
        double yt0 = -(e.m * Rec.y)/Met;
        double zt0 = -(e.m * Rec.z)/Met;

        double vtx0 = -(e.m * Vec.x)/Met;
        double vty0 = -(e.m * Vec.y)/Met;
        double vtz0 = -(e.m * Vec.z)/Met;
        // >> Projectile
        double xp0 = 0;
        double yp0 = b;
        double zp0 = -50;

        double vpx0 = 0;
        double vpy0 = 0;
        double vpz0 = Math.sqrt((2*KE)/p.m);

        //================================== SETTING VALUES ===========================
        p.R.set(xp0, yp0, zp0);
        e.R.set(xe0, ye0, ze0);
        t.R.set(xt0, yt0, zt0);

        p.V.set(vpx0, vpy0, vpz0);
        e.V.set(vex0, vey0, vez0);
        t.V.set(vtx0, vty0, vtz0);

        updateABC();
        updateDiffABC();
    }

    public Vector calcRelativeForceB(Vector _A, Vector _B, Vector _C) {
        //
        // * Optimization
        double A = _A.length();
        double B = _B.length();
        double C = _C.length();

        double termA = ((e.Z/(A*A)) * (-t.dZc(A) + t.Vc(A)));
        double termB = ((1.0/(B*B)) * (-t.dZc(B)*p.Zc(B) - t.Zc(B)*p.dZc(B) + p.Zc(B)*t.Zc(B)/B));
        double termC = ((e.Z/(C*C)) * (-p.dZc(C) + p.Vc(C)));

        Vector force;
        // FB
        force = Vector.sum(Vector.times(_B, termB * (1 / p.m + 1 / t.m)), Vector.subtract(Vector.times(_C, termC / p.m), Vector.times(_A, termA / t.m)));

        return force;
    }

    public Vector[] calcRelativeForces(Vector _A, Vector _B, Vector _C) {
        //
        // * Optimization
        double A = _A.length();
        double B = _B.length();
        double C = _C.length();

        double termA = ((e.Z/(A*A)) * (-t.dZc(A) + t.Vc(A)));
        double termB = ((1.0/(B*B)) * (-t.dZc(B)*p.Zc(B) - t.Zc(B)*p.dZc(B) + p.Zc(B)*t.Zc(B)/B));
        double termC = ((e.Z/(C*C)) * (-p.dZc(C) + p.Vc(C)));

        Vector[] forces = new Vector[2];
        // FA
        forces[0] = Vector.sum(Vector.times(_A, termA * (1 / e.m + 1 / t.m)), Vector.subtract(Vector.times(_C, termC / e.m), Vector.times(_B, termB / t.m)));
        // FB
        forces[1] = Vector.sum(Vector.times(_B, termB * (1 / p.m + 1 / t.m)), Vector.subtract(Vector.times(_C, termC / p.m), Vector.times(_A, termA / t.m)));

        return forces;
    }

    public void calcOneStepRK4B(double h)
    {
        Vector force;
        Vector[] KB = new Vector[4];
        Vector[] LB = new Vector[4];

        Vector _A = new Vector(this.A);
        Vector _B = new Vector(this.B);
        Vector _C = new Vector(this.C);

        //**************** 1ST ORDER ***************
        KB[0] = Vector.times(dB, h);
        force = this.calcRelativeForceB(_A, _B, _C);
        LB[0] = Vector.times(force, h);

        _B = Vector.sum(B, Vector.times(KB[0], 0.5));
        _C = Vector.subtract(_A, _B);

        //**************** 2ND ORDER ***************
        KB[1] = Vector.times(Vector.sum(dB, Vector.times(LB[0], 0.5)), h);
        force = this.calcRelativeForceB(_A, _B, _C);
        LB[1] = Vector.times(force, h);

        _B = Vector.sum(B, Vector.times(KB[1], 0.5));
        _C = Vector.subtract(_A, _B);

        //**************** 3RD ORDER ***************
        KB[2] = Vector.times(Vector.sum(dB, Vector.times(LB[1], 0.5)), h);
        force = this.calcRelativeForceB(_A, _B, _C);
        LB[2] = Vector.times(force, h);

        _B = Vector.sum(B, KB[2]);
        _C = Vector.subtract(_A, _B);

        //**************** 4TH ORDER ***************
        KB[3] = Vector.times(Vector.sum(dB, LB[2]), h);
        force = this.calcRelativeForceB(_A, _B, _C);
        LB[3] = Vector.times(force, h);

        B = Vector.sum(B, Vector.times(Vector.sum(Vector.sum(KB[0], KB[3]), Vector.times(Vector.sum(KB[1], KB[2]), 2.0)), 1.0 / 6.0));
        dB = Vector.sum(dB, Vector.times(Vector.sum(Vector.sum(LB[0], LB[3]), Vector.times(Vector.sum(LB[1], LB[2]), 2.0)), 1.0 / 6.0));

        C = Vector.subtract(A, B);
        dC = Vector.subtract(dA, dB);
    }

    public void calcOneStepRK4ABC(double h)
    {
        Vector[] forces;
        Vector[] KA = new Vector[4];
        Vector[] LA = new Vector[4];
        Vector[] KB = new Vector[4];
        Vector[] LB = new Vector[4];

        Vector _A = new Vector(this.A);
        Vector _B = new Vector(this.B);
        Vector _C = new Vector(this.C);

        //**************** 1ST ORDER ***************
        KA[0] = Vector.times(dA, h);
        KB[0] = Vector.times(dB, h);

        forces = this.calcRelativeForces(_A, _B, _C);
        LA[0] = Vector.times(forces[0], h);
        LB[0] = Vector.times(forces[1], h);

        _A = Vector.sum(A, Vector.times(KA[0], 0.5));
        _B = Vector.sum(B, Vector.times(KB[0], 0.5));
        _C = Vector.subtract(_A, _B);

        //**************** 2ND ORDER ***************
        KA[1] = Vector.times(Vector.sum(dA, Vector.times(LA[0], 0.5)), h);
        KB[1] = Vector.times(Vector.sum(dB, Vector.times(LB[0], 0.5)), h);

        forces = this.calcRelativeForces(_A, _B, _C);
        LA[1] = Vector.times(forces[0], h);
        LB[1] = Vector.times(forces[1], h);

        _A = Vector.sum(A, Vector.times(KA[1], 0.5));
        _B = Vector.sum(B, Vector.times(KB[1], 0.5));
        _C = Vector.subtract(_A, _B);

        //**************** 3RD ORDER ***************
        KA[2] = Vector.times(Vector.sum(dA, Vector.times(LA[1], 0.5)), h);
        KB[2] = Vector.times(Vector.sum(dB, Vector.times(LB[1], 0.5)), h);

        forces = this.calcRelativeForces(_A, _B, _C);
        LA[2] = Vector.times(forces[0], h);
        LB[2] = Vector.times(forces[1], h);

        _A = Vector.sum(A, KA[2]);
        _B = Vector.sum(B, KB[2]);
        _C = Vector.subtract(_A, _B);

        //**************** 4TH ORDER ***************
        KA[3] = Vector.times(Vector.sum(dA, LA[2]), h);
        KB[3] = Vector.times(Vector.sum(dB, LB[2]), h);

        forces = this.calcRelativeForces(_A, _B, _C);
        LA[3] = Vector.times(forces[0], h);
        LB[3] = Vector.times(forces[1], h);

        A = Vector.sum(A, Vector.times(Vector.sum(Vector.sum(KA[0], KA[3]), Vector.times(Vector.sum(KA[1], KA[2]), 2.0)), 1.0 / 6.0));
        dA = Vector.sum(dA, Vector.times(Vector.sum(Vector.sum(LA[0], LA[3]), Vector.times(Vector.sum(LA[1], LA[2]), 2.0)), 1.0 / 6.0));

        B = Vector.sum(B, Vector.times(Vector.sum(Vector.sum(KB[0], KB[3]), Vector.times(Vector.sum(KB[1], KB[2]), 2.0)), 1.0 / 6.0));
        dB = Vector.sum(dB, Vector.times(Vector.sum(Vector.sum(LB[0], LB[3]), Vector.times(Vector.sum(LB[1], LB[2]), 2.0)), 1.0 / 6.0));

        C = Vector.subtract(A, B);
        dC = Vector.subtract(dA, dB);
    }

    public double RungeKuttaABC() {
        // 10-Jul-15 Fix NaN by moving stopping condition and a new one
        // 30-Jun-15 Optimized
        // Test OK same result with matlab(very slightly different) 23-Jun-15
        // 24-Jun-15 Add a break condition when this functions run too long, specifically over 1 minute.
        double Rcritical;
        double limitDis = 1.0e-3;
        double limitDisPrj = 3.0e-1;
        double vB;
        double vA;
        boolean I = true;

        double time = 0;
        double h = 0.001;

        Rcritical = Math.pow(Math.sqrt(p.Zc(C.length())) + Math.sqrt(t.Zc(A.length())), 2)/Eb;

        while (time < 1000) {
            // Stopping conditions
            if (B.length() > 51 || A.length() > 51) {
                break;
            }

            time += h;
            if (I) {

                if (B.length() > Rcritical) {
                    I = true;
                } else {
                    I = false;
                    continue;
                }

                this.calcOneStepRK4B(h);

                vB = dB.length();
                if (vB*h < limitDisPrj) {
                    h = limitDisPrj/vB;
                }
            } else {
                vA = dA.length();
                h = limitDis/vA;
                this.calcOneStepRK4ABC(h);

                vA = dA.length();
                if (vA <= 15)
                    limitDis = 1.0e-2;
                else
                    limitDis = 1.0e-3;
            }
        }

        System.out.println("time " + time);
        if (time > 500) {
            System.out.println("b " + this.b + " Eb " + this.Eb + " KE " + this.KE);
        }

        return time;
    }

    public Outcome getOutcome() {
        double Eet = this.Eet();
        double Eep = this.Eep();
        int type = 5;
        double angle = 0;
        //Type of interation
        if (Eet > 0 && Eep > 0) {
            // Ionization
            type = 1;
            Eet *= 27.2116; // -> SI unit
            //angle = this.theta();
            angle = this.calcAzimuthalAngle();
        } else if (Eet > 0 && Eep < 0) {
            // Electron capture
            type = 2;
            Eet = 0;
        }
        else if (Eet <= 0 && Eep > 0) // Remained bound to target
        {
            type = 3;
            Eet = 0;
        }
        else // Non-physical
        {
            System.out.println("b " + this.b + " Eb " + this.Eb + " KE " + this.KE + " Eet " + Eet + " Eep " + Eep);
            type = 5;
            Eet = 0;
        }

        Outcome outcome;
        outcome = new Outcome(type, Eet, angle);
        return outcome;
    }
}
