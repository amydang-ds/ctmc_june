package ctmc.scheme.ionization;

import java.text.DecimalFormat;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 * Particle model:
 *                  mass m
 *                  atomic number Z
 *                  remained electrons N
 *                  position vector R
 *                  velocity vector V
 *                  effective charge parameters: etat, netat
 *
 */
public class Particle {

    public enum ParticleName {
        CARBON6P, CARBON5P, CARBON4P, CARBON3P, CARBON2P, CARBON1P, CARBON, WATER, PROTON, HELIUM2P
    }

    private ParticleName pName;

    private double etat;
    private double netat;

    private double Z;
    private double N;
    public double m;
    public Vector R;
    public Vector V;

    // Initialize basic parameters
    public Particle(ParticleName pName) {
        // Tested OK 21-jun-15

        double etat0 = 0;
        double netat0 = 0;
        double etat1 = 0;
        double netat1 = 0;

        this.pName = pName;
        this.R = new Vector();
        this.V = new Vector();

        switch (this.pName) {
            case CARBON6P:
                m = 12.0*1836;
                Z = 6;
                N = 0;
                etat0 = 1;
                etat1 = 0;
                netat0 = 0;
                netat1 = 0;
                break;
            case CARBON5P:
                m = 12.0*1836;
                Z = 6;
                N = 1;
                etat0 = 2.625;
                etat1 = 1.2996;
                netat0 = 1.770;
                netat1 = 1.1402;
                break;
            case CARBON4P:
                m = 12.0*1836;
                Z = 6;
                N = 2;
                etat0 = 2.164;
                etat1 = 0.9764;
                netat0 = 1.750;
                netat1 = 0.6821;
                break;
            case CARBON3P:
                m = 12.0*1836;
                Z = 6;
                N = 3;
                etat0 = 1.30;
                etat1 = 0.6465;
                netat0 = 1.880;
                netat1 = 0.5547;
                break;
            case CARBON2P:
                m = 12.0*1836;
                Z = 6;
                N = 4;
                etat0 = 1.031;
                etat1 = 0.4924;
                netat0 = 2.0;
                netat1 = 0.4939;
                break;
            case CARBON1P:
                m = 12.0*1836;
                Z = 6;
                N = 5;
                etat0 = 1.065;
                etat1 = 0.48;
                netat0 = 2.13;
                netat1 = 0.4434;
                break;
            case CARBON:
                m = 12.0*1836;
                Z = 6;
                N = 6;
                etat0 = 1.179;
                etat1 = 0.4677;
                netat0 = 2.270;
                netat1 = 0.4143;
                break;
            case WATER:
                m = 18.0*1836;
                Z = 10;
                N = 9;
                etat0 = 1.792;
                etat1 = 0.4515;
                netat0 = 2.710;
                netat1 = 0.3671;
                break;
            case PROTON:
                m = 1.0*1836;
                Z = 1;
                N = 0;
                etat0 = 1;
                etat1 = 0;
                netat0 = 0;
                netat1 = 0;
                break;
            case HELIUM2P:
                m = 4.0*1836;
                Z = 2;
                N = 0;
                etat0 = 1;
                etat1 = 0;
                netat0 = 0;
                netat1 = 0;
                break;
        }
        netat = netat0 + netat1*(Z-N-1);
        etat = etat0 + etat1*(Z-N-1);
    }

    public void reset() {
        this.V.reset();
        this.R.reset();
    }

    // Eq. 7 Art. Liamsuwan
    // Omega = 1/{(neta/etat)[exp(r.etat)-1]+1}
    public double Omega(double r) {
        return 1 / (((this.netat * (Math.exp(r * this.etat) - 1)) / this.etat) + 1);
    }

    // Eq. 6 art. Liamsuwan
    // Zc = Z - N(1 - Omega(r))
    public double Zc(double r) {
        return this.Z - (this.N * (1 - this.Omega(r)));
    }

    public double Vc(double r) {
        return this.Zc(r)/r;
    }

    // dZ/dr = -N.neta.exp(r.etat).Omega^2
    public double dZc(double r) {
        return -this.N*this.netat*Math.exp(r * this.etat)*this.Omega(r)*this.Omega(r);
    }

    public String toString() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(1);
        df.setMinimumFractionDigits(1);
        String str = "\nParticle " + pName + " m = " + df.format(m) + " Z = " + df.format(Z) + " N = " + df.format(N);
        str += "\nPosition: ";
        str += R.toString();
        str += "\nVelocity: ";
        str += V.toString();
        str += "\n";
        return str;
    }
}
