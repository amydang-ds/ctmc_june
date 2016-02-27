package ctmc.scheme.ionization;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.Random;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 *
 * CODE IS FIXED FOR WATER MOLECULE
 *
 * Generate a sequence of Rec based on a classical distribution function
 * by applying Monte Carlo method - Rejection Sampling
 */
public class RecGenerator {

    // Electron binding energy enums
    public enum Orbital {
        EB1, EB2, EB3, EB4, EB5
    }

    // atomic Constant.ha system
    private static final double Eb1 = 12.62/Constant.ha;
    private static final double Eb2 = 14.75/Constant.ha;
    private static final double Eb3 = 18.51/Constant.ha;
    private static final double Eb4 = 32.4/Constant.ha;
    private static final double Eb5 = 539.7/Constant.ha;
    private double Eb;

    // electron mass is 1 a.u
    private static final double me = 1;
    // a nucleon = 1836 a.u -> Water molecule is 18x
    private static final double m = 18.0*1836;
    private static final double mec = me*m/(me+m);
    // 2 Hydrogen and 1 Oxygen -> charge: 10
    private static final double Z = 10;
    // 1 electron is active, the remain is 9
    private static final double N = 9;


    // Table 2. art. Thiamsuwan parameters for calculating effective charge
    private static final double etat0 = 1.792;
    private static final double etat1 = 0.4515;
    private static final double netat0 = 2.710;
    private static final double netat1 = 0.3671;

    private static final double netat = netat0 + netat1*(Z-N-1);
    private static final double etat = etat0 + etat1*(Z-N-1);

    /**
     *  Parameters of Rejection Sampling
     */
    private double r1, r2, rmax;
    private double ymax;

    private Random rand;

    // Contructor will instantiate the random object
    // and set the selected value to Eb
    public RecGenerator(Orbital orb) {
        rand = new Random();

        switch (orb) {
            case EB1:
                Eb = Eb1;
                break;
            case EB2:
                Eb = Eb2;
                break;
            case EB3:
                Eb = Eb3;
                break;
            case EB4:
                Eb = Eb4;
                break;
            case EB5:
                Eb = Eb5;
                break;
        }

        // Find values of all parameters
        double a = 0.01;
        double b = 10000.0;
        rmax = FindRmax(a, b);
        ymax = FindYmax();
        setRange();
    }

    public static double getValue(Orbital orb) {
        switch (orb) {
            case EB1:
                return Eb1;
            case EB2:
                return Eb2;
            case EB3:
                return Eb3;
            case EB4:
                return Eb4;
            case EB5:
                return Eb5;
        }
        return Eb1;
    }

    public static double getEb1() {
        return Eb1;
    }

    public static double getEb2() {
        return Eb2;
    }

    public static double getEb3() {
        return Eb3;
    }

    public static double getEb4() {
        return Eb4;
    }

    public static double getEb5() {
        return Eb5;
    }

    // Eq. 7 Art. Liamsuwan
    // Omega = 1/{(neta/etat)[exp(r.etat)-1]+1}
    public static double Omega(double r) {
        return 1 / (((netat * (Math.exp(r * etat) - 1)) / etat) + 1);
    }

    // Eq. 6 art. Liamsuwan
    // Zc = Z - N(1 - Omega(r))
    public static double Zc(double r) {
        return Z - N*(1-Omega(r));
    }

    // Vc = Zc/r
    public static double Vc(double r) {
        return Zc(r)/r;
    }

    // Eq. (3) art.Liamsuwan
    // Pec = 2mec.sqrt(-|Eb|+V(r))
    public double Pec(double r) {
        return Math.sqrt(2*mec*(-Math.abs(this.Eb)+Vc(r)));
    }

    // Golden search for rmax
    // rmax is the solution of eq. -|Eb| + V(r) = 0
    private double FindRmax(double x, double y) {
        double a;
        double b;
        double c;
        double tolerant = 1.0e-14;
        a = x;
        b = y;
        while (Math.abs(b-a) >= tolerant) {
            c = (a + b)/2.0;
            if ((Vc(a) - this.Eb)*(Vc(c) - this.Eb) > 0)
                a = c;
            else
                b = c;
        }
        c = a;
        return c;
    }

    private void setRange() {
        // Setting range of Rec
        r1 = 0.03;
        r2 = 0.9*rmax;
        if (Eb == Eb5)
            r2 = 0.87*rmax;
    }

    private double FindYmax() {
        // linear space vector of r values
        double[] r = new double[10000];
        r[0] = 1.0e-7;

        double delta = (rmax - r[0])/10000.0;
        for (int i = 0; i < 10000-1; i++) {
            r[i+1] = r[i] + delta;
        }

        // ----------------- find ymax ---------------
        double ytemp;
        double ymax = 1.0e-6;

        // Eq. (5) art. Liamsuwan
        // y = mec.r^2.sqrt(2mec(-|Eb|+V(r)))
        for (int i = 0; i < 10000; i++){
            ytemp = mec*r[i]*r[i]*Math.sqrt(2 * mec * (-Math.abs(this.Eb) + Vc(r[i])));
            // finding ymax value
            if (ytemp > ymax)
                ymax = ytemp;
        }
        return ymax;
    }

    // return the next value of rec in the sequence
    // Sampling rejection method
    public double nextRec() {
        double rec_temp;
        double Ptemp;

        while (true){
            rec_temp = rand.nextDouble()*rmax;
            // rec is in the range
            if (rec_temp >= r1 && rec_temp <=r2){
                Ptemp = mec*rec_temp*rec_temp*Math.sqrt(2 * mec * (-Math.abs(Eb) + Vc(rec_temp)));
                // generate a random number (0,1) and compare to w(re)/w(rmax)
                // Ptemp is the value on the distribution curve
                // Accept when rand*ymax under the curve
                if (Ptemp/ymax > rand.nextDouble())
                    break;
            }
        }
        return rec_temp;
    }

    public double probabilityFunctionRec(double r) {
        return mec*r*r*Math.sqrt(2 * mec * (-Math.abs(Eb) + Vc(r)));
    }

    public double probabilityFunctionPec(double p) {
        double p0 = 2*mec*Eb;
        double K3 = 8*Math.pow(p0, 5)/Math.pow(Constant.Pi, 2);
        return K3/Math.pow(p0*p0 + p*p,4);
    }

    public static void getProbCurves()  throws IOException, NullPointerException{
        for (RecGenerator.Orbital orb : RecGenerator.Orbital.values()) {
            RecGenerator gen = new RecGenerator(orb);

            double Rec = 0;
            double Pec = 0;

            double max = 25;
            double min = 0.01;
            double step = 0.02;

            double var = min;
            DecimalFormat df = new DecimalFormat();
            df.setMaximumFractionDigits(10);
            df.setMinimumFractionDigits(10);

            String strFilePath = "curveRec_" + String.valueOf(orb) + ".dat";
            String strFilePath1 = "curvePec_" + String.valueOf(orb) + ".dat";

            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
            BufferedWriter bw1 = new BufferedWriter(new FileWriter(new File(strFilePath1)));

            while (var < max) {
                Rec = gen.probabilityFunctionRec(var);
                Pec = gen.probabilityFunctionPec(var);

                bw.write(df.format(var) +"\t"+ df.format(Rec));
                bw.newLine();

                bw1.write(df.format(var) +"\t"+ df.format(Pec));
                bw1.newLine();

                var += step;
            }

            bw.close();
            bw1.close();
        }
    }
}
