package ctmc.scheme.loss;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.Locale;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 * CrossSection Object contains
 *                          DDCS
 *                          SDCS
 *                          TCS
 */
public class CrossSection {
    public static final int binsEnergy = 20;
    public static final int binsAngle = 17;

    // ENERGY
    private static final double logEmin = 0.7;
    private static final double logEmax = 3.2578947368421;//3.13;
    private static final double Emin = Math.pow(10.0,logEmin);
    private static final double Emax = Math.pow(10.0,logEmax);
    private static final double deltalogE = 0.127894736842105;//(logEmax - logEmin)/binsEnergy;
    // ANGLE
    private static final double angleMin = 5.0* Constant.Pi/180.0;
    private static final double angleMax = 175.0* Constant.Pi/180.0;
    private static final double deltaAngle = 10.0* Constant.Pi/180.0;

    private double Psi;
    private double Psc;
    private double Pti;
    private double Pdi;
    private double Pdc;

    private double[] SDCS;
    private double[][] DDCS;

    public CrossSection() {
        Psi = 0;
        Psc = 0;
        Pti = 0;
        Pdi = 0;
        Pdc = 0;

        SDCS = new double[binsEnergy];
        DDCS = new double[binsEnergy][binsAngle];
    }

    public static boolean limitEnergy(double energy) {
        return (energy >= Emin) && (energy <= Emax);
    }

    public static boolean limitAngle(double angle) {
        return (angle >= angleMin) && (angle <= angleMax);
    }

    public static int getEnergyBin(double energy) {
        return (int) ((Math.log10(energy)-logEmin)/deltalogE);
    }

    public static int getAngleBin(double angle) {
        return (int) ((angle - angleMin)/deltaAngle);
    }

    public static double getEnergyRange(int binEnergy) {
        if (binEnergy == 0) {
            return Math.pow(10.0, logEmin + deltalogE*binEnergy);
        } else {
            double energyLeft = Math.pow(10.0, logEmin + deltalogE * (binEnergy - 1));
            double energyRight = Math.pow(10.0, logEmin + deltalogE * binEnergy);
            return (energyRight - energyLeft);
        }
    }

    public static double getDeltaAngle() {
        return deltaAngle;
    }

    public static double getEnergy(int binEnergy) {
        return Math.pow(10.0, logEmin + deltalogE*binEnergy);
    }

    public static double getAngleMid(int binAngle) {
        double angleLeft = angleMin + deltaAngle*binAngle;
        double angleRight = angleMin + deltaAngle*(binAngle + 1);
        return (angleRight + angleLeft)/2.0;
    }

    public static double getEmin() {
        return Emin;
    }

    public void addPsi(double Psi) {
        this.Psi += Psi;
    }

    public void addPsc(double Psc) {
        this.Psc += Psc;
    }

    public void addPti(double Pti) {
        this.Pti += Pti;
    }

    public void addPdi(double Pdi) {
        this.Pdi += Pdi;
    }

    public void addPdc(double Pdc) {
        this.Pdc += Pdc;
    }

    public void addSDCS(int binEnergy, double value) {
        if (binEnergy < binsEnergy && binEnergy >= 0) {
            this.SDCS[binEnergy] += value;
        } else {
            System.out.println("bin Energy does not exist! " + binEnergy);
        }
    }

    public void addDDCS(int binEnergy, int binAngle, double value) {
        if (binEnergy < binsEnergy && binEnergy >= 0) {
            if (binAngle < binsAngle && binAngle >= 0) {
                this.DDCS[binEnergy][binAngle] += value;
            } else {
                System.out.println("bin Angle does not exist! " + binAngle);
            }
        } else {
            System.out.println("bin Energy does not exist! " + binEnergy);
        }
    }

    public String toStringPsi(double KE) {
        return Double.toString(KE) + "\t" + Double.toString(this.Psi) + "\n";
    }

    public String toStringPdi(double KE) {
        return Double.toString(KE) + "\t" + Double.toString(this.Pdi) + "\n";
    }

    public String toStringPsc(double KE) {
        return Double.toString(KE) + "\t" + Double.toString(this.Psc) + "\n";
    }

    public String toStringPdc(double KE) {
        return Double.toString(KE) + "\t" + Double.toString(this.Pdc) + "\n";
    }

    public String toStringPti(double KE) {
        return Double.toString(KE) + "\t" + Double.toString(this.Pti) + "\n";
    }

    public String toStringSDCS(double KE) {
        StringBuilder sdcs = new StringBuilder();
        for (int i = 0; i < binsEnergy; ++i) {
            sdcs.append(Double.toString(KE) + "\t" + getEnergy(i) + "\t" + this.SDCS[i] + "\n");
        }
        return sdcs.toString();
    }

    public String toStringDDCS(double KE) {
        StringBuilder ddcs = new StringBuilder();
        for (int i = 0; i < binsEnergy; ++i) {
            for (int j = 0; j < binsAngle; ++j) {
                ddcs.append(Double.toString(KE) + "\t" + getEnergy(i) + "\t" + getAngleMid(j) + "\t" + this.DDCS[i][j] + "\n");
            }
        }
        return ddcs.toString();
    }

    public String toString() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(10);
        df.setMinimumFractionDigits(10);
        StringBuilder sb = new StringBuilder();
        sb.append("-----------------------------------------------------------\n" +
                "                     TOTAL CROSS SECTION                   \n" +
                "-----------------------------------------------------------\n");
        sb.append("Single ionization : ");
        sb.append(Psi);
        sb.append("\nSingle capture : ");
        sb.append(Psc);
        sb.append("\nTransfer ionization : ");
        sb.append(Pti);
        sb.append("\nDouble ionization : ");
        sb.append(Pdi);
        sb.append("\nDouble capture : ");
        sb.append(Pdc);

        sb.append("\n-----------------------------------------------------------\n" +
                "               SINGLY DIFFERENTIAL CROSS SECTION           \n" +
                "-----------------------------------------------------------\n");
        for (int i = 0; i < binsEnergy; ++i) {
            sb.append(this.SDCS[i]);
            sb.append("\n");
        }

        sb.append("-----------------------------------------------------------\n" +
                "               DOUBLY DIFFERENTIAL CROSS SECTION           \n" +
                "-----------------------------------------------------------\n");
        sb.append("Energy range\n");
        for (int j = 0; j < binsEnergy; ++j) {
            sb.append(getEnergy(j));
            sb.append("\n");
        }
        sb.append("\n");

        for (int i = 0; i < binsAngle; ++i) {
            sb.append("" + (int) (getAngleMid(i)*180/ Constant.Pi) + ")-------------------------------");
            for (int j = 0; j < binsEnergy; ++j) {
                sb.append("\n");
                sb.append(this.DDCS[j][i]);
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    public String toStringBackUp() {
        StringBuilder sb = new StringBuilder();
        sb.append(Psi);
        sb.append("\t");
        sb.append(Psc);
        sb.append("\t");
        sb.append(Pti);
        sb.append("\t");
        sb.append(Pdi);
        sb.append("\t");
        sb.append(Pdc);
        sb.append("\n");

        for (int i = 0; i < binsEnergy; ++i) {
            sb.append(this.SDCS[i]);
            sb.append("\n");
        }

        for (int i = 0; i < binsAngle; ++i) {
            for (int j = 0; j < binsEnergy; ++j) {
                sb.append(this.DDCS[j][i]);
                sb.append("\n");
            }
        }
        return sb.toString();
    }

    public boolean readBackUp(String fileName) throws IOException, ParseException{
        boolean readSuccess = true;
        CrossSection readCrossSection = new CrossSection();
        NumberFormat format = NumberFormat.getInstance(Locale.US);
        String strFilePath = "" + fileName;
        BufferedReader br = null;

        br = new BufferedReader(new FileReader(strFilePath));


        String line = null;

        // Read TCS
        line = br.readLine();
        String[] columns = line.split("\t");
        if (columns.length == 5) {
            this.Psi = format.parse(columns[0]).doubleValue();
            this.Psc = format.parse(columns[1]).doubleValue();
            this.Pti = format.parse(columns[2]).doubleValue();
            this.Pdi = format.parse(columns[3]).doubleValue();
            this.Pdc = format.parse(columns[4]).doubleValue();
        } else {
            System.out.println("There is an error on the first line. Please check your file!");
            readSuccess = false;
        }

        // Read SDCS
        line = br.readLine();
        int countLines = 0;
        while (line != null && countLines < binsEnergy) {
            this.SDCS[countLines] = format.parse(line).doubleValue();
            ++countLines;
            line = br.readLine();
        }

        System.out.println("count sdcs " + countLines);

        if (countLines < binsEnergy) {
            System.out.println("Cannot read SDCS correctly. Please check your file!");
            readSuccess = false;
        }

        // Read DDCS
        countLines = 0;
        int countLinesOfOneAngle = 0;
        for (int i = 0; i < binsAngle; ++i) {
            countLinesOfOneAngle = 0;
            while (line != null && countLinesOfOneAngle < binsEnergy) {
                this.DDCS[countLinesOfOneAngle][i] = format.parse(line).doubleValue();
                ++countLinesOfOneAngle;
                ++countLines;
                line = br.readLine();
            }
            System.out.println("count of one angle " + countLinesOfOneAngle);
            if (countLinesOfOneAngle < binsEnergy) {
                System.out.println("Cannot read DDCS correctly. Please check your file!");
                readSuccess = false;
            }
        }
        System.out.println("count Lines " + countLines);
        System.out.println("binsExbinsA " + binsEnergy*binsAngle);

        br.close();

        return readSuccess;
    }

    public void multiplyUnits() {
        this.Psi *= 2.0* Constant.Pi* Constant.a0* Constant.a0;
        this.Psc *= 2.0* Constant.Pi* Constant.a0* Constant.a0;
        this.Pti *= 2.0* Constant.Pi* Constant.a0* Constant.a0;
        this.Pdi *= 2.0* Constant.Pi* Constant.a0* Constant.a0;
        this.Pdc *= 2.0* Constant.Pi* Constant.a0* Constant.a0;

        for (int i = 0; i < binsEnergy; ++i) {
            this.SDCS[i] *= 2.0* Constant.Pi* Constant.a0* Constant.a0;
            for (int j = 0; j < binsAngle; ++j) {
                this.DDCS[i][j] *= Constant.a0* Constant.a0;
            }
        }
    }
}
