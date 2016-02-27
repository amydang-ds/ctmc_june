package ctmc.scheme.ionization;

import java.io.*;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Future;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 * This Class will handle all simulation cases
 * and the concurrency of CTMC program
 */
public class SimulationHandler {
    public static final Particle.ParticleName carbonIon = Particle.ParticleName.CARBON6P;
    //public static final double[] Barray = {0.03, 0.23, 0.43, 0.63, 0.83, 1.03, 1.23, 1.43, 1.63, 1.83, 2.03, 2.23, 2.43, 2.63, 2.83, 3.03, 3.23, 3.43, 3.63, 3.83, 4.03, 4.23, 4.43, 4.63, 4.83, 5.03, 5.23, 5.43, 5.63, 5.83, 6.03, 6.23, 6.43, 6.63, 6.83, 7.03};

    private static final double minB = 0.03;
    private static final double deltaB = 0.2;
    private static final double deltaBorb5 = 0.02;
    private static final int binsB = 35;
    private static final int binsEb = 4;
    private static final int simulationNumber = 10000;
    private double KE = 0;
    private int startIndex = 0;

    public SimulationHandler(double KE, int startIndex) {
        this.setKE(KE);
        this.setStartIndex(startIndex);
    }

    public int getStartIndex() {
        return this.startIndex;
    }

    public void setStartIndex(int startIndex) {
        this.startIndex = startIndex;
    }

    public double getKE() {
        return this.KE;
    }

    public static int getBinsB() {
        return binsB;
    }

    public static int getBinsEb() {
        return binsEb;
    }

    public static int getMaxIndex() {
        return binsB*binsEb;
    }

    public void setKE(double KE) {
        this.KE = KE;
    }

    public double getKineticEnergy() {
        switch (carbonIon) {
            case PROTON:
                return this.KE * 1 * 1e3 / Constant.ha;
            case HELIUM2P:
                return this.KE * 4 * 1e3 / Constant.ha;
        }
        return this.KE*12*1e3/Constant.ha;
    }

    public static int getIndexB(int index) {
        return index/binsEb;
    }

    public static double getImpactParameter(int index) {
        int indexB = getIndexB(index);
        return getIndexEb(index) < 4 ? (deltaB*indexB + minB) : (deltaBorb5*indexB + minB);
    }

    public static int getIndexEb(int index) {
        //Tested 25-Jun-15
        return index%binsEb;
    }

    public static RecGenerator.Orbital getOrbital(int index) {
        int indexEb = getIndexEb(index);

        switch (indexEb) {
            case 0:
                return RecGenerator.Orbital.EB1;
            case 1:
                return RecGenerator.Orbital.EB2;
            case 2:
                return RecGenerator.Orbital.EB3;
            case 3:
                return RecGenerator.Orbital.EB4;
            case 4:
                return RecGenerator.Orbital.EB5;
        }
        return RecGenerator.Orbital.EB1;
    }

    public void writeStatusFile(int index) {
        String strFilePath = null;
        BufferedWriter bw = null;
        try {
            strFilePath = "running_status.txt";
            bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
            bw.write("" + this.KE + "\t" + index);
            bw.close();
            System.out.println("Running status has been saved " + strFilePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] readStatusFile() throws IOException, ParseException{
        NumberFormat format = NumberFormat.getInstance(Locale.US);
        String strFilePath = "running_status.txt";
        BufferedReader br = null;

        double[] status = null;

        br = new BufferedReader(new FileReader(strFilePath));

        String line = null;
        line = br.readLine();
        while (line != null) {
            String[] columns = line.split("\t");
            if (columns.length == 2) {
                status = new double[2];
                status[0] = format.parse(columns[0]).doubleValue();
                status[1] = format.parse(columns[1]).intValue();
            }
            line = br.readLine();
        }
        br.close();

        return status;
    }

    public static void clearStatusFile() {
        String strFilePath = null;
        BufferedWriter bw = null;
        try {
            strFilePath = "running_status.txt";
            bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
            bw.close();
            System.out.println("Running status has been cleared " + strFilePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static int getSimulationNumber() {
        return simulationNumber;
    }

    public int getMaxJob() {
        return binsB*binsEb;
    }

    public void writeCrossSection(CrossSection cs, boolean backUp) {
        String strFilePath = null;
        BufferedWriter bw = null;
        try {
            if (backUp) {
                strFilePath = Double.toString(KE);
                bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
                bw.write(cs.toStringBackUp());
            } else {
                strFilePath = Double.toString(KE) + "keV_for_viewing.txt";
                bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
                bw.write(cs.toString());
            }
            bw.close();
            System.out.println("New file created! filename " + strFilePath);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void writeCrossSection(CrossSection cs) {
        this.writeCrossSection(cs, false);
    }

    public void analyzeOneJob(int index, CrossSection crossSection, ArrayList<Outcome> outcomes) {
        int countNonPhys = 0;
        int countIon = 0;
        int countCap = 0;
        int countBound = 0;
        int[] countSDCS = new int[CrossSection.binsEnergy];
        int[][] countDDCS = new int[CrossSection.binsEnergy][CrossSection.binsAngle];

        double Pion;
        double Psi;
        double Pti;
        double Pdi;
        double Pcap;
        double Psc;
        double Pdc;
        double Pnon;

        double impactParameter = getImpactParameter(index);
        double delta_b;

        if (getIndexEb(index) < 4) {
            delta_b = deltaB;
        } else {
            delta_b = deltaBorb5;
        }

        int countIgnore = 0;
        int countIgnore2 = 0;
        int binEnergy;
        int binAngle;
        for (Outcome outcome : outcomes) {
            // Ionization
            if (outcome.getType() == 1) {
                ++countIon;
                if (CrossSection.limitEnergy(outcome.getEnergy())) {
                    binEnergy = CrossSection.getEnergyBin(outcome.getEnergy());
                    ++countSDCS[binEnergy];
                    if (CrossSection.limitAngle(outcome.getAngle())) {
                        binAngle = CrossSection.getAngleBin(outcome.getAngle());
                        ++countDDCS[binEnergy][binAngle];
                    } else {
                        ++countIgnore2;
                    }
                } else {
                    System.out.println("Ignore " + outcome.getEnergy());
                    ++countIgnore;
                }
            } else if (outcome.getType() == 2) {
                // Capture
                ++countCap;
            } else if (outcome.getType() == 3) {
                // Remained bound
                ++countBound;
            } else {
                // Non-physical
                ++countNonPhys;
            }
        }
        int numberOfSims = getSimulationNumber();
        double validSimulations;
        validSimulations = numberOfSims - countNonPhys;
        System.out.println("count cap = " + countCap);
        System.out.println("count ignore = " + countIgnore);
        System.out.println("count ignore2 = " + countIgnore2);
        // each b-value
        Pion = ((double)countIon)/validSimulations;
        Pcap = ((double)countCap)/validSimulations;
        Pnon = 1.0 - Pion - Pcap;

        Psi = impactParameter*(2.0*Pion*Pnon)*delta_b;
        Psc = impactParameter*(2.0*Pcap*Pnon)*delta_b;

        Pti = impactParameter*(2.0*Pion*Pcap)*delta_b;
        Pdi = impactParameter*(2.0*Pion*Pion)*delta_b;
        Pdc = impactParameter*(2.0*Pcap * Pcap)*delta_b;

        crossSection.addPsi(Psi);
        crossSection.addPsc(Psc);
        crossSection.addPti(Pti);
        crossSection.addPdi(Pdi);
        crossSection.addPdc(Pdc);

        double sdcs;
        double ddcs;
        for (int i = 0; i < CrossSection.binsEnergy; ++i) {
            sdcs = (((double) countSDCS[i])/validSimulations)*impactParameter*delta_b/CrossSection.getEnergyRange(i);
            crossSection.addSDCS(i, sdcs);
            for (int j = 0; j < CrossSection.binsAngle; ++j) {
                ddcs = (((double) countDDCS[i][j])/validSimulations)*impactParameter*delta_b/(CrossSection.getEnergyRange(i)*Math.sin(CrossSection.getAngleMid(j))*CrossSection.getDeltaAngle());
                crossSection.addDDCS(i, j, ddcs);
            }
        }
    }

    public Future<ArrayList<Outcome>> startNextJob(ExecutorService executor, int index) {
        JobCallable job = new JobCallable(this.getKineticEnergy(), index);
        Future<ArrayList<Outcome>> future = executor.submit(job);
        System.out.println("Job " + index);
        return future;
    }

    public ArrayList<Outcome> readJobData(int index) throws IOException, ParseException {
        ArrayList<Outcome> outcomes = new ArrayList<Outcome>();
        NumberFormat format = NumberFormat.getInstance(Locale.US);
        String strFilePath = "" + index + ".txt";
        BufferedReader br = null;
        int type;
        double energy;
        double angle;

        br = new BufferedReader(new FileReader(strFilePath));

        String line = null;
        line = br.readLine();
        while (line != null) {
            String[] columns = line.split("\t");
            if (columns.length == 3) {
                Number number1 = format.parse(columns[1]);
                energy = number1.doubleValue();
                Number number2 = format.parse(columns[2]);
                angle = number2.doubleValue();
                Outcome outcome = new Outcome(Integer.parseInt(columns[0]), energy, angle);
                outcomes.add(outcome);
            }
            line = br.readLine();
        }
        br.close();

        return outcomes;
    }

    public void writeJob(int index, ArrayList<Outcome> outcomes) throws IOException {
        String strFilePath = "" + index + ".txt";
        BufferedWriter bw = null;

        bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
        System.out.println("New file created! index " + index);
        for (Iterator<Outcome> iterator = outcomes.iterator(); iterator.hasNext(); ) {
            Outcome outcome = iterator.next();
            System.out.println(outcome);
            bw.write(outcome.toString());
            bw.newLine();
        }
        bw.close();

    }

    public static void testRec() throws IOException, NullPointerException{

        for (RecGenerator.Orbital orb : RecGenerator.Orbital.values()) {
            RecGenerator r_gen = new RecGenerator(orb);
            Random rand = new Random();
            double r;
            double p;
            int lines = 10000000;
            DecimalFormat df = new DecimalFormat();
            df.setMaximumFractionDigits(10);
            df.setMinimumFractionDigits(10);

            String strFilePath = "Rec_" + String.valueOf(orb) + ".dat";
            String strFilePath1 = "Pec_" + String.valueOf(orb) + ".dat";

            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(strFilePath)));
            BufferedWriter bw1 = new BufferedWriter(new FileWriter(new File(strFilePath1)));

            for (int i = 0; i < lines; ++i) {
                r = r_gen.nextRec();
                System.out.println(i + "   " + r);
                p = r_gen.Pec(r);

                bw.write(df.format(r));
                bw.newLine();

                bw1.write(df.format(p));
                bw1.newLine();

            }
            bw.close();
            bw1.close();
        }
    }

    public void testOneSimulation() {
        // Test one simulation
        RecGenerator r_gen = new RecGenerator(RecGenerator.Orbital.EB1);
        Random rand1 = new Random();
        Random rand2 = new Random();
        Random rand3 = new Random();

        double KE = 12.*4.00e+5/27.2114;
        double b = 0.7;
        IonizationScheme ionizationScheme = new IonizationScheme(RecGenerator.Orbital.EB1, KE, b, carbonIon);
        ionizationScheme.initialize(r_gen, rand1, rand2, rand3);
        System.out.println(ionizationScheme.toString());
        double time = ionizationScheme.RungeKuttaABC();
        System.out.println(ionizationScheme.toString());
        System.out.println("time = " + time);
        System.out.println("Eet " + ionizationScheme.Eet());
        System.out.println("Eep " + ionizationScheme.Eep());
        System.out.println("theta " + ionizationScheme.theta());
        System.out.println("theta " + ionizationScheme.calcAzimuthalAngle());
        System.out.println("outcome " + ionizationScheme.getOutcome());
    }

    public static void testRec3DDistribution() throws IOException {
        RecGenerator r_gen = new RecGenerator(RecGenerator.Orbital.EB1);
        // These Random objects have different seeds and give different sequence of number (No need to worry about that!)
        Random rand1 = new Random();
        //Random rand2 = new Random();
        //Random rand3 = new Random();

        double rec;
        double Pi = 3.1415926535897D;
        //Random spherical coordinates
        double teta;
        double phi;
        double neta;

        int lines = 10000;
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(10);
        df.setMinimumFractionDigits(10);

        String strFilePath = "Rec_3D.txt";
        BufferedWriter bw = null;

        bw = new BufferedWriter(new FileWriter(new File(strFilePath)));


        for (int i = 0; i < lines; ++i) {
            rec = r_gen.nextRec();
            teta = rand1.nextDouble() * Pi;
            phi = rand1.nextDouble() * 2.0 * Pi;
            neta = rand1.nextDouble() * 2.0 * Pi;

            Vector Rec = new Vector(Math.sin(teta)*Math.cos(phi), Math.sin(teta)*Math.sin(phi), Math.cos(teta));
            Rec.times(rec);

            bw.write(df.format(Rec.x) +"\t"+ df.format(Rec.y)+"\t"+ df.format(Rec.z));
            bw.newLine();

        }

        bw.close();

    }

    public static void getPotentialDistribution (Particle.ParticleName particleName) throws IOException {
        Particle particle = new Particle(particleName);

        double distanceMax = 7;
        double distanceMin = 0.01;
        double step = 0.02;
        double distance = distanceMin;
        double potential = 0;

        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(10);
        df.setMinimumFractionDigits(10);

        String strFilePath = "Zc_" + particleName + ".dat";
        BufferedWriter bw = null;

        bw = new BufferedWriter(new FileWriter(new File(strFilePath)));

        while (distance < distanceMax) {
            potential = particle.Zc(distance);

            bw.write(df.format(distance) +"\t"+ df.format(potential));
            bw.newLine();

            distance += step;
        }
        bw.close();
    }
}
