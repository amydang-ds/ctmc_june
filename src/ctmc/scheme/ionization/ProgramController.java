package ctmc.scheme.ionization;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 *
 Carbon 6
 10  50  100  200  300 500 1000 6000 10000
other
 10
 50
 100
 500
 1000
 10000

 */
public class ProgramController {
    // KINETIC ENERGY
    // private static double[] KEarray = {10, 50, 100, 200, 300, 500, 1000, 6000, 10000, 50000, 100000};
    private static double[] KEarray = {50000, 100000};
    //private static double[] KEarray = {10, 50, 100, 200, 300, 500, 1000, 6000, 10000};
    //private static double[] KEarray = {6000};
    // private static double[] KEarray = {10, 50, 100, 500, 1000, 10000};

    public static int numberOfKE = KEarray.length;

    public static double getKE(int i) {
        return KEarray[i];
    }

    public static int getIndexKE(double KE) {
        int resultIndex = 0;
        for (int i = 0; i < KEarray.length; ++i) {
            if (getKE(i) == KE) {
                resultIndex = i;
                break;
            }
        }
        System.out.println("index KE " + resultIndex);
        System.out.println("KE " + getKE(resultIndex));
        return resultIndex;
    }

    public static void simulateOneKineticEnergy(double KE, int startIndex) throws IOException, ParseException, InterruptedException, ExecutionException {
        ExecutorService executor = Executors.newFixedThreadPool(200);
        ArrayList<Future<ArrayList<Outcome>>> outcomeSets = new ArrayList<Future<ArrayList<Outcome>>>();
        SimulationHandler simulationHandler = new SimulationHandler(KE, startIndex);

        for (int i = simulationHandler.getStartIndex(); i < simulationHandler.getMaxJob(); ++i) {
            outcomeSets.add(simulationHandler.startNextJob(executor, i));
        }

        System.out.println("outcome sets: " + outcomeSets.size());
        for (int i = 0; i < outcomeSets.size(); i++) {
            Future<ArrayList<Outcome>> fut = outcomeSets.get(i);
            ArrayList<Outcome> outcomes = fut.get();
            simulationHandler.writeJob(i, outcomes);
        }
        executor.shutdown();

        CrossSection crossSection = new CrossSection();

        for (int i = 0; i < simulationHandler.getMaxJob(); ++i) {
            ArrayList<Outcome> outcomes = simulationHandler.readJobData(i);
            for (Outcome outcome : outcomes) {
                System.out.println(outcome);
            }
            simulationHandler.analyzeOneJob(i, crossSection, outcomes);
        }
        crossSection.multiplyUnits(); // must do to get final cross sections
        simulationHandler.writeCrossSection(crossSection, true); // back up for power-cut
        simulationHandler.writeCrossSection(crossSection, false); // for viewing
    }

    public static void simulateOneKineticEnergy(double KE)  throws IOException, ParseException, InterruptedException, ExecutionException {
        int startIndex = 0;
        simulateOneKineticEnergy(KE, startIndex);
    }

    public static void writeCrossSectionOneKE(SimulationHandler simulationHandler) throws IOException, ParseException{
        CrossSection crossSection = new CrossSection();
        for (int i = 0; i < simulationHandler.getMaxJob(); ++i) {
            ArrayList<Outcome> outcomes = simulationHandler.readJobData(i);
            for (Outcome outcome : outcomes) {
                System.out.println(outcome);
            }
            simulationHandler.analyzeOneJob(i, crossSection, outcomes);
        }
        crossSection.multiplyUnits(); // must do to get final cross sections
        simulationHandler.writeCrossSection(crossSection, true); // back up for power-cut
        simulationHandler.writeCrossSection(crossSection, false); // for viewing
    }

    public static void exportAllCrossSections(CrossSection[] crossSections) {
        // NEED TO BE REFACTORED -> ARRAY ATTR.
        StringBuilder si = new StringBuilder();
        StringBuilder di = new StringBuilder();
        StringBuilder sc = new StringBuilder();
        StringBuilder dc = new StringBuilder();
        StringBuilder ti = new StringBuilder();
        StringBuilder sdcs = new StringBuilder();
        StringBuilder ddcs = new StringBuilder();
        double KE;
        // write cross-section data
        for (int i = 0; i < KEarray.length; ++i) {
            KE = getKE(i);
            si.append(crossSections[i].toStringPsi(KE));
            di.append(crossSections[i].toStringPdi(KE));

            sc.append(crossSections[i].toStringPsc(KE));
            dc.append(crossSections[i].toStringPdc(KE));

            ti.append(crossSections[i].toStringPti(KE));

            sdcs.append(crossSections[i].toStringSDCS(KE));
            ddcs.append(crossSections[i].toStringDDCS(KE));
        }

        String strFilePathSI = "SI.dat";
        BufferedWriter bw_si = null;
        String strFilePathDI = "DI.dat";
        BufferedWriter bw_di = null;

        String strFilePathSC = "SC.dat";
        BufferedWriter bw_sc = null;
        String strFilePathDC = "DC.dat";
        BufferedWriter bw_dc = null;

        String strFilePathTI = "TI.dat";
        BufferedWriter bw_ti = null;

        String strFilePathSDCS = "SDCS.dat";
        BufferedWriter bw_sdcs = null;
        String strFilePathDDCS = "DDCS.dat";
        BufferedWriter bw_ddcs = null;
        try {
            bw_si = new BufferedWriter(new FileWriter(new File(strFilePathSI)));
            System.out.println("New file created! filename " + strFilePathSI);
            bw_si.write(si.toString());
            bw_si.close();

            bw_di = new BufferedWriter(new FileWriter(new File(strFilePathDI)));
            System.out.println("New file created! filename " + strFilePathDI);
            bw_di.write(di.toString());
            bw_di.close();

            bw_sc = new BufferedWriter(new FileWriter(new File(strFilePathSC)));
            System.out.println("New file created! filename " + strFilePathSC);
            bw_sc.write(sc.toString());
            bw_sc.close();

            bw_dc = new BufferedWriter(new FileWriter(new File(strFilePathDC)));
            System.out.println("New file created! filename " + strFilePathDC);
            bw_dc.write(dc.toString());
            bw_dc.close();

            bw_ti = new BufferedWriter(new FileWriter(new File(strFilePathTI)));
            System.out.println("New file created! filename " + strFilePathTI);
            bw_ti.write(ti.toString());
            bw_ti.close();

            bw_sdcs = new BufferedWriter(new FileWriter(new File(strFilePathSDCS)));
            System.out.println("New file created! filename " + strFilePathSDCS);
            bw_sdcs.write(sdcs.toString());
            bw_sdcs.close();

            bw_ddcs = new BufferedWriter(new FileWriter(new File(strFilePathDDCS)));
            System.out.println("New file created! filename " + strFilePathDDCS);
            bw_ddcs.write(ddcs.toString());
            bw_ddcs.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static CrossSection[] readAllCrossSection() throws IOException, ParseException{
        CrossSection[] results = new CrossSection[KEarray.length];

        for (int i = 0; i < KEarray.length; ++i) {
            System.out.println(Double.toString(getKE(i)));
            results[i] = new CrossSection();
            if (results[i].readBackUp(Double.toString(getKE(i))) ) {
                System.out.println("Read file " + getKE(i) + "successful!");
            } else {
                System.out.println("Failed to read file " + getKE(i));
            }
        }

        return results;
    }

    public static void runCTMC()  throws Exception{
        double[] status = SimulationHandler.readStatusFile();

        if (status == null) { // A fresh start, status is empty
            System.out.println("A FRESH START!!");
            // inverse loop for efficiency
            for (int i = ProgramController.numberOfKE - 1; i >= 0; --i) {
                simulateOneKineticEnergy(getKE(i));
            }
        } else { // Restart after power-cut happened
            // loop starts from backUp KE down to 0
            int backedUpIndexKE = getIndexKE(status[0]);
            int startIndex = (int) status[1];

            if ((int) status[1] == SimulationHandler.getMaxIndex() - 1) { // Finished backed-up KE, move on to next KE
                // Unfortunately, the cross-section file has not been saved so calculate cross-section before moving on
                SimulationHandler simulationHandler = new SimulationHandler(getKE(backedUpIndexKE), startIndex);
                writeCrossSectionOneKE(simulationHandler);

                ++backedUpIndexKE; // next KE
                startIndex = -1;
            } // Haizz what if the cross section file has not been writen =((
            ++startIndex;

            System.out.println("Restart at KE " + getKE(backedUpIndexKE) + " index " + startIndex);
            // restart the interrupted loop
            for (int i = backedUpIndexKE; i >= 0; --i) {
                simulateOneKineticEnergy(ProgramController.getKE(i), startIndex);
            }
        }
        SimulationHandler.clearStatusFile();

        // read cross-sections of all projectile's energy before combine to one file
        CrossSection[] dataSet = readAllCrossSection();

        exportAllCrossSections(dataSet);
    }

}
