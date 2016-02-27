import ctmc.scheme.ionization.*;

import java.util.Random;

public class Main {

    public static void main(String[] args) throws Exception{
        long start = System.currentTimeMillis();


        //ProgramController.runCTMC();

        //SimulationHandler.getPotentialDistribution(Particle.ParticleName.WATER);
        //RecGenerator.getProbCurves();
        SimulationHandler.testRec();

        long end = System.currentTimeMillis();
        long elapse_time = (end - start)/(1000*60);
        System.out.println("Computing time (minutes) " + elapse_time);

    }
}
