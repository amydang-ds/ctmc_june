package ctmc.scheme.loss;

import java.util.ArrayList;
import java.util.Random;
import java.util.concurrent.Callable;

/**
 * Created by Mymy on 24-Jun-15.
 */
public class JobCallable implements Callable<ArrayList<Outcome>> {
    private final int index;
    private final double KE;

    public JobCallable(double KE, int index) {
        this.index = index;
        this.KE = KE;
    }

    @Override
    public ArrayList<Outcome> call() throws Exception {
        // 10-Jul-15 r_gen and ionScheme must have the same orb
        ArrayList<Outcome> resultSet = new ArrayList<Outcome>();
        Outcome outcome;

        double b = SimulationHandler.getImpactParameter(this.index);
        int numberOfSimulations = SimulationHandler.getSimulationNumber();
        RecGenerator.Orbital orb = SimulationHandler.getOrbital(this.index);
        System.out.println("KE " + this.KE + " b " + b + " orb " + orb);

        RecGenerator r_gen = new RecGenerator(orb, SimulationHandler.carbonIon);
        // These Random objects have different seeds and give different sequence of number (No need to worry about that!)
        Random rand1 = new Random();
        Random rand2 = new Random();
        Random rand3 = new Random();

        LossScheme lossScheme = new LossScheme(orb, this.KE, b, SimulationHandler.carbonIon);

        for (int i = 0; i < numberOfSimulations; ++i) {
            lossScheme.initialize(r_gen, rand1, rand2, rand3);
            lossScheme.RungeKuttaABC();
            outcome = lossScheme.getOutcome();
            resultSet.add(outcome);
        }
        return resultSet;
    }

    public int getIndex() {
        return index;
    }
}
