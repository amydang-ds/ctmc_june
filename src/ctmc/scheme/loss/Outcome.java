package ctmc.scheme.loss;

import java.text.DecimalFormat;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 * Outcome Object is the result of a collision event
 * including:
 *          the interaction type
 *          energy of the escaped electron
 *          scattering angle of the electron
 */
public class Outcome {
    private int type;
    private double energy;
    private double angle;

    public Outcome(int type, double energy, double angle) {
        this.type = type;
        this.energy = energy;
        this.angle = angle;
    }

    public Outcome() {
    }

    public int getType() {
        return type;
    }

    public void setType(int type) {
        this.type = type;
    }

    public double getEnergy() {
        return energy;
    }

    public void setEnergy(double energy) {
        this.energy = energy;
    }

    public double getAngle() {
        return angle;
    }

    public void setAngle(double angle) {
        this.angle = angle;
    }

    public String toString() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(15);
        df.setMinimumFractionDigits(15);
        String str = "";
        str += type + "\t" + df.format(energy) + "\t" + df.format(angle);
        return str;
    }
}
