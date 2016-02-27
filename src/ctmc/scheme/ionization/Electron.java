package ctmc.scheme.ionization;

import java.text.DecimalFormat;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 * Electron model:
 *                  position vector R
 *                  velocity vector V
 *                  mass m
 *                  charge Z
 */
public class Electron {
    public double m;
    public double Z;
    public Vector R;
    public Vector V;

    public Electron() {
        m = 1;
        Z = -1;

        R = new Vector();
        V = new Vector();
    }

    public void reset() {
        m = 1;
        Z = -1;

        R.reset();
        V.reset();
    }

    public String toString() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(1);
        df.setMinimumFractionDigits(1);
        String str = "\nElectron m = " + df.format(m) + " Z = " + df.format(Z);
        str += "\nPosition: ";
        str += R.toString();
        str += "\nVelocity: ";
        str += V.toString();
        str += "\n";
        return str;
    }
}
