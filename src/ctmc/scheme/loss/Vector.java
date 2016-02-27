package ctmc.scheme.loss;

import java.text.DecimalFormat;

/**
 * Created by Mymy Huntress on 31-Jan-16.
 * 3-Dimension Vector
 */
public class Vector {
    public double x;
    public double y;
    public double z;

    public Vector() {
    }

    public Vector(Vector that) {
        this.x = that.x;
        this.y = that.y;
        this.z = that.z;
    }

    public Vector(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public void set(double x, double y, double z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public double length() {
        return Math.sqrt(x*x + y*y + z*z);
    }

    public double square() {
        return x*x+y*y+z*z;
    }

    public String toString() {
        DecimalFormat df = new DecimalFormat();
        df.setMaximumFractionDigits(25);
        df.setMinimumFractionDigits(25);
        String str = "(";
        str += df.format(this.x);
        str += "\t";
        str += df.format(this.y);
        str += "\t";
        str += df.format(this.z);
        str += ")";
        return str;
    }

    public void reset() {
        this.x = 0;
        this.y = 0;
        this.z = 0;
    }

    /* * * * *  return a vector  * * */
    public static Vector sum(Vector me, Vector that) {
        Vector sum = new Vector(me.x + that.x, me.y + that.y, me.z + that.z);
        return sum;
    }


    public static Vector subtract(Vector me, Vector that) {
        Vector rs = new Vector(me.x - that.x, me.y - that.y, me.z - that.z);
        return rs;
    }

    public static Vector times(Vector me, Vector that) {
        Vector rs = new Vector(me.x * that.x, me.y * that.y, me.z * that.z);
        return rs;
    }

    public static Vector times(Vector me, double constant) {
        Vector rs = new Vector(me.x * constant, me.y * constant, me.z * constant);
        return rs;
    }

    /* * * * * resulting new values on this vector * * */
    public void plus(Vector that) {
        this.x += that.x;
        this.y += that.y;
        this.z += that.z;
    }

    public void minus(Vector that) {
        this.x -= that.x;
        this.y -= that.y;
        this.z -= that.z;
    }

    // multiply a number
    public void times(double constant) {
        this.x *= constant;
        this.y *= constant;
        this.z *= constant;
    }
}

