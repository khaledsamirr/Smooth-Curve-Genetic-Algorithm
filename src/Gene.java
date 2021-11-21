public class Gene {
    private double coefficient;

    public Gene(double coefficient) {
        this.coefficient = coefficient;
    }

    public void setCoefficient(double coefficient) {
        this.coefficient = coefficient;
    }

    public double getCoefficient() {
        return coefficient;
    }

    public String toString() {
        return String.valueOf(coefficient);
    }
}
