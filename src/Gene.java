public class Gene {
    private double coefficient;

    public Gene(Point point, double coefficient) {
        this.coefficient=coefficient;
    }
    public void setCoefficient(double coefficient) {
        this.coefficient = coefficient;
    }

    public double getCoefficient() {
        return coefficient;
    }
}
