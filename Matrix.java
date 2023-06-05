import java.util.concurrent.ThreadLocalRandom;

public class Matrix {

    private final int d;
    private int[][] m;

    public Matrix(int dim) {
        d = dim;
        m = new int[d][d];
    }

    /*
    public void randLoad(int low, int high) {
        for (int i = 0; i++ i<d){
            for (int j = 0; j++; j < d) {

            }
        }
    }
*/

    public void randLoad(int excHi) {
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                m[i][j] = ThreadLocalRandom.current().nextInt(0, excHi + 1);
            }
        }
    }

    public int getD() {
        return d;
    }

    public void load(int i, int j, int val) {
        m[i][j] = val;
    }

    public int getVal(int i, int j) {
        return m[i][j];
    }
}
