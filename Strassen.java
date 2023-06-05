import java.io.*;
import java.lang.*;

class Strassen {
    static Matrix standard(Matrix A, Matrix B) {
        int d = A.getD();
        Matrix C = new Matrix(d);
        for (int i = 0; i < d; i++) {
            for (int k = 0; k < d; k++) {
                int sum = 0;
                for (int j = 0; j < d; j++) {
                    sum += A.getVal(i, j) * B.getVal(j, k);
                }
                C.load(i, k, sum);
            }
        }
        return C;
    }

    /*Matrix strassen(Matrix X, Matrix Y) {
        int div = d / 2;
        int n0 = 0;
        Matrix A, B, C, D, E, F, G, H, P1, P2, P3, P4, P5, P6, P6 = new Matrix(div);
        for (int i = 0; i < div; i++) {
            for (int j = 0; j < div; j++) {
                rLoad(X, A, i, j);
                rLoad(Y, E, i, j)
            }
            for (int j = div; j < d; j++) {
                rLoad(X, B, i, j);
                rLoad(Y, F, i, j);
            }
        }
        for (int i = div; i < d; i++) {
            for (int j = 0; j < div; j++) {
                rLoad(X, C, i, j);
                rLoad(Y, G, i, j)
            }
            for (int j = div; j < d; j++) {
                rLoad(X, D, i, j);
                rLoad(Y, H, i, j);
            }
        }
    }*/

    // Maybe? Better memory
    // mults: [A, FmH, ApB, H, CpD, E, D, GmE, ApD, EpH, BmD, GpH, CmA, EpF]
    // prods: [P1, P2, P3, P4, P5, P6, P6, P7]
    // Done: A, E, ApB, EpF, CmA, GmE, CpD, ApD, D, BmD, H, FmH, EpH, GpH
    static Matrix strassen(Matrix X, Matrix Y) {
        int n0 = 273;
        int d = X.getD();
        if (d <= n0) {
            return standard(X, Y);
        }
        int div;
        boolean oddFlag = false;
        // d is odd
        if (d % 2 != 0) {
            oddFlag = true;
            d++;
            div = d / 2;
            Matrix Xtmp = new Matrix(d);
            Matrix Ytmp = new Matrix(d);
            for (int i = 0; i < d - 1; i++) {
                for (int j = 0; j < d - 1; j++) {
                    Xtmp.load(i, j, X.getVal(i, j));
                    Ytmp.load(i, j, Y.getVal(i, j));
                }
                Xtmp.load(i, d - 1, 0);
                Ytmp.load(i, d - 1, 0);
            }
            for (int j = 0; j < d; j++) {
                Xtmp.load(d - 1, j, 0);
                Ytmp.load(d - 1, j, 0);
            }
            X = Xtmp;
            Y = Ytmp;
        } else {
            div = d / 2;
        }
        Matrix[] mults = new Matrix[14];
        for (int i = 0; i < 14; i++) {
            mults[i] = new Matrix(div);
        }
        for (int i = 0; i < div; i++) {
            for (int j = 0; j < div; j++) {
                // A
                rLoad(X, mults[0], i, j, i, j);
                // E
                rLoad(Y, mults[5], i, j, i, j);
                // ApB first
                rLoad(X, mults[2], i, j, i, j);
                // ApD first
                rLoad(X, mults[8], i, j, i, j);
                // EpF first
                rLoad(Y, mults[13], i, j, i, j);
                // EpH first
                rLoad(Y, mults[9], i, j, i, j);
            }
            for (int j = div; j < d; j++) {
                int adjj = j - div;
                // ApB second
                addLoad(X, mults[2], i, j, i, adjj, mults[2].getVal(i, adjj));
                // BmD first
                rLoad(X, mults[10], i, j, i, adjj);
                // FmH first
                rLoad(Y, mults[1], i, j, i, adjj);
                // EpF second
                addLoad(Y, mults[13], i, j, i, adjj, mults[13].getVal(i, adjj));
            }
        }
        for (int i = div; i < d; i++) {
            for (int j = 0; j < div; j++) {
                int adji = i - div;
                // CpD first
                rLoad(X, mults[4], i, j, adji, j);
                // CmA
                subToLoad(X, mults[12], i, j, adji, j, mults[0].getVal(adji, j));
                // GmE
                subToLoad(Y, mults[7], i, j, adji, j, mults[5].getVal(adji, j));
                // GpH first
                rLoad(Y, mults[11], i, j, adji, j);
            }
            for (int j = div; j < d; j++) {
                int adji = i - div;
                int adjj = j - div;
                // D
                rLoad(X, mults[6], i, j, adji, adjj);
                // ApD second
                addLoad(X, mults[8], i, j, adji, adjj, mults[8].getVal(adji, adjj));
                // CpD second
                addLoad(X, mults[4], i, j, adji, adjj, mults[4].getVal(adji, adjj));
                // BmD second
                subFromLoad(X, mults[10], i, j, adji, adjj, mults[10].getVal(adji, adjj));
                // H
                rLoad(Y, mults[3], i, j, adji, adjj);
                // FmH second
                subFromLoad(Y, mults[1], i, j, adji, adjj, mults[1].getVal(adji, adjj));
                // EpH second
                addLoad(Y, mults[9], i, j, adji, adjj, mults[9].getVal(adji, adjj));
                // GpH second
                addLoad(Y, mults[11], i, j, adji, adjj, mults[11].getVal(adji, adjj));
            }
        }
        // A becomes P1
        mults[0] = strassen(mults[0], mults[1]);
        // ApB becomes P2
        mults[2] = strassen(mults[2], mults[3]);
        // CpD becomes P3
        mults[4] = strassen(mults[4], mults[5]);
        // D becomes P4
        mults[6] = strassen(mults[6], mults[7]);
        // ApD becomes P5
        mults[8] = strassen(mults[8], mults[9]);
        // BmD becomes P6
        mults[10] = strassen(mults[10], mults[11]);
        // CmA becomes P7
        mults[12] = strassen(mults[12], mults[13]);

        for (int i = 0; i < div; i++) {
            for (int j = 0; j < div; j++) {
                // P6 becomes P5 + P6
                addLoad(mults[8], mults[10], i, j, i, j, mults[10].getVal(i, j));
                // P5 + P6 becomes P5 + P6 + P4
                addLoad(mults[6], mults[10], i, j, i, j, mults[10].getVal(i, j));
                // P4 + P5 + P6 becomes AE + BG (P4 + P5 + P6 - P2) in mults[10]
                subFromLoad(mults[2], mults[10], i, j, i, j, mults[10].getVal(i, j));

                // P7 becomes P7 + P5
                addLoad(mults[8], mults[12], i, j, i, j, mults[12].getVal(i, j));
                // P7 + P5 becomes P7 + P5 + P1
                addLoad(mults[0], mults[12], i, j, i, j, mults[12].getVal(i, j));
                // P7 + P5 + P1 becomes CF + DH (P7 + P5 + P1 - P3)
                subFromLoad(mults[4], mults[12], i, j, i, j, mults[12].getVal(i, j));
            }
        }
        // Must be separate to not taint other additions.
        for (int i = 0; i < div; i++) {
            for (int j = 0; j < div; j++) {
                // P3 becomes CE + DG (P3 + P4)
                addLoad(mults[6], mults[4], i, j, i, j, mults[4].getVal(i, j));
                // P1 becomes AF + BH (P1 + P2)
                addLoad(mults[2], mults[0], i, j, i, j, mults[0].getVal(i, j));
            }
        }
        // Now, mults = [AF + BH, FmH, P2, H, CE + DG, E, P4, GmE, P5, EpH, AE + BG, GpH, CF + DH, EpF]
        Matrix Prod = new Matrix(d);
        for (int i = 0; i < div; i++) {
            for (int j = 0; j < div; j++) {
                rLoad(mults[10], Prod, i, j, i, j);
            }
            for (int j = div; j < d; j++) {
                int adjj = j - div;
                rLoad(mults[0], Prod, i, adjj, i, j);
            }
        }
        for (int i = div; i < d; i++) {
            for (int j = 0; j < div; j++) {
                int adji = i - div;
                rLoad(mults[4], Prod, adji, j, i, j);
            }
            for (int j = div; j < d; j++) {
                int adji = i - div;
                int adjj = j - div;
                rLoad(mults[12], Prod, adji, adjj, i, j);
            }
        }
        if (oddFlag) {
            Matrix prodTmp = new Matrix(d - 1);
            for (int i = 0; i < d - 1; i++) {
                for (int j = 0; j < d - 1; j++) {
                    prodTmp.load(i, j, Prod.getVal(i, j));
                }
            }
            Prod = prodTmp;
        }
        return Prod;
    }

    static void rLoad(Matrix in, Matrix out, int ini, int inj, int outi, int outj) {
        out.load(outi, outj, in.getVal(ini, inj));
    }

    static void addLoad(Matrix in, Matrix out, int ini, int inj, int outi, int outj, int addTo) {
        out.load(outi, outj, in.getVal(ini, inj) + addTo);
    }

    static void subFromLoad(Matrix in, Matrix out, int ini, int inj, int outi, int outj, int subFrom) {
        out.load(outi, outj, subFrom - in.getVal(ini, inj));
    }

    static void subToLoad(Matrix in, Matrix out, int ini, int inj, int outi, int outj, int sub) {
        out.load(outi, outj, in.getVal(ini, inj) - sub);
    }

    // Figure out how to read data from input ascii file in args[2]
    public static void main(String[] args) {
        int flag = Integer.parseInt(args[0]);
        int d = Integer.parseInt(args[1]);
        if (flag == 0) {
            File inputfile = new File(args[2]);
            try {
                BufferedReader r = new BufferedReader(new FileReader(inputfile));
                int val;
                // Check if d is power of 2
                /*double base = Math.log(d) / Math.log(2);
                if (base != (int) base) {
                    int nextPow = d + 1;
                    double nextBase = Math.log(nextPow) / Math.log(2);
                    while (nextBase != (int) nextBase) {
                        nextPow++;
                        nextBase = Math.log(nextPow) / Math.log(2);
                    }
                    Matrix A = new Matrix(nextPow);
                    Matrix B = new Matrix(nextPow);
                    for (int i = 0; i < d; i++) {
                        for (int j = 0; j < d; j++) {
                            val = Integer.parseInt(r.readLine());
                            A.load(i, j, val);
                        }
                        for (int j = d; j < nextPow; j++) {
                            A.load(i, j, 0);
                        }
                    }
                    for (int i = d; i < nextPow; i++) {
                        for (int j = 0; j < nextPow; j++) {
                            A.load(i, j, 0);
                        }
                    }
                    for (int i = 0; i < d; i++) {
                        for (int j = 0; j < d; j++) {
                            val = Integer.parseInt(r.readLine());
                            B.load(i, j, val);
                        }
                        for (int j = d; j < nextPow; j++) {
                            B.load(i, j, 0);
                        }
                    }
                    for (int i = d; i < nextPow; i++) {
                        for (int j = 0; j < nextPow; j++) {
                            B.load(i, j, 0);
                        }
                    }
                    Matrix C = strassen(A, B);
                    *//*long start = System.nanoTime();
                    Matrix C = strassen(A, B);
                    long finish = System.nanoTime();
                    long elapsed = finish - start;
                    System.out.println(elapsed);*//*
                    for (int i = 0; i < d; i++) {
                        System.out.print(C.getVal(i, i) + "\n");
                    }
                } else {*/
                Matrix A = new Matrix(d);
                Matrix B = new Matrix(d);
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        val = Integer.parseInt(r.readLine());
                        A.load(i, j, val);
                    }
                }
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        val = Integer.parseInt(r.readLine());
                        B.load(i, j, val);
                    }
                }
                Matrix C = strassen(A, B);
                    /*long start = System.nanoTime();
                    Matrix C = strassen(A, B);
                    long finish = System.nanoTime();
                    long elapsed = finish - start;
                    // System.out.println(elapsed);*/
                for (int i = 0; i < d; i++) {
                    System.out.println(C.getVal(i, i));
                }
            } catch (FileNotFoundException e) {
                System.out.println("Third arg must be file.");
            } catch (IOException e) {
                System.out.println("File must contain 2d^2 lines of integers.");
            }
        } else if (flag == 1) {
            long elapsed = 0;
            long start;
            long finish;
            // Check if d is even
            /*if ((d & 1) != 0) {
                int nextEven = d + 1;
                Matrix A = new Matrix(nextEven);
                Matrix B = new Matrix(nextEven);
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        if (Math.random() < 0.5) {
                            A.load(i, j, 0);
                        } else {
                            A.load(i, j, 1);
                        }
                    }
                    for (int j = d; j < nextEven; j++) {
                        A.load(i, j, 0);
                    }
                }
                for (int i = d; i < nextEven; i++) {
                    for (int j = 0; j < nextEven; j++) {
                        A.load(i, j, 0);
                    }
                }
                for (int i = 0; i < d; i++) {
                    for (int j = 0; j < d; j++) {
                        if (Math.random() < 0.5) {
                            B.load(i, j, 0);
                        } else {
                            B.load(i, j, 1);
                        }
                    }
                    for (int j = d; j < nextEven; j++) {
                        B.load(i, j, 0);
                    }
                }
                for (int i = d; i < nextEven; i++) {
                    for (int j = 0; j < nextEven; j++) {
                        B.load(i, j, 0);
                    }
                }
                Matrix C = new Matrix(nextEven);
                start = System.nanoTime();
                C = strassen(A, B);
                finish = System.nanoTime();
                elapsed = finish - start;
                for (int i = 0; i < d; i++) {
                    System.out.println(C.getVal(i, i));
                }
            } else {*/
            Matrix A = new Matrix(d);
            Matrix B = new Matrix(d);
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < d; j++) {
                    if (Math.random() < 0.5) {
                        A.load(i, j, 0);
                        // System.out.print(" " + 0);
                    } else {
                        A.load(i, j, 1);
                        // System.out.print(" " + 1);
                    }
                }
                // System.out.println();
            }
            // System.out.println();
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < d; j++) {
                    if (Math.random() < 0.5) {
                        B.load(i, j, 0);
                        // System.out.print(" " + 0);
                    } else {
                        B.load(i, j, 1);
                        // System.out.print(" " + 1);
                    }
                }
                // System.out.println();
            }
            start = System.nanoTime();
            Matrix C = strassen(A, B);
            finish = System.nanoTime();
            elapsed = finish - start;
            // System.out.println();
            /*for (int i = 0; i < d; i++) {
                System.out.println(C.getVal(i, i));
            }*/
            System.out.println((double) elapsed / 1000000000);
        } else if (flag == 2) {
            double p = 0.05;
            Matrix A = new Matrix(d);
            for (int i = 0; i < d; i++) {
                for (int j = 0; j < d; j++) {
                    if (Math.random() < p) {
                        A.load(i, j, 1);
                    } else {
                        A.load(i, j, 0);
                    }
                }
            }
            int triangles = 0;
            Matrix C = strassen(A, A);
            Matrix D = strassen(A, C);
            for (int i = 0; i < d; i++) {
                triangles += D.getVal(i, i);
            }
            System.out.println(triangles / 6);
        }
    }
}
