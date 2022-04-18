import java.io.*;
import java.util.Random;

public class Runner {
    public static void main(String[] args) {
        int N = 1000;
        int R = 350;
        double sigma = 0.7;
        double mu = 1;
        double u_th = 0.98;
        double dt = 0.01;

        double[] u = new double[N];
        double[] circles = new double[N];

        Random rand = new Random();

        for (int i = 0; i < N; i++) {
            double r = rand.nextInt(3000) % 10;
            u[i] = r / 10.0;
            circles[i] = 0;
        }
        
        double t = 0;
        double time_max = 300;

        while (t < time_max) {
            for (int i=0; i<N; i++) {
                u[i] += dt*(mu - u[i]);
                double coupling = 0;
                if (i-R<0){
                    for (int j=((N + ((i-R)%N)) % N); j<N; j++) {
                        coupling += sigma*(u[i]-u[j])/(2*R);                  
                    }
                    for (int j=0; j<i+R+1; j++) {
                        coupling += sigma*(u[i]-u[j])/(2*R); 
                    }
                    u[i] += dt*coupling;
                }
                else if (i+R > N-1) {
                    for (int j=i-R; j<N; j++) {
                        coupling += sigma*(u[i]-u[j])/(2*R);  
                    }
                    for (int j=0; j<((N + ((i+R)%N)) % N + 1); j++) {
                        coupling += sigma*(u[i]-u[j])/(2*R);  
                    }
                    u[i] += dt*coupling;
                }
                else {
                    for (int j=i-R; j<i+R+1; j++){
                        coupling += sigma*(u[i]-u[j])/(2*R); 
                    }
                    u[i] += dt*coupling;
                }

                if (u[i] > u_th) {
                    u[i] = 0;
                    circles[i] += 1;
                }
            }
            System.out.println(t);
            t += dt;

            try {
                PrintWriter writer = new PrintWriter("u.txt");
                for(int k=0; k<N; k++)
                {
                    writer.println(u[k]); // output to txt
                }
                writer.close();
            } catch (FileNotFoundException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
            

        }
    }
}
