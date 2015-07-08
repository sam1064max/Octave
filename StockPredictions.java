import java.util.Scanner;

class StockPredictions {
public static void main(String[] args) {
        Scanner in = new Scanner(System.in);
        
        double _m;
        _m = in.nextDouble();
        
        int _k;
        _k = in.nextInt();
        
        int _d;
        _d = in.nextInt();
        
        String[] _name = new String[_k];
        int[] _owned = new int[_k];
        double[][] _prices = new double[_k][5];
        
        String _name_item;
        int _owned_item;
        double _prices_item_item;
        
        for(int _i = 0; _i < _k; _i++) {
            _name_item = in.next();
            _name[_i] = _name_item;
            
            _owned_item = in.nextInt();
            _owned[_i] = _owned_item;
            
            for(int _j = 0; _j<5; _j++) {
                _prices_item_item = in.nextDouble();
                _prices[_i][_j] = _prices_item_item;                
            }
        }
        
        printTransactions(_m, _k, _d, _name, _owned, _prices);
        
    }
	static void printTransactions(double m, int k, int d, String[] name, int[] owned, double[][] prices) {
		
		for(int index=0;index<name.length;index++)
		{
			
			int m1 = prices[index].length;
			int n = 1;
			double x[][] = new double[m1][n + 1];
			double y[][] = new double[m1][1];
			for (int i = 0; i < m1; i++) {
				x[i][0] = 1;
				for (int j = 1; j <= n; j++)
					x[i][j] = Math.pow(i,j);
					y[i][0] = prices[index][i];
			}

			double ans[][] = regression(x, y);
			int p = 1;
			double z[][] = new double[p][n + 1];
			z[0][0] = 1;
			for (int j = 1; j <= n; j++)
				z[0][j] = Math.pow(m1,j);

			double finalans[][] = mul(z, ans);
			
			int tra = 0;
			int sell[] = new int[m1];
			int buy = Integer.MAX_VALUE;
			double big=0;
			for (int i=0;i<m1;i++) {
				if(owned[index]>0 && prices[index][m1-1]-finalans[0][0]<0)
				{
					sell[tra++]=i;
				}
				else if(big<prices[index][m1-1]-finalans[0][0])
				{
					big = prices[index][m1-1]-finalans[0][0];
					buy = i;
				}
			}
			System.out.println(tra);
			for (int i=0;i<tra;i++)
				if(sell[i]<buy)
					System.out.println(name[sell[i]]+" "+"SELL "+owned[sell[i]]);
				else
				{
					System.out.println(name[buy]+" "+"BUY "+owned[buy]);
					buy = Integer.MAX_VALUE;
				}
			print(finalans);
		}		
	}
	

	public static void print(double a[][]) {
		for (int i = 0; i < a.length; ++i) {
			for (int j = 0; j < a[0].length; ++j) {
					System.out.print(a[i][j] + "  ");
			}
			System.out.println();
		}
	}
	
	public static double[][] regression(double x[][],double y[][]) {
		double t[][] = transpose(x);
		double mul[][] = mul(t, x);
		double theta[][] = mul(invert(mul), t);
		return mul(theta, y);
	}
	public static double[][] invert(double a[][]) {
		int n = a.length;
		double x[][] = new double[n][n];
		double b[][] = new double[n][n];
		int index[] = new int[n];
		for (int i = 0; i < n; ++i)
			b[i][i] = 1;

		// Transform the matrix into an upper triangle
		gaussian(a, index);

		// Update the matrix b[i][j] with the ratios stored
		for (int i = 0; i < n - 1; ++i)
			for (int j = i + 1; j < n; ++j)
				for (int k = 0; k < n; ++k)
					b[index[j]][k] -= a[index[j]][i] * b[index[i]][k];

		// Perform backward substitutions
		for (int i = 0; i < n; ++i) {
			x[n - 1][i] = b[index[n - 1]][i] / a[index[n - 1]][n - 1];
			for (int j = n - 2; j >= 0; --j) {
				x[j][i] = b[index[j]][i];
				for (int k = j + 1; k < n; ++k) {
					x[j][i] -= a[index[j]][k] * x[k][i];
				}
				x[j][i] /= a[index[j]][j];
			}
		}
		return x;
	}

	// Method to carry out the partial-pivoting Gaussian
	// elimination. Here index[] stores pivoting order.

	public static void gaussian(double a[][], int index[]) {
		int n = index.length;
		double c[] = new double[n];

		// Initialize the index
		for (int i = 0; i < n; ++i)
			index[i] = i;

		// Find the rescaling factors, one from each row
		for (int i = 0; i < n; ++i) {
			double c1 = 0;
			for (int j = 0; j < n; ++j) {
				double c0 = Math.abs(a[i][j]);
				if (c0 > c1)
					c1 = c0;
			}
			c[i] = c1;
		}

		// Search the pivoting element from each column
		int k = 0;
		for (int j = 0; j < n - 1; ++j) {
			double pi1 = 0;
			for (int i = j; i < n; ++i) {
				double pi0 = Math.abs(a[index[i]][j]);
				pi0 /= c[index[i]];
				if (pi0 > pi1) {
					pi1 = pi0;
					k = i;
				}
			}

			// Interchange rows according to the pivoting order
			int itmp = index[j];
			index[j] = index[k];
			index[k] = itmp;
			for (int i = j + 1; i < n; ++i) {
				double pj = a[index[i]][j] / a[index[j]][j];

				// Record pivoting ratios below the diagonal
				a[index[i]][j] = pj;

				// Modify other elements accordingly
				for (int l = j + 1; l < n; ++l)
					a[index[i]][l] -= pj * a[index[j]][l];
			}
		}
	}

	public static double[][] transpose(double a[][]) {
		double transpose[][] = new double[a[0].length][a.length];

		for (int c = 0; c < a.length; c++) {
			for (int d = 0; d < a[0].length; d++)
				transpose[d][c] = a[c][d];
		}
		return transpose;
	}

	public static double[][] mul(double first[][], double second[][]) {
		double sum = 0;
		double multiply[][] = new double[first.length][second[0].length];
		for (int c = 0; c < first.length; c++) {
			for (int d = 0; d < second[0].length; d++) {
				for (int k = 0; k < first[0].length; k++) {
					sum = sum + first[c][k] * second[k][d];
				}

				multiply[c][d] = sum;
				sum = 0;
			}
		}
		return multiply;
	}
}
