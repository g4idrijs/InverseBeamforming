using System;
using System.Collections.Generic;

namespace InverseBeamforming
{
	public partial class Filter
	{
		public static class NewIIRFilter
		{
			
			static double t_DS_Gold, xi_DS_Gold;
			static double IC_DS_Gold00 = 0, IC_DS_Gold01 = 0;
			static double IC_DS_Gold10 = 0, IC_DS_Gold11 = 0;
			static double IC_DS_Gold20 = 0, IC_DS_Gold21 = 0;
			static double IC_DS_Gold30 = 0, IC_DS_Gold31 = 0;
			static double IC_DS_Gold40 = 0, IC_DS_Gold41 = 0;
			static double IC_DS_Gold50 = 0, IC_DS_Gold51 = 0;
			static double IC_DS_Gold60 = 0, IC_DS_Gold61 = 0;

			/// <summary>
			/// Filters the input array with the DS Gold filter
			/// </summary>
			/// <param name="x">Array to filter</param>
			/// <returns>Filtered signal.</returns>
			/// <remarks>This function first filters the input array forwards, then filters the result backwards. This produces a zero phase filter</remarks>
			public static double[] DS_GoldFilter(double[] x)
			{
				int i = 0;
				bool reverse = false;
				int N = x.Length;
				double[] y = new double[N];
				while(i > -1)
				{
					xi_DS_Gold = x[i];
					// Stage 0
					t_DS_Gold = (0.00070413321580150461 * xi_DS_Gold) - (-1.7310923922594577 * IC_DS_Gold00) - (0.9996862587110722 * IC_DS_Gold01);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold00) + (-1 * IC_DS_Gold01);
					IC_DS_Gold01 = IC_DS_Gold00;
					IC_DS_Gold00 = t_DS_Gold;
					// Stage 1
					t_DS_Gold = (0.00070413321580150461 * xi_DS_Gold) - (-1.7324659955166632 * IC_DS_Gold10) - (0.99968700391354992 * IC_DS_Gold11);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold10) + (-1 * IC_DS_Gold11);
					IC_DS_Gold11 = IC_DS_Gold10;
					IC_DS_Gold10 = t_DS_Gold;
					// Stage 2
					t_DS_Gold = (0.00070393452932227517 * xi_DS_Gold) - (-1.7307395977902 * IC_DS_Gold20) - (0.99912137137990831 * IC_DS_Gold21);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold20) + (-1 * IC_DS_Gold21);
					IC_DS_Gold21 = IC_DS_Gold20;
					IC_DS_Gold20 = t_DS_Gold;
					// Stage 3
					t_DS_Gold = (0.00070393452932227517 * xi_DS_Gold) - (-1.7318417632443899 * IC_DS_Gold30) - (0.99912304488606463 * IC_DS_Gold31);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold30) + (-1 * IC_DS_Gold31);
					IC_DS_Gold31 = IC_DS_Gold30;
					IC_DS_Gold30 = t_DS_Gold;
					// Stage 4
					t_DS_Gold = (0.00070379713244000128 * xi_DS_Gold) - (-1.7306469944097309 * IC_DS_Gold40) - (0.9987311306864024 * IC_DS_Gold41);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold40) + (-1 * IC_DS_Gold41);
					IC_DS_Gold41 = IC_DS_Gold40;
					IC_DS_Gold40 = t_DS_Gold;
					// Stage 5
					t_DS_Gold = (0.00070379713244000128 * xi_DS_Gold) - (-1.7312588884159359 * IC_DS_Gold50) - (0.99873247220986539 * IC_DS_Gold51);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold50) + (-1 * IC_DS_Gold51);
					IC_DS_Gold51 = IC_DS_Gold50;
					IC_DS_Gold50 = t_DS_Gold;
					// Stage 6
					t_DS_Gold = (0.00070374811553745875 * xi_DS_Gold) - (-1.7308324675544136 * IC_DS_Gold60) - (0.99859250376892517 * IC_DS_Gold61);
					xi_DS_Gold = (1 * t_DS_Gold) + (0 * IC_DS_Gold60) + (-1 * IC_DS_Gold61);
					IC_DS_Gold61 = IC_DS_Gold60;
					IC_DS_Gold60 = t_DS_Gold;
					y[i] = xi_DS_Gold;
					if (reverse)
					{
						i--;
					}
					else
					{
						i++;
					}
					if (i == N-1)
					{
						reverse=true;
					}
				}
				return y;
			}

			static double t_RF_Gold, xi_RF_Gold;
			static double IC_RF_Gold00 = 0, IC_RF_Gold01 = 0;
			static double IC_RF_Gold10 = 0, IC_RF_Gold11 = 0;
			static double IC_RF_Gold20 = 0, IC_RF_Gold21 = 0;
			static double IC_RF_Gold30 = 0, IC_RF_Gold31 = 0;
			static double IC_RF_Gold40 = 0, IC_RF_Gold41 = 0;
			static double IC_RF_Gold50 = 0, IC_RF_Gold51 = 0;
			static double IC_RF_Gold60 = 0, IC_RF_Gold61 = 0;

			/// <summary>
			/// Filters the input array with the RF Gold filter
			/// </summary>
			/// <param name="x">Array to filter</param>
			/// <returns>Filtered signal.</returns>
			/// <remarks>This function first filters the input array forwards, then filters the result backwards. This produces a zero phase filter</remarks>
			public static double[] RF_GoldFilter(double[] x)
			{
				int i = 0;
				bool reverse = false;
				int N = x.Length;
				double[] y = new double[N];
				while(i > -1)
				{
					xi_RF_Gold = x[i];
					// Stage 0
					t_RF_Gold = (0.021725833739074118 * xi_RF_Gold) - (-1.7019811484262308 * IC_RF_Gold00) - (0.98997846909011589 * IC_RF_Gold01);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold00) + (-1 * IC_RF_Gold01);
					IC_RF_Gold01 = IC_RF_Gold00;
					IC_RF_Gold00 = t_RF_Gold;
					// Stage 1
					t_RF_Gold = (0.021725833739074118 * xi_RF_Gold) - (-1.7449473174581205 * IC_RF_Gold10) - (0.99068822751544738 * IC_RF_Gold11);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold10) + (-1 * IC_RF_Gold11);
					IC_RF_Gold11 = IC_RF_Gold10;
					IC_RF_Gold10 = t_RF_Gold;
					// Stage 2
					t_RF_Gold = (0.021539942426096515 * xi_RF_Gold) - (-1.6913419189584451 * IC_RF_Gold20) - (0.97236003479660005 * IC_RF_Gold21);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold20) + (-1 * IC_RF_Gold21);
					IC_RF_Gold21 = IC_RF_Gold20;
					IC_RF_Gold20 = t_RF_Gold;
					// Stage 3
					t_RF_Gold = (0.021539942426096515 * xi_RF_Gold) - (-1.7263759324434287 * IC_RF_Gold30) - (0.97392868802171106 * IC_RF_Gold31);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold30) + (-1 * IC_RF_Gold31);
					IC_RF_Gold31 = IC_RF_Gold30;
					IC_RF_Gold30 = t_RF_Gold;
					// Stage 4
					t_RF_Gold = (0.021414060129567619 * xi_RF_Gold) - (-1.6892102394364403 * IC_RF_Gold40) - (0.96079297678572728 * IC_RF_Gold41);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold40) + (-1 * IC_RF_Gold41);
					IC_RF_Gold41 = IC_RF_Gold40;
					IC_RF_Gold40 = t_RF_Gold;
					// Stage 5
					t_RF_Gold = (0.021414060129567619 * xi_RF_Gold) - (-1.7088692439197495 * IC_RF_Gold50) - (0.96203728461035953 * IC_RF_Gold51);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold50) + (-1 * IC_RF_Gold51);
					IC_RF_Gold51 = IC_RF_Gold50;
					IC_RF_Gold50 = t_RF_Gold;
					// Stage 6
					t_RF_Gold = (0.021369664771290865 * xi_RF_Gold) - (-1.6955905028879716 * IC_RF_Gold60) - (0.95726067045741825 * IC_RF_Gold61);
					xi_RF_Gold = (1 * t_RF_Gold) + (0 * IC_RF_Gold60) + (-1 * IC_RF_Gold61);
					IC_RF_Gold61 = IC_RF_Gold60;
					IC_RF_Gold60 = t_RF_Gold;
					y[i] = xi_RF_Gold;
					if (reverse)
					{
						i--;
					}
					else
					{
						i++;
					}
					if (i == N-1)
					{
						reverse=true;
					}
				}
				return y;
			}

			static double t_DS_Hadamard, xi_DS_Hadamard;
			static double IC_DS_Hadamard00 = 0, IC_DS_Hadamard01 = 0;
			static double IC_DS_Hadamard10 = 0, IC_DS_Hadamard11 = 0;
			static double IC_DS_Hadamard20 = 0, IC_DS_Hadamard21 = 0;
			static double IC_DS_Hadamard30 = 0, IC_DS_Hadamard31 = 0;
			static double IC_DS_Hadamard40 = 0, IC_DS_Hadamard41 = 0;
			static double IC_DS_Hadamard50 = 0, IC_DS_Hadamard51 = 0;
			static double IC_DS_Hadamard60 = 0, IC_DS_Hadamard61 = 0;

			/// <summary>
			/// Filters the input array with the DS Hadamard filter
			/// </summary>
			/// <param name="x">Array to filter</param>
			/// <returns>Filtered signal.</returns>
			/// <remarks>This function first filters the input array forwards, then filters the result backwards. This produces a zero phase filter</remarks>
			public static double[] DS_HadamardFilter(double[] x)
			{
				int i = 0;
				bool reverse = false;
				int N = x.Length;
				double[] y = new double[N];
				while(i > -1)
				{
					xi_DS_Hadamard = x[i];
					// Stage 0
					t_DS_Hadamard = (0.00069918304950063398 * xi_DS_Hadamard) - (-1.7310991339141522 * IC_DS_Hadamard00) - (0.99968846696256808 * IC_DS_Hadamard01);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard00) + (-1 * IC_DS_Hadamard01);
					IC_DS_Hadamard01 = IC_DS_Hadamard00;
					IC_DS_Hadamard00 = t_DS_Hadamard;
					// Stage 1
					t_DS_Hadamard = (0.00069918304950063398 * xi_DS_Hadamard) - (-1.7324630760212818 * IC_DS_Hadamard10) - (0.99968920172147402 * IC_DS_Hadamard11);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard10) + (-1 * IC_DS_Hadamard11);
					IC_DS_Hadamard11 = IC_DS_Hadamard10;
					IC_DS_Hadamard10 = t_DS_Hadamard;
					// Stage 2
					t_DS_Hadamard = (0.00069898714599631474 * xi_DS_Hadamard) - (-1.7307488178286681 * IC_DS_Hadamard20) - (0.99912755238269946 * IC_DS_Hadamard21);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard20) + (-1 * IC_DS_Hadamard21);
					IC_DS_Hadamard21 = IC_DS_Hadamard20;
					IC_DS_Hadamard20 = t_DS_Hadamard;
					// Stage 3
					t_DS_Hadamard = (0.00069898714599631474 * xi_DS_Hadamard) - (-1.731843226949493 * IC_DS_Hadamard30) - (0.99912920244988612 * IC_DS_Hadamard31);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard30) + (-1 * IC_DS_Hadamard31);
					IC_DS_Hadamard31 = IC_DS_Hadamard30;
					IC_DS_Hadamard30 = t_DS_Hadamard;
					// Stage 4
					t_DS_Hadamard = (0.00069885167293218199 * xi_DS_Hadamard) - (-1.7306568598311589 * IC_DS_Hadamard40) - (0.99874005148952527 * IC_DS_Hadamard41);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard40) + (-1 * IC_DS_Hadamard41);
					IC_DS_Hadamard41 = IC_DS_Hadamard40;
					IC_DS_Hadamard40 = t_DS_Hadamard;
					// Stage 5
					t_DS_Hadamard = (0.00069885167293218199 * xi_DS_Hadamard) - (-1.7312644460704316 * IC_DS_Hadamard50) - (0.99874137422446563 * IC_DS_Hadamard51);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard50) + (-1 * IC_DS_Hadamard51);
					IC_DS_Hadamard51 = IC_DS_Hadamard50;
					IC_DS_Hadamard50 = t_DS_Hadamard;
					// Stage 6
					t_DS_Hadamard = (0.00069880334222482043 * xi_DS_Hadamard) - (-1.7308410239235303 * IC_DS_Hadamard60) - (0.99860239331555034 * IC_DS_Hadamard61);
					xi_DS_Hadamard = (1 * t_DS_Hadamard) + (0 * IC_DS_Hadamard60) + (-1 * IC_DS_Hadamard61);
					IC_DS_Hadamard61 = IC_DS_Hadamard60;
					IC_DS_Hadamard60 = t_DS_Hadamard;
					y[i] = xi_DS_Hadamard;
					if (reverse)
					{
						i--;
					}
					else
					{
						i++;
					}
					if (i == N-1)
					{
						reverse=true;
					}
				}
				return y;
			}

			static double t_RF_Hadamard, xi_RF_Hadamard;
			static double IC_RF_Hadamard00 = 0, IC_RF_Hadamard01 = 0;
			static double IC_RF_Hadamard10 = 0, IC_RF_Hadamard11 = 0;
			static double IC_RF_Hadamard20 = 0, IC_RF_Hadamard21 = 0;
			static double IC_RF_Hadamard30 = 0, IC_RF_Hadamard31 = 0;
			static double IC_RF_Hadamard40 = 0, IC_RF_Hadamard41 = 0;
			static double IC_RF_Hadamard50 = 0, IC_RF_Hadamard51 = 0;
			static double IC_RF_Hadamard60 = 0, IC_RF_Hadamard61 = 0;

			/// <summary>
			/// Filters the input array with the RF Hadamard filter
			/// </summary>
			/// <param name="x">Array to filter</param>
			/// <returns>Filtered signal.</returns>
			/// <remarks>This function first filters the input array forwards, then filters the result backwards. This produces a zero phase filter</remarks>
			public static double[] RF_HadamardFilter(double[] x)
			{
				int i = 0;
				bool reverse = false;
				int N = x.Length;
				double[] y = new double[N];
				while(i > -1)
				{
					xi_RF_Hadamard = x[i];
					// Stage 0
					t_RF_Hadamard = (0.022266280124460171 * xi_RF_Hadamard) - (-1.7012202640216549 * IC_RF_Hadamard00) - (0.98972023876610726 * IC_RF_Hadamard01);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard00) + (-1 * IC_RF_Hadamard01);
					IC_RF_Hadamard01 = IC_RF_Hadamard00;
					IC_RF_Hadamard00 = t_RF_Hadamard;
					// Stage 1
					t_RF_Hadamard = (0.022266280124460171 * xi_RF_Hadamard) - (-1.7452703745995208 * IC_RF_Hadamard10) - (0.99046576476670678 * IC_RF_Hadamard11);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard10) + (-1 * IC_RF_Hadamard11);
					IC_RF_Hadamard11 = IC_RF_Hadamard10;
					IC_RF_Hadamard10 = t_RF_Hadamard;
					// Stage 2
					t_RF_Hadamard = (0.022071113020399803 * xi_RF_Hadamard) - (-1.6903230217699421 * IC_RF_Hadamard20) - (0.97165880295864016 * IC_RF_Hadamard21);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard20) + (-1 * IC_RF_Hadamard21);
					IC_RF_Hadamard21 = IC_RF_Hadamard20;
					IC_RF_Hadamard20 = t_RF_Hadamard;
					// Stage 3
					t_RF_Hadamard = (0.022071113020399803 * xi_RF_Hadamard) - (-1.7262550793803035 * IC_RF_Hadamard30) - (0.97330586662077301 * IC_RF_Hadamard31);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard30) + (-1 * IC_RF_Hadamard31);
					IC_RF_Hadamard31 = IC_RF_Hadamard30;
					IC_RF_Hadamard30 = t_RF_Hadamard;
					// Stage 4
					t_RF_Hadamard = (0.021939019125336238 * xi_RF_Hadamard) - (-1.6881569639409681 * IC_RF_Hadamard40) - (0.95981624579774649 * IC_RF_Hadamard41);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard40) + (-1 * IC_RF_Hadamard41);
					IC_RF_Hadamard41 = IC_RF_Hadamard40;
					IC_RF_Hadamard40 = t_RF_Hadamard;
					// Stage 5
					t_RF_Hadamard = (0.021939019125336238 * xi_RF_Hadamard) - (-1.7083250317154333 * IC_RF_Hadamard50) - (0.96112243501656613 * IC_RF_Hadamard51);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard50) + (-1 * IC_RF_Hadamard51);
					IC_RF_Hadamard51 = IC_RF_Hadamard50;
					IC_RF_Hadamard50 = t_RF_Hadamard;
					// Stage 6
					t_RF_Hadamard = (0.021892446462511397 * xi_RF_Hadamard) - (-1.6947127137144216 * IC_RF_Hadamard60) - (0.95621510707497714 * IC_RF_Hadamard61);
					xi_RF_Hadamard = (1 * t_RF_Hadamard) + (0 * IC_RF_Hadamard60) + (-1 * IC_RF_Hadamard61);
					IC_RF_Hadamard61 = IC_RF_Hadamard60;
					IC_RF_Hadamard60 = t_RF_Hadamard;
					y[i] = xi_RF_Hadamard;
					if (reverse)
					{
						i--;
					}
					else
					{
						i++;
					}
					if (i == N-1)
					{
						reverse=true;
					}
				}
				return y;
			}

		}

	}
}
