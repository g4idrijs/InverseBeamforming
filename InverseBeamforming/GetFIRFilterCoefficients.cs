using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public partial class Filter
	{
		/// <summary>
		///Gets the FIR coefficients for the DS filter.
		/// </summary>
		/// <param name="coefs">Array to hold the coefficients of the filter</param>
		public static void getFIRFilterCoefficients_Gold_DS(out double[] coefs)
		{
			coefs = new double[]{-0.015337964032003045,0,0.10353139777807874,0.28733616882443702,0.38344979493380493,0.28733616882443702,0.10353139777807874,0,-0.015337964032003045};
		}

		/// <summary>
		///Gets the FIR coefficients for the RF filter.
		/// </summary>
		/// <param name="coefs">Array to hold the coefficients of the filter</param>
		public static void getFIRFilterCoefficients_Gold_RF(out double[] coefs)
		{
			coefs = new double[]{-0.015313247593339584,0,0.10349943035323542,0.28734107073054704,0.38349798844687483,0.28734107073054704,0.10349943035323542,0,-0.015313247593339584};
		}

		/// <summary>
		///Gets the FIR coefficients for the DS filter.
		/// </summary>
		/// <param name="coefs">Array to hold the coefficients of the filter</param>
		public static void getFIRFilterCoefficients_Hadamard_DS(out double[] coefs)
		{
			coefs = new double[]{-0.01533796439293844,4.387877368209604e-15,0.10353139824470713,0.28733616875281504,0.38344979423029402,0.28733616875281504,0.10353139824470713,4.387877368209604e-15,-0.01533796439293844};
		}

		/// <summary>
		///Gets the FIR coefficients for the RF filter.
		/// </summary>
		/// <param name="coefs">Array to hold the coefficients of the filter</param>
		public static void getFIRFilterCoefficients_Hadamard_RF(out double[] coefs)
		{
			coefs = new double[]{-0.015311995576924479,6.8569637641755706e-17,0.1034978102962285,0.28734131890582154,0.38350043066811251,0.28734131890582154,0.1034978102962285,6.8569637641755706e-17,-0.015311995576924479};
		}
	}
}
