using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using static System.Math;
using static InverseBeamforming.Filter.LowPassPrototypes;
using static InverseBeamforming.Filter.QuadRootsCode;

namespace InverseBeamforming
{
	public partial class Filter
	{
		private static readonly int NUM_SAMPLES = 2048*8;

		/// <summary>
		/// Do an example filter operation
		/// </summary>
		/// <param name="CoeffFilename"></param>
		public static void ExampleIIRCall(string CoeffFilename,string FFTfilename, string IIRFilename)
		{
			int j, N;
			IIR_Filter.TIIRFilterParams IIRFilt;  // Defined in IIRFilterCode.h
			IIR_Filter.TIIRCoeff IIRCoeff=new IIR_Filter.TIIRCoeff(true);

			// This structure must be filled before calling CalcIIRFilterCoeff().
			IIRFilt.IIRPassType = IIR_Filter.TIIRPassTypes.BPF;        // iirLPF, iirHPF, iirBPF, iirNOTCH, iirALLPASS  (defined in IIRFilterCode.h)
			IIRFilt.OmegaC = 0.2;                // 0.0 < OmegaC < 1.0        3 dB freq for low and high pass filters, center freq for band pass and notch filters.
			IIRFilt.BW = 0.1;                    // 0.0 < BandWidth < 1.0     3 dB bandwidth for bandpass and notch filters
			IIRFilt.dBGain = 0.0;                // -60.0 < dBGain < 60.0     All filters

			// Define the low pass prototype. These are range checked in LowPassPrototypes.cpp
			IIRFilt.ProtoType = TFilterPoly.BUTTERWORTH;   // BUTTERWORTH, CHEBYSHEV, GAUSSIAN, BESSEL, ADJUSTABLE, INVERSE_CHEBY, PAPOULIS, ELLIPTIC  (defined in LowPassPrototypes.h)
			IIRFilt.NumPoles = 6;                // 1 <= NumPoles <= 12, 15, 20 Depending on the filter.
			IIRFilt.Ripple = 0.1;                // 0.0 <= Ripple <= 1.0 dB     Chebyshev and Elliptic (less for high order Chebyshev).
			IIRFilt.StopBanddB = 60.0;           // 20 <= StopBand <= 120 dB    Inv Cheby and Elliptic
			IIRFilt.Gamma = 0.0;                 // -1.0 <= Gamma <= 1.0        Adjustable Gauss  Controls the transition BW.


			// This will fill the IIRCoeff struct with the 2nd order IIR coefficients.
			//IIRCoeff = IIR_Filter.CalcIIRFilterCoeff(IIRFilt);
			getIIRCoefficients_Hadamard_RF(ref IIRCoeff);

			// If desired, this will create an Nth order poly from the 2nd order polys in IIRCoeff.
			double[] DenomCoeff = new double[25];
			double[] NumerCoeff = new double[25];
			N = RebuildPoly(IIRCoeff.NumSections, ref DenomCoeff, ref IIRCoeff.a2, ref IIRCoeff.a1, ref IIRCoeff.a0);
			N = RebuildPoly(IIRCoeff.NumSections, ref NumerCoeff, ref IIRCoeff.b2, ref IIRCoeff.b1, ref IIRCoeff.b0);

			// This calculates the frequency response of the filter by doing a DFT of the IIR coefficients.
			double[] RealHofZ = new double[NUM_SAMPLES];   // Real and imag parts of H(z). Used with the function IIRFreqResponse.
			double[] ImagHofZ = new double[NUM_SAMPLES];
			IIR_Filter.IIRFreqResponse(IIRCoeff, IIRCoeff.NumSections, ref RealHofZ, ref ImagHofZ, NUM_SAMPLES);

			// This is an alternative way to calculate the filter's frequency response using the FFT.
			// We send an impulse through the filter, and calc the FFT of the filters output.
			// Since the FFT scales the output of a forward transform by 1/N, we use N = NUM_SAMPLES instead of 1 for the impulse.
			double[] Samples = new double[NUM_SAMPLES];
			for (j = 0; j < NUM_SAMPLES; j++)
				Samples[j] = RealHofZ[j] = ImagHofZ[j] = 0.0;

			Samples[0] = NUM_SAMPLES;                                 // The impulse.
			IIR_Filter.FilterWithIIR(IIRCoeff, ref Samples, ref RealHofZ, NUM_SAMPLES);  // Filter the impulse. RealHofZ is used to store the filtered output.

			using (StreamWriter OutputFile = new StreamWriter(IIRFilename))
			{
				// Print the IIR output to a file
				for (j = 0; j < RealHofZ.Length; j++)
				{
					OutputFile.WriteLine(RealHofZ[j]);
				}
			}

			FFTCode.FFT(ref RealHofZ, ref ImagHofZ, NUM_SAMPLES, FFTCode.TTransFormType.FORWARD);            // The FFT's results are returned in input arrays, RealHofZ and ImagHofZ.

			using (StreamWriter OutputFile = new StreamWriter(FFTfilename))
			{
				// Print the fft results
				for (j = 0; j < RealHofZ.Length; j++)
				{
					OutputFile.WriteLine(String.Format("{0}, {1}", RealHofZ[j],ImagHofZ[j]));
				}
			}

			using (StreamWriter OutputFile = new StreamWriter(CoeffFilename))
			{

				// Print the IIR coefficients to a text file in 3 formats.
				for (j = 0; j < IIRCoeff.NumSections; j++)
				{
					OutputFile.WriteLine(String.Format("\n Section {0}", j));
					OutputFile.WriteLine(String.Format("\n a0= {0}  a1= {1}  a2= {2}", IIRCoeff.a0[j], IIRCoeff.a1[j], IIRCoeff.a2[j]));
					OutputFile.WriteLine(String.Format("\n b0= {0}  b1= {1}  b2= {2} ", IIRCoeff.b0[j], IIRCoeff.b1[j], IIRCoeff.b2[j]));
					OutputFile.WriteLine("");
				}
				for (j = 0; j < IIRCoeff.NumSections; j++)
				{
					OutputFile.WriteLine(String.Format("\n  {0} \n  {1} \n  {2}", IIRCoeff.a0[j], IIRCoeff.a1[j], IIRCoeff.a2[j]));
					OutputFile.WriteLine(String.Format("\n  {0} \n  {1} \n  {2} ", IIRCoeff.b0[j], IIRCoeff.b1[j], IIRCoeff.b2[j]));
					OutputFile.WriteLine( "\n ");
				}

				OutputFile.WriteLine( "\n Nth Order Coeff. \n b's are numerator, a's are denominator.");
				for (j = N; j >= 0; j--)
				{
					OutputFile.WriteLine(String.Format("\n b{0} {1} ", j, NumerCoeff[j]));
				}
				OutputFile.WriteLine( "\n ");
				for (j = N; j >= 0; j--)
				{
					OutputFile.WriteLine(String.Format("\n a{0} {1} ", j, DenomCoeff[j]));
				}

				OutputFile.WriteLine( "\n ");
			} 

		}


		public static class IIR_Filter
		{
			private static readonly double OVERFLOW_LIMIT = 1.0E20;
			private static double MaxRegVal;
			private static double[] RegX1, RegX2, RegY1, RegY2;
			private static bool MessageShown = false;

			public enum TIIRPassTypes { LPF, HPF, BPF, Notch, AllPass };

			public struct TIIRCoeff
			{
				public double[] a0, a1, a2, a3, a4;
				public double[] b0, b1, b2, b3, b4;

				public int NumSections;

				public TIIRCoeff(bool init)
				{
					a0 = new double[ARRAY_DIM];
					a1= new double[ARRAY_DIM];
					a2= new double[ARRAY_DIM];
					a3= new double[ARRAY_DIM];
					a4= new double[ARRAY_DIM];
					b0= new double[ARRAY_DIM];
					b1= new double[ARRAY_DIM];
					b2= new double[ARRAY_DIM];
					b3= new double[ARRAY_DIM];
					b4= new double[ARRAY_DIM];
					NumSections = 0;
				}
			}

			public struct TIIRFilterParams
			{
				public TIIRPassTypes IIRPassType;
				public double OmegaC;
				public double BW;
				public double dBGain;

				public TFilterPoly ProtoType;
				public int NumPoles;
				public double Ripple;
				public double StopBanddB;
				public double Gamma;
			}

			public static TIIRCoeff CalcIIRFilterCoeff(TIIRFilterParams IIRFilt)
			{
				int j, k;
				double Scalar, SectionGain;
				double[] Coeff = new double[5];
				double A, B, C, D, E, F, T, Q, Arg;
				double[] a2 = new double[ARRAY_DIM];
				double[] a1 = new double[ARRAY_DIM];
				double[] a0 = new double[ARRAY_DIM];
				double[] b2 = new double[ARRAY_DIM];
				double[] b1 = new double[ARRAY_DIM];
				double[] b0 = new double[ARRAY_DIM];
				Complex[] Roots = new Complex[5];

				TIIRCoeff IIR=new TIIRCoeff(true);                // Gets returned by this function.
				TLowPassParams LowPassFilt;   // Passed to the CalcLowPassProtoCoeff() function.
				TSPlaneCoeff SPlaneCoeff;     // Filled by the CalcLowPassProtoCoeff() function.


				// We can set the TLowPassParams variables directly from the TIIRFilterParams variables.
				LowPassFilt.ProtoType = IIRFilt.ProtoType;
				LowPassFilt.NumPoles = IIRFilt.NumPoles;
				LowPassFilt.Ripple = IIRFilt.Ripple;
				LowPassFilt.Gamma = IIRFilt.Gamma;
				LowPassFilt.StopBanddB = IIRFilt.StopBanddB;

				// Get the low pass prototype 2nd order s plane coefficients.
				SPlaneCoeff = CalcLowPassProtoCoeff(LowPassFilt);


				// Init the IIR structure.
				for (j = 0; j < ARRAY_DIM; j++)
				{
					IIR.a0[j] = 0.0;
					IIR.b0[j] = 0.0;
					IIR.a1[j] = 0.0;
					IIR.b1[j] = 0.0;
					IIR.a2[j] = 0.0;
					IIR.b2[j] = 0.0;
					IIR.a3[j] = 0.0;
					IIR.b3[j] = 0.0;
					IIR.a4[j] = 0.0;
					IIR.b4[j] = 0.0;
				}

				// Set the number of IIR filter sections we will be generating.
				IIR.NumSections = (IIRFilt.NumPoles + 1) / 2;
				if (IIRFilt.IIRPassType == TIIRPassTypes.BPF || IIRFilt.IIRPassType == TIIRPassTypes.Notch)
					IIR.NumSections = IIRFilt.NumPoles;


				// For All Pass filters, the numerator is set to the denominator values as shown here.
				// If the prototype was an Inv Cheby or Elliptic, the S plane numerator is discarded.
				// Use the Gauss as the prototype for the best all pass results (most linear phase).
				// The all pass H(s) = ( As^2 - Bs + C ) / ( As^2 + Bs + C )
				if (IIRFilt.IIRPassType == TIIRPassTypes.AllPass)
				{
					for (j = 0; j < SPlaneCoeff.NumSections; j++)
					{
						SPlaneCoeff.N2[j] = SPlaneCoeff.D2[j];
						SPlaneCoeff.N1[j] = -SPlaneCoeff.D1[j];
						SPlaneCoeff.N0[j] = SPlaneCoeff.D0[j];
					}
				}

				// T sets the IIR filter's corner frequency, or center freqency.
				// The Bilinear transform is defined as:  s = 2/T * tan(Omega/2) = 2/T * (1 - z)/(1 + z)
				T = 2.0 * Tan(IIRFilt.OmegaC * PI/2);
				Q = 1.0 + IIRFilt.OmegaC;             // Q is used for band pass and notch filters.
				if (Q > 1.95)
					Q = 1.95;
				Q = 0.8 * Tan(Q * PI/4);            // This is a correction factor for Q.
				Q = IIRFilt.OmegaC / IIRFilt.BW / Q;  // This is the corrected Q.


				// Calc the IIR coefficients.
				// SPlaneCoeff.NumSections is the number of 1st and 2nd order s plane factors.
				k = 0;
				for (j = 0; j < SPlaneCoeff.NumSections; j++)
				{
					A = SPlaneCoeff.D2[j]; // We use A - F to make the code easier to read.
					B = SPlaneCoeff.D1[j];
					C = SPlaneCoeff.D0[j];
					D = SPlaneCoeff.N2[j];
					E = SPlaneCoeff.N1[j]; // N1 is always zero, except for the all pass. Consequently, the equations below can be simplified a bit by removing E.
					F = SPlaneCoeff.N0[j];

					// b's are the numerator  a's are the denominator
					if (IIRFilt.IIRPassType == TIIRPassTypes.LPF || IIRFilt.IIRPassType == TIIRPassTypes.AllPass) // Low Pass and All Pass
					{
						if (A == 0.0 && D == 0.0) // 1 pole case
						{
							Arg = (2.0 * B + C * T);
							IIR.a2[j] = 0.0;
							IIR.a1[j] = (-2.0 * B + C * T) / Arg;
							IIR.a0[j] = 1.0;

							IIR.b2[j] = 0.0;
							IIR.b1[j] = (-2.0 * E + F * T) / Arg * C / F;
							IIR.b0[j] = (2.0 * E + F * T) / Arg * C / F;
						}
						else // 2 poles
						{
							Arg = (4.0 * A + 2.0 * B * T + C * T * T);
							IIR.a2[j] = (4.0 * A - 2.0 * B * T + C * T * T) / Arg;
							IIR.a1[j] = (2.0 * C * T * T - 8.0 * A) / Arg;
							IIR.a0[j] = 1.0;

							// With all pole filters, our LPF numerator is (z+1)^2, so all our Z Plane zeros are at -1
							IIR.b2[j] = (4.0 * D - 2.0 * E * T + F * T * T) / Arg * C / F;
							IIR.b1[j] = (2.0 * F * T * T - 8.0 * D) / Arg * C / F;
							IIR.b0[j] = (4 * D + F * T * T + 2.0 * E * T) / Arg * C / F;
						}
					}

					if (IIRFilt.IIRPassType == TIIRPassTypes.HPF) // High Pass
					{
						if (A == 0.0 && D == 0.0) // 1 pole
						{
							Arg = 2.0 * C + B * T;
							IIR.a2[j] = 0.0;
							IIR.a1[j] = (B * T - 2.0 * C) / Arg;
							IIR.a0[j] = 1.0;

							IIR.b2[j] = 0.0;
							IIR.b1[j] = (E * T - 2.0 * F) / Arg * C / F;
							IIR.b0[j] = (E * T + 2.0 * F) / Arg * C / F;
						}
						else  // 2 poles
						{
							Arg = A * T * T + 4.0 * C + 2.0 * B * T;
							IIR.a2[j] = (A * T * T + 4.0 * C - 2.0 * B * T) / Arg;
							IIR.a1[j] = (2.0 * A * T * T - 8.0 * C) / Arg;
							IIR.a0[j] = 1.0;

							// With all pole filters, our HPF numerator is (z-1)^2, so all our Z Plane zeros are at 1
							IIR.b2[j] = (D * T * T - 2.0 * E * T + 4.0 * F) / Arg * C / F;
							IIR.b1[j] = (2.0 * D * T * T - 8.0 * F) / Arg * C / F;
							IIR.b0[j] = (D * T * T + 4.0 * F + 2.0 * E * T) / Arg * C / F;
						}
					}

					if (IIRFilt.IIRPassType == TIIRPassTypes.BPF) // Band Pass
					{
						if (A == 0.0 && D == 0.0) // 1 pole
						{
							Arg = 4.0 * B * Q + 2.0 * C * T + B * Q * T * T;
							a2[k] = (B * Q * T * T + 4.0 * B * Q - 2.0 * C * T) / Arg;
							a1[k] = (2.0 * B * Q * T * T - 8.0 * B * Q) / Arg;
							a0[k] = 1.0;

							b2[k] = (E * Q * T * T + 4.0 * E * Q - 2.0 * F * T) / Arg * C / F;
							b1[k] = (2.0 * E * Q * T * T - 8.0 * E * Q) / Arg * C / F;
							b0[k] = (4.0 * E * Q + 2.0 * F * T + E * Q * T * T) / Arg * C / F;
							k++;
						}
						else //2 Poles
						{
							IIR.a4[j] = (16.0 * A * Q * Q + A * Q * Q * T * T * T * T + 8.0 * A * Q * Q * T * T - 2.0 * B * Q * T * T * T - 8.0 * B * Q * T + 4.0 * C * T * T) * F;
							IIR.a3[j] = (4.0 * T * T * T * T * A * Q * Q - 4.0 * Q * T * T * T * B + 16.0 * Q * B * T - 64.0 * A * Q * Q) * F;
							IIR.a2[j] = (96.0 * A * Q * Q - 16.0 * A * Q * Q * T * T + 6.0 * A * Q * Q * T * T * T * T - 8.0 * C * T * T) * F;
							IIR.a1[j] = (4.0 * T * T * T * T * A * Q * Q + 4.0 * Q * T * T * T * B - 16.0 * Q * B * T - 64.0 * A * Q * Q) * F;
							IIR.a0[j] = (16.0 * A * Q * Q + A * Q * Q * T * T * T * T + 8.0 * A * Q * Q * T * T + 2.0 * B * Q * T * T * T + 8.0 * B * Q * T + 4.0 * C * T * T) * F;

							// With all pole filters, our BPF numerator is (z-1)^2 * (z+1)^2 so the zeros come back as +/- 1 pairs
							IIR.b4[j] = (8.0 * D * Q * Q * T * T - 8.0 * E * Q * T + 16.0 * D * Q * Q - 2.0 * E * Q * T * T * T + D * Q * Q * T * T * T * T + 4.0 * F * T * T) * C;
							IIR.b3[j] = (16.0 * E * Q * T - 4.0 * E * Q * T * T * T - 64.0 * D * Q * Q + 4.0 * D * Q * Q * T * T * T * T) * C;
							IIR.b2[j] = (96.0 * D * Q * Q - 8.0 * F * T * T + 6.0 * D * Q * Q * T * T * T * T - 16.0 * D * Q * Q * T * T) * C;
							IIR.b1[j] = (4.0 * D * Q * Q * T * T * T * T - 64.0 * D * Q * Q + 4.0 * E * Q * T * T * T - 16.0 * E * Q * T) * C;
							IIR.b0[j] = (16.0 * D * Q * Q + 8.0 * E * Q * T + 8.0 * D * Q * Q * T * T + 2.0 * E * Q * T * T * T + 4.0 * F * T * T + D * Q * Q * T * T * T * T) * C;

							// T = 2 makes these values approach 0.0 (~ 1.0E-12) The root solver needs 0.0 for numerical reasons.
							if (Abs(T - 2.0) < 0.0005)
							{
								IIR.a3[j] = 0.0;
								IIR.a1[j] = 0.0;
								IIR.b3[j] = 0.0;
								IIR.b1[j] = 0.0;
							}

							// We now have a 4th order poly in the form a4*s^4 + a3*s^3 + a2*s^2 + a2*s + a0
							// We find the roots of this so we can form two 2nd order polys.
							Coeff[0] = IIR.a4[j];
							Coeff[1] = IIR.a3[j];
							Coeff[2] = IIR.a2[j];
							Coeff[3] = IIR.a1[j];
							Coeff[4] = IIR.a0[j];
							P51OneRevC.FindRoots(4, ref Coeff, ref Roots);

							// In effect, the root finder scales the poly by 1/a4 so we have to apply this factor back into
							// the two 2nd order equations we are forming.
							Scalar = Sqrt(Abs(IIR.a4[j]));

							// Form the two 2nd order polys from the roots.
							a2[k] = Scalar;
							a1[k] = -(Roots[0] + Roots[1]).Real * Scalar;
							a0[k] = (Roots[0] * Roots[1]).Real * Scalar;
							k++;
							a2[k] = Scalar;
							a1[k] = -(Roots[2] + Roots[3]).Real * Scalar;
							a0[k] = (Roots[2] * Roots[3]).Real * Scalar;
							k--;

							// Now do the same with the numerator.
							Coeff[0] = IIR.b4[j];
							Coeff[1] = IIR.b3[j];
							Coeff[2] = IIR.b2[j];
							Coeff[3] = IIR.b1[j];
							Coeff[4] = IIR.b0[j];


							if (IIRFilt.ProtoType == TFilterPoly.INVERSE_CHEBY || IIRFilt.ProtoType == TFilterPoly.ELLIPTIC)
							{
								P51OneRevC.FindRoots(4, ref Coeff, ref Roots);
							}
							else // With all pole filters (Butter, Cheb, etc), we know we have these 4 real roots. The root finder won't necessarily pair real roots the way we need, so rather than compute these, we simply set them.
							{
								Roots[0] = new Complex(-1.0, 0.0);
								Roots[1] = new Complex(1.0, 0.0);
								Roots[2] = new Complex(-1.0, 0.0);
								Roots[3] = new Complex(1.0, 0.0);
							}

							Scalar = Sqrt(Abs(IIR.b4[j]));

							b2[k] = Scalar;
							if (IIRFilt.ProtoType == TFilterPoly.INVERSE_CHEBY || IIRFilt.ProtoType == TFilterPoly.ELLIPTIC)
							{
								b1[k] = -(Roots[0] + Roots[1]).Real * Scalar; // = 0.0
							}
							else  // else the prototype is an all pole filter
							{
								b1[k] = 0.0;  // b1 = 0 for all pole filters, but the addition above won't always equal zero exactly.
							}
							b0[k] = (Roots[0] * Roots[1]).Real * Scalar;

							k++;

							b2[k] = Scalar;
							if (IIRFilt.ProtoType == TFilterPoly.INVERSE_CHEBY || IIRFilt.ProtoType == TFilterPoly.ELLIPTIC)
							{
								b1[k] = -(Roots[2] + Roots[3]).Real * Scalar;
							}
							else // All pole
							{
								b1[k] = 0.0;
							}
							b0[k] = (Roots[2] * Roots[3]).Real * Scalar;
							k++;
							// Go below to see where we store these 2nd order polys back into IIR
						}
					}

					if (IIRFilt.IIRPassType ==TIIRPassTypes.Notch) // Notch
					{
						if (A == 0.0 && D == 0.0) // 1 pole
						{
							Arg = 2.0 * B * T + C * Q * T * T + 4.0 * C * Q;
							a2[k] = (4.0 * C * Q - 2.0 * B * T + C * Q * T * T) / Arg;
							a1[k] = (2.0 * C * Q * T * T - 8.0 * C * Q) / Arg;
							a0[k] = 1.0;

							b2[k] = (4.0 * F * Q - 2.0 * E * T + F * Q * T * T) / Arg * C / F;
							b1[k] = (2.0 * F * Q * T * T - 8.0 * F * Q) / Arg * C / F;
							b0[k] = (2.0 * E * T + F * Q * T * T + 4.0 * F * Q) / Arg * C / F;
							k++;
						}
						else
						{
							IIR.a4[j] = (4.0 * A * T * T - 2.0 * B * T * T * T * Q + 8.0 * C * Q * Q * T * T - 8.0 * B * T * Q + C * Q * Q * T * T * T * T + 16.0 * C * Q * Q) * -F;
							IIR.a3[j] = (16.0 * B * T * Q + 4.0 * C * Q * Q * T * T * T * T - 64.0 * C * Q * Q - 4.0 * B * T * T * T * Q) * -F;
							IIR.a2[j] = (96.0 * C * Q * Q - 8.0 * A * T * T - 16.0 * C * Q * Q * T * T + 6.0 * C * Q * Q * T * T * T * T) * -F;
							IIR.a1[j] = (4.0 * B * T * T * T * Q - 16.0 * B * T * Q - 64.0 * C * Q * Q + 4.0 * C * Q * Q * T * T * T * T) * -F;
							IIR.a0[j] = (4.0 * A * T * T + 2.0 * B * T * T * T * Q + 8.0 * C * Q * Q * T * T + 8.0 * B * T * Q + C * Q * Q * T * T * T * T + 16.0 * C * Q * Q) * -F;

							// Our Notch Numerator isn't simple. [ (4+T^2)*z^2 - 2*(4-T^2)*z + (4+T^2) ]^2
							IIR.b4[j] = (2.0 * E * T * T * T * Q - 4.0 * D * T * T - 8.0 * F * Q * Q * T * T + 8.0 * E * T * Q - 16.0 * F * Q * Q - F * Q * Q * T * T * T * T) * C;
							IIR.b3[j] = (64.0 * F * Q * Q + 4.0 * E * T * T * T * Q - 16.0 * E * T * Q - 4.0 * F * Q * Q * T * T * T * T) * C;
							IIR.b2[j] = (8.0 * D * T * T - 96.0 * F * Q * Q + 16.0 * F * Q * Q * T * T - 6.0 * F * Q * Q * T * T * T * T) * C;
							IIR.b1[j] = (16.0 * E * T * Q - 4.0 * E * T * T * T * Q + 64.0 * F * Q * Q - 4.0 * F * Q * Q * T * T * T * T) * C;
							IIR.b0[j] = (-4.0 * D * T * T - 2.0 * E * T * T * T * Q - 8.0 * E * T * Q - 8.0 * F * Q * Q * T * T - F * Q * Q * T * T * T * T - 16.0 * F * Q * Q) * C;

							// T = 2 (OmegaC = 0.5) makes these values approach 0.0 (~ 1.0E-12). The root solver wants 0.0 for numerical reasons.
							if (Abs(T - 2.0) < 0.0005)
							{
								IIR.a3[j] = 0.0;
								IIR.a1[j] = 0.0;
								IIR.b3[j] = 0.0;
								IIR.b1[j] = 0.0;
							}

							// We now have a 4th order poly in the form a4*s^4 + a3*s^3 + a2*s^2 + a2*s + a0
							// We find the roots of this so we can form two 2nd order polys.
							Coeff[0] = IIR.a4[j];
							Coeff[1] = IIR.a3[j];
							Coeff[2] = IIR.a2[j];
							Coeff[3] = IIR.a1[j];
							Coeff[4] = IIR.a0[j];


							// In effect, the root finder scales the poly by 1/a4 so we have to apply this factor back into
							// the two 2nd order equations we are forming.
							P51OneRevC.FindRoots(4, ref Coeff, ref Roots);
							Scalar = Sqrt(Abs(IIR.a4[j]));
							a2[k] = Scalar;
							a1[k] = -(Roots[0] + Roots[1]).Real * Scalar;
							a0[k] = (Roots[0] * Roots[1]).Real * Scalar;

							k++;
							a2[k] = Scalar;
							a1[k] = -(Roots[2] + Roots[3]).Real * Scalar;
							a0[k] = (Roots[2] * Roots[3]).Real * Scalar;
							k--;

							// Now do the same with the numerator.
							Coeff[0] = IIR.b4[j];
							Coeff[1] = IIR.b3[j];
							Coeff[2] = IIR.b2[j];
							Coeff[3] = IIR.b1[j];
							Coeff[4] = IIR.b0[j];
							P51OneRevC.FindRoots(4, ref Coeff, ref Roots);

							Scalar = Sqrt(Abs(IIR.b4[j]));
							b2[k] = Scalar;
							b1[k] = -(Roots[0] + Roots[1]).Real * Scalar;
							b0[k] = (Roots[0] * Roots[1]).Real * Scalar;

							k++;
							b2[k] = Scalar;
							b1[k] = -(Roots[2] + Roots[3]).Real * Scalar;
							b0[k] = (Roots[2] * Roots[3]).Real * Scalar;
							k++;
						}
					}
				}

				if (IIRFilt.IIRPassType == TIIRPassTypes.BPF || IIRFilt.IIRPassType == TIIRPassTypes.Notch)
				{
					// In the calcs above for the BPF and Notch, we didn't set a0=1, so we do it here.
					for (j = 0; j < IIR.NumSections; j++)
					{
						b2[j] /= a0[j];
						b1[j] /= a0[j];
						b0[j] /= a0[j];
						a2[j] /= a0[j];
						a1[j] /= a0[j];
						a0[j] = 1.0;
					}

					for (j = 0; j < IIR.NumSections; j++)
					{
						IIR.a0[j] = a0[j];
						IIR.a1[j] = a1[j];
						IIR.a2[j] = a2[j];
						IIR.b0[j] = b0[j];
						IIR.b1[j] = b1[j];
						IIR.b2[j] = b2[j];
					}
				}

				// Adjust the b's or a0 for the desired Gain.
				SectionGain = Pow(10.0, IIRFilt.dBGain / 20.0);
				SectionGain = Pow(SectionGain, 1.0 / (double)IIR.NumSections);
				for (j = 0; j < IIR.NumSections; j++)
				{
					IIR.b0[j] *= SectionGain;
					IIR.b1[j] *= SectionGain;
					IIR.b2[j] *= SectionGain;
					// This is an alternative to adjusting the b's
					// IIR.a0[j] = SectionGain;
				}

				return (IIR);
			}

			public static void FilterWithIIR(TIIRCoeff IIRCoeff, ref double[] Signal, ref double[] FilteredSignal, int NumSigPts)
			{
				double y;
				int j, k;

				for (j = 0; j < NumSigPts; j++)
				{
					k = 0;
					y = SectCalc(j, k, Signal[j], IIRCoeff);
					for (k = 1; k < IIRCoeff.NumSections; k++)
					{
						y = SectCalc(j, k, y, IIRCoeff);
					}
					FilteredSignal[j] = y;
				}

			}

			public static double SectCalc(int j, int k, double x, TIIRCoeff IIRCoeff)
			{
				double y, CenterTap;

				if (RegX1 == null)
					RegX1 = new double[ARRAY_DIM];
				if (RegX2 == null)
					RegX2 = new double[ARRAY_DIM];
				if (RegY1 == null)
					RegY1 = new double[ARRAY_DIM];
				if (RegY2 == null)
					RegY2 = new double[ARRAY_DIM];

				// Zero the regiisters on the 1st call or on an overflow condition. The overflow limit used
				// here is small for double variables, but a filter that reaches this threshold is broken.
				if ((j == 0 && k == 0) || MaxRegVal > OVERFLOW_LIMIT)
				{
					if (MaxRegVal > OVERFLOW_LIMIT && !MessageShown)
					{
						// ShowMessage("ERROR: Math Over Flow in IIR Section Calc. \nThe register values exceeded 1.0E20 \n");
						MessageShown = true; // So this message doesn't get shown thousands of times.
					}

					MaxRegVal = 1.0E-12;
					for (int i = 0; i < ARRAY_DIM; i++)
					{
						RegX1[i] = 0.0;
						RegX2[i] = 0.0;
						RegY1[i] = 0.0;
						RegY2[i] = 0.0;
					}
				}

				CenterTap = x * IIRCoeff.b0[k] + IIRCoeff.b1[k] * RegX1[k] + IIRCoeff.b2[k] * RegX2[k];
				y = IIRCoeff.a0[k] * CenterTap - IIRCoeff.a1[k] * RegY1[k] - IIRCoeff.a2[k] * RegY2[k];

				RegX2[k] = RegX1[k];
				RegX1[k] = x;
				RegY2[k] = RegY1[k];
				RegY1[k] = y;

				// MaxRegVal is used to prevent overflow.  Overflow seldom occurs, but will
				// if the filter has faulty coefficients. MaxRegVal is usually less than 100.0
				if (Abs(CenterTap) > MaxRegVal)
					MaxRegVal = Abs(CenterTap);
				if (Abs(y) > MaxRegVal)
					MaxRegVal = Abs(y);
				return (y);
			}

			public static void IIRFreqResponse(TIIRCoeff IIR, int NumSections, ref double[] RealHofZ, ref double[] ImagHofZ, int NumPts)
			{
				int j, n;
				double Arg;
				Complex z1, z2, HofZ, Denom;
				for (j = 0; j < NumPts; j++)
				{
					Arg = PI * (double)j / (double)NumPts;
					z1 = new Complex(Cos(Arg), -Sin(Arg));  // z = e^(j*omega)
					z2 = z1 * z1;                     // z squared

					HofZ = new Complex(1.0, 0.0);
					for (n = 0; n < NumSections; n++)
					{
						HofZ *= IIR.a0[n];  // This can be in the denominator, but only if a0=1. a0 can be other than 1.0 to adjust the filter's gain. See the bottom of the CalcIIRFilterCoeff() function.
						HofZ *= IIR.b0[n] + IIR.b1[n] * z1 + IIR.b2[n] * z2;  // Numerator
						Denom = 1.0 + IIR.a1[n] * z1 + IIR.a2[n] * z2;        // Denominator
						if (Denom.Magnitude < 1.0E-12)
							Denom = 1.0E-12;             // A pole on the unit circle would cause this to be zero, so this should never happen. It could happen however if the filter also has a zero at this frequency. Then H(z) needs to be determined by L'Hopitals rule at this freq.
						HofZ /= Denom;
					}
					RealHofZ[j] = HofZ.Real;
					ImagHofZ[j] = HofZ.Imaginary;
				}
			}
		}

		public static class QuadRootsCode
		{
			private static readonly double LDBL_EPSILON = 2.2204460492503131E-16;
			private static readonly double M_Sqrt3_2 = 0.8660254037844386467637231;
			private static readonly double PI_2 = PI / 2;
			private static readonly double ZERO_MINUS = -8.88178419700125232E-16;
			private static readonly double ZERO_PLUS = 8.88178419700125232E-16;

			/// <summary>
			/// This function is the quadratic formula with P[0] = 1
			/// y=P[0]x^2 + P[1]x+P[2]
			/// </summary>
			/// <param name="P">Coefficients of the quadratic</param>
			/// <param name="RealRoot">Real roots of the quadratic</param>
			/// <param name="ImagRoot">Imaginary roots of the quadratic</param>
			public static void QuadRoots(double[] P, ref double[] RealRoot, ref double[] ImagRoot, int offset = 0)
			{
				double D;
				D = P[1] * P[1] - 4 * P[2];



				if (D >= 0) //1 or 2 real roots
				{
					RealRoot[0 + offset] = (-P[1] - D.Sqrt()) * .5;
					RealRoot[0 + offset] = (-P[1] + D.Sqrt()) * .5;
					ImagRoot[0 + offset] = 0;
					ImagRoot[1 + offset] = 0;
				}
				else //2 Complex roots
				{
					RealRoot[0 + offset] = -P[1] * .5;
					RealRoot[1 + offset] = -P[1] * .5;
					ImagRoot[0 + offset] = (-D).Sqrt() * .5;
					ImagRoot[1 + offset] = -ImagRoot[0];
				}

			}

			/// <summary>
			/// Finds the roots of a cubic equation
			/// y = P0x^3 + P1x^2 + P2x+ P3, P[0] = 1
			/// </summary>
			/// <param name="P">Coefficients of the cubic equation</param>
			/// <param name="RealRoot">Real roots of the cubic</param>
			/// <param name="ImagRoot">Imaginary roots of the cubic</param>
			public static void CubicRoots(double[] P, ref double[] RealRoot, ref double[] ImagRoot, int offset = 0)
			{
				int j;
				double s, t, b, c, d, Scalar;
				bool CubicPolyReversed = false;

				// Scale the polynomial so that P[N] = +/-1. This moves the roots toward unit circle.
				Scalar = Pow(Abs(P[3]), 1.0 / 3.0);
				for (j = 1; j <= 3; j++)
					P[j] /= Pow(Scalar, (double)j);

				if (Abs(P[3]) < Abs(P[2]) && P[2] > 0.0)
				{
					ReversePoly(ref P,3);
					CubicPolyReversed = true;
				}

				s = P[1] / 3.0;
				b = (6.0 * P[1] * P[1] * P[1] - 27.0 * P[1] * P[2] + 81.0 * P[3]) / 162.0;
				t = (P[1] * P[1] - 3.0 * P[2]) / 9.0;
				c = t * t * t;
				d = 2.0 * P[1] * P[1] * P[1] - 9.0 * P[1] * P[2] + 27.0 * P[3];
				d = d * d / 2916.0 - c;

				// if(d > 0) 1 complex and 1 real root. We use LDBL_EPSILON to account for round off err.
				if (d > LDBL_EPSILON)
				{
					d = Pow((d.Sqrt() + Abs(b)), 1.0 / 3.0);
					if (d != 0.0)
					{
						if (b > 0)
							b = -d;
						else
							b = d;
						c = t / b;
					}
					d = M_Sqrt3_2 * (b - c);
					b = b + c;
					c = -b / 2.0 - s;

					RealRoot[0 + offset] = (b - s);
					ImagRoot[0 + offset] = 0.0;
					RealRoot[1 + offset] = RealRoot[2] = c;
					ImagRoot[1 + offset] = d;
					ImagRoot[2 + offset] = -ImagRoot[1];
				}

				else // d < 0.0 3 real roots
				{
					if (b == 0.0)
						d = PI_2 / 3.0; // b can be as small as 1.0E-25
					else
						d = Atan(Sqrt(Abs(d)) / Abs(b)) / 3.0;

					if (b < 0.0)
						b = 2.0 * Sqrt(Abs(t));
					else
						b = -2.0 * Sqrt(Abs(t));

					c = Cos(d) * b;
					t = -M_Sqrt3_2 * Sin(d) * b - 0.5 * c;

					RealRoot[0 + offset] = (t - s);
					RealRoot[1 + offset] = -(t + c + s);
					RealRoot[2 + offset] = (c - s);
					ImagRoot[0 + offset] = 0.0;
					ImagRoot[1 + offset] = 0.0;
					ImagRoot[2 + offset] = 0.0;
				}

				// If we reversed the poly, the roots need to be inverted.
				if (CubicPolyReversed)
					InvertRoots(3, ref RealRoot, ref ImagRoot);

				// Apply the Scalar to the roots.
				for (j = 0; j < 3; j++)
					RealRoot[j] *= Scalar;
				for (j = 0; j < 3; j++)
					ImagRoot[j] *= Scalar;
			}

			/// <summary>
			/// This finds the roots of y = P0x^4 + P1x^3 + P2x^2 + P3x + P4 P[0] = 1
			/// </summary>
			/// <param name="P">Quartic equation coefficients</param>
			/// <param name="RealRoot">Real roots of the cubic</param>
			/// <param name="ImagRoot">Imaginary roots of the cubic</param>
			public static void BiQuadRoots(double[] P, ref double[] RealRoot, ref double[] ImagRoot, int offset = 0)
			{
				int j;
				double a, b, c, d, e, Q3Limit, Scalar, MinRoot;
				double[] Q = new double[ARRAY_DIM];
				bool QuadPolyReversed = false;

				RealRoot = new double[ARRAY_DIM];
				ImagRoot = new double[ARRAY_DIM];

				// Scale the polynomial so that P[N] = +/- 1. This moves the roots toward unit circle.
				Scalar = Pow(Abs(P[4]), 0.25);
				for (j = 1; j <= 4; j++)
					P[j] /= Pow(Scalar, (double)j);

				// Having P[1] < P[3] helps with the Q[3] calculation and test.
				if (Abs(P[1]) > Abs(P[3]))
				{
					ReversePoly(ref P,4);
					QuadPolyReversed = true;
				}

				a = P[2] - P[1] * P[1] * 0.375;
				b = P[3] + P[1] * P[1] * P[1] * 0.125 - P[1] * P[2] * 0.5;
				c = P[4] + 0.0625 * P[1] * P[1] * P[2] - 0.01171875 * P[1] * P[1] * P[1] * P[1] - 0.25 * P[1] * P[3];
				e = P[1] * 0.25;

				Q[0] = 1.0;
				Q[1] = P[2] * 0.5 - P[1] * P[1] * 0.1875;
				Q[2] = (P[2] * P[2] - P[1] * P[1] * P[2] + 0.1875 * P[1] * P[1] * P[1] * P[1] - 4.0 * P[4] + P[1] * P[3]) * 0.0625;
				Q[3] = -b * b * 0.015625;


				/* The value of Q[3] can cause problems when it should have calculated to zero (just above) but
				is instead ~ -1.0E-17 because of round off errors. Consequently, we need to determine whether
				a tiny Q[3] is due to roundoff, or if it is legitimately small. It can legitimately have values
				of ~ -1E-28. When this happens, we assume Q[2] should also be small. Q[3] can also be tiny with
				2 sets of equal real roots. Then P[1] and P[3], are approx equal. */

				Q3Limit = ZERO_MINUS;
				if (Abs(Abs(P[1]) - Abs(P[3])) >= ZERO_PLUS &&
				Q[3] > ZERO_MINUS && Abs(Q[2]) < 1.0E-5)
					Q3Limit = 0.0;

				if (Q[3] < Q3Limit && Abs(Q[2]) < 1.0E20 * Abs(Q[3]))
				{
					CubicRoots(Q, ref RealRoot, ref ImagRoot);

					// Find the smallest positive real root. One of the real roots is always positive.
					MinRoot = 1.0E100;
					for (j = 0; j < 3; j++)
					{
						if (ImagRoot[j] == 0.0 && RealRoot[j] > 0 && RealRoot[j] < MinRoot)
							MinRoot = RealRoot[j];
					}

					d = 4.0 * MinRoot;
					a += d;
					if (a * b < 0.0)
						Q[1] = -Sqrt(d);
					else
						Q[1] = Sqrt(d);
					b = 0.5 * (a + b / Q[1]);
				}
				else
				{
					if (Q[2] < 0.0) // 2 sets of equal imag roots
					{
						b = Sqrt(Abs(c));
						d = b + b - a;
						if (d > 0.0)
							Q[1] = Sqrt(Abs(d));
						else
							Q[1] = 0.0;
					}
					else
					{
						if (Q[1] > 0.0)
							b = 2.0 * Sqrt(Abs(Q[2])) + Q[1];
						else
							b = -2.0 * Sqrt(Abs(Q[2])) + Q[1];
						Q[1] = 0.0;
					}
				}

				// Calc the roots from two 2nd order polys and subtract e from the real part.
				if (Abs(b) > 1.0E-8)
				{
					Q[2] = c / b;
					QuadRoots(Q, ref RealRoot, ref ImagRoot);

					Q[1] = -Q[1];
					Q[2] = b;
					QuadRoots(Q, ref RealRoot, ref ImagRoot, 2);

					for (j = 0; j < 4; j++)
						RealRoot[j] -= e;
				}
				else // b==0 with 4 equal real roots
				{
					for (j = 0; j < 4; j++)
						RealRoot[j] = -e;
					for (j = 0; j < 4; j++)
						ImagRoot[j] = 0.0;
				}

				// If we reversed the poly, the roots need to be inverted.
				if (QuadPolyReversed)
					InvertRoots(4, ref RealRoot, ref ImagRoot);

				// Apply the Scalar to the roots.
				for (j = 0; j < 4; j++)
					RealRoot[j] *= Scalar;
				for (j = 0; j < 4; j++)
					ImagRoot[j] *= Scalar;
			}

			/// <summary>
			/// A reversed polynomial has its roots at the same angle, but reflected about the unit circle.
			/// </summary>
			/// <param name="P">polynomial to reverse</param>
			public static void ReversePoly(ref double[] P, int N)
			{
				double temp;
				for (int i = 0; i < N; i++)
				{
					temp = P[i];
					P[i] = P[N - i];
					P[N - i] = temp;
				}

				if (P[0] != 0)
				{
					for (int i = N; i >= 0; i--)
					{
						P[i] /= P[0];
					}
				}
			}

			/// <summary>
			/// This is used in conjunction with ReversePoly
			/// </summary>
			/// <param name="RealRoot">Real roots</param>
			/// <param name="ImagRoot">imaginary roots</param>
			public static void InvertRoots(int N, ref double[] RealRoot, ref double[] ImagRoot)
			{
				int j;
				double Mag;
				for (j = 0; j < N; j++)
				{
					// complex math for 1/x
					Mag = RealRoot[j] * RealRoot[j] + ImagRoot[j] * ImagRoot[j];
					if (Mag != 0.0)
					{
						RealRoot[j] /= Mag;
						ImagRoot[j] /= -Mag;
					}
				}
			}
		}

		public static class LowPassPrototypes
		{
			public static readonly int MAX_POLE_COUNT = 20;
			public static readonly int ARRAY_DIM = 50;

			public enum TOurSortTypes { stMax, stMin };

			/// <summary>
			/// These are the available filter polynomials. NOT_IIR is for code testing.
			/// </summary>
			public enum TFilterPoly
			{
				BUTTERWORTH,
				GAUSSIAN,
				BESSEL,
				ADJUSTABLE,
				CHEBYSHEV,
				INVERSE_CHEBY,
				PAPOULIS,
				ELLIPTIC,
				NOT_IIR
			};

			/// <summary>
			/// These coeff form H(s) = (N2*s^2 + N1*s + N0) / (D2*s^2 + D1*s + D0)
			/// NumSections is the number of 1st and 2nd order polynomial factors .
			/// </summary>
			public struct TSPlaneCoeff
			{
				public double[] N2;
				public double[] N1;
				public double[] N0;

				public double[] D2;
				public double[] D1;
				public double[] D0;
				public int NumSections;

				public TSPlaneCoeff(bool initalize = true)
				{
					N2=new double[ARRAY_DIM];
					N1=new double[ARRAY_DIM];
					N0=new double[ARRAY_DIM];

					D2=new double[ARRAY_DIM];
					D1=new double[ARRAY_DIM];
					D0=new double[ARRAY_DIM];
					NumSections = 0;
				}
			};

			/// <summary>
			/// This structure defines the low pass filter prototype.
			/// The 3 dB corner frequency is 1 rad/sec for all filters.
			/// </summary>
			public struct TLowPassParams
			{
				public TFilterPoly ProtoType; // Butterworth, Cheby, etc.
				public int NumPoles; // Pole count
				public double Ripple; // Passband Ripple for the Elliptic and Chebyshev
				public double StopBanddB; // Stop Band Attenuation in dB for the Elliptic and Inverse Cheby
				public double Gamma; // Controls the transition bandwidth on the Adjustable Gauss. -1 <= Gamma <= 1
			};

			/// <summary>
			///
			/// </summary>
			/// <param name="Filt">low pass prototype (NumPoles, Ripple, etc.)</param>
			/// <returns>SPlaneCoeff filled with the 2nd order S plane coefficients.</returns>
			public static TSPlaneCoeff CalcLowPassProtoCoeff(TLowPassParams Filt)
			{
				int j, DenomCount=0, NumerCount, NumRoots, ZeroCount=1;
				Complex[] Poles = new Complex[ARRAY_DIM];
				Complex[] Zeros = new Complex[ARRAY_DIM];

				TSPlaneCoeff Coeff = new TSPlaneCoeff(true); // The return value.

				// Init the S Plane Coeff. H(s) = (N2*s^2 + N1*s + N0) / (D2*s^2 + D1*s + D0)
				for (j = 0; j < ARRAY_DIM; j++)
				{
					Coeff.N2[j] = 0.0;
					Coeff.N1[j] = 0.0;
					Coeff.N0[j] = 1.0;
					Coeff.D2[j] = 0.0;
					Coeff.D1[j] = 0.0;
					Coeff.D0[j] = 1.0;
				}
				Coeff.NumSections = 0;


				// We need to range check the various argument values here.
				// These are the practical limits the max number of poles.
				if (Filt.NumPoles < 1)
					Filt.NumPoles = 1;
				if (Filt.NumPoles > MAX_POLE_COUNT)
					Filt.NumPoles = MAX_POLE_COUNT;
				if (Filt.ProtoType == TFilterPoly.ELLIPTIC || Filt.ProtoType == TFilterPoly.INVERSE_CHEBY)
				{
					if (Filt.NumPoles > 15)
						Filt.NumPoles = 15;
				}
				if (Filt.ProtoType == TFilterPoly.GAUSSIAN || Filt.ProtoType == TFilterPoly.BESSEL)
				{
					if (Filt.NumPoles > 12)
						Filt.NumPoles = 12;
				}

				// Gamma is used by the Adjustable Gauss.
				if (Filt.Gamma < -1.0)
					Filt.Gamma = -1.0; // -1 gives ~ Gauss response
				if (Filt.Gamma > 1.0)
					Filt.Gamma = 1.0; // +1 gives ~ Butterworth response.

				// Ripple is used by the Chebyshev and Elliptic
				if (Filt.Ripple < 0.0001)
					Filt.Ripple = 0.0001;
				if (Filt.Ripple > 1.0)
					Filt.Ripple = 1.0;

				// With the Chebyshev we need to use less ripple for large pole counts to keep the poles out of the RHP.
				if (Filt.ProtoType == TFilterPoly.CHEBYSHEV && Filt.NumPoles > 15)
				{
					double MaxRipple = 1.0;
					if (Filt.NumPoles == 16)
						MaxRipple = 0.5;
					if (Filt.NumPoles == 17)
						MaxRipple = 0.4;
					if (Filt.NumPoles == 18)
						MaxRipple = 0.25;
					if (Filt.NumPoles == 19)
						MaxRipple = 0.125;
					if (Filt.NumPoles >= 20)
						MaxRipple = 0.10;
					if (Filt.Ripple > MaxRipple)
						Filt.Ripple = MaxRipple;
				}

				// StopBanddB is used by the Inverse Chebyshev and the Elliptic
				// It is given in positive dB values.
				if (Filt.StopBanddB < 20.0)
					Filt.StopBanddB = 20.0;
				if (Filt.StopBanddB > 120.0)
					Filt.StopBanddB = 120.0;


				// There isn't such a thing as a 1 pole Chebyshev, or 1 pole Bessel, etc.
				// A one pole filter is simply 1/(s+1).
				NumerCount = 0; // init
				if (Filt.NumPoles == 1)
				{
					Coeff.D1[0] = 1.0;
					DenomCount = 1; // DenomCount is the number of denominator factors (1st or 2nd order).
				}
				else if (Filt.ProtoType == TFilterPoly.BUTTERWORTH)
				{
					NumRoots =LowPassRoots.ButterworthPoly(Filt.NumPoles, ref Poles);
					DenomCount =GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					// A Butterworth doesn't require frequncy scaling with SetCornerFreq().
				}

				else if (Filt.ProtoType == TFilterPoly.ADJUSTABLE) // Adjustable Gauss
				{
					NumRoots = LowPassRoots.AdjustablePoly(Filt.NumPoles, ref Poles, Filt.Gamma);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				else if (Filt.ProtoType == TFilterPoly.CHEBYSHEV)
				{
					NumRoots = LowPassRoots.ChebyshevPoly(Filt.NumPoles, Filt.Ripple, ref Poles);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				else if (Filt.ProtoType == TFilterPoly.INVERSE_CHEBY)
				{
					NumRoots = LowPassRoots.InvChebyPoly(Filt.NumPoles, Filt.StopBanddB, ref Poles, ref Zeros, ref ZeroCount);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					NumerCount = GetFilterCoeff(ZeroCount, ref Zeros, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				else if (Filt.ProtoType == TFilterPoly.ELLIPTIC)
				{
					NumRoots = LowPassRoots.EllipticPoly(Filt.NumPoles, Filt.Ripple, Filt.StopBanddB, ref Poles, ref Zeros, ref ZeroCount);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					NumerCount = GetFilterCoeff(ZeroCount, ref Zeros, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				// Papoulis works OK, but it doesn't accomplish anything the Chebyshev can't.
				else if (Filt.ProtoType == TFilterPoly.PAPOULIS)
				{
					NumRoots = LowPassRoots.PapoulisPoly(Filt.NumPoles, ref Poles);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				else if (Filt.ProtoType == TFilterPoly.BESSEL)
				{
					NumRoots = LowPassRoots.BesselPoly(Filt.NumPoles, ref Poles);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				else if (Filt.ProtoType == TFilterPoly.GAUSSIAN)
				{
					NumRoots = LowPassRoots.GaussianPoly(Filt.NumPoles, ref Poles);
					DenomCount = GetFilterCoeff(NumRoots, ref Poles, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0);
					SetCornerFreq(DenomCount, ref Coeff.D2, ref Coeff.D1, ref Coeff.D0, ref Coeff.N2, ref Coeff.N1, ref Coeff.N0);
				}

				Coeff.NumSections = DenomCount;

				// If we have an odd pole count, there will be 1 less zero than poles, so we need to shift the
				// zeros down in the arrays so the 1st zero (which is zero) and aligns with the real pole.
				if (NumerCount != 0 && Filt.NumPoles % 2 == 1)
				{
					for (j = NumerCount; j >= 0; j--)
					{
						Coeff.N2[j + 1] = Coeff.N2[j]; // Coeff.N1's are always zero
						Coeff.N0[j + 1] = Coeff.N0[j];
					}
					Coeff.N2[0] = 0.0; // Set the 1st zero to zero for odd pole counts.
					Coeff.N0[0] = 1.0;
				}

				return (Coeff);

			}

			public static void SetCornerFreq(int PolyCount, ref double[] D2, ref double[] D1, ref double[] D0, ref double[] N2, ref double[] N1, ref double[] N0)
			{
				int j, n;
				double Omega=0, FreqScalar, Zeta, Gain;
				Complex s, H;

				Gain = 1.0;
				for (j = 0; j < PolyCount; j++)
					Gain *= D0[j] / N0[j];

				// Evaluate H(s) by increasing Omega until |H(s)| < -3 dB
				for (j = 1; j < 6000; j++)
				{
					Omega = (double)j / 512.0; // The step size for Omega is 1/512 radians.
					s = new Complex(0.0, Omega);

					H = new Complex(1.0, 0.0);
					for (n = 0; n < PolyCount; n++)
					{
						H = H * (N2[n] * s * s + N1[n] * s + N0[n]) / (D2[n] * s * s + D1[n] * s + D0[n]);
					}
					H *= Gain;
					if (H.Magnitude < 0.7071)
						break;  // -3 dB
				}

				FreqScalar = 1.0 / Omega;

				// Freq scale the denominator. We hold the damping factor Zeta constant.
				for (j = 0; j < PolyCount; j++)
				{
					Omega = Sqrt(D0[j]);
					if (Omega == 0.0)
						continue;   // should never happen
					Zeta = D1[j] / Omega / 2.0;
					if (D2[j] != 0.0)           // 2nd degree poly
					{
						D0[j] = Omega * Omega * FreqScalar * FreqScalar;
						D1[j] = 2.0 * Zeta * Omega * FreqScalar;
					}
					else  // 1st degree poly
					{
						D0[j] *= FreqScalar;
					}
				}

				// Scale the numerator.   H(s) = (N2*s^2 + N1*s + N0) / (D2*s^2 + D1*s + D0)
				// N1 is always zero. N2 is either 1 or 0. If N2 = 0, then N0 = 1 and there isn't a zero to scale.
				// For all pole filters (Butter, Cheby, etc) N2 = 0 and N0 = 1.
				for (j = 0; j < PolyCount; j++)
				{
					if (N2[j] == 0.0)
						continue;
					N0[j] *= FreqScalar * FreqScalar;
				}

			}

			public static int GetFilterCoeff(int RootCount, ref Complex[] Roots, ref double[] A2, ref double[] A1, ref double[] A0)
			{
				int PolyCount, j, k;

				SortRootsByZeta(ref Roots, RootCount, TOurSortTypes.stMin);   // stMin places the most negative real part 1st.

				// Check for duplicate roots. The Inv Cheby generates duplcate imag roots, and the
				// Elliptic generates duplicate real roots. We set duplicates to a RHP value.
				for (j = 0; j < RootCount - 1; j++)
				{
					for (k = j + 1; k < RootCount; k++)
					{
						if (Abs(Roots[j].Real - Roots[k].Real) < 1.0E-3 &&
							Abs(Roots[j].Imaginary - Roots[k].Imaginary) < 1.0E-3)
						{
							Roots[k] = new Complex((double)k, 0.0); // RHP roots are ignored below, Use k is to prevent duplicate checks for matches.
						}
					}
				}

				// This forms the 2nd order coefficients from the root value.
				// We ignore roots in the Right Hand Plane.
				PolyCount = 0;
				for (j = 0; j < RootCount; j++)
				{
					if (Roots[j].Real > 0.0)
						continue; // Right Hand Plane
					if (Roots[j].Real == 0.0 && Roots[j].Imaginary == 0.0)
						continue; // At the origin.  This should never happen.

					if (Roots[j].Real == 0.0) // Imag Root (A poly zero)
					{
						A2[PolyCount] = 1.0;
						A1[PolyCount] = 0.0;
						A0[PolyCount] = Roots[j].Imaginary * Roots[j].Imaginary;
						j++;
						PolyCount++;
					}
					else if (Roots[j].Imaginary == 0.0) // Real Pole
					{
						A2[PolyCount] = 0.0;
						A1[PolyCount] = 1.0;
						A0[PolyCount] = -Roots[j].Real;
						PolyCount++;
					}
					else // Complex Pole
					{
						A2[PolyCount] = 1.0;
						A1[PolyCount] = -2.0 * Roots[j].Real;
						A0[PolyCount] = Roots[j].Real * Roots[j].Real + Roots[j].Imaginary * Roots[j].Imaginary;
						j++;
						PolyCount++;
					}
				}

				return (PolyCount);

			}

			public static int RebuildPoly(int PolyCount, ref double[] PolyCoeff, ref double[] A2, ref double[] A1, ref double[] A0)
			{
				int j, k, n;
				double[] Sum = new double[P51OneRevC.P51_ARRAY_SIZE];
				for (j = 0; j <= 2 * PolyCount; j++)
					PolyCoeff[j] = 0.0;
				for (j = 0; j < P51OneRevC.P51_ARRAY_SIZE; j++)
					Sum[j] = 0.0;

				PolyCoeff[2] = A2[0];
				PolyCoeff[1] = A1[0];
				PolyCoeff[0] = A0[0];

				for (j = 1; j < PolyCount; j++)
				{
					for (n = 0; n <= 2 * j; n++)
					{
						Sum[n + 2] += PolyCoeff[n] * A2[j];
						Sum[n + 1] += PolyCoeff[n] * A1[j];
						Sum[n + 0] += PolyCoeff[n] * A0[j];
					}
					for (k = 0; k <= 2 * j + 2; k++)
						PolyCoeff[k] = Sum[k];
					for (k = 0; k < P51OneRevC.P51_ARRAY_SIZE; k++)
						Sum[k] = 0.0;
				}

				// Want to return the poly order. This will be 2 * PolyCount if there aren't any 1st order polys.
				// 1st order Polys create leading zeros. N 1st order polys Gives N leading zeros.
				for (j = 2 * PolyCount; j >= 0; j--)
					if (PolyCoeff[j] != 0.0)
						break;
				return (j);
			}

			public static void SortRootsByZeta(ref Complex[] Roots, int Count, TOurSortTypes SortType)
			{
				if (Count >=P51OneRevC.P51_MAXDEGREE)
				{
					//ShowMessage("Count > P51_MAXDEGREE in TPolyForm::SortRootsByZeta()");
					return;
				}

				int j, k;
				int[] RootJ= new int[P51OneRevC.P51_ARRAY_SIZE];
				double[] SortValue=new double[P51OneRevC.P51_ARRAY_SIZE];
				Complex[] TempRoots=new Complex[P51OneRevC.P51_ARRAY_SIZE];

				// Set an inconsequential real or imag part to zero.
				for (j = 0; j < Count; j++)
				{
					if (Abs(Roots[j].Real) * 1.0E3 < Abs(Roots[j].Imaginary))
						Roots[j] =new Complex( 0.0, Roots[j].Imaginary);
					if (Abs(Roots[j].Imaginary) * 1.0E3 < Abs(Roots[j].Real))
						Roots[j] =new Complex(Roots[j].Real, 0.0);
				}

				// Sort the roots.
				for (j = 0; j < Count; j++)
					RootJ[j] = j;  // Needed for HeapIndexSort
				if (Roots[0].Real != 0.0) // Cplx roots
				{
					for (j = 0; j < Count; j++)
						SortValue[j] = Roots[j].Real;
				}
				else  // Imag roots, so we sort on imag part.
				{
					for (j = 0; j < Count; j++)
						SortValue[j] = Abs(Roots[j].Imaginary);
				}
				HeapIndexSort(ref SortValue, ref RootJ, Count, SortType);  // stMin gives the most negative root on top

				for (j = 0; j < Count; j++)
				{
					k = RootJ[j];   // RootJ is the sort index
					TempRoots[j] = Roots[k];
				}
				for (j = 0; j < Count; j++)
				{
					Roots[j] = TempRoots[j];
				}

			}

			public static bool HeapIndexSort(ref double[] Data, ref int[] Index, int N, TOurSortTypes SortType)
			{
				int i, j, k, m, IndexTemp;
				long FailSafe, NSquared; // need this for big sorts

				NSquared = (long)N * (long)N;
				m = N / 2;
				k = N - 1;
				for (FailSafe = 0; FailSafe < NSquared; FailSafe++) // typical FailSafe value on return is N*log2(N)
				{
					if (m > 0)
						IndexTemp = Index[--m];
					else
					{
						IndexTemp = Index[k];
						Index[k] = Index[0];
						if (--k == 0)
						{
							Index[0] = IndexTemp;
							return (true);
						}
					}

					i = m + 1;
					j = 2 * i;

					if (SortType == TOurSortTypes.stMax)
						while (j < k + 2)
						{
							FailSafe++;
							if (j <= k && Data[Index[j - 1]] > Data[Index[j]])
								j++;
							if (Data[IndexTemp] > Data[Index[j - 1]])
							{
								Index[i - 1] = Index[j - 1];
								i = j;
								j += i;
							}
							else
								break;
						}

					else // SortType == stMin
						while (j < k + 2)
						{
							FailSafe++;
							if (j <= k && Data[Index[j - 1]] < Data[Index[j]])
								j++;
							if (Data[IndexTemp] < Data[Index[j - 1]])
							{
								Index[i - 1] = Index[j - 1];
								i = j;
								j += i;
							}
							else
								break;
						}

					Index[i - 1] = IndexTemp;
				}
				return (false);
			}
		}

		public static class LowPassRoots
		{
			private static readonly int MAX_ELLIP_ITER = 15;
			private static readonly int ELLIPARRAYSIZE = 20;
			private static readonly double PI_2 = PI / 2;

			private static double Factorial(int N)
			{
				int j;
				double Fact = 1.0;
				for (j = 1; j <= N; j++)
					Fact *= (double)j;
				return (Fact);
			}

			public static void ReverseCoeff(ref double[] P, int N)
			{
				int j;
				double Temp;
				for (j = 0; j <= N / 2; j++)
				{
					Temp = P[j];
					P[j] = P[N - j];
					P[N - j] = Temp;
				}

				for (j = N; j >= 1; j--)
				{
					if (P[0] != 0.0)
						P[j] /= P[0];
				}
				P[0] = 1.0;

			}

			public static int ButterworthPoly(int NumPoles, ref Complex[] Roots)
			{
				int j, n, N;
				double Theta;

				N = NumPoles;
				n = 0;
				for (j = 0; j < N / 2; j++)
				{
					Theta = PI * (double)(2 * j + N + 1) / (double)(2 * N);
					Roots[n++] = new Complex(Cos(Theta), Sin(Theta));
					Roots[n++] = new Complex(Cos(Theta), -Sin(Theta));
				}
				if (N % 2 == 1)
					Roots[n++] = new Complex(-1.0, 0.0); // The real root for odd pole counts.
				return (N);
			}

			public static int ChebyshevPoly(int NumPoles, double Ripple, ref Complex[] Roots)
			{
				int j, n, N;
				double Sigma, Omega;
				double Arg, Theta, Epsilon;

				N = NumPoles;
				Epsilon = Pow(10.0, Ripple / 10.0) - 1.0;
				Epsilon = Sqrt(Epsilon);
				if (Epsilon < 0.00001)
					Epsilon = 0.00001;
				if (Epsilon > 0.996)
					Epsilon = 0.996;
				Epsilon = 1.0 / Epsilon;
				Arg = Log(Epsilon + Sqrt(Epsilon * Epsilon + 1.0)) / (double)N; // = asinh(Epsilon) / (double)N;
				n = 0;
				for (j = 0; j < N / 2; j++)
				{
					Theta = (2 * j + 1) * PI_2 / (double)N;
					Sigma = -Sinh(Arg) * Sin(Theta);
					Omega = Cosh(Arg) * Cos(Theta);
					Roots[n++] = new Complex(Sigma, Omega);
					Roots[n++] = new Complex(Sigma, -Omega);
				}
				if (N % 2 == 1)
					Roots[n++] = new Complex(-Sinh(Arg), 0.0); // The real root for odd pole counts.
				return (N);
			}

			public static int GaussianPoly(int NumPoles, ref Complex[] Roots)
			{
				int j, N, RootsCount;
				double[] GaussCoeff = new double[102];

				N = NumPoles;
				GaussCoeff[0] = 1.0;
				GaussCoeff[1] = 0.0;
				for (j = 2; j <= 2 * N; j += 2)
				{
					GaussCoeff[j] = 1.0 / Factorial(j / 2);
					GaussCoeff[j + 1] = 0.0;
					if ((j / 2) % 2 == 1)
						GaussCoeff[j] *= -1.0;
				}

				// The coefficients are generated in reverse order needed for P51.
				ReverseCoeff(ref GaussCoeff, N * 2);
				RootsCount = P51OneRevC.FindRoots(N * 2, ref GaussCoeff, ref Roots);
				return (RootsCount);
			}

			public static int AdjustablePoly(int NumPoles, ref Complex[] Roots, double Gamma)
			{
				int j, N, RootsCount;
				double[] GaussCoeff = new double[P51OneRevC.P51_ARRAY_SIZE];

				N = NumPoles;
				if (Gamma > 0.0)
					Gamma *= 2.0; // Gamma < 0 is the orig Gauss and Bessel responses. Gamma > 0 has an asymptotic response, so we double it, which also makes the user interface a bit nicer. i.e. -1 <= Gamma <= 1

				GaussCoeff[0] = 1.0;
				GaussCoeff[1] = 0.0;
				for (j = 2; j <= 2 * N; j += 2)
				{
					GaussCoeff[j] = Pow(Factorial(j / 2), Gamma); // Gamma = -1 is orig Gauss poly, Gamma = 1 approaches a Butterworth response.
					GaussCoeff[j + 1] = 0.0;
					if ((j / 2) % 2 == 1)
						GaussCoeff[j] *= -1.0;
				}

				// The coefficients are generated in reverse order needed for P51.
				ReverseCoeff(ref GaussCoeff, N * 2);
				RootsCount = P51OneRevC.FindRoots(N * 2, ref GaussCoeff, ref Roots);

				// Scale the imag part of the root by 1.1 to get a response closer to a Butterworth when Gamma = -2
				for (j = 0; j < N * 2; j++)
					Roots[j] = new Complex(Roots[j].Real, Roots[j].Imaginary * 1.10);
				return (RootsCount);
			}

			public static int BesselPoly(int NumPoles, ref Complex[] Roots)
			{
				int k, N, RootsCount;
				double b;
				double[] PolyCoeff = new double[P51OneRevC.P51_ARRAY_SIZE];

				N = NumPoles;
				for (k = N - 1; k >= 0; k--)
				{
					// b is calc'd as a double because of all the division, but the result is essentially a large int.
					b = Factorial(2 * N - k) / Factorial(k) / Factorial(N - k) / Pow(2.0, (double)(N - k));
					PolyCoeff[k] = b;
				}
				PolyCoeff[N] = 1.0;

				// The coefficients are generated in reverse order needed for P51.
				ReverseCoeff(ref PolyCoeff, N);
				RootsCount = P51OneRevC.FindRoots(N, ref PolyCoeff, ref Roots);
				return (RootsCount);
			}

			public static int InvChebyPoly(int NumPoles, double StopBanddB, ref Complex[] ChebyPoles, ref Complex[] ChebyZeros, ref int ZeroCount)
			{
				int j, k, N, PolesCount;
				double Arg, Epsilon;
				double[] ChebPolyCoeff = new double[P51OneRevC.P51_ARRAY_SIZE];
				double[] PolyCoeff = new double[P51OneRevC.P51_ARRAY_SIZE];
				Complex[] SquaredPolyCoeff = new Complex[P51OneRevC.P51_ARRAY_SIZE];
				Complex A, B;

				N = NumPoles;
				Epsilon = 1.0 / (Pow(10.0, StopBanddB / 10.0) - 1.0);  // actually Epsilon Squared

				// This algorithm is from the paper by Richard J Mathar. It generates the coefficients for the Cheb poly.
				// It stores the Nth order coefficient in ChebPolyCoeff[N], and so on. Every other Cheb coeff is 0. See Wikipedia for a table that this code will generate.
				for (j = 0; j <= N / 2; j++)
				{
					Arg = Factorial(N - j - 1) / Factorial(j) / Factorial(N - 2 * j);
					if (j % 2 == 1)
						Arg *= -1.0;
					Arg *= Pow(2.0, (double)(N - 2 * j)) * (double)N / 2.0;
					ChebPolyCoeff[N - 2 * j] = Arg;
					ChebPolyCoeff[N - (2 * j + 1)] = 0.0;
				}


				// Now square the Chebshev polynomial where we assume s = jw.  To get the signs correct,
				// we need to take j to the Power. Then its a simple matter of adding Powers and
				// multiplying coefficients. j and k represent the exponents. That is, j=3 is the x^3 coeff, and so on.
				for (j = 0; j <= 2 * N; j++)
					SquaredPolyCoeff[j] = new Complex(0.0, 0.0);

				for (j = 0; j <= N; j++)
					for (k = 0; k <= N; k++)
					{
						A = (new Complex(0.0, 1.0)).Pow(((double)j) * ChebPolyCoeff[j]);
						B = (new Complex(0.0, 1.0)).Pow(((double)k) * ChebPolyCoeff[k]);
						SquaredPolyCoeff[j + k] = SquaredPolyCoeff[j + k] + A * B; // these end up entirely real.
					}

				// Denominator
				// Now we multiply the coefficients by Epsilon and add 1 to the denominator poly.
				k = 0;
				for (j = 0; j <= 2 * N; j++)
					ChebPolyCoeff[j] = SquaredPolyCoeff[j].Real * Epsilon;
				ChebPolyCoeff[0] += 1.0;
				for (j = 0; j <= 2 * N; j++)
					PolyCoeff[k++] = ChebPolyCoeff[j];  // Note this order is reversed from the Chebyshev routine.
				k--;
				PolesCount = P51OneRevC.FindRoots(k, ref PolyCoeff, ref ChebyPoles);


				// Numerator
				k = 0;
				for (j = 0; j <= 2 * N; j++)
					ChebPolyCoeff[j] = SquaredPolyCoeff[j].Real;  // Not using Epsilon here so the check for 0 on the next line is easier. Since the root finder normalizes the poly, it gets factored out anyway.
				for (j = 0; j <= 2 * N; j++)
					if (Abs(ChebPolyCoeff[j]) > 0.01)
						break;    // Get rid of the high order zeros. There will be eithe 0ne or two zeros to delete.
				for (; j <= 2 * N; j++)
					PolyCoeff[k++] = ChebPolyCoeff[j];
				k--;
				ZeroCount = P51OneRevC.FindRoots(k, ref PolyCoeff, ref ChebyZeros);

				return (PolesCount);

			}

			public static int PapoulisPoly(int NumPoles, ref Complex[] Roots)
			{
				int j, N, RootsCount;
				double Epsilon;
				double[] PolyCoeff = new double[P51OneRevC.P51_ARRAY_SIZE];

				N = NumPoles;
				for (j = 0; j < 2 * N; j++)
					PolyCoeff[j] = 0.0; // so we don't have to fill all the zero's.

				switch (N)
				{
					case 1:  // 1 pole
						PolyCoeff[2] = 1.0;
						break;

					case 2:  // 2 pole
						PolyCoeff[4] = 1.0;
						break;

					case 3:  // 3 pole
						PolyCoeff[6] = 3.0;
						PolyCoeff[4] = -3.0;
						PolyCoeff[2] = 1.0;
						break;

					case 4:
						PolyCoeff[8] = 6.0;
						PolyCoeff[6] = -8.0;
						PolyCoeff[4] = 3.0;
						break;

					case 5:
						PolyCoeff[10] = 20.0;
						PolyCoeff[8] = -40.0;
						PolyCoeff[6] = 28.0;
						PolyCoeff[4] = -8.0;
						PolyCoeff[2] = 1.0;
						break;

					case 6:
						PolyCoeff[12] = 50.0;
						PolyCoeff[10] = -120.0;
						PolyCoeff[8] = 105.0;
						PolyCoeff[6] = -40.0;
						PolyCoeff[4] = 6.0;
						break;

					case 7:
						PolyCoeff[14] = 175.0;
						PolyCoeff[12] = -525.0;
						PolyCoeff[10] = 615.0;
						PolyCoeff[8] = -355.0;
						PolyCoeff[6] = 105.0;
						PolyCoeff[4] = -15.0;
						PolyCoeff[2] = 1.0;
						break;

					case 8:
						PolyCoeff[16] = 490.0;
						PolyCoeff[14] = -1680.0;
						PolyCoeff[12] = 2310.0;
						PolyCoeff[10] = -1624.0;
						PolyCoeff[8] = 615.0;
						PolyCoeff[6] = -120.0;
						PolyCoeff[4] = 10.0;
						break;

					case 9:
						PolyCoeff[18] = 1764.0;
						PolyCoeff[16] = -7056.0;
						PolyCoeff[14] = 11704.0;
						PolyCoeff[12] = -10416.0;
						PolyCoeff[10] = 5376.0;
						PolyCoeff[8] = -1624.0;
						PolyCoeff[6] = 276.0;
						PolyCoeff[4] = -24.0;
						PolyCoeff[2] = 1.0;
						break;

					case 10:
						PolyCoeff[20] = 5292.0;
						PolyCoeff[18] = -23520.0;
						PolyCoeff[16] = 44100.0;
						PolyCoeff[14] = -45360.0;
						PolyCoeff[12] = 27860.0;
						PolyCoeff[10] = -10416.0;
						PolyCoeff[8] = 2310.0;
						PolyCoeff[6] = -280.0;
						PolyCoeff[4] = 15.0;
						break;

					case 11:
						PolyCoeff[22] = 19404;
						PolyCoeff[20] = -97020.0;
						PolyCoeff[18] = 208740.0;
						PolyCoeff[16] = -252840.0;
						PolyCoeff[14] = 189420.0;
						PolyCoeff[12] = -90804.0;
						PolyCoeff[10] = 27860.0;
						PolyCoeff[8] = -5320.0;
						PolyCoeff[6] = 595.0;
						PolyCoeff[4] = -35.0;
						PolyCoeff[2] = 1.0;
						break;

					case 12:
						PolyCoeff[24] = 60984.0;
						PolyCoeff[22] = -332640.0;
						PolyCoeff[20] = 790020.0;
						PolyCoeff[18] = -1071840.0;
						PolyCoeff[16] = 916020.0;
						PolyCoeff[14] = -512784.0;
						PolyCoeff[12] = 189420.0;
						PolyCoeff[10] = -45360.0;
						PolyCoeff[8] = 6720.0;
						PolyCoeff[6] = -560.0;
						PolyCoeff[4] = 21.0;
						break;

					case 13:
						PolyCoeff[26] = 226512.0;
						PolyCoeff[24] = -1359072.0;
						PolyCoeff[22] = 3597264.0;
						PolyCoeff[20] = -5528160.0;
						PolyCoeff[18] = 5462820.0;
						PolyCoeff[16] = -3632112.0;
						PolyCoeff[14] = 1652232.0;
						PolyCoeff[12] = -512784.0;
						PolyCoeff[10] = 106380.0;
						PolyCoeff[8] = -14160.0;
						PolyCoeff[6] = 1128.0;
						PolyCoeff[4] = -48.0;
						PolyCoeff[2] = 1.0;
						break;


					case 14:
						PolyCoeff[28] = 736164.0;
						PolyCoeff[26] = -4756752.0;
						PolyCoeff[24] = 13675662.0;
						PolyCoeff[22] = -23063040.0;
						PolyCoeff[20] = 25322220.0;
						PolyCoeff[18] = -18993744.0;
						PolyCoeff[16] = 9934617.0;
						PolyCoeff[14] = -3632112.0;
						PolyCoeff[12] = 916020.0;
						PolyCoeff[10] = -154560.0;
						PolyCoeff[8] = 16506.0;
						PolyCoeff[6] = -1008.0;
						PolyCoeff[4] = 28.0;
						break;

					case 15:
						PolyCoeff[30] = 2760615.0;
						PolyCoeff[28] = -19324305.0;
						PolyCoeff[26] = 60747687.0;
						PolyCoeff[24] = -113270157.0;
						PolyCoeff[22] = 139378239.0;
						PolyCoeff[20] = -119144025.0;
						PolyCoeff[18] = 72539775.0;
						PolyCoeff[16] = -31730787.0;
						PolyCoeff[14] = 9934617.0;
						PolyCoeff[12] = -2191959.0;
						PolyCoeff[10] = 331065.0;
						PolyCoeff[8] = -32655.0;
						PolyCoeff[6] = 1953.0;
						PolyCoeff[4] = -63.0;
						PolyCoeff[2] = 1.0;
						break;

					case 16:
						PolyCoeff[32] = 9202050.0;
						PolyCoeff[30] = -68708640.0;
						PolyCoeff[28] = 231891660.0;
						PolyCoeff[26] = -467747280.0;
						PolyCoeff[24] = 628221594.0;
						PolyCoeff[22] = -592431840.0;
						PolyCoeff[20] = 403062660.0;
						PolyCoeff[18] = -200142800.0;
						PolyCoeff[16] = 72539775.0;
						PolyCoeff[14] = -18993744.0;
						PolyCoeff[12] = 3515820.0;
						PolyCoeff[10] = -443520.0;
						PolyCoeff[8] = 35910.0;
						PolyCoeff[6] = -1680.0;
						PolyCoeff[4] = 36.0;
						break;

					case 17:
						PolyCoeff[34] = 34763300.0;
						PolyCoeff[32] = -278106400.0;
						PolyCoeff[30] = 1012634480.0;
						PolyCoeff[28] = -2221579360.0;
						PolyCoeff[26] = 3276433160.0;
						PolyCoeff[24] = -3431908480.0;
						PolyCoeff[22] = 2629731104.0;
						PolyCoeff[20] = -1496123200.0;
						PolyCoeff[18] = 634862800.0;
						PolyCoeff[16] = -200142800.0;
						PolyCoeff[14] = 46307800.0;
						PolyCoeff[12] = -7696304.0;
						PolyCoeff[10] = 888580.0;
						PolyCoeff[8] = -67760.0;
						PolyCoeff[6] = 3160.0;
						PolyCoeff[4] = -80.0;
						PolyCoeff[2] = 1.0;
						break;

					case 18:
						PolyCoeff[36] = 118195220.0;
						PolyCoeff[34] = -1001183040.0;
						PolyCoeff[32] = 3879584280.0;
						PolyCoeff[30] = -9110765664.0;
						PolyCoeff[28] = 14480345880.0;
						PolyCoeff[26] = -16474217760.0;
						PolyCoeff[24] = 13838184360.0;
						PolyCoeff[22] = -8725654080.0;
						PolyCoeff[20] = 4158224928.0;
						PolyCoeff[18] = -1496123200.0;
						PolyCoeff[16] = 403062660.0;
						PolyCoeff[14] = -79999920.0;
						PolyCoeff[12] = 11397540.0;
						PolyCoeff[10] = -1119888.0;
						PolyCoeff[8] = 71280.0;
						PolyCoeff[6] = -2640.0;
						PolyCoeff[4] = 45.0;
						break;

					case 19:
						PolyCoeff[38] = 449141836.0;
						PolyCoeff[36] = -4042276524.0;
						PolyCoeff[34] = 16732271556.0;
						PolyCoeff[32] = -42233237904.0;
						PolyCoeff[30] = 72660859128.0;
						PolyCoeff[28] = -90231621480.0;
						PolyCoeff[26] = 83545742280.0;
						PolyCoeff[24] = -58751550000.0;
						PolyCoeff[22] = 31671113760.0;
						PolyCoeff[20] = -13117232128.0;
						PolyCoeff[18] = 4158224928.0;
						PolyCoeff[16] = -999092952.0;
						PolyCoeff[14] = 178966788.0;
						PolyCoeff[12] = -23315292.0;
						PolyCoeff[10] = 2130876.0;
						PolyCoeff[8] = -129624.0;
						PolyCoeff[6] = 4851.0;
						PolyCoeff[4] = -99.0;
						PolyCoeff[2] = 1.0;
						break;

					case 20:
						PolyCoeff[40] = 1551580888.0;
						PolyCoeff[38] = -14699187360.0;
						PolyCoeff[36] = 64308944700.0;
						PolyCoeff[34] = -172355177280.0;
						PolyCoeff[32] = 316521742680.0;
						PolyCoeff[30] = -422089668000.0;
						PolyCoeff[28] = 422594051880.0;
						PolyCoeff[26] = -323945724960.0;
						PolyCoeff[24] = 192167478360.0;
						PolyCoeff[22] = -88572527680.0;
						PolyCoeff[20] = 31671113760.0;
						PolyCoeff[18] = -8725654080.0;
						PolyCoeff[16] = 1829127300.0;
						PolyCoeff[14] = -286125840.0;
						PolyCoeff[12] = 32458140.0;
						PolyCoeff[10] = -2560272.0;
						PolyCoeff[8] = 131670.0;
						PolyCoeff[6] = -3960.0;
						PolyCoeff[4] = 55.0;
						break;
				}
				Epsilon = 0.1; // This controls the amount of pass band roll off.  0.01 < Epsilon < 0.250

				// The poly is in terms of omega, but we need it in term of s = jw. So we need to
				// multiply the approp coeff by neg 1 to account for j. Then mult by epsilon.
				for (j = 0; j <= 2 * N; j++)
				{
					if ((j / 2) % 2 == 1)
						PolyCoeff[j] *= -1.0;
					PolyCoeff[j] *= Epsilon;
				}

				// Now add 1 to the poly.
				PolyCoeff[0] = 1.0;

				// The coefficients are in reverse order needed for P51.
				ReverseCoeff(ref PolyCoeff, N * 2);
				RootsCount = P51OneRevC.FindRoots(N * 2, ref PolyCoeff, ref Roots);

				return (RootsCount);
			}

			public static int EllipticPoly(int FiltOrder, double Ripple, double DesiredSBdB, ref Complex[] EllipPoles, ref Complex[] EllipZeros, ref int ZeroCount)
			{
				int j, k, n, LastK=0;
				double[] K = new double[ELLIPARRAYSIZE];
				double[] G = new double[ELLIPARRAYSIZE];
				double[] Epsilon = new double[ELLIPARRAYSIZE];
				double A, D, SBdB, dBErr, RealPart, ImagPart;
				double DeltaK=.99, PrevErr=Double.MaxValue, Deriv;
				Complex C;

				for (j = 0; j < ELLIPARRAYSIZE; j++)
					K[j] = G[j] = Epsilon[j] = 0.0;
				if (Ripple < 0.001)
					Ripple = 0.001;
				if (Ripple > 1.0)
					Ripple = 1.0;
				Epsilon[0] = Sqrt(Pow(10.0, Ripple / 10.0) - 1.0);

				// Estimate K[0] to get the algorithm started.
				K[0] = (double)(FiltOrder - 2) * 0.1605 + 0.016;
				if (K[0] < 0.01)
					K[0] = 0.01;
				if (K[0] > 0.7)
					K[0] = 0.7;

				// This loop calculates K[0] for the desired stopband attenuation. It typically loops < 5 times.
				for (j = 0; j < MAX_ELLIP_ITER; j++)
				{
					// Compute K with a forward Landen Transformation.
					for (k = 1; k < 10; k++)
					{
						K[k] = Pow(K[k - 1] / (1.0 + Sqrt(1.0 - K[k - 1] * K[k - 1])), 2.0);   // eq. 10
						if (K[k] <= 1.0E-6)
							break;
					}
					LastK = k;

					// Compute G with a backwards Landen Transformation.
					G[LastK] = 4.0 * Pow(K[LastK] / 4.0, (double)FiltOrder);
					for (k = LastK; k >= 1; k--)
					{
						G[k - 1] = 2.0 * Sqrt(G[k]) / (1.0 + G[k]);  // eq. 9
					}

					if (G[0] <= 0.0)
						G[0] = 1.0E-10;
					SBdB = 10.0 * Log10(1.0 + Pow(Epsilon[0] / G[0], 2.0)); // Current stopband attenuation dB
					dBErr = DesiredSBdB - SBdB;

					if (Abs(dBErr) < 0.1)
						break;
					if (j == 0) // Do this on the 1st loop so we can calc a derivative.
					{
						if (dBErr > 0)
							DeltaK = 0.005;
						else
							DeltaK = -0.005;
						PrevErr = dBErr;
					}
					else
					{
						// Use Newtons Method to adjust K[0].
						Deriv = (PrevErr - dBErr) / DeltaK;
						PrevErr = dBErr;
						if (Deriv == 0.0)
							break; // This happens when K[0] hits one of the limits set below.
						DeltaK = dBErr / Deriv;
						if (DeltaK > 0.1)
							DeltaK = 0.1;
						if (DeltaK < -0.1)
							DeltaK = -0.1;
					}
					K[0] -= DeltaK;
					if (K[0] < 0.001)
						K[0] = 0.001;  // must not be < 0.0
					if (K[0] > 0.990)
						K[0] = 0.990;  // if > 0.990 we get a pole in the RHP. This means we were unable to set the stop band atten to the desired level (the Ripple is too large for the Pole Count).
				}


				// Epsilon[0] was calulated above, now calculate Epsilon[LastK] from G
				for (j = 1; j <= LastK; j++)
				{
					A = (1.0 + G[j]) * Epsilon[j - 1] / 2.0;  // eq. 37
					Epsilon[j] = A + Sqrt(A * A + G[j]);
				}

				// Calulate the poles and zeros.
				ImagPart = Log((1.0 + Sqrt(1.0 + Epsilon[LastK] * Epsilon[LastK])) / Epsilon[LastK]) / (double)FiltOrder;  // eq. 22
				n = 0;
				for (j = 1; j <= FiltOrder / 2; j++)
				{
					RealPart = (double)(2 * j - 1) * PI / 2 / (double)FiltOrder;   // eq. 19
					C =(new Complex(0.0, -1.0) / (new Complex(-RealPart, ImagPart)).Cos());      // eq. 20
					D = 1.0 / Cos(RealPart);
					for (k = LastK; k >= 1; k--)
					{
						C = (C - K[k] / C) / (1.0 + K[k]);  // eq. 36
						D = (D + K[k] / D) / (1.0 + K[k]);
					}

					EllipPoles[n] = 1.0 / C;
					EllipPoles[n + 1] = EllipPoles[n].Conjugate();
					EllipZeros[n] = new Complex(0.0, D / K[0]);
					EllipZeros[n + 1] = EllipZeros[n].Conjugate();
					n += 2;
				}
				ZeroCount = n; // n is the num zeros

				if (FiltOrder % 2 == 1)   // The real pole for odd pole counts
				{
					A = 1.0 / Sinh(ImagPart);
					for (k = LastK; k >= 1; k--)
					{
						A = (A - K[k] / A) / (1.0 + K[k]);      // eq. 38
					}
					EllipPoles[n] = new Complex(-1.0 / A, 0.0);
					n++;
				}

				return (n); // n is the num poles. There will be 1 more pole than zeros for odd pole counts.

			}
		}

		public static class P51OneRevC
		{
			public enum TUpdateStatus { UPDATED, BAD_ANGLE, ZERO_DEL, DAMPER_ON, DAMPER_OFF };

			public static readonly int P51_ARRAY_SIZE = 102;
			public static readonly int P51_MAXDEGREE = 100;
			public static readonly double TINY_VALUE = 1.0E-30;
			public static readonly int MAX_NUM_K = 4;
			public static readonly int MAX_NUM_ANGLES = 180;
			public static readonly int REAL_ITER_MAX = 20;
			public static readonly int QUAD_ITER_MAX = 20;
			public static readonly double LDBL_EPSILON = 1.084202172485504434E-19;
			public static readonly double HUGE_VALUE = 1.0E200;
			private static int FirstDamperlIter;
			private static double PrevQPN;

			public static int FindRoots(int N, ref double[] Coeff, ref Complex[] Roots)
			{
				int j;
				double[] P = new double[P51_ARRAY_SIZE];
				double[] RealRoot = new double[P51_ARRAY_SIZE];
				double[] ImagRoot = new double[P51_ARRAY_SIZE];

				for (j = 0; j <= N; j++)
					P[j] = Coeff[j]; // double to long double
				N = PFiftyOne(ref P, N, ref RealRoot, ref ImagRoot);
				for (j = 0; j < N; j++)
					Roots[j] = new Complex(RealRoot[j], ImagRoot[j]); // long double to double
				return (N);
			}

			public static int PFiftyOne(ref double[] Coeff, int Degree, ref double[] RealRoot, ref double[] ImagRoot)
			{
				if (Degree > P51_MAXDEGREE || Degree < 0)
				{
					//ShowMessage("Poly Degree is out of range in the P51 Root Finder.");
					return (0);
				}

				TUpdateStatus UpdateStatus = TUpdateStatus.UPDATED;
				int N, NZ, j, Iter, AngleNumber, TypeOfK;
				double RealZero = 0, QuadX;
				double[] TUV = new double[P51_ARRAY_SIZE];
				double[] P, QuadQP, RealQP, QuadK, RealK, QK;

				N = Degree; // N is decremented as roots are found.

				P = new double[N + 2];
				QuadQP = new double[N + 2];
				RealQP = new double[N + 2];
				QuadK = new double[N + 2];
				RealK = new double[N + 2];
				QK = new double[N + 2];
				if (P == null || QuadQP == null || RealQP == null || QuadK == null || RealK == null || QK == null)
				{
					//ShowMessage("Memory not Allocated in PFiftyOne root finder.");
					return (0);
				}

				for (j = 0; j <= N; j++)
					P[j] = Coeff[j]; // Copy the coeff. P gets modified.
				for (j = 0; j < N; j++)
					RealRoot[j] = ImagRoot[j] = 0.0; // Init to zero, in case any leading or trailing zeros are removed.

				// Remove trailing zeros. A tiny P[N] relative to P[N-1] is not allowed.
				while (Abs(P[N]) <= TINY_VALUE * Abs(P[N - 1]) && N > 0)
				{
					N--;
				}

				// Remove leading zeros.
				while (P[0] == 0.0 && N > 0)
				{
					for (j = 0; j < N; j++)
						P[j] = P[j + 1];
					N--;
				}

				// P[0] must = 1
				if (P[0] != 1.0)
				{
					for (j = 1; j <= N; j++)
						P[j] /= P[0];
					P[0] = 1.0;
				}

				TypeOfK = 0;
				while (N > 4 && TypeOfK < MAX_NUM_K)
				{
					NZ = 0; // Num Zeros found. (Used in the loop controls below.)
					QuadX = Pow(Abs(P[N]), 1.0 / (double)N) / 2.0; // QuadX is used to init TUV

					for (TypeOfK = 0; TypeOfK < MAX_NUM_K && NZ == 0; TypeOfK++) // Iterate on the different possible QuadK inits.
					{
						for (AngleNumber = 0; AngleNumber < MAX_NUM_ANGLES && NZ == 0; AngleNumber++) // Iterate on the angle used to init TUV.
						{
							SetTUVandK(ref P, N, ref TUV, ref RealK, ref QuadK, QuadX, AngleNumber, TypeOfK); // Init TUV and both K's
							for (Iter = 0; Iter < N && NZ == 0; Iter++) // Allow N calls to QuadIterate for N*QUAD_ITER_MAX iterations, then try a different angle.
							{
								NZ = QuadIterate(Iter, ref P, ref QuadQP, ref QuadK, ref QK, N, ref TUV, ref UpdateStatus); // NZ = 2 for a pair of complex roots or 2 real roots.

								if (NZ == 0) // Try for a real root.
								{
									if (Abs(QuadK[N - 2]) > TINY_VALUE * Abs(P[N]))
										RealZero = -P[N] / QuadK[N - 2]; // This value gets refined by QuadIterate.
									else
										RealZero = 0.0;
									NZ = RealIterate(Iter, ref P, ref RealQP, ref RealK, ref QK, N, ref RealZero); // NZ = 1 for a single real root.
								}

								if (NZ == 0 && UpdateStatus == TUpdateStatus.BAD_ANGLE)
									break; // If RealIterate failed and UpdateTUV called this a bad angle, it's pointless to iterate further on this angle.

							} // Iter loop Note the use of NZ in the loop controls.
						} // AngleNumber loop
					} // TypeOfK loop



					// Done iterating. If NZ==0 at this point, we failed and will exit below.
					// Decrement N, and set P to the quotient QP. QP = P/TUV or QP = P/(x-RealZero)
					if (NZ == 2) // Store a pair of complex roots or 2 real roots.
					{
						j = Degree - N;
						QuadRoots(TUV, ref RealRoot, ref ImagRoot, j);
						N -= NZ;
						for (j = 0; j <= N; j++)
							P[j] = QuadQP[j];
						TypeOfK = 0;
					}

					if (NZ == 1) // Store a single real root
					{
						j = Degree - N;
						RealRoot[j] = RealZero;
						ImagRoot[j] = 0.0;
						N -= NZ;
						for (j = 0; j <= N; j++)
							P[j] = RealQP[j];
						TypeOfK = 0;
					}



					// Remove any trailing zeros on P. P[N] should never equal zero, but can approach zero
					// because of roundoff errors. If P[N] is zero or tiny relative to P[N-1], we take the hit,
					// and place a root at the origin. This needs to be checked, but virtually never happens.
					while (Abs(P[N]) <= TINY_VALUE * Abs(P[N - 1]) && N > 0)
					{
						j = Degree - N;
						RealRoot[j] = 0.0;
						ImagRoot[j] = 0.0;
						N--;
						//ShowMessage("Removed a zero at the origin.");
					}

				} // The outermost loop while(N > 2)

				// Done, except for the last 1 or 2 roots. If N isn't 1 or 2 at this point, we failed.
				if (N == 1)
				{
					j = Degree - N;
					RealRoot[j] = -P[1] / P[0];
					ImagRoot[j] = 0.0;
					return (Degree);
				}

				if (N == 2)
				{
					j = Degree - N;
					QuadRoots(P, ref RealRoot, ref ImagRoot, j);
					return (Degree);
				}

				if (N == 3)
				{
					j = Degree - N;
					CubicRoots(P, ref RealRoot, ref ImagRoot, j);
					return (Degree);
				}

				if (N == 4)
				{
					j = Degree - N;
					BiQuadRoots(P, ref RealRoot, ref ImagRoot, j);
					return (Degree);
				}

				// ShowMessage("The P51 root finder failed to converge on a solution.");
				return (0);

			}

			public static int QuadIterate(int P51_Iter, ref double[] P, ref double[] QP, ref double[] K, ref double[] QK, int N, ref double[] TUV, ref TUpdateStatus UpdateStatus)
			{
				int Iter;
				double Err, MinErr, ErrScalar, QKCheck;

				ErrScalar = 1.0 / (16.0 * Pow((double)N, 3.0) * Abs(P[N]));

				P51_Iter *= QUAD_ITER_MAX;
				Err = MinErr = 1.0E100;
				UpdateStatus = TUpdateStatus.UPDATED;
				QuadSynDiv(ref P, N, ref TUV, ref QP); // Init QP
				QuadSynDiv(ref K, N - 1, ref TUV, ref QK); // Init QK

				for (Iter = 0; Iter < QUAD_ITER_MAX; Iter++)
				{
					UpdateTUV(P51_Iter + Iter, ref P, N, ref QP, ref K, ref QK, ref TUV, ref UpdateStatus);
					if (UpdateStatus == TUpdateStatus.BAD_ANGLE)
					{
						return (0); // Failure, need a different angle.
					}

					Err = Abs(QP[N - 1]) + Abs(QP[N + 1]); // QP[N-1] & QP[N+1] are the remainder terms of P/TUV.
					Err *= ErrScalar; // Normalize the error.

					if (Err < LDBL_EPSILON)
					{
						return (2); // Success!! 2 roots have been found.
					}

					// ZERO_DEL means both DelU and DelV were ~ 0 in UpdateTUV which means the algorithm has stalled.
					// It might be stalled in a dead zone with large error, or stalled because it can't adjust u and v with a resolution fine enough to meet our convergence criteria.
					if (UpdateStatus == TUpdateStatus.ZERO_DEL)
					{
						if (Err < 4.0 * (double)N * LDBL_EPSILON) // Small error, this is the best we can do.
						{
							UpdateStatus = TUpdateStatus.UPDATED;
							return (2);
						}
						else // Large error, get a different angle
						{
							UpdateStatus = TUpdateStatus.BAD_ANGLE;
							return (0);
						}
					}

					QKCheck = Abs(QK[N - 2]) + Abs(QK[N]); // QK[N-2] & QK[N] are the remainder terms of K/TUV.
					QKCheck *= ErrScalar;

					// Huge values indicate the algorithm is diverging and overflow is imminent. This can indicate
					// a single real root, or that we simply need a different angle on TUV. This happens frequently.
					if (Err > HUGE_VALUE || QKCheck > HUGE_VALUE)
					{
						UpdateStatus = TUpdateStatus.BAD_ANGLE;
						return (0);
					}

					// Record our best result thus far. We turn on the damper in UpdateTUV if the errs increase.
					if (Err < MinErr)
					{
						UpdateStatus = TUpdateStatus.DAMPER_OFF;
						MinErr = Err;
					}
					else if (Iter > 2)
					{

						UpdateStatus = TUpdateStatus.DAMPER_ON;
					}
				}

				// If we get here, we didn't converge, but TUV is getting updated.
				// If RealIterate can't find a real zero, this function will get called again.
				UpdateStatus = TUpdateStatus.UPDATED;
				return (0);
			}

			public static void UpdateTUV(int Iter, ref double[] P, int N, ref double[] QP, ref double[] K, ref double[] QK, ref double[] TUV, ref TUpdateStatus UpdateStatus)
			{
				int j;
				double DelU, DelV, Denom;
				double E2, E3, E4, E5;

				if (Iter == 0)
					FirstDamperlIter = 0;  // Reset this static var.
				if (UpdateStatus == TUpdateStatus.DAMPER_ON)
					FirstDamperlIter = Iter;

				// Update K, unless E3 is tiny relative to E2. The algorithm will work its way out of a tiny E3.
				// These equations are from the Jenkins and Traub paper "A Three Stage Algorithm for Real Polynomials Using Quadratic Iteration" Equation 9.8
				E2 = QP[N] * QP[N] + TUV[1] * QP[N] * QP[N - 1] + TUV[2] * QP[N - 1] * QP[N - 1];
				E3 = QP[N] * QK[N - 1] + TUV[1] * QP[N] * QK[N - 2] + TUV[2] * QP[N - 1] * QK[N - 2];

				if (Abs(E3) * HUGE_VALUE > Abs(E2))
				{
					E2 /= -E3;
					for (j = 1; j <= N - 2; j++)
						K[j] = E2 * QK[j - 1] + QP[j]; // At covergence, E2 ~ 0, so K ~ QP.
				}
				else
				{
					for (j = 1; j <= N - 2; j++)
						K[j] = QK[j - 1];
				}
				K[0] = QP[0]; // QP[0] = 1.0 always
				K[N - 1] = 0.0;


				QuadSynDiv(ref K, N - 1, ref TUV, ref QK); // Update QK   QK = K/TUV

				// These equations are modified versions of those used in the original Jenkins Traub Fortran algorithm RealPoly, found at  www.netlib.org/toms/493
				E3 = QP[N] * QK[N - 1] + TUV[1] * QP[N] * QK[N - 2] + TUV[2] * QP[N - 1] * QK[N - 2];
				E4 = QK[N - 1] * QK[N - 1] + TUV[1] * QK[N - 1] * QK[N - 2] + TUV[2] * QK[N - 2] * QK[N - 2];
				E5 = QP[N - 1] * QK[N - 1] - QP[N] * QK[N - 2];

				Denom = E5 * K[N - 2] * TUV[2] + E4 * P[N];
				if (Abs(Denom) * HUGE_VALUE < Abs(P[N]))
				{
					UpdateStatus = TUpdateStatus.BAD_ANGLE;  // Denom is tiny, overflow is imminent, get a new angle.
					return;
				}

				// Calc DelU and DelV. If they are effectively zero, bump them by epsilon.
				DelU = E3 * K[N - 2] * TUV[2] / Denom;
				if (Abs(DelU) < LDBL_EPSILON * Abs(TUV[1]))
				{
					if (DelU < 0.0)
						DelU = -Abs(TUV[1]) * LDBL_EPSILON;
					else
						DelU = Abs(TUV[1]) * LDBL_EPSILON;
				}

				DelV = -E5 * K[N - 2] * TUV[2] * TUV[2] / Denom;
				if (Abs(DelV) < LDBL_EPSILON * Abs(TUV[2]))
				{
					if (DelV < 0.0)
						DelV = -Abs(TUV[2]) * LDBL_EPSILON;
					else
						DelV = Abs(TUV[2]) * LDBL_EPSILON;
				}

				// If we haven't converged by QUAD_ITER_MAX iters, we need to test DelU and DelV for effectiveness.
				if (Iter >= QUAD_ITER_MAX - 1)
				{
					// We can't improve u and v any further because both DelU and DelV ~ 0 This usually happens when we are near a solution, but we don't have the precision needed to ine u and v enough to meet our convergence criteria. This can also happen in a dead zone where the errors are large, which means we need a different angle on TUV. We test for this in the QuadIterate function.
					if (Abs(DelU) < 8.0 * LDBL_EPSILON * Abs(TUV[1]) && Abs(DelV) < 8.0 * LDBL_EPSILON * Abs(TUV[2]))
					{
						UpdateStatus = TUpdateStatus.ZERO_DEL;
						return;
					}
					// A big change after this many iterations means we are wasting our time on this angle.
					if (Abs(DelU) > 10.0 * Abs(TUV[1]) || Abs(DelV) > 10.0 * Abs(TUV[2]))
					{
						UpdateStatus = TUpdateStatus.BAD_ANGLE;
						return;
					}
				}

				// Dampen the changes for 3 iterations after Damper was set in QuadIterate.
				if (Iter - FirstDamperlIter < 3)
				{
					DelU *= 0.75;
					DelV *= 0.75;
				}

				// Update U and V
				TUV[1] += DelU;
				if (Abs(TUV[2] + DelV) < TINY_VALUE)
					DelV *= 0.9;  // If this, we would set TUV[2] = 0 which we can't allow, so we use 90% of DelV.
				TUV[2] += DelV;

				if (Abs(TUV[2]) < Abs(TUV[1]) * TINY_VALUE)
				{
					UpdateStatus = TUpdateStatus.BAD_ANGLE; // TUV[2] is effectively 0, which is never allowed.
					return;
				}

				UpdateStatus = TUpdateStatus.UPDATED;    // TUV was updated.
				QuadSynDiv(ref P, N, ref TUV, ref QP);  // Update QP  QP = P/TUV
			}

			public static int RealIterate(int P51_Iter, ref double[] P, ref double[] QP, ref double[] K, ref double[] QK, int N, ref double RealZero)
			{
				int Iter, k;
				double X, DelX, Damper, Err, ErrScalar;

				ErrScalar = 1.0 / (16.0 * Pow((double)N, 2.0) * Abs(P[N]));

				X = RealZero;        // Init with our best guess for X.
				if (P51_Iter == 0)
					PrevQPN = 0.0;
				QK[0] = K[0];
				for (k = 1; k <= N - 1; k++)
				{
					QK[k] = QK[k - 1] * X + K[k];
				}

				for (Iter = 0; Iter < REAL_ITER_MAX; Iter++)
				{
					// Calculate a new QP.  This is poly division  QP = P/(x+X)  QP[0] to QP[N-1] is the quotient.
					// The remainder is QP[N], which is P(X), the error term.
					QP[0] = P[0];
					for (k = 1; k <= N; k++)
					{
						QP[k] = QP[k - 1] * X + P[k];
					}
					Err = Abs(QP[N]) * ErrScalar; // QP[N] is the error. ErrScalar accounts for the wide variations in N and P[N].

					if (Err < LDBL_EPSILON)
					{
						RealZero = X;
						return (1);      // Success!!
					}
					else if (Err > HUGE_VALUE)
						return (0);  // Overflow is imminent.

					// Calculate a new K.  QK[N-1] is the remainder of K /(x-X).
					// QK[N-1] is approximately P'(X) when P(X) = QP[N] ~ 0
					if (Abs(QK[N - 1]) > Abs(P[N] * TINY_VALUE))
					{
						DelX = -QP[N] / QK[N - 1];
						K[0] = QP[0];
						for (k = 1; k <= N - 1; k++)
						{
							K[k] = DelX * QK[k - 1] + QP[k];
						}
					}
					else  // Use this if QK[N-1] ~ 0
					{
						K[0] = 0.0;
						for (k = 1; k <= N - 1; k++)
							K[k] = QK[k - 1];
					}

					if (Abs(K[N - 1]) > HUGE_VALUE)
						return (0); // Overflow is imminent.

					// Calculate a new QK.  This is poly division QK = K /(x+X).  QK[0] to QK[N-2] is the quotient.
					QK[0] = K[0];
					for (k = 1; k <= N - 1; k++)
					{
						QK[k] = QK[k - 1] * X + K[k];
					}
					if (Abs(QK[N - 1]) <= TINY_VALUE * Abs(P[N]))
						return (0);  // QK[N-1] ~ 0 will cause overflow below.

					// This algorithm routinely oscillates back and forth about a zero with nearly equal pos and neg error magnitudes.
					// If the current and previous error add to give a value less than the current error, we dampen the change.
					Damper = 1.0;
					if (Abs(QP[N] + PrevQPN) < Abs(QP[N]))
						Damper = 0.5;
					PrevQPN = QP[N];

					// QP[N] is P(X) and at convergence QK[N-1] ~ P'(X), so this is ~ Newtons Method.
					DelX = QP[N] / QK[N - 1] * Damper;
					if (X != 0.0)
					{
						if (Abs(DelX) < LDBL_EPSILON * Abs(X))  // If true, the algorithm is stalled, so bump X by 2 epsilons.
						{
							if (DelX < 0.0)
								DelX = -2.0 * X * LDBL_EPSILON;
							else
								DelX = 2.0 * X * LDBL_EPSILON;
						}
					}
					else  // X = 0
					{
						if (DelX == 0.0)
							return (0); // Stalled at the origin, so exit.
					}
					X -= DelX;  // Update X

				} // end of loop

				// If we get here, we failed to converge.
				return (0);
			}

			public static void DerivOfP(double[] P, int N, ref double[] dP)
			{
				int j;
				double Power;
				for (j = 0; j < N; j++)
				{
					Power = N - j;
					dP[j] = Power * P[j];
				}
				dP[N] = 0.0;

			}

			public static void QuadSynDiv(ref double[] P, int N, ref double[] TUV, ref double[] Q)
			{
				int j;
				Q[0] = P[0];
				Q[1] = P[1] - TUV[1] * Q[0];
				for (j = 2; j <= N; j++)
				{
					Q[j] = P[j] - TUV[1] * Q[j - 1] - TUV[2] * Q[j - 2];
				}

				// Here we calculate the final remainder term used to calculate the convergence error.
				// This and Q[N-1] are the remainder values you get if you do this poly division manually.
				Q[N + 1] = Q[N - 1] * TUV[1] + Q[N];  // = b*u + a
			}

			public static void SetTUVandK(ref double[] P, int N, ref double[] TUV, ref double[] RealK, ref double[] QuadK, double X, int AngleNumber, int TypeOfQuadK)
			{
				int j;
				double a, Theta;

				// These angles define our search pattern in the complex plane. We start in the 1st quadrant,
				// go to the 2nd, then the real axis, etc. The roots are conjugates, so there isn't a need
				// to go to the 3rd or 4th quadrants. The first 2 angles find about 99% of all roots.
				double[] Angle = new double[]
				 {45.0, 135.0, 0.0, 90.0, 15.0, 30.0, 60.0, 75.0, 105.0, 120.0, 150.0, 165.0,  // 12 angles
				  6.0, 51.0, 96.0, 141.0, 12.0, 57.0, 102.0, 147.0, 21.0, 66.0, 111.0, 156.0,
				  27.0, 72.0, 117.0, 162.0, 36.0, 81.0, 126.0, 171.0, 42.0, 87.0, 132.0, 177.0,
				  3.0, 48.0, 93.0, 138.0, 9.0, 54.0, 99.0, 144.0, 18.0, 63.0, 108.0, 153.0,
				  24.0, 69.0, 114.0, 159.0, 33.0, 78.0, 123.0, 168.0, 39.0, 84.0, 129.0, 174.0,  // 60 angles
				  46.0, 136.0, 91.0, 1.0, 16.0, 31.0, 61.0, 76.0, 106.0, 121.0, 151.0, 166.0,
				  7.0, 52.0, 97.0, 142.0, 13.0, 58.0, 103.0, 148.0, 22.0, 67.0, 112.0, 157.0,
				  28.0, 73.0, 118.0, 163.0, 37.0, 82.0, 127.0, 172.0, 43.0, 88.0, 133.0, 178.0,
				  4.0, 49.0, 94.0, 139.0, 10.0, 55.0, 100.0, 145.0, 19.0, 64.0, 109.0, 154.0,
				  25.0, 70.0, 115.0, 160.0, 34.0, 79.0, 124.0, 169.0, 40.0, 85.0, 130.0, 175.0,
				  47.0, 137.0, 92.0, 2.0, 17.0, 32.0, 62.0, 77.0, 107.0, 122.0, 152.0, 167.0,
				  8.0, 53.0, 98.0, 143.0, 14.0, 59.0, 104.0, 149.0, 23.0, 68.0, 113.0, 158.0,
				  29.0, 74.0, 119.0, 164.0, 38.0, 83.0, 128.0, 173.0, 44.0, 89.0, 134.0, 179.0,
				  5.0, 50.0, 95.0, 140.0, 11.0, 56.0, 101.0, 146.0, 20.0, 65.0, 110.0, 155.0,
				  26.0, 71.0, 116.0, 161.0, 35.0, 80.0, 125.0, 170.0, 41.0, 86.0, 131.0, 176.0 }; // 180 angles


				// Initialize TUV to form  (x - (a + jb)) * (x - (a - jb)) =  x^2 - 2ax + a^2 + b^2
				// We init TUV for complex roots, except at angle 0, where we use real roots at +/- X
				if (AngleNumber == 2) // At 0 degrees we int to 2 real roots at +/- X.
				{
					TUV[0] = 1.0;   // t
					TUV[1] = 0.0;   // u
					TUV[2] = -(X * X);  // v
				}
				else  // We init to a complex root at  a +/- jb
				{
					Theta = Angle[AngleNumber] / 180.0 * PI;
					a = X * Cos(Theta);
					//b = X * sinl(Theta);
					TUV[0] = 1.0;
					TUV[1] = -2.0 * a;
					TUV[2] = X * X;   // = a*a + b*b because Cos^2 + sin^2 = 1
				}

				// The code below initializes the K polys used by RealIterate and QuadIterate.

				// Initialize the K poly used in RealIterate to P'.
				DerivOfP(P, N, ref RealK);
				RealK[N] = 0.0;


				// Initialize the K poly used in QuadIterate. Initializing QuadK to P" works virtually
				// 100% of the time, but if P51 stalls on a difficult root, these alternate inits give
				// us a way to work out of the stall. All these inits work almost as well as P".
				if (TypeOfQuadK == 0)  // Init QuadK 2nd derivative of P
				{
					double[] Temp = new double[N + 2];
					if (Temp == null)
					{
						//ShowMessage("Memory not Allocated in PFiftyOne SetTUVandK.");
						return;
					}

					DerivOfP(P, N, ref Temp);
					DerivOfP(Temp, N - 1, ref QuadK);
					QuadK[N] = QuadK[N - 1] = 0.0;
				}

				else if (TypeOfQuadK == 1) // Set QuadK to QP, because K = QP at convergence.
				{
					QuadSynDiv(ref P, N, ref TUV, ref QuadK);
					QuadK[N] = QuadK[N - 1] = 0.0;
				}

				else if (TypeOfQuadK == 2) // Set QuadK to the 1st derivative of P
				{
					for (j = 0; j <= N - 2; j++)
						QuadK[j] = RealK[j + 1];
					QuadK[N] = QuadK[N - 1] = 0.0;
				}

				else  // Set QuadK to zero, except QuadK[0].
				{
					for (j = 1; j <= N; j++)
						QuadK[j] = 0.0;
					QuadK[0] = 1.0;
				}

				if (QuadK[0] == 0.0)
					QuadK[0] = 1.0; // This can happen if TypeOfQuadK == 2 and P[1] == 0.0
				for (j = N - 2; j > 0; j--)
					QuadK[j] /= QuadK[0];
				QuadK[0] = 1.0;
			}
		}

		public static class FFTCode
		{
			private static readonly int MINIMUM_FFT_SIZE = 8;
			private static readonly int MAXIMUM_FFT_SIZE = 1048576;
			private static readonly double M_Sqrt_2 = Math.Sqrt(2)/2;

			public enum TWindowType
			{
				wtFIRSTWINDOW, wtNONE, wtKAISER, wtSINC, wtHANNING,
				wtHAMMING, wtBLACKMAN, wtFLATTOP, wtBLACKMAN_HARRIS,
				wtBLACKMAN_NUTTALL, wtNUTTALL, wtKAISER_BESSEL, wtTRAPEZOID,
				wtGAUSS, wtSINE, wtTEST
			};

			public enum TTransFormType { FORWARD, INVERSE };
			//---------------------------------------------------------------------------
			// This calculates the required FFT size for a given number of points.
			public static int RequiredFFTSize(int NumPts)
			{
				int N = MINIMUM_FFT_SIZE;
				while (N < NumPts && N < MAXIMUM_FFT_SIZE)
				{
					N *= 2;
				}
				return N;
			}

			//---------------------------------------------------------------------------

			// This verifies that the FFT Size N = 2^M.   M is returned
			// N must be >= 8 for the Twiddle calculations
			public static int IsValidFFTSize(int N)
			{
				if (N < MINIMUM_FFT_SIZE || N > MAXIMUM_FFT_SIZE || (N & (N - 1)) != 0)
					return (0);   // N & (N - 1) ensures a Power of 2
				return ((int)(Log((double)N) / Log(2) + 0.5));         // return M where N = 2^M
			}

			//---------------------------------------------------------------------------

			// Fast Fourier Transform
			// This code puts DC in bin 0 and scales the output of a forward transform by 1/N.
			// InputR and InputI are the real and imaginary input arrays of length N.
			// The output values are returned in the Input arrays.
			// TTransFormType is either FORWARD or INVERSE (defined in the header file)
			// 256 pts in 50 us
			public static void FFT(ref double[] InputR, ref double[] InputI, int N, TTransFormType Type)
			{
				int j, LogTwoOfN;
				int[] RevBits;
				double[] BufferR, BufferI, TwiddleR, TwiddleI;
				double OneOverN;

				// Verify the FFT size and type.
				LogTwoOfN = IsValidFFTSize(N);
				if (LogTwoOfN == 0 || (Type != TTransFormType.FORWARD && Type != TTransFormType.INVERSE))
				{
					// ShowMessage("Invalid FFT type or size.");
					return;
				}

				// Memory allocation for all the arrays.
				BufferR = new double[N];
				BufferI = new double[N];
				TwiddleR = new double[N / 2];
				TwiddleI = new double[N / 2];
				RevBits = new int[N];

				if (BufferR == null || BufferI == null ||
				   TwiddleR == null || TwiddleI == null || RevBits == null)
				{
					// ShowMessage("FFT Memory Allocation Error");
					return;
				}

				ReArrangeInput(ref InputR, ref InputI, ref BufferR, ref BufferI, ref RevBits, N);
				FillTwiddleArray(ref TwiddleR, ref TwiddleI, N, Type);
				Transform(ref InputR, ref InputI, ref BufferR, ref BufferI, ref TwiddleR, ref TwiddleI, N);


				// The ReArrangeInput function swapped Input[] and Buffer[]. Then Transform()
				// swapped them again, LogTwoOfN times. Ultimately, these swaps must be done
				// an even number of times, or the pointer to Buffer gets returned.
				// So we must do one more swap here, for N = 16, 64, 256, 1024, ...
				OneOverN = 1.0;
				if (Type == TTransFormType.FORWARD)
					OneOverN = 1.0 / (double)N;

				if (LogTwoOfN % 2 == 1)
				{
					for (j = 0; j < N; j++)
						InputR[j] = InputR[j] * OneOverN;
					for (j = 0; j < N; j++)
						InputI[j] = InputI[j] * OneOverN;
				}
				else // if(LogTwoOfN % 2 == 0) then the results are still in Buffer.
				{
					for (j = 0; j < N; j++)
						InputR[j] = BufferR[j] * OneOverN;
					for (j = 0; j < N; j++)
						InputI[j] = BufferI[j] * OneOverN;
				}
			}
			//---------------------------------------------------------------------------

			// This puts the input arrays in bit reversed order.
			// The while loop generates an array of bit reversed numbers. e.g.
			// e.g. N=8: RevBits = 0,4,2,6,1,5,3,7   N=16: RevBits = 0,8,4,12,2,10,6,14,1,9,5,13,3,11,7,15
			public static void ReArrangeInput(ref double[] InputR, ref double[] InputI, ref double[] BufferR, ref double[] BufferI, ref int[] RevBits, int N)
			{
				int j, k, J, K;

				J = N / 2;
				K = 1;
				RevBits[0] = 0;
				while (J >= 1)
				{
					for (k = 0; k < K; k++)
					{
						RevBits[k + K] = RevBits[k] + J;
					}
					K *= 2;
					J /= 2;
				}


				// Move the rearranged input values to Buffer.
				// Take note of the pointer swaps at the top of the transform algorithm.
				for (j = 0; j < N; j++)
				{
					BufferR[j] = InputR[RevBits[j]];
					BufferI[j] = InputI[RevBits[j]];
				}

			}

			//---------------------------------------------------------------------------

			/*
			 The Pentium takes a surprising amount of time to calculate the sine and Cosine.
			 You may want to make the twiddle arrays static if doing repeated FFTs of the same size.
			 This uses 4 fold symmetry to calculate the twiddle factors. As a result, this function
			 requires a minimum FFT size of 8.
			*/
			public static void FillTwiddleArray(ref double[] TwiddleR, ref double[] TwiddleI, int N, TTransFormType Type)
			{
				int j;
				double Theta, TwoPiOverN;

				TwoPiOverN = 2 * PI / (double)N;

				if (Type == TTransFormType.FORWARD)
				{
					TwiddleR[0] = 1.0;
					TwiddleI[0] = 0.0;
					TwiddleR[N / 4] = 0.0;
					TwiddleI[N / 4] = -1.0;
					TwiddleR[N / 8] = M_Sqrt_2;
					TwiddleI[N / 8] = -M_Sqrt_2;
					TwiddleR[3 * N / 8] = -M_Sqrt_2;
					TwiddleI[3 * N / 8] = -M_Sqrt_2;
					for (j = 1; j < N / 8; j++)
					{
						Theta = (double)j * -TwoPiOverN;
						TwiddleR[j] = Cos(Theta);
						TwiddleI[j] = Sin(Theta);
						TwiddleR[N / 4 - j] = -TwiddleI[j];
						TwiddleI[N / 4 - j] = -TwiddleR[j];
						TwiddleR[N / 4 + j] = TwiddleI[j];
						TwiddleI[N / 4 + j] = -TwiddleR[j];
						TwiddleR[N / 2 - j] = -TwiddleR[j];
						TwiddleI[N / 2 - j] = TwiddleI[j];
					}
				}

				else
				{
					TwiddleR[0] = 1.0;
					TwiddleI[0] = 0.0;
					TwiddleR[N / 4] = 0.0;
					TwiddleI[N / 4] = 1.0;
					TwiddleR[N / 8] = M_Sqrt_2;
					TwiddleI[N / 8] = M_Sqrt_2;
					TwiddleR[3 * N / 8] = -M_Sqrt_2;
					TwiddleI[3 * N / 8] = M_Sqrt_2;
					for (j = 1; j < N / 8; j++)
					{
						Theta = (double)j * TwoPiOverN;
						TwiddleR[j] = Cos(Theta);
						TwiddleI[j] = Sin(Theta);
						TwiddleR[N / 4 - j] = TwiddleI[j];
						TwiddleI[N / 4 - j] = TwiddleR[j];
						TwiddleR[N / 4 + j] = -TwiddleI[j];
						TwiddleI[N / 4 + j] = TwiddleR[j];
						TwiddleR[N / 2 - j] = -TwiddleR[j];
						TwiddleI[N / 2 - j] = TwiddleI[j];
					}
				}

			}

			//---------------------------------------------------------------------------

			// The Fast Fourier Transform.
			public static void Transform(ref double[] InputR, ref double[] InputI, ref double[] BufferR, ref double[] BufferI, ref double[] TwiddleR, ref double[] TwiddleI, int N)
			{
				int j, k, J, K, I, T;
				double[] TempPointer;
				double TempR, TempI;

				J = N / 2;     // J increments down to 1
				K = 1;       // K increments up to N/2
				while (J > 0) // Loops Log2(N) times.
				{
					// Swap pointers, instead doing this: for(j=0; j<N; j++) Input[j] = Buffer[j];
					// We start with a swap because of the swap in ReArrangeInput.
					TempPointer = InputR;
					InputR = BufferR;
					BufferR = TempPointer;
					TempPointer = InputI;
					InputI = BufferI;
					BufferI = TempPointer;

					I = 0;
					for (j = 0; j < J; j++)
					{
						T = 0;
						for (k = 0; k < K; k++) // Loops N/2 times for every value of J and K
						{
							TempR = InputR[K + I] * TwiddleR[T] - InputI[K + I] * TwiddleI[T];
							TempI = InputR[K + I] * TwiddleI[T] + InputI[K + I] * TwiddleR[T];
							BufferR[I] = InputR[I] + TempR;
							BufferI[I] = InputI[I] + TempI;
							BufferR[I + K] = InputR[I] - TempR;
							BufferI[I + K] = InputI[I] - TempI;
							I++;
							T += J;
						}
						I += K;
					}
					K *= 2;
					J /= 2;
				}

			}

			//-----------------------------------------------------------------------------------------------

			/*
			 The only difficulty in writing an FFT is figuring out how to calculate the various array indexes.
			 This shows how the index values change when doing a 16 pt FFT.
			 This print only has value if you compare it to a butterfly chart. Then you can
			 see how an FFT works. Use a 16 point decimation in time butterfly chart. We have one here:
			 http://www.iowahills.com/FFTCode.html

			 Note: The code above uses real variables. This print out came from code using complex variables as shown here.
			 Buffer[I]   = Input[I] + Input[K+I] * Twiddle[T];
			 Buffer[I+K] = Input[I] - Input[K+I] * Twiddle[T];

			 N = 16
			 J = 8    K = 1
			 Buffer[0]  = Input[0]  + Input[1]  * Twiddle[0]   I = 0
			 Buffer[1]  = Input[0]  - Input[1]  * Twiddle[0]   I = 0
			 Buffer[2]  = Input[2]  + Input[3]  * Twiddle[0]   I = 2
			 Buffer[3]  = Input[2]  - Input[3]  * Twiddle[0]   I = 2
			 Buffer[4]  = Input[4]  + Input[5]  * Twiddle[0]   etc.
			 Buffer[5]  = Input[4]  - Input[5]  * Twiddle[0]
			 Buffer[6]  = Input[6]  + Input[7]  * Twiddle[0]
			 Buffer[7]  = Input[6]  - Input[7]  * Twiddle[0]
			 Buffer[8]  = Input[8]  + Input[9]  * Twiddle[0]
			 Buffer[9]  = Input[8]  - Input[9]  * Twiddle[0]
			 Buffer[10] = Input[10] + Input[11] * Twiddle[0]
			 Buffer[11] = Input[10] - Input[11] * Twiddle[0]
			 Buffer[12] = Input[12] + Input[13] * Twiddle[0]
			 Buffer[13] = Input[12] - Input[13] * Twiddle[0]
			 Buffer[14] = Input[14] + Input[15] * Twiddle[0]
			 Buffer[15] = Input[14] - Input[15] * Twiddle[0]

			 J = 4    K = 2
			 Buffer[0]  = Input[0]  + Input[2]  * Twiddle[0]
			 Buffer[2]  = Input[0]  - Input[2]  * Twiddle[0]
			 Buffer[1]  = Input[1]  + Input[3]  * Twiddle[4]
			 Buffer[3]  = Input[1]  - Input[3]  * Twiddle[4]
			 Buffer[4]  = Input[4]  + Input[6]  * Twiddle[0]
			 Buffer[6]  = Input[4]  - Input[6]  * Twiddle[0]
			 Buffer[5]  = Input[5]  + Input[7]  * Twiddle[4]
			 Buffer[7]  = Input[5]  - Input[7]  * Twiddle[4]
			 Buffer[8]  = Input[8]  + Input[10] * Twiddle[0]
			 Buffer[10] = Input[8]  - Input[10] * Twiddle[0]
			 Buffer[9]  = Input[9]  + Input[11] * Twiddle[4]
			 Buffer[11] = Input[9]  - Input[11] * Twiddle[4]
			 Buffer[12] = Input[12] + Input[14] * Twiddle[0]
			 Buffer[14] = Input[12] - Input[14] * Twiddle[0]
			 Buffer[13] = Input[13] + Input[15] * Twiddle[4]
			 Buffer[15] = Input[13] - Input[15] * Twiddle[4]

			 J = 2    K = 4
			 Buffer[0]  = Input[0]  + Input[4]  * Twiddle[0]
			 Buffer[4]  = Input[0]  - Input[4]  * Twiddle[0]
			 Buffer[1]  = Input[1]  + Input[5]  * Twiddle[2]
			 Buffer[5]  = Input[1]  - Input[5]  * Twiddle[2]
			 Buffer[2]  = Input[2]  + Input[6]  * Twiddle[4]
			 Buffer[6]  = Input[2]  - Input[6]  * Twiddle[4]
			 Buffer[3]  = Input[3]  + Input[7]  * Twiddle[6]
			 Buffer[7]  = Input[3]  - Input[7]  * Twiddle[6]
			 Buffer[8]  = Input[8]  + Input[12] * Twiddle[0]
			 Buffer[12] = Input[8]  - Input[12] * Twiddle[0]
			 Buffer[9]  = Input[9]  + Input[13] * Twiddle[2]
			 Buffer[13] = Input[9]  - Input[13] * Twiddle[2]
			 Buffer[10] = Input[10] + Input[14] * Twiddle[4]
			 Buffer[14] = Input[10] - Input[14] * Twiddle[4]
			 Buffer[11] = Input[11] + Input[15] * Twiddle[6]
			 Buffer[15] = Input[11] - Input[15] * Twiddle[6]

			 J = 1    K = 8
			 Buffer[0]  = Input[0]  + Input[8]  * Twiddle[0]
			 Buffer[8]  = Input[0]  - Input[8]  * Twiddle[0]
			 Buffer[1]  = Input[1]  + Input[9]  * Twiddle[1]
			 Buffer[9]  = Input[1]  - Input[9]  * Twiddle[1]
			 Buffer[2]  = Input[2]  + Input[10] * Twiddle[2]
			 Buffer[10] = Input[2]  - Input[10] * Twiddle[2]
			 Buffer[3]  = Input[3]  + Input[11] * Twiddle[3]
			 Buffer[11] = Input[3]  - Input[11] * Twiddle[3]
			 Buffer[4]  = Input[4]  + Input[12] * Twiddle[4]
			 Buffer[12] = Input[4]  - Input[12] * Twiddle[4]
			 Buffer[5]  = Input[5]  + Input[13] * Twiddle[5]
			 Buffer[13] = Input[5]  - Input[13] * Twiddle[5]
			 Buffer[6]  = Input[6]  + Input[14] * Twiddle[6]
			 Buffer[14] = Input[6]  - Input[14] * Twiddle[6]
			 Buffer[7]  = Input[7]  + Input[15] * Twiddle[7]
			 Buffer[15] = Input[7]  - Input[15] * Twiddle[7]

			*/

			//-----------------------------------------------------------------------------------------------


			// Discrete Fourier Transform ( textbook code )
			// This takes the same arguments as the FFT function.
			// 256 pts in 1.720 ms
			public static void DFT(ref double[] InputR, ref double[] InputI, int N, TTransFormType Type)
			{
				int j, k, n;
				double[] SumR, SumI;
				double Sign, Arg;
				double[] TwiddleReal, TwiddleImag;

				SumR = new double[N];
				SumI = new double[N];
				TwiddleReal = new double[N];
				TwiddleImag = new double[N];

				if (SumR == null || SumI == null ||
				   TwiddleReal == null || TwiddleImag == null ||
				   (Type != TTransFormType.FORWARD && Type != TTransFormType.INVERSE))
				{
					// ShowMessage("Incorrect DFT Type or unable to allocate memory");
					return;
				}

				// Calculate the twiddle factors and initialize the Sum arrays.
				if (Type == TTransFormType.FORWARD)
					Sign = -1.0;
				else
					Sign = 1.0;

				for (j = 0; j < N; j++)
				{
					Arg = 2 * PI * (double)j / (double)N;
					TwiddleReal[j] = Cos(Arg);
					TwiddleImag[j] = Sign * Sin(Arg);
					SumR[j] = SumI[j] = 0.0;
				}


				//Calculate the DFT
				for (j = 0; j < N; j++) // Sum index
					for (k = 0; k < N; k++) // Input index
					{
						n = (j * k) % N;
						SumR[j] += TwiddleReal[n] * InputR[k] - TwiddleImag[n] * InputI[k];
						SumI[j] += TwiddleReal[n] * InputI[k] + TwiddleImag[n] * InputR[k];
					}

				// Scale the result if doing a forward DFT, and move the result to the input arrays.
				if (Type == TTransFormType.FORWARD)
				{
					for (j = 0; j < N; j++)
						InputR[j] = SumR[j] / (double)N;
					for (j = 0; j < N; j++)
						InputI[j] = SumI[j] / (double)N;
				}
				else  // Inverse DFT
				{
					for (j = 0; j < N; j++)
						InputR[j] = SumR[j];
					for (j = 0; j < N; j++)
						InputI[j] = SumI[j];
				}
			}

			//-----------------------------------------------------------------------------------------------

			// This is a DFT for real valued samples. Since Samples is real, it can only do a forward DFT.
			// The results are returned in OutputR  OutputI
			// 256 pts in 700 us    30 us to calc the twiddles
			public static void RealSigDFT(ref double[] Samples, ref double[] OutputR, ref double[] OutputI, int N)
			{
				int j, k;
				double Arg;
				double[] TwiddleReal,TwiddleImag;

				TwiddleReal = new double[N];
				TwiddleImag = new double[N];
				if (TwiddleReal == null || TwiddleImag == null)
				{
					// ShowMessage("Failed to allocate memory in RealSigDFT");
					return;
				}

				for (j = 0; j < N; j++)
				{
					Arg = 2 * PI * (double)j / (double)N;
					TwiddleReal[j] = Cos(Arg);
					TwiddleImag[j] = -Sin(Arg);
				}


				// Compute the DFT.
				// We have a real input, so only do the pos frequencies. i.e. j<N/2
				for (j = 0; j <= N / 2; j++)
				{
					OutputR[j] = 0.0;
					OutputI[j] = 0.0;
					for (k = 0; k < N; k++)
					{
						OutputR[j] += Samples[k] * TwiddleReal[(j * k) % N];
						OutputI[j] += Samples[k] * TwiddleImag[(j * k) % N];
					}

					// Scale the result
					OutputR[j] /= (double)N;
					OutputI[j] /= (double)N;
				}

				// The neg freq components are the conj of the pos components because the input signal is real.
				for (j = 1; j < N / 2; j++)
				{
					OutputR[N - j] = OutputR[j];
					OutputI[N - j] = -OutputI[j];
				}
			}

			//---------------------------------------------------------------------------

			// This is a single frequency DFT.
			// This code uses iteration to calculate the Twiddle factors.
			// To evaluate the frequency response of an FIR filter at Omega, set
			// Samples[] = FirCoeff[]   N = NumTaps  0.0 <= Omega <= 1.0
			// 256 pts in 15.6 us
			public static double SingleFreqDFT(ref double[] Samples, int N, double Omega)
			{
				int k;
				double SumR, SumI, zR, zI, TwiddleR, TwiddleI, Temp;

				TwiddleR = Cos(Omega * PI);
				TwiddleI = -Sin(Omega * PI);
				zR = 1.0;    // z, as in e^(j*omega)
				zI = 0.0;
				SumR = 0.0;
				SumI = 0.0;

				for (k = 0; k < N; k++)
				{
					SumR += Samples[k] * zR;
					SumI += Samples[k] * zI;

					// Calculate the complex exponential z by taking it to the kth Power.
					Temp = zR * TwiddleR - zI * TwiddleI;
					zI = zR * TwiddleI + zI * TwiddleR;
					zR = Temp;
				}

				/*
				// This is the more conventional implementation of the loop above.
				// It is a bit more accurate, but slower.
				for(k=0; k<N; k++)
				 {
				  SumR += Samples[k] *  Cos((double)k * Omega * PI);
				  SumI += Samples[k] * -Sin((double)k * Omega * PI);
				 }
				*/

				return (Sqrt(SumR * SumR + SumI * SumI));
				// return( ComplexD(SumR, SumI) );// if phase is needed.
			}

			//---------------------------------------------------------------------------

			// Goertzel is essentially a single frequency DFT, but without phase information.
			// Its simplicity allows it to run about 3 times faster than a single frequency DFT.
			// It is typically used to find a tone embedded in a signal. A DTMF tone for example.
			// 256 pts in 6 us
			public static double Goertzel(ref double[] Samples, int N, double Omega)
			{
				int j;
				double Reg0, Reg1, Reg2;        // 3 shift registers
				double CosVal, Mag;
				Reg1 = Reg2 = 0.0;

				CosVal = 2.0 * Cos(PI * Omega);
				for (j = 0; j < N; j++)
				{
					Reg0 = Samples[j] + CosVal * Reg1 - Reg2;
					Reg2 = Reg1;  // Shift the values.
					Reg1 = Reg0;
				}
				Mag = Reg2 * Reg2 + Reg1 * Reg1 - CosVal * Reg1 * Reg2;

				if (Mag > 0.0)
					Mag = Sqrt(Mag);
				else
					Mag = 1.0E-12;

				return (Mag);
			}

			//---------------------------------------------------------------------------

			/*
			 These are the window definitions. These windows can be used for either
			 FIR filter design or with an FFT for spectral analysis.
			 For definitions, see this article:  http://en.wikipedia.org/wiki/Window_function

			 This function has 6 inputs
			 Data is the array, of length N, containing the data to to be windowed.
			 This data is either an FIR filter sinc pulse, or the data to be analyzed by an fft.

			 WindowType is an enum defined in the header file.
			 e.g. wtKAISER, wtSINC, wtHANNING, wtHAMMING, wtBLACKMAN, ...

			 Alpha sets the width of the flat top.
			 Windows such as the Tukey and Trapezoid are defined to have a variably wide flat top.
			 As can be seen by its definition, the Tukey is just a Hanning window with a flat top.
			 Alpha can be used to give any of these windows a partial flat top, except the Flattop and Kaiser.
			 Alpha = 0 gives the original window. (i.e. no flat top)
			 To generate a Tukey window, use a Hanning with 0 < Alpha < 1
			 To generate a Bartlett window (triangular), use a Trapezoid window with Alpha = 0.
			 Alpha = 1 generates a rectangular window in all cases. (except the Flattop and Kaiser)


			 Beta is used with the Kaiser, Sinc, and Sine windows only.
			 These three windows are used primarily for FIR filter design. Then
			 Beta controls the filter's transition bandwidth and the sidelobe levels.
			 All other windows ignore Beta.

			 UnityGain controls whether the gain of these windows is set to unity.
			 Only the Flattop window has unity gain by design. The Hanning window, for example, has a gain
			 of 1/2.  UnityGain = true  sets the gain to 1, which preserves the signal's energy level
			 when these windows are used for spectral analysis.

			 Don't use this with FIR filter design however. Since most of the enegy in an FIR sinc pulse
			 is in the middle of the window, the window needs a peak amplitude of one, not unity gain.
			 Setting UnityGain = true will simply cause the resulting FIR filter to have excess gain.

			 If using these windows for FIR filters, start with the Kaiser, Sinc, or Sine windows and
			 adjust Beta for the desired transition BW and sidelobe levels (set Alpha = 0).
			 While the FlatTop is an excellent window for spectral analysis, don't use it for FIR filter design.
			 It has a peak amplitude of ~ 4.7 which causes the resulting FIR filter to have about this much gain.
			 It works poorly for FIR filters even if you adjust its peak amplitude.
			 The Trapezoid also works poorly for FIR filter design.

			 If using these windows with an fft for spectral analysis, start with the Hanning, Gauss, or Flattop.
			 When choosing a window for spectral analysis, you must trade off between resolution and amplitude
			 accuracy. The Hanning has the best resolution while the Flatop has the best amplitude accuracy.
			 The Gauss is midway between these two for both accuracy and resolution. These three were
			 the only windows available in the HP 89410A Vector Signal Analyzer. Which is to say, these three
			 are the probably the best windows for general purpose signal analysis.
			*/

			public static void WindowData(ref double[] Data, int N, TWindowType WindowType, double Alpha, double Beta, bool UnityGain)
			{
				if (WindowType == TWindowType.wtNONE)
					return;

				int j, M, TopWidth;
				double dM;
				double[] WinCoeff;

				if (WindowType == TWindowType.wtKAISER || WindowType == TWindowType.wtFLATTOP)
					Alpha = 0.0;

				if (Alpha < 0.0)
					Alpha = 0.0;
				if (Alpha > 1.0)
					Alpha = 1.0;

				if (Beta < 0.0)
					Beta = 0.0;
				if (Beta > 10.0)
					Beta = 10.0;

				WinCoeff = new double[N + 2];
				if (WinCoeff == null)
				{
					// ShowMessage("Failed to allocate memory in WindowData() ");
					return;
				}

				TopWidth = (int)(Alpha * (double)N);
				if (TopWidth % 2 != 0)
					TopWidth++;
				if (TopWidth > N)
					TopWidth = N;
				M = N - TopWidth;
				dM = M + 1;


				// Calculate the window for N/2 points, then fold the window over (at the bottom).
				// TopWidth points will be set to 1.
				if (WindowType == TWindowType.wtKAISER)
				{
					double Arg;
					for (j = 0; j < M; j++)
					{
						Arg = Beta * Sqrt(1.0 - Pow(((double)(2 * j + 2) - dM) / dM, 2.0));
						WinCoeff[j] = Bessel(Arg) / Bessel(Beta);
					}
				}

				else if (WindowType == TWindowType.wtSINC)  // Lanczos
				{
					for (j = 0; j < M; j++)
						WinCoeff[j] = Sinc((double)(2 * j + 1 - M) / dM * PI);
					for (j = 0; j < M; j++)
						WinCoeff[j] = Pow(WinCoeff[j], Beta);
				}

				else if (WindowType == TWindowType.wtSINE)  // Hanning if Beta = 2
				{
					for (j = 0; j < M / 2; j++)
						WinCoeff[j] = Sin((double)(j + 1) * PI / dM);
					for (j = 0; j < M / 2; j++)
						WinCoeff[j] = Pow(WinCoeff[j], Beta);
				}

				else if (WindowType == TWindowType.wtHANNING)
				{
					for (j = 0; j < M / 2; j++)
						WinCoeff[j] = 0.5 - 0.5 * Cos((double)(j + 1) * 2 * PI / dM);
				}

				else if (WindowType == TWindowType.wtHAMMING)
				{
					for (j = 0; j < M / 2; j++)
						WinCoeff[j] = 0.54 - 0.46 * Cos((double)(j + 1) * 2 * PI / dM);
				}

				else if (WindowType == TWindowType.wtBLACKMAN)
				{
					for (j = 0; j < M / 2; j++)
					{
						WinCoeff[j] = 0.42
						- 0.50 * Cos((double)(j + 1) * 2 * PI / dM)
						+ 0.08 * Cos((double)(j + 1) * 2 * PI * 2.0 / dM);
					}
				}


				// Defined at: http://www.bth.se/fou/forskinfo.nsf/0/130c0940c5e7ffcdc1256f7f0065ac60/$file/ICOTA_2004_ttr_icl_mdh.pdf
				else if (WindowType == TWindowType.wtFLATTOP)
				{
					for (j = 0; j <= M / 2; j++)
					{
						WinCoeff[j] = 1.0
						- 1.93293488969227 * Cos((double)(j + 1) * 2 * PI / dM)
						+ 1.28349769674027 * Cos((double)(j + 1) * 2 * PI * 2.0 / dM)
						- 0.38130801681619 * Cos((double)(j + 1) * 2 * PI * 3.0 / dM)
						+ 0.02929730258511 * Cos((double)(j + 1) * 2 * PI * 4.0 / dM);
					}
				}


				else if (WindowType == TWindowType.wtBLACKMAN_HARRIS)
				{
					for (j = 0; j < M / 2; j++)
					{
						WinCoeff[j] = 0.35875
						- 0.48829 * Cos((double)(j + 1) * 2 * PI / dM)
						+ 0.14128 * Cos((double)(j + 1) * 2 * PI * 2.0 / dM)
						- 0.01168 * Cos((double)(j + 1) * 2 * PI * 3.0 / dM);
					}
				}

				else if (WindowType == TWindowType.wtBLACKMAN_NUTTALL)
				{
					for (j = 0; j < M / 2; j++)
					{
						WinCoeff[j] = 0.3535819
						- 0.4891775 * Cos((double)(j + 1) * 2 * PI / dM)
						+ 0.1365995 * Cos((double)(j + 1) * 2 * PI * 2.0 / dM)
						- 0.0106411 * Cos((double)(j + 1) * 2 * PI * 3.0 / dM);
					}
				}

				else if (WindowType == TWindowType.wtNUTTALL)
				{
					for (j = 0; j < M / 2; j++)
					{
						WinCoeff[j] = 0.355768
						- 0.487396 * Cos((double)(j + 1) * 2 * PI / dM)
						+ 0.144232 * Cos((double)(j + 1) * 2 * PI * 2.0 / dM)
						- 0.012604 * Cos((double)(j + 1) * 2 * PI * 3.0 / dM);
					}
				}

				else if (WindowType == TWindowType.wtKAISER_BESSEL)
				{
					for (j = 0; j <= M / 2; j++)
					{
						WinCoeff[j] = 0.402
						- 0.498 * Cos(2 * PI * (double)(j + 1) / dM)
						+ 0.098 * Cos(2.0 * 2 * PI * (double)(j + 1) / dM)
						+ 0.001 * Cos(3.0 * 2 * PI * (double)(j + 1) / dM);
					}
				}

				else if (WindowType == TWindowType.wtTRAPEZOID) // Rectangle for Alpha = 1  Triangle for Alpha = 0
				{
					int K = M / 2;
					if ((M % 2)==1)
						K++;
					for (j = 0; j < K; j++)
						WinCoeff[j] = (double)(j + 1) / (double)K;
				}


				// This definition is from http://en.wikipedia.org/wiki/Window_function (Gauss Generalized normal window)
				// We set their p = 2, and use Alpha in the numerator, instead of Sigma in the denominator, as most others do.
				// Alpha = 2.718 puts the Gauss window response midway between the Hanning and the Flattop (basically what we want).
				// It also gives the same BW as the Gauss window used in the HP 89410A Vector Signal Analyzer.
				else if (WindowType == TWindowType.wtGAUSS)
				{
					for (j = 0; j < M / 2; j++)
					{
						WinCoeff[j] = ((double)(j + 1) - dM / 2.0) / (dM / 2.0) * 2.7183;
						WinCoeff[j] *= WinCoeff[j];
						WinCoeff[j] = Exp(-WinCoeff[j]);
					}
				}

				else // Error.
				{
					// ShowMessage("Incorrect window type in WindowFFTData");
					return;
				}

				// Fold the coefficients over.
				for (j = 0; j < M / 2; j++)
					WinCoeff[N - j - 1] = WinCoeff[j];

				// This is the flat top if Alpha > 0. Cannot be applied to a Kaiser or Flat Top.
				if (WindowType != TWindowType.wtKAISER && WindowType != TWindowType.wtFLATTOP)
				{
					for (j = M / 2; j < N - M / 2; j++)
						WinCoeff[j] = 1.0;
				}


				// UnityGain = true will set the gain of these windows to 1. Don't use this with FIR filter design.
				if (UnityGain)
				{
					double Sum = 0.0;
					for (j = 0; j < N; j++)
						Sum += WinCoeff[j];
					Sum /= (double)N;
					if (Sum != 0.0)
						for (j = 0; j < N; j++)
							WinCoeff[j] /= Sum;
				}

				// Apply the window to the data.
				for (j = 0; j < N; j++)
					Data[j] *= WinCoeff[j];
				
			}

			//---------------------------------------------------------------------------

			// This gets used with the Kaiser window.
			public static double Bessel(double x)
			{
				double Sum = 0.0, XtoIPower;
				int i, j, Factorial;
				for (i = 1; i < 10; i++)
				{
					XtoIPower = Pow(x / 2.0, (double)i);
					Factorial = 1;
					for (j = 1; j <= i; j++)
						Factorial *= j;
					Sum += Pow(XtoIPower / (double)Factorial, 2.0);
				}
				return (1.0 + Sum);
			}

			//-----------------------------------------------------------------------------
			
			public static double Sinc(double x)
			{
				if (x > -1.0E-5 && x < 1.0E-5)
					return (1.0);
				return (Sin(x) / x);
			}
		}


		/// <summary>
		/// Infinite impulse response filter (old style analog filters)
		/// </summary>
		class IIRFilter
		{
			/// <summary>
			/// The type of filter
			/// </summary>
			public enum FilterType
			{
				None = 0,
				LP,
				HP,
				BP
			}

			/// <summary>
			/// The filter prototype
			/// </summary>
			public enum ProtoType
			{
				None = 0,
				Butterworth,
				Chebyshev,
			}

			const int kHistMask = 31;
			const int kHistSize = 32;

			private int m_order;
			private ProtoType m_protoType;
			private FilterType m_filterType;

			private float m_fp1;
			private float m_fp2;
			private float m_fN;
			private float m_ripple;
			private float m_sampleRate;
			private double[] m_real;
			private double[] m_imag;
			private double[] m_z;
			private double[] m_aCoeff;
			private double[] m_bCoeff;
			private double[] m_inHistory;
			private double[] m_outHistory;
			private int m_histIdx;
			private bool m_invertDenormal;

			public IIRFilter()
			{
			}

			/// <summary>
			/// Returns true if all the filter parameters are valid
			/// </summary>
			public bool FilterValid
			{
				get
				{
					if (m_order < 1 || m_order > 16 ||
					m_protoType == ProtoType.None ||
					m_filterType == FilterType.None ||
					m_sampleRate <= 0.0f ||
					m_fN <= 0.0f)
						return false;

					switch (m_filterType)
					{
						case FilterType.LP:
							if (m_fp2 <= 0.0f)
								return false;
							break;

						case FilterType.BP:
							if (m_fp1 <= 0.0f || m_fp2 <= 0.0f || m_fp1 >= m_fp2)
								return false;
							break;

						case FilterType.HP:
							if (m_fp1 <= 0.0f)
								return false;
							break;
					}

					// For bandpass, the order must be even
					if (m_filterType == FilterType.BP && (m_order & 1) != 0)
						return false;

					return true;
				}
			}

			/// <summary>
			/// Set the filter prototype
			/// </summary>
			public ProtoType Proto
			{
				get { return m_protoType; }

				set
				{
					m_protoType = value;
					Design();
				}
			}

			/// <summary>
			/// Set the filter type
			/// </summary>
			public FilterType Type
			{
				get { return m_filterType; }

				set
				{
					m_filterType = value;
					Design();
				}
			}

			public int Order
			{
				get { return m_order; }

				set
				{
					m_order = Min(16, Max(1, Abs(value)));

					if (m_filterType == FilterType.BP && Odd(m_order))
						m_order++;

					Design();
				}
			}

			public float SampleRate
			{
				get { return m_sampleRate; }

				set
				{
					m_sampleRate = value;
					m_fN = 0.5f * m_sampleRate;
					Design();
				}
			}

			public float FreqLow
			{
				get { return m_fp1; }

				set
				{
					m_fp1 = value;
					Design();
				}
			}

			public float FreqHigh
			{
				get { return m_fp2; }

				set
				{
					m_fp2 = value;
					Design();
				}
			}

			public float Ripple
			{
				get { return m_ripple; }

				set
				{
					m_ripple = value;
					Design();
				}
			}

			/// <summary>
			/// Returns true if n is odd
			/// </summary>
			/// <param name="n"></param>
			/// <returns></returns>
			private bool Odd(int n)
			{
				return (n & 1) == 1;
			}

			/// <summary>
			/// Square
			/// </summary>
			/// <param name="f"></param>
			/// <returns></returns>
			private float Sqr(float value)
			{
				return value * value;
			}

			/// <summary>
			/// Square
			/// </summary>
			/// <param name="f"></param>
			/// <returns></returns>
			private double Sqr(double value)
			{
				return value * value;
			}

			/// <summary>
			/// Determines poles and zeros of IIR filter
			/// based on bilinear transform method
			/// </summary>
			private void LocatePolesAndZeros()
			{
				m_real = new double[m_order + 1];
				m_imag = new double[m_order + 1];
				m_z = new double[m_order + 1];
				double ln10 = Log(10.0);

				// Butterworth, Chebyshev parameters
				int n = m_order;

				if (m_filterType == FilterType.BP)
					n = n / 2;

				int ir = n % 2;
				int n1 = n + ir;
				int n2 = (3 * n + ir) / 2 - 1;
				double f1;

				switch (m_filterType)
				{
					case FilterType.LP:
						f1 = m_fp2;
						break;

					case FilterType.HP:
						f1 = m_fN - m_fp1;
						break;

					case FilterType.BP:
						f1 = m_fp2 - m_fp1;
						break;

					default:
						f1 = 0.0;
						break;
				}

				double tanw1 = Tan(0.5 * PI * f1 / m_fN);
				double tansqw1 = Sqr(tanw1);

				// Real and Imaginary parts of low-pass poles
				double t, a = 1.0, r = 1.0, i = 1.0;

				for (int k = n1; k <= n2; k++)
				{
					t = 0.5 * (2 * k + 1 - ir) * PI / (double)n;

					switch (m_protoType)
					{
						case ProtoType.Butterworth:
							double b3 = 1.0 - 2.0 * tanw1 * Cos(t) + tansqw1;
							r = (1.0 - tansqw1) / b3;
							i = 2.0 * tanw1 * Sin(t) / b3;
							break;

						case ProtoType.Chebyshev:
							double d = 1.0 - Exp(-0.05 * m_ripple * ln10);
							double e = 1.0 / Sqrt(1.0 / Sqr(1.0 - d) - 1.0);
							double x = Pow(Sqrt(e * e + 1.0) + e, 1.0 / (double)n);
							a = 0.5 * (x - 1.0 / x);
							double b = 0.5 * (x + 1.0 / x);
							double c3 = a * tanw1 * Cos(t);
							double c4 = b * tanw1 * Sin(t);
							double c5 = Sqr(1.0 - c3) + Sqr(c4);
							r = 2.0 * (1.0 - c3) / c5 - 1.0;
							i = 2.0 * c4 / c5;
							break;
					}

					int m = 2 * (n2 - k) + 1;
					m_real[m + ir] = r;
					m_imag[m + ir] = Abs(i);
					m_real[m + ir + 1] = r;
					m_imag[m + ir + 1] = -Abs(i);
				}

				if (Odd(n))
				{
					if (m_protoType == ProtoType.Butterworth)
						r = (1.0 - tansqw1) / (1.0 + 2.0 * tanw1 + tansqw1);

					if (m_protoType == ProtoType.Chebyshev)
						r = 2.0 / (1.0 + a * tanw1) - 1.0;

					m_real[1] = r;
					m_imag[1] = 0.0;
				}

				switch (m_filterType)
				{
					case FilterType.LP:
						for (int m = 1; m <= n; m++)
							m_z[m] = -1.0;
						break;

					case FilterType.HP:
						// Low-pass to high-pass transformation
						for (int m = 1; m <= n; m++)
						{
							m_real[m] = -m_real[m];
							m_z[m] = 1.0;
						}
						break;

					case FilterType.BP:
						// Low-pass to bandpass transformation
						for (int m = 1; m <= n; m++)
						{
							m_z[m] = 1.0;
							m_z[m + n] = -1.0;
						}

						double f4 = 0.5 * PI * m_fp1 / m_fN;
						double f5 = 0.5 * PI * m_fp2 / m_fN;
						double aa = Cos(f4 + f5) / Cos(f5 - f4);
						double aR, aI, h1, h2, p1R, p2R, p1I, p2I;

						for (int m1 = 0; m1 <= (m_order - 1) / 2; m1++)
						{
							int m = 1 + 2 * m1;
							aR = m_real[m];
							aI = m_imag[m];

							if (Abs(aI) < 0.0001)
							{
								h1 = 0.5 * aa * (1.0 + aR);
								h2 = Sqr(h1) - aR;
								if (h2 > 0.0)
								{
									p1R = h1 + Sqrt(h2);
									p2R = h1 - Sqrt(h2);
									p1I = 0.0;
									p2I = 0.0;
								}
								else
								{
									p1R = h1;
									p2R = h1;
									p1I = Sqrt(Abs(h2));
									p2I = -p1I;
								}
							}
							else
							{
								double fR = aa * 0.5 * (1.0 + aR);
								double fI = aa * 0.5 * aI;
								double gR = Sqr(fR) - Sqr(fI) - aR;
								double gI = 2 * fR * fI - aI;
								double sR = Sqrt(0.5 * Abs(gR + Sqrt(Sqr(gR) + Sqr(gI))));
								double sI = gI / (2.0 * sR);
								p1R = fR + sR;
								p1I = fI + sI;
								p2R = fR - sR;
								p2I = fI - sI;
							}

							m_real[m] = p1R;
							m_real[m + 1] = p2R;
							m_imag[m] = p1I;
							m_imag[m + 1] = p2I;
						}

						if (Odd(n))
						{
							m_real[2] = m_real[n + 1];
							m_imag[2] = m_imag[n + 1];
						}

						for (int k = n; k >= 1; k--)
						{
							int m = 2 * k - 1;
							m_real[m] = m_real[k];
							m_real[m + 1] = m_real[k];
							m_imag[m] = Abs(m_imag[k]);
							m_imag[m + 1] = -Abs(m_imag[k]);
						}

						break;
				}
			}

			/// <summary>
			/// Calculate all the values
			/// </summary>
			public void Design()
			{
				if (!this.FilterValid)
					return;

				m_aCoeff = new double[m_order + 1];
				m_bCoeff = new double[m_order + 1];
				m_inHistory = new double[kHistSize];
				m_outHistory = new double[kHistSize];

				double[] newA = new double[m_order + 1];
				double[] newB = new double[m_order + 1];

				// Find filter poles and zeros
				LocatePolesAndZeros();

				// Compute filter coefficients from pole/zero values
				m_aCoeff[0] = 1.0;
				m_bCoeff[0] = 1.0;

				for (int i = 1; i <= m_order; i++)
				{
					m_aCoeff[i] = 0.0;
					m_bCoeff[i] = 0.0;
				}

				int k = 0;
				int n = m_order;
				int pairs = n / 2;

				if (Odd(m_order))
				{
					// First subfilter is first order
					m_aCoeff[1] = -m_z[1];
					m_bCoeff[1] = -m_real[1];
					k = 1;
				}

				for (int p = 1; p <= pairs; p++)
				{
					int m = 2 * p - 1 + k;
					double alpha1 = -(m_z[m] + m_z[m + 1]);
					double alpha2 = m_z[m] * m_z[m + 1];
					double beta1 = -2.0 * m_real[m];
					double beta2 = Sqr(m_real[m]) + Sqr(m_imag[m]);

					newA[1] = m_aCoeff[1] + alpha1 * m_aCoeff[0];
					newB[1] = m_bCoeff[1] + beta1 * m_bCoeff[0];

					for (int i = 2; i <= n; i++)
					{
						newA[i] = m_aCoeff[i] + alpha1 * m_aCoeff[i - 1] + alpha2 * m_aCoeff[i - 2];
						newB[i] = m_bCoeff[i] + beta1 * m_bCoeff[i - 1] + beta2 * m_bCoeff[i - 2];
					}

					for (int i = 1; i <= n; i++)
					{
						m_aCoeff[i] = newA[i];
						m_bCoeff[i] = newB[i];
					}
				}

				// Ensure the filter is normalized
				FilterGain(1000);
			}

			/// <summary>
			/// Reset the history buffers
			/// </summary>
			public void Reset()
			{
				if (m_inHistory != null)
					m_inHistory.Clear();

				if (m_outHistory != null)
					m_outHistory.Clear();

				m_histIdx = 0;
			}

			/// <summary>
			/// Reset the filter, and fill the appropriate history buffers with the specified value
			/// </summary>
			/// <param name="startValue"></param>
			public void Reset(double startValue)
			{
				m_histIdx = 0;

				if (m_inHistory == null || m_outHistory == null)
					return;

				m_inHistory.Fill(startValue);

				if (m_inHistory != null)
				{
					switch (m_filterType)
					{
						case FilterType.LP:
							m_outHistory.Fill(startValue);
							break;

						default:
							m_outHistory.Clear();
							break;
					}
				}
			}

			/// <summary>
			/// Apply the filter to the buffer
			/// </summary>
			/// <param name="bufIn"></param>
			public void FilterBuffer(float[] srcBuf, long srcPos, float[] dstBuf, long dstPos, long nLen)
			{
				const double kDenormal = 0.000000000000001;
				double denormal = m_invertDenormal ? -kDenormal : kDenormal;
				m_invertDenormal = !m_invertDenormal;

				for (int sampleIdx = 0; sampleIdx < nLen; sampleIdx++)
				{
					double sum = 0.0f;

					m_inHistory[m_histIdx] = srcBuf[srcPos + sampleIdx] + denormal;

					for (int idx = 0; idx < m_aCoeff.Length; idx++)
						sum += m_aCoeff[idx] * m_inHistory[(m_histIdx - idx) & kHistMask];

					for (int idx = 1; idx < m_bCoeff.Length; idx++)
						sum -= m_bCoeff[idx] * m_outHistory[(m_histIdx - idx) & kHistMask];

					m_outHistory[m_histIdx] = sum;
					m_histIdx = (m_histIdx + 1) & kHistMask;
					dstBuf[dstPos + sampleIdx] = (float)sum;
				}
			}

			public float FilterSample(float inVal)
			{
				double sum = 0.0f;

				m_inHistory[m_histIdx] = inVal;

				for (int idx = 0; idx < m_aCoeff.Length; idx++)
					sum += m_aCoeff[idx] * m_inHistory[(m_histIdx - idx) & kHistMask];

				for (int idx = 1; idx < m_bCoeff.Length; idx++)
					sum -= m_bCoeff[idx] * m_outHistory[(m_histIdx - idx) & kHistMask];

				m_outHistory[m_histIdx] = sum;
				m_histIdx = (m_histIdx + 1) & kHistMask;

				return (float)sum;
			}

			/// <summary>
			/// Get the gain at the specified number of frequency points
			/// </summary>
			/// <param name="freqPoints"></param>
			/// <returns></returns>
			public float[] FilterGain(int freqPoints)
			{
				// Filter gain at uniform frequency intervals
				float[] g = new float[freqPoints];
				double theta, s, c, sac, sas, sbc, sbs;
				float gMax = -100.0f;
				float sc = 10.0f / (float)Log(10.0f);
				double t = PI / (freqPoints - 1);

				for (int i = 0; i < freqPoints; i++)
				{
					theta = i * t;

					if (i == 0)
						theta = PI * 0.0001;

					if (i == freqPoints - 1)
						theta = PI * 0.9999;

					sac = 0.0f;
					sas = 0.0f;
					sbc = 0.0f;
					sbs = 0.0f;

					for (int k = 0; k <= m_order; k++)
					{
						c = Cos(k * theta);
						s = Sin(k * theta);
						sac += c * m_aCoeff[k];
						sas += s * m_aCoeff[k];
						sbc += c * m_bCoeff[k];
						sbs += s * m_bCoeff[k];
					}

					g[i] = sc * (float)Log((Sqr(sac) + Sqr(sas)) / (Sqr(sbc) + Sqr(sbs)));
					gMax = Max(gMax, g[i]);
				}

				// Normalize to 0 dB maximum gain
				for (int i = 0; i < freqPoints; i++)
					g[i] -= gMax;

				// Normalize numerator (a) coefficients
				float normFactor = (float)Pow(10.0, -0.05 * gMax);

				for (int i = 0; i <= m_order; i++)
					m_aCoeff[i] *= normFactor;

				return g;
			}
		}
	}
}
