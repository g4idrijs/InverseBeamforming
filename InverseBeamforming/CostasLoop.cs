using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public class CostasLoop
	{

		private int _order;
		private double _loopBW;
		private double _damping;
		private double _alpha;
		private double _beta;
		private delegate double PhaseDetector(Complex sample);
		private PhaseDetector _phaseDetector;

		private double _phase;
		private double _freq;
		private double _max_freq;
		private double _min_freq;
		private double _noise = 1;
		private double _error = 0;
		private static readonly double M_TWOPI = 2 * Math.PI;

		public CostasLoop(double loopBW, int order, bool useSNR)
		{
			_order = order;
			_loopBW = loopBW;

			// Set up the phase detector to use based on the constellation order
			switch (_order)
			{
				case 2:
					if (useSNR)
						_phaseDetector = phase_detector_snr_2;
					else
						_phaseDetector = phase_detector_2;
					break;

				case 4:
					if (useSNR)
						_phaseDetector = phase_detector_snr_4;
					else
						_phaseDetector = phase_detector_4;
					break;

				case 8:
					if (useSNR)
						_phaseDetector = phase_detector_snr_8;
					else
						_phaseDetector = phase_detector_8;
					break;

				default:
					_phaseDetector = phase_detector_2;
					break;
			}
		}

		/// <summary>
		/// Run the costas loop
		/// </summary>
		/// <param name="input">input waveform</param>
		/// <returns>Waveform reduced to baseband</returns>
		public Complex[] runCostasLoop(Complex[] input)
		{
			Complex[] output = new Complex[input.Length];

			Complex nco_out;

			for (int i = 0; i < input.Length; i++)
			{
				nco_out = new Complex(Math.Cos(_phase), Math.Sin(_phase));
				output[i] = input[i] * nco_out;

				_error = _phaseDetector(output[i]);
				_error = branchlessClip(_error, 1);

				advanceLoop(_error);
				phaseWrap();
				frequencyLimit();
			}
			return output;
		}

		/// <summary>
		/// Wraps the phase to be between -2pi and 2pi
		/// </summary>
		private void phaseWrap()
		{
			while (_phase > M_TWOPI)
			{
				_phase -= M_TWOPI;
			}
			while (_phase < -M_TWOPI)
			{
				_phase += M_TWOPI;
			}
		}

		/// <summary>
		/// Advance the control loop based on the current gain settings and the inputted error signal. 
		/// </summary>
		/// <param name="error">Error signal from the feedback loop</param>
		private void advanceLoop(double error)
		{
			_freq = _freq + _beta * error;
			_phase = _phase + _freq + _alpha * error;
		}

		/// <summary>
		/// Limits the frequency to the range given
		/// </summary>
		void frequencyLimit()
		{
			if (_freq > _max_freq)
				_freq = _max_freq;
			else if (_freq < _min_freq)
				_freq = _min_freq;
		}

		/// <summary>
		/// Update the system gains from the loop bandwidth and damping factor. 
		/// </summary>
		private void updateGains()
		{
			double denom = (1 + 2 * _damping * _loopBW + _loopBW * _loopBW);
			_alpha = (4 * _damping * _loopBW) / denom;
			_beta = (4 * _loopBW * _loopBW) / denom;
		}

		/// <summary>
		/// This bounds x by +/- clip without a branch 
		/// </summary>
		/// <param name="x">variable to clip</param>
		/// <param name="clip">Limit on how far away from 0 the output can be</param>
		/// <returns>Returns x clipped to clip without any branches</returns>
		static double branchlessClip(double x, double clip)
		{
			double x1 = Math.Abs(x + clip);
			double x2 = Math.Abs(x - clip);
			x1 -= x2;
			return 0.5 * x1;
		}







		double phase_detector_8(Complex sample)
		{
			/* This technique splits the 8PSK ellation into 2 squashed
		   QPSK ellations, one when I is larger than Q and one
		   where Q is larger than I. The error is then calculated
		   proportionally to these squashed ellations by the 
		   K = Math.Sqrt(2)-1.

		   The signal magnitude must be > 1 or K will incorrectly bias
		   the error value.

		   Ref: Z. Huang, Z. Yi, M. Zhang, K. Wang, "8PSK demodulation for
		   new generation DVB-S2", IEEE Proc. Int. Conf. Communications,
		   Circuits and Systems, Vol. 2, pp. 1447 - 1450, 2004.
			*/

			double K = (Math.Sqrt(2.0) - 1);
			if (Math.Abs(sample.Real) >= Math.Abs(sample.Imaginary))
			{
				return ((sample.Real > 0 ? 1.0 : -1.0) * sample.Imaginary -
					(sample.Imaginary > 0 ? 1.0 : -1.0) * sample.Real * K);
			}
			else
			{
				return ((sample.Real > 0 ? 1.0 : -1.0) * sample.Imaginary * K -
					(sample.Imaginary > 0 ? 1.0 : -1.0) * sample.Real);
			}
		}

		double
		phase_detector_4(Complex sample)
		{
			return ((sample.Real > 0 ? 1.0 : -1.0) * sample.Imaginary -
				(sample.Imaginary > 0 ? 1.0 : -1.0) * sample.Real);
		}

		double
		phase_detector_2(Complex sample)
		{
			return (sample.Real * sample.Imaginary);
		}

		double
		phase_detector_snr_8(Complex sample)
		{
			double K = (Math.Sqrt(2.0) - 1);
			double snr = sample.Magnitude * sample.Magnitude / _noise;
			if (Math.Abs(sample.Real) >= Math.Abs(sample.Imaginary))
			{
				return ((tanhf_lut(snr * sample.Real) * sample.Imaginary) -
					  (tanhf_lut(snr * sample.Imaginary) * sample.Real * K));
			}
			else
			{
				return ((tanhf_lut(snr * sample.Real) * sample.Imaginary * K) -
					  (tanhf_lut(snr * sample.Imaginary) * sample.Real));
			}
		}

		double
		phase_detector_snr_4(Complex sample)
		{
			double snr = sample.Magnitude * sample.Magnitude / _noise;
			return ((tanhf_lut(snr * sample.Real) * sample.Imaginary) -
					(tanhf_lut(snr * sample.Imaginary) * sample.Real));
		}

		double
		phase_detector_snr_2(Complex sample)
		{
			double snr = sample.Magnitude * sample.Magnitude / _noise;
			return tanhf_lut(snr * sample.Real) * sample.Imaginary;
		}


		// This is a table of tanh(x) for x in [-2, 2] used in tanh_lut. 
		private static double[] tanh_lut_table = { -0.96402758, -0.96290241, -0.96174273, -0.96054753, -0.95931576,
													-0.95804636, -0.95673822, -0.95539023, -0.95400122, -0.95257001,
													-0.95109539, -0.9495761 , -0.94801087, -0.94639839, -0.94473732,
													-0.94302627, -0.94126385, -0.93944862, -0.93757908, -0.93565374,
													-0.93367104, -0.93162941, -0.92952723, -0.92736284, -0.92513456,
													-0.92284066, -0.92047938, -0.91804891, -0.91554743, -0.91297305,
													-0.91032388, -0.90759795, -0.9047933 , -0.90190789, -0.89893968,
													-0.89588656, -0.89274642, -0.88951709, -0.88619637, -0.88278203,
													-0.87927182, -0.87566342, -0.87195453, -0.86814278, -0.86422579,
													-0.86020115, -0.85606642, -0.85181914, -0.84745683, -0.84297699,
													-0.83837709, -0.83365461, -0.82880699, -0.82383167, -0.81872609,
													-0.81348767, -0.80811385, -0.80260204, -0.7969497 , -0.79115425,
													-0.78521317, -0.77912392, -0.772884  , -0.76649093, -0.75994227,
													-0.75323562, -0.74636859, -0.73933889, -0.73214422, -0.7247824 ,
													-0.71725127, -0.70954876, -0.70167287, -0.6936217 , -0.68539341,
													-0.67698629, -0.66839871, -0.65962916, -0.65067625, -0.64153871,
													-0.6322154 , -0.62270534, -0.61300768, -0.60312171, -0.59304692,
													-0.58278295, -0.57232959, -0.56168685, -0.55085493, -0.53983419,
													-0.52862523, -0.51722883, -0.50564601, -0.49387799, -0.48192623,
													-0.46979241, -0.45747844, -0.44498647, -0.4323189 , -0.41947836,
													-0.40646773, -0.39329014, -0.37994896, -0.36644782, -0.35279057,
													-0.33898135, -0.32502449, -0.31092459, -0.2966865 , -0.28231527,
													-0.26781621, -0.25319481, -0.23845682, -0.22360817, -0.208655  ,
													-0.19360362, -0.17846056, -0.16323249, -0.14792623, -0.13254879,
													-0.11710727, -0.10160892, -0.08606109, -0.07047123, -0.05484686,
													-0.0391956 , -0.02352507, -0.00784298,  0.00784298,  0.02352507,
													0.0391956 ,  0.05484686,  0.07047123,  0.08606109,  0.10160892,
													0.11710727,  0.13254879,  0.14792623,  0.16323249,  0.17846056,
													0.19360362,  0.208655  ,  0.22360817,  0.23845682,  0.25319481,
													0.26781621,  0.28231527,  0.2966865 ,  0.31092459,  0.32502449,
													0.33898135,  0.35279057,  0.36644782,  0.37994896,  0.39329014,
													0.40646773,  0.41947836,  0.4323189 ,  0.44498647,  0.45747844,
													0.46979241,  0.48192623,  0.49387799,  0.50564601,  0.51722883,
													0.52862523,  0.53983419,  0.55085493,  0.56168685,  0.57232959,
													0.58278295,  0.59304692,  0.60312171,  0.61300768,  0.62270534,
													0.6322154 ,  0.64153871,  0.65067625,  0.65962916,  0.66839871,
													0.67698629,  0.68539341,  0.6936217 ,  0.70167287,  0.70954876,
													0.71725127,  0.7247824 ,  0.73214422,  0.73933889,  0.74636859,
													0.75323562,  0.75994227,  0.76649093,  0.772884  ,  0.77912392,
													0.78521317,  0.79115425,  0.7969497 ,  0.80260204,  0.80811385,
													0.81348767,  0.81872609,  0.82383167,  0.82880699,  0.83365461,
													0.83837709,  0.84297699,  0.84745683,  0.85181914,  0.85606642,
													0.86020115,  0.86422579,  0.86814278,  0.87195453,  0.87566342,
													0.87927182,  0.88278203,  0.88619637,  0.88951709,  0.89274642,
													0.89588656,  0.89893968,  0.90190789,  0.9047933 ,  0.90759795,
													0.91032388,  0.91297305,  0.91554743,  0.91804891,  0.92047938,
													0.92284066,  0.92513456,  0.92736284,  0.92952723,  0.93162941,
													0.93367104,  0.93565374,  0.93757908,  0.93944862,  0.94126385,
													0.94302627,  0.94473732,  0.94639839,  0.94801087,  0.9495761 ,
													0.95109539,  0.95257001,  0.95400122,  0.95539023,  0.95673822,
													0.95804636,  0.95931576,  0.96054753,  0.96174273,  0.96290241,
													0.96402758 };


		/// <summary>
		/// A look-up table (LUT) tanh calcuation. This function returns an 
		/// estimate to tanh(x) based on a 256-point LUT between -2 and
		/// 2. If x is less than -2, it returns -1; if greater than 2, it returns 1. 
		/// 
		/// This LUT form of the tanh is "hidden" in this code because it
		/// is likely too coarse an estimate for any real uses of a
		/// tanh. It is useful, however, in certain control loop
		/// applications where the input is expected to be within these
		/// bounds and the noise will be greater than the quanitzation of 
		/// this small LUT. For more accurate forms of tanh, see
		/// volk_32f_tanh_32f. 
		/// </summary>
		/// <param name="x">Value to find the hyperbolic tangent of</param>
		/// <returns>Coarse tanh of x</returns>
		static double tanhf_lut(double x)
		{
			if (x > 2)
				return 1;
			else if (x <= -2)
				return -1;
			else
			{
				int index = (int)(128 + 64 * x);
				return tanh_lut_table[index];
			}
		}
	}
}
