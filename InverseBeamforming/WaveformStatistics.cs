using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	public static class WaveformStatistics
	{

		/// <summary>
		/// Get the power in the waveform
		/// </summary>
		/// <param name="waveform">Waveform to find the power of</param>
		/// <returns>The power in the waveform (Watts)</returns>
		public static double PowerInWaveform(double[] waveform)
		{
			//The power is just the mean value squared plus the variance
			return Math.Pow(waveform.Mean(), 2) + waveform.Variance();
		}

		/// <summary>
		/// Get the power in the waveform
		/// </summary>
		/// <param name="waveform">Waveform to find the power of</param>
		/// <param name="mean">Mean value of the waveform</param>
		/// <returns>The power in the waveform (Watts)</returns>
		public static double PowerInWaveform(double[] waveform, double mean)
		{
			//The power is just the mean value squared plus the variance
			return Math.Pow(mean, 2) + waveform.Variance(mean);
		}

		/// <summary>
		/// Convert the ratio to decibels
		/// </summary>
		/// <param name="signalPower">Power in the signal</param>
		/// <param name="noisePower">Power in the noise</param>
		/// <returns>Signal to noise ratio in dB</returns>
		public static double Calc_dB(double signalPower, double noisePower)
		{
			return 10 * Math.Log10(signalPower / noisePower);
		}
	}
}
