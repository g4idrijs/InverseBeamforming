using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Collects many different waveform modulations into 1 class
	/// </summary>
	public partial class Modulations
	{
		/// <summary>
		/// Provides methods to perform BPSK modulation on a bit stream
		/// </summary>
		public class BPSK_Modulation : ModulationType
		{
			/// <summary>
			/// Create an instance with the given parameters
			/// </summary>
			/// <param name="carrierFrequency">Carrier frequency of any resulting modulated waveforms</param>
			/// <param name="samplingRate">Sampling rates to generate the waveforms at</param>
			/// <param name="seed">Random seed to start a random number generator at</param>
			/// <param name="samplesPerSymbol">Number of samples in one symbol</param>
			/// <param name="signalPower">Average power of the resulting modulated waveforms</param>
			/// <param name="firCoefficients">Coefficients of an FIR filter to filter with before bit estimation.</param>
			/// <param name="numberSymbolsPerWaveform">Number of symbols created in each waveform</param>
			public BPSK_Modulation(double carrierFrequency, int samplingRate, int seed, int samplesPerSymbol, double signalPower, double[] firCoefficients, int numberSymbolsPerWaveform)
				: base(seed, carrierFrequency, samplingRate, samplesPerSymbol, signalPower, firCoefficients, numberSymbolsPerWaveform, 2)
			{
				var pi2 = 2 * Math.PI * carrierFrequency / samplingRate;

				this._bitsToCommunicationsSymbols = new byte[] { 0, 1 };

				for (int i = 0; i < samplesPerSymbol; i++)
				{
					this._reference[1, i] = Math.Cos(pi2 * i);
					this._reference[0, i] = -_reference[1, i];
				}
			}

			/// <summary>
			/// Provides BPSK modulation for a given set of bits
			/// </summary>
			/// <param name="bitsToModulate">The bit array to modulate</param>
			/// <returns>A modulated waveform generated from the given bit array</returns>
			public override double[] ModulateBits(byte[] bitsToModulate)
			{
				//Generate the (full length) time vector
				var vecLength = bitsToModulate.Length * this._samplesPerSymbol;
				var time = getTimeArray(vecLength);

				//Generate the (bit length) phase vector
				double phase;
				int index;
				var waveform = new double[vecLength];

				for (int i = 0; i < bitsToModulate.Length; i++)
				{
					phase = Math.PI - Math.PI * bitsToModulate[i];

					for (int k = 0; k < this._samplesPerSymbol; k++)
					{
						index = i * this._samplesPerSymbol + k;
						waveform[index] = Math.Sqrt(2 * this._signalPower) * Math.Cos(2 * Math.PI * this._carrierFrequency * time[index] - phase);
					}
				}

				return waveform;
			}
		}
	}
}
