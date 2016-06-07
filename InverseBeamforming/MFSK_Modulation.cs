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
		public class MFSK_Modulation : ModulationType
		{
			/// <summary>
			/// Set of frequencies used by the communications symbols
			/// </summary>
			public double[] Frequencies
			{
				get { return this._frequencies; }
				set
				{
					this._frequencies = value;
				}
			}
			private double[] _frequencies;

			/// <summary>
			/// Create a new instance of a MFSK modulator
			/// </summary>
			/// <param name="seed">Random seed to start a random number generator at</param>
			/// <param name="carrierFrequency">Carrier frequency of any resulting modulated waveforms</param>
			/// <param name="samplingRate">Sampling rates to generate the waveforms at</param>
			/// <param name="samplesPerSymbol">Number of samples in one symbol</param>
			/// <param name="signalPower">Average power of the resulting modulated waveforms</param>
			/// <param name="firCoefficients">Coefficients of an FIR filter to filter with before bit estimation.</param>
			/// <param name="numberSymbolsPerWaveform">Number of symbols created in each waveform</param>
			/// <param name="M">Number of unique communications symbols</param>
			/// <param name="frequencies">Frequencies used for each communications symbol</param>
			public MFSK_Modulation(int seed, double carrierFrequency, int samplingRate, int samplesPerSymbol, double signalPower, double[] firCoefficients, int numberSymbolsPerWaveform, int M, double[] frequencies, byte[] bitsToCommunicationsSymbols)
				: base(seed, carrierFrequency, samplingRate, samplesPerSymbol, signalPower, firCoefficients, numberSymbolsPerWaveform, M)
			{
				//Set up the frequencies to be used to generate the reference waveforms
				this.Frequencies = frequencies;

				//Get a scale for each sample
				double pi2 = Math.PI * 2 / samplingRate;

				//Loop through each symbol
				for (int i = 0; i < M; i++)
				{
					//Loop through the samples in the symbol
					for (int k = 0; k < samplesPerSymbol; k++)
					{
						this._reference[i, k] = Math.Sin(pi2 * frequencies[i] * k);
					}
				}

				//Set the bit to communication symbol matching
				this._bitsToCommunicationsSymbols = bitsToCommunicationsSymbols;
			}

			/// <summary>
			/// Provides MFSK modulation for a given set of bits
			/// </summary>
			/// <param name="bitsToModulate">The bit array to modulate</param>
			/// <returns>A modulated waveform generated from the given bit array</returns>
			public override double[] ModulateBits(byte[] bitsToModulate)
			{
				if (!checkParams())
					throw new ArgumentException("One of the parameters is not the correct length. (Length of the frequency array, number of rows in the reference array, number of rows in the bits->communications symbols array, or bits per symbol>log_2(M).");

				//Generate the (full length) time vector
				var vecLength = bitsToModulate.Length * SamplesPerSymbol;
				var time = getTimeArray(SamplesPerSymbol);
				double[] waveform = new double[vecLength];
				double freq = this._frequencies[0];

				//Loop through each group of bits
				for (int i = 0; i < bitsToModulate.Length; i++)
				{
					//Determine which symbol to use for the current waveform
					freq = _frequencies[(int)bitsToModulate[i]];

					//Get the samaple values for the current symbol
					for (int k = 0; k < this._samplesPerSymbol; k++)
					{
						waveform[i * this._samplesPerSymbol + k] = Math.Sqrt(2 * this._signalPower) * Math.Sin(2 * Math.PI * freq * time[k]);
					}
				}

				return waveform;
			}

			/// <summary>
			/// Check the parameters to make sure that all are the right length
			/// </summary>
			/// <returns>True if they are okay, false if there is a problem.</returns>
			private bool checkParams()
			{
				if (this._frequencies.Length != this.M ||
					this._reference.GetLength(0) != this.M ||
					this._bitsToCommunicationsSymbols.Length != this.M)
					return false;
				else
					return true;
			}
		}
	}
}
