using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

/// <summary>
/// Contains all the components necesary to perform simulations regaurding digital communications systems
/// and inverse beamforming
/// </summary>
namespace InverseBeamforming
{
	/// <summary>
	/// Collects many different waveform modulations into 1 class
	/// </summary>
	public partial class Modulations
	{
		/// <summary>
		/// Provide methods to complete MPSK modulation
		/// </summary>
		public class MPSK_Modulation : ModulationType
		{
			/// <summary>
			/// Phases used for the modulation of the bits
			/// </summary>
			public double[] Phases
			{
				get { return this._phases; }
				set
				{
					foreach (var phase in value)
					{
						if (phase > 2 * Math.PI)
							throw new ArgumentOutOfRangeException("The phases must be in radians between 0 and 2pi.");
					}
					this._phases = value;
				}
			}
			private double[] _phases;

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
			/// <param name="M">Number of unique communications symbols</param>
			public MPSK_Modulation(double carrierFrequency, int samplingRate, int seed, int samplesPerSymbol, double signalPower, double[] firCoefficients, int numberSymbolsPerWaveform, int M, double[] phases, byte[] bitsToCommunicationsSymbols)
				: base(carrierFrequency, samplingRate, samplesPerSymbol, signalPower, firCoefficients, numberSymbolsPerWaveform, M)
			{
				this.Phases = phases;
				var pi2 = 2 * Math.PI * carrierFrequency / samplingRate;

				signalPower = Math.Sqrt(2 * signalPower);

				//Loop through each symbol
				for (int i = 0; i < M; i++)
				{
					//Loop through the samples in the symbol
					for (int k = 0; k < samplesPerSymbol; k++)
					{
						this._reference[i][k] = signalPower * Math.Cos(pi2 * k + phases[i]);
					}
				}

				//Set the bit to communication symbol matching
				this._bitsToCommunicationsSymbols = bitsToCommunicationsSymbols;
			}

			/// <summary>
			/// Copy Constructor
			/// </summary>
			/// <param name="old">Old class to copy</param>
			public MPSK_Modulation(ModulationType old) : base(old)
			{}

			/// <summary>
			/// Copy Constructor
			/// </summary>
			/// <param name="old">Old class to copy</param>
			public MPSK_Modulation(MPSK_Modulation old) : base(old)
			{
				this._phases = old._phases;
				this._bitsToCommunicationsSymbols = old._bitsToCommunicationsSymbols;
				this._reference = (double[][])old._reference.Clone();
			}

			/// <summary>
			/// Modulates bits using the parameters already set
			/// </summary>
			/// <param name="bitsToModulate">Array of bits to modulate</param>
			/// <returns>Waveform representing the modulated bits</returns>
			public override double[] ModulateBits(byte[] bitsToModulate)
			{
				//Generate the (full length) time vector
				var vecLength = bitsToModulate.Length * this._samplesPerSymbol;
				var time = getTimeArray(this._samplesPerSymbol);
				double phase = this._phases[0];

				//Generate the (bit length) phase vector
				var waveform = new double[vecLength];

				for (int i = 0; i < bitsToModulate.Length; i++)
				{
					phase = this._phases[(int)bitsToModulate[i]];

					for (int k = 0; k < this._samplesPerSymbol; k++)
					{
						waveform[i * this._samplesPerSymbol + k] = Math.Sqrt(2 * this._signalPower) * Math.Cos(2 * Math.PI * this._carrierFrequency * time[k] + phase);
					}
				}

				return waveform;
			}
		}
	}
}
