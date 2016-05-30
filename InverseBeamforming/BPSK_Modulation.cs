using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Provides methods to perform BPSK modulation on a bit stream
	/// </summary>
	public class BPSK_Modulation : ModulationType
	{
		private double[] _reference1;
		private double[] _reference0;

		/// <summary>
		/// Create an instance with the given parameters
		/// </summary>
		/// <param name="carrierFrequency">Carrier frequency of any resulting modulated waveforms</param>
		/// <param name="samplingRate">Sampling rates to generate the waveforms at</param>
		/// <param name="seed">Random seed to start a random number generator at</param>
		/// <param name="samplesPerSymbol">Number of samples in one symbol</param>
		/// <param name="signalPower">Average power of the resulting modulated waveforms</param>
		/// <param name="firCoefficients">Coefficients of an FIR filter to filter with before bit estimation.</param>
		public BPSK_Modulation(double carrierFrequency, int samplingRate, int seed, int samplesPerSymbol, double signalPower, double[] firCoefficients, int numberSymbolsPerWaveform)
			: base(seed, carrierFrequency, samplingRate, samplesPerSymbol, signalPower, firCoefficients, numberSymbolsPerWaveform)
		{
			var pi2 = 2 * Math.PI *carrierFrequency / samplingRate;
			this._reference0 = new double[samplesPerSymbol];
			this._reference1 = new double[samplesPerSymbol];

			for (int i=0; i<samplesPerSymbol; i++)
			{
				this._reference1[i] = Math.Cos(pi2 * i);
				this._reference0[i] = -_reference1[i];
			}
		}

		/// <summary>
		/// Demodulates a waveform by estimating the bits
		/// </summary>
		/// <param name="waveform">Waveform to estimate bits from</param>
		/// <returns>Estimated bits from the waveform</returns>
		/// <remarks>This method assumes perfect synchronization with the received waveform, and perfect symbol boundaries.</remarks>
		public override byte[] DemodulateWaveform(double[] waveform)
		{
			//Create array to hold the demodulated bits
			byte[] bits = new byte[waveform.Length / this._samplesPerSymbol];

			//Create an array to hold all of the z-Scores
			double z1, z2;
			
			//Loop through each symbol
			for(int i=0; i< bits.Length; i++)
			{
				//Reset the z scores
				z1 = 0;
				z2 = 0;
				//Loop through each sample i the signal
				for(int k=0; k<_samplesPerSymbol; k++)
				{
					//Correlate the symbol to the reference signals
					z1+=waveform[i * _samplesPerSymbol + k] * _reference1[k];
					z2+=waveform[i * _samplesPerSymbol + k] * _reference0[k];
				}
				//If the z score from the first reference is higher than the second reference, then call the bit a 1
				if(z1>=z2)
				{
					bits[i] = 1;
				}
				else	//Else call it a 0
				{
					bits[i] = 0;
				}
			}

			//Return the bit array
			return bits;
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
			var time = new double[vecLength];
			for(int i=0; i< vecLength; i++)
			{
				time[i] = (double)i / this._sampleRate;
			}

			//Generate the (bit length) phase vector
			double phase;
			int index;
			var waveform = new double[vecLength];

			for (int i = 0; i < bitsToModulate.Length; i++)
			{
				phase = Math.PI - Math.PI * bitsToModulate[i];

				for(int k=0; k< this._samplesPerSymbol; k++)
				{
					index = i * this._samplesPerSymbol + k;
					waveform[index]= Math.Sqrt(2*this._signalPower) * Math.Cos(2 * Math.PI * this._carrierFrequency * time[index] - phase);
				}
			}

			return waveform;
		}
	}
}
