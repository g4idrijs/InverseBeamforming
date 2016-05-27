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
		public BPSK_Modulation(double carrierFrequency, int samplingRate, int seed, int samplesPerSymbol, double signalPower)
			: base(seed, carrierFrequency, samplingRate, samplesPerSymbol, signalPower){ }

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
				phase = Math.PI*(double)bitsToModulate[i];

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
