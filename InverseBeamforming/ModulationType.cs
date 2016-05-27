using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace InverseBeamforming
{
	/// <summary>
	/// Defines an interface that allows multiple modulation types to be created
	/// </summary>
	public abstract class ModulationType
	{
		/// <summary>
		/// String containing the type of modulation used
		/// </summary>
		public string Modulation { get; }

		/// <summary>
		/// Dobule containing the carrier frequency of output waveform (in Hz)
		/// </summary>
		public double CarrierFrequency
		{
			get
			{
				return this._carrierFrequency;
			}
			set
			{
				if (value < 0)
					throw new ArgumentOutOfRangeException("CarrierFrequency", "The Carrier Frequency must be positive.");
				else
					this._carrierFrequency = value;
			}
		}
		protected double _carrierFrequency;

		/// <summary>
		/// Int containing the sampling rate of the output waveform (in Samples/second)
		/// </summary>
		public int SampleRate
		{
			get { return this._sampleRate; }
			set
			{
				if (value <= 0)
					throw new ArgumentOutOfRangeException("SampleRate", "The sample rate must be greater than 0.");

				this._sampleRate = value;
			}
		}
		protected int _sampleRate;

		/// <summary>
		/// Random number generator to generate random bits to modulate
		/// </summary>
		public Random rng
		{
			set
			{
				this._rng = value;
			}
		}
		protected Random _rng;

		/// <summary>
		/// The number of samples in one symbol period (Samples/symbol)
		/// </summary>
		public int SamplesPerSymbol
		{
			get { return this._samplesPerSymbol; }
			set
			{
				if (value <= 0)
					throw new ArgumentOutOfRangeException("SamplesPerSymbol", "Samples per symbol cannot be less than or equal to zero.");

				this._samplesPerSymbol = value;
			}
		}
		protected int _samplesPerSymbol;

		/// <summary>
		/// Coefficients of an FIR filter
		/// </summary>
		public double[] FIR_Coefficients
		{
			get { return this._firCoefficients; }
			set
			{
				this._firCoefficients = value;
			}
		}
		protected double[] _firCoefficients;

		/// <summary>
		/// Power to be contained in the signal waveform after modulation
		/// </summary>
		public double SignalPower
		{
			get { return this._signalPower; }
			set
			{
				if (value <= 0)
					throw new ArgumentOutOfRangeException("SignalPower", "Signal power cannot be 0 or less.");

				this._signalPower = value;
			}
		}
		protected double _signalPower;

		/// <summary>
		/// Creates a new instance of the ModulationType class with seed as the seed for the RNG
		/// </summary>
		/// <param name="seed">Seed to initialize the RNG</param>
		public ModulationType(int seed, double carrierFrequency, int samplingRate, int samplesPerSymbol, double signalPower)
		{
			this._rng = new Random(seed);
			this.CarrierFrequency = carrierFrequency;
			this.SampleRate = samplingRate;
			this.SamplesPerSymbol = samplesPerSymbol;
			this.SignalPower = signalPower;
		}

		/// <summary>
		/// (Needs to be implemented in inheriting class)
		/// Modulates numberBits and produces a waveform at the carrier frequency and with the given sampling rate
		/// </summary>
		/// <param name="numberBits">Number of bits to modulate</param>
		/// <returns>Waveform containg the modulated bits</returns>
		public double[] ModulateBits(int numberBits)
		{
			var bits = this.GenerateRandomBits(numberBits);

			return this.ModulateBits(bits);
		}

		/// <summary>
		/// (Needs to be implemented in inheriting class)
		/// Modulates the given bits at the sampling rate and carrier frequency given in the constructor
		/// </summary>
		/// <param name="bitsToModulate">The collection of bits to modulate.
		/// (Does not produce a random set of bits to modulate)</param>
		/// <returns>Waveform containing the modulated bits</returns>
		public abstract double[] ModulateBits(byte[] bitsToModulate);

		/// <summary>
		/// Generates numberBits of random bits.
		/// </summary>
		/// <param name="numberBits">The number of bits to generate</param>
		/// <returns>Byte array containg the generated bits</returns>
		public byte[] GenerateRandomBits(int numberBits)
		{
			//Create the array
			var bits = new byte[numberBits];
			
			//Fill each bit with a random value
			for(int i=0; i<numberBits; i++)
			{
				bits[i] = (byte)this._rng.Next(0, 2);
			}
			//Return the array
			return bits;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Demodulation
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// <summary>
		/// Demodulates the waveform based on the internal parameters passed when the class was constructed
		/// </summary>
		/// <param name="waveform">Waveform containing the direct samples</param>
		/// <returns>Estimated bits corresponding to the modulated waveform</returns>
		public abstract byte[] DemodulateWaveform(double[] waveform);

		public double[] FIR_Filter(ref double[] signal)
		{

			int numTaps = this._firCoefficients.Length;
			int numSigPoints = signal.Length;
			int top=0,n;
			double[] reg = new double[numTaps];
			double[] filteredSignal = new double[numTaps];
			double y;

			for(int j=0; j< numSigPoints; j++)
			{
				reg[top] = signal[j];
				y = 0;
				n = 0;

				for(int k= top; k>=0; k--)
				{
					y += this._firCoefficients[n++] * reg[k];
				}

			}

			return filteredSignal;
		}
	}
}
