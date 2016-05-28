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
		public ModulationType(int seed, double carrierFrequency, int samplingRate, int samplesPerSymbol, double signalPower, double[] firCoefficients)
		{
			this._rng = new Random(seed);
			this.CarrierFrequency = carrierFrequency;
			this.SampleRate = samplingRate;
			this.SamplesPerSymbol = samplesPerSymbol;
			this.SignalPower = signalPower;
			this.FIR_Coefficients = firCoefficients;
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

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Demodulation
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// <summary>
		/// Demodulates the waveform based on the internal parameters passed when the class was constructed
		/// </summary>
		/// <param name="waveform">Waveform containing the direct samples</param>
		/// <returns>Estimated bits corresponding to the modulated waveform</returns>
		public abstract byte[] DemodulateWaveform(double[] waveform);

		/// <summary>
		/// Implements FIR filtering on the signal using the coefficients given when the instance was constructed
		/// </summary>
		/// <param name="signal">Signal to filter</param>
		/// <returns>Filtered signal</returns>
		public double[] FIR_Filter(double[] signal)
		{

			int numTaps = this._firCoefficients.Length;
			int numSigPoints = signal.Length;
			int top=0,n,k;
			double[] reg = new double[numTaps];
			double[] filteredSignal = new double[numSigPoints];
			double y;

			for(int j=0; j< numSigPoints; j++)
			{
				reg[top] = signal[j];
				y = 0;
				n = 0;

				for(k= top; k>=0; k--)
				{
					y += this._firCoefficients[n++] * reg[k];
				}
				for (k = numTaps-1; k > top; k--)
				{
					y += this._firCoefficients[n++] * reg[k];
				}
				filteredSignal[j] = y;
				top++;
				if(top>=numTaps)
				{
					top = 0;
				}
			}

			return filteredSignal;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Bit functions
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
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

		/// <summary>
		/// Gets the number of bits that are different in two equal length bit streams
		/// </summary>
		/// <param name="inBits">First bit stream</param>
		/// <param name="outBits">Second bit stream</param>
		/// <returns>Number of bit positions that are different between the two streams</returns>
		public int numberDifferentBits(byte[] inBits, byte[] outBits)
		{
			int numberDifferent = 0;

			//Loop and check each bit
			for(int i=0; i<inBits.Length; i++)
			{
				//If they are different, increment the count
				if (inBits[i] != outBits[i])
					numberDifferent++;
			}

			//Return the number that are different
			return numberDifferent;
		}

		/// <summary>
		/// Returns the number of differences between two bit arrays with two categories
		/// </summary>
		/// <param name="inBits">Input bit array</param>
		/// <param name="outBits">Output bit array</param>
		/// <returns>Number of differences where the input bit is a 1 and the output is a 0, and vice versa.
		/// The first element is the 0->1 case, and the second element is the 1->0 case.</returns>
		public int[] numberDifferentBitsEachType(byte[] inBits, byte[] outBits)
		{
			int[] numberDifferent = new int[2];

			//Loop and check each bit
			for (int i = 0; i < inBits.Length; i++)
			{
				//If the input bit is a zero and the output bit is a one, increment the first counter
				if (inBits[i] == 0 && outBits[i] == 1)
					numberDifferent[0]++;
				//Else if the inpput bit is a one and the output bit is a zero, increment the second counter
				else if (inBits[i] == 1 && outBits[i] == 0)
					numberDifferent[1]++;
			}

			//Return the two counts
			return numberDifferent;
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Noise
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// <summary>
		/// Creates an array of randomly generated white gaussian noise
		/// </summary>
		/// <param name="numberBits">Number of bits long the array should be</param>
		/// <param name="coefficient">Coefficient that determines the power of the noise</param>
		/// <returns>Noise array with specified power and length</returns>
		public double[] AdditiveWhiteGaussianNoise(int numberBits, double coefficient)
		{
			//Create an array to hold the noise sample
			double[] awgn = new double[numberBits * this._samplesPerSymbol];
			double x, y;

			//Loop through each sample
			for(int i=0; i< awgn.Length; i++)
			{
				//Uses the Box-Muller transform to get a Gaussian distributed RV
				x = this._rng.NextDouble();
				y = this._rng.NextDouble();
				awgn[i] = Math.Sqrt(coefficient) * (Math.Sqrt(-2 * Math.Log(x)) * Math.Sin(2 * Math.PI * y));
			}

			//Return the noise array
			return awgn;
		}

		/// <summary>
		/// Adds additive white gaussian noise to a waveform
		/// </summary>
		/// <param name="waveform">Waveform to add the noise to.</param>
		/// <param name="coefficient">Power of the noise</param>
		/// <returns>Actaul power of the noise</returns>
		public double AdditiveWhiteGaussianNoise(ref double[] waveform, double coefficient)
		{
			//Create an array to hold the noise sample
			double x, y, mean = 0;
			double[] awgn = new double[waveform.Length];

			if (coefficient != 0)
			{
				//Loop through each sample
				for (int i = 0; i < waveform.Length; i++)
				{
					//Uses the Box-Muller transform to get a Gaussian distributed RV
					x = this._rng.NextDouble();
					y = this._rng.NextDouble();
					awgn[i] += Math.Sqrt(coefficient) * (Math.Sqrt(-2 * Math.Log(x)) * Math.Sin(2 * Math.PI * y));
					waveform[i] += awgn[i];
					mean += awgn[i];
				}
			} 
			return PowerInWaveform(awgn, mean / waveform.Length);
		}

		/// <summary>
		/// Adds additive white gaussian noise to a wavefor
		/// (Does not return anything)
		/// </summary>
		/// <param name="waveform">Waveform to add the noise to.</param>
		/// <param name="coefficient">Power of the noise</param>
		public void AdditiveWhiteGaussianNoiseNR(ref double[] waveform, double coefficient)
		{
			//Create an array to hold the noise sample
			double x, y;

			if (coefficient != 0)
			{
				//Loop through each sample
				for (int i = 0; i < waveform.Length; i++)
				{
					//Uses the Box-Muller transform to get a Gaussian distributed RV
					x = this._rng.NextDouble();
					y = this._rng.NextDouble();
					waveform[i] += Math.Sqrt(coefficient) * (Math.Sqrt(-2 * Math.Log(x)) * Math.Sin(2 * Math.PI * y));
				}
			}
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//Stats
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		/// <summary>
		/// Get the power in the waveform
		/// </summary>
		/// <param name="waveform">Waveform to find the power of</param>
		/// <returns>The power in the waveform (Watts)</returns>
		public double PowerInWaveform(double[] waveform)
		{
			//The power is just the mean value squared plus the variance
			return Math.Pow(waveform.Mean(),2) + waveform.Variance();
		}

		/// <summary>
		/// Get the power in the waveform
		/// </summary>
		/// <param name="waveform">Waveform to find the power of</param>
		/// <param name="mean">Mean value of the waveform</param>
		/// <returns>The power in the waveform (Watts)</returns>
		public double PowerInWaveform(double[] waveform, double mean)
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
		public double Calc_dB(double signalPower, double noisePower)
		{
			return 10 * Math.Log10(signalPower / noisePower);
		}

		/// <summary>
		/// Runs a simulation for the given parameters and outputs a bit error rate for that simulation
		/// </summary>
		/// <param name="numberToGetWrongEventually">Number of bits to eventually mis-estimate</param>
		/// <param name="numberBitsPerIteration">Number of bits simulated in each iteration of the loop</param>
		/// <returns>Overall bit error rate of the simulation</returns>
		public double RunSimulationOneNoisePower(int numberToGetWrongEventually, int numberBitsPerIteration, double noisePower)
		{
			byte[] inbits;
			double[] waveform;
			byte[] outbits;

			int totalNumWrong = 0, totalBitsSimulated = 0;
			while (totalNumWrong < numberToGetWrongEventually)
			{
				inbits = GenerateRandomBits(numberBitsPerIteration);
				waveform = ModulateBits(inbits);
				AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);
				outbits = DemodulateWaveform(waveform);
				totalNumWrong += numberDifferentBits(inbits, outbits);
				totalBitsSimulated += numberBitsPerIteration;
			}
			return totalNumWrong / (double)totalBitsSimulated;
		}
	}
}
