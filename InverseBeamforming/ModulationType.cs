using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reactive.Linq;
using System.Security.Cryptography;
using System.Text;
using System.Threading;
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
		/// Defines an interface that allows multiple modulation types to be created
		/// </summary>
		public abstract class ModulationType
		{
			#region Properties

			/// <summary>
			/// Describes the types of filters available
			/// </summary>
			public enum EFilterTypeToUse { IIR, FIR};

			/// <summary>
			/// Describes the filter choices
			/// </summary>
			public enum EFilterToUse { RF, DS };
			
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
			protected static RNGCryptoServiceProvider _rng = new RNGCryptoServiceProvider();

			/// <summary>
			/// Thread local random number generator (Used to generate the random gaussian noise
			/// </summary>
			protected ThreadLocal<Random> _rand;

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
			/// Coefficients of an FIR filter for the RF filter stage
			/// </summary>
			public double[] FIR_Coefficients_RF
			{
				get { return this._firCoefficients_RF; }
				set
				{
					this._firCoefficients_RF = value;
				}
			}
			protected double[] _firCoefficients_RF;

			/// <summary>
			/// Coefficients of an FIR filter for the DS filter stage
			/// </summary>
			public double[] FIR_Coefficients_DS
			{
				get { return this._firCoefficients_DS; }
				set
				{
					this._firCoefficients_DS = value;
				}
			}
			protected double[] _firCoefficients_DS;

			/// <summary>
			/// Coefficients of an FIR filter for the DS filter stage
			/// </summary>
			//public Filter.IIR_Filter.TIIRCoeff IIR_Coefficients_DS { get { return this._iirCoefficients_DS; } }
			//protected Filter.IIR_Filter.TIIRCoeff _iirCoefficients_DS;

			/// <summary>
			/// Coefficients of an FIR filter for the RF filter stage
			/// </summary>
			//public Filter.IIR_Filter.TIIRCoeff IIR_Coefficients_RF { get { return this._iirCoefficients_RF; } }
			//protected Filter.IIR_Filter.TIIRCoeff _iirCoefficients_RF;

			/// <summary>
			/// Denominator array of coefficients for the DS filter
			/// </summary>
			public double[] FilterCoefficients_A_DS { get { return this._filterCoefficients_A_DS; } }
			protected double[] _filterCoefficients_A_DS;

			/// <summary>
			/// Numerator array of coefficients for the DS filter
			/// </summary>
			public double[] FilterCoefficients_B_DS { get { return this._filterCoefficients_B_DS; } }
			protected double[] _filterCoefficients_B_DS;

			/// <summary>
			/// Denominator array of coefficients for the DS filter
			/// </summary>
			public double[] FilterCoefficients_A_RF { get { return this._filterCoefficients_A_RF; } }
			protected double[] _filterCoefficients_A_RF;

			/// <summary>
			/// Numerator array of coefficients for the DS filter
			/// </summary>
			public double[] FilterCoefficients_B_RF { get { return this._filterCoefficients_B_RF; } }
			protected double[] _filterCoefficients_B_RF;

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
			/// From I_SpreadingCode, Code matrix to be used to spread signals
			/// </summary>
			public byte[,] CodeMatrix
			{
				get
				{
					return this._codeMatrix;
				}
				set
				{
					this._codeMatrix = value;
				}
			}
			protected byte[,] _codeMatrix;

			/// <summary>
			/// Number of chips per symbol
			/// </summary>
			public int NumChips
			{
				get
				{
					return this._numChips;
				}

				set
				{
					if (value > 0 && _samplesPerSymbol % value == 0)
						this._numChips = value;
					else
						throw new ArgumentOutOfRangeException("NumChips", "Value must be postive and evenly divide into SamplesPerSymbol");
				}
			}
			protected int _numChips;

			/// <summary>
			/// Number of symbols per waveform
			/// </summary>
			public int NumberSymbolsPerWaveform
			{
				get
				{
					return this._numberSymbolsPerWaveform;
				}

				set
				{
					if (value >= 0)
						this._numberSymbolsPerWaveform = value;
					else
						throw new ArgumentOutOfRangeException("NumberSymbolsPerWaveform", "Must be positive.");
				}
			}
			protected int _numberSymbolsPerWaveform;

			/// <summary>
			/// Number of unique communications symbols available for use
			/// </summary>
			public int M
			{
				get { return this._M; }
				set
				{
					if (value > 1)
						this._M = value;
					else
						throw new ArgumentOutOfRangeException("M", "The number of unique symbols M must be 2 or greater.");
				}
			}
			protected int _M;

			/// <summary>
			/// Holds the reference waveforms for correlation implelmentation of a matched filter
			/// </summary>
			public double[][] Reference { get { return this._reference; } }
			protected double[][] _reference;

			/// <summary>
			/// Bit arrays corresponding to each communications symbol
			/// </summary>
			public byte[] BitsToCommunicationsSymbols
			{
				get { return this._bitsToCommunicationsSymbols; }
				set { this._bitsToCommunicationsSymbols = value; }
			}
			protected byte[] _bitsToCommunicationsSymbols;

			/// <summary>
			/// Supplies methods to spread waveforms
			/// </summary>
			public SpreadingCode CodeDivision;
			protected SpreadingCode _codeDivision;


			#endregion

			#region Constructors

			/// <summary>
			/// Creates a new instance of the ModulationType class with seed as the seed for the RNG
			/// </summary>
			/// <param name="seed">Seed to initialize the RNG</param>
			public ModulationType(double carrierFrequency, int samplingRate, int samplesPerSymbol, double signalPower, double[] firCoefficients, int numberSymbolsPerWaveform, int M)
			{
				this._rand = new ThreadLocal<Random>(() => new Random(Guid.NewGuid().GetHashCode()));
				this.CarrierFrequency = carrierFrequency;
				this.SampleRate = samplingRate;
				this.SamplesPerSymbol = samplesPerSymbol;
				this.SignalPower = signalPower;
				this.FIR_Coefficients_RF = firCoefficients;
				this.NumberSymbolsPerWaveform = numberSymbolsPerWaveform;
				this.M = M;

				//Initalize the reference array
				this._reference = new double[M][];
				for (int i = 0; i < this._reference.Length; i++)
				{
					this._reference[i] = new double[samplesPerSymbol];
				}
			}

			/// <summary>
			/// Copy constructor. Makes a new copy of a previous instance
			/// </summary>
			/// <param name="old">Previous instance to make a new copy of</param>
			public ModulationType(ModulationType old)
				: this(old._carrierFrequency, old._sampleRate, old._samplesPerSymbol, old._signalPower, old._firCoefficients_RF.NullCloneSafely(), old._numberSymbolsPerWaveform, old._M)
			{
				this._reference = (double[][])old._reference.Clone();
			}

			#endregion

			#region Modulation

			/// <summary>
			/// (Needs to be implemented in inheriting class)
			/// Modulates numberBits and produces a waveform at the carrier frequency and with the given sampling rate
			/// </summary>
			/// <param name="numberBits">Number of bits to modulate</param>
			/// <returns>Waveform containg the modulated bits</returns>
			public double[] ModulateBits()
			{
				byte[] bits = new byte[this._numberSymbolsPerWaveform];
				//Generate and modulate bits
				this.GenerateRandomBits(ref bits);
				return this.ModulateBits(bits);
			}

			/// <summary>
			/// Generates random bits in an existing array, then modulates and produces a waveform at the carrier frequency and with the given sampling rate
			/// </summary>
			/// <param name="numberBits">Number of bits to modulate</param>
			/// <param name="bits">Byte array to overwrite with new random bits</param>
			/// <returns>Waveform containg the modulated bits</returns>
			public double[] ModulateRandomBits(ref byte[] bits)
			{
				//Generate and modulate bits
				this.GenerateRandomBits(ref bits);
				return this.ModulateBits(bits);
			}
			
			/// <summary>
			/// Modulates bits based on the reference signals generated according to each bit
			/// </summary>
			/// <param name="bitsToModulate">The collection of bits to modulate.
			/// (Does not produce a random set of bits to modulate)</param>
			/// <returns>Waveform containing the modulated bits</returns>
			public virtual double[] ModulateBits(byte[] bitsToModulate)
			{
				double[] waveform = new double[_numberSymbolsPerWaveform * _samplesPerSymbol];

				//Loop through each symbol
				for (int i = 0; i < _numberSymbolsPerWaveform; i++)
				{
					//Loop through each chip
					for (int k = 0; k < _numChips; k++)
					{
						_reference[bitsToModulate[i]].CopyTo(waveform, i * _samplesPerSymbol + k * _samplesPerSymbol / _numChips);
					}
				}
				return waveform;
			}

			#endregion

			#region Demodulation

			/// <summary>
			/// Demodulates a waveform by estimating the bits using a correlation receiver structure
			/// </summary>
			/// <param name="waveform">Waveform to estimate bits from</param>
			/// <returns>Estimated bits from the waveform</returns>
			/// <remarks>This method assumes perfect synchronization with the received waveform, and perfect symbol boundaries.</remarks>
			public byte[] CorrelationReceiver(double[] waveform)
			{
				//Create array to hold the demodulated bits
				byte[] bits = new byte[waveform.Length / this._samplesPerSymbol];

				//Create an array to hold all of the z-Scores
				double[] z = new double[M];
				int bigZ = 0;

				//Loop through each symbol
				for (int i = 0; i < bits.Length; i++)
				{
					//Reset big Z
					bigZ = 0;

					//Reset the z scores
					for (int j = 0; j < M; j++)
					{
						z[j] = 0;
					}

					//Loop through each sample i the signal
					for (int k = 0; k < _samplesPerSymbol; k++)
					{
						//Correlate the symbol to the reference signals
						for (int j = 0; j < M; j++)
						{
							z[j] += waveform[i * _samplesPerSymbol + k] * _reference[j][k];
						}
					}

					for (int j = 0; j < M; j++)
					{
						//If there is a new biggest z-score, then save its bit
						if (z[j] > z[bigZ])
							bigZ = j;
					}

					bits[i] = this._bitsToCommunicationsSymbols[bigZ];
				}

				//Return the bit array
				return bits;
			}

			#endregion

			#region Filters

			/// <summary>
			/// Implements FIR filtering on the signal using the coefficients given when the instance was constructed
			/// </summary>
			/// <param name="signal">Signal to filter</param>
			/// <param name="filt">Use the RF or DS filter (RF by default)</param>
			/// <returns>Filtered signal</returns>
			public double[] FIR_Filter(double[] signal, EFilterToUse filt = EFilterToUse.RF)
			{
				double[] filter;

				if (filt == EFilterToUse.DS)
					filter = this._firCoefficients_DS;
				else
					filter = this._firCoefficients_RF;

				int numTaps = filter.Length;
				int numSigPoints = signal.Length;
				int top = 0, n, k;

				double[] reg = new double[numTaps];
				double[] filteredSignal = new double[numSigPoints];
				double y;

				for (int j = 0; j < numSigPoints; j++)
				{
					reg[top] = signal[j];
					y = 0;
					n = 0;

					for (k = top; k >= 0; k--)
					{
						y += filter[n++] * reg[k];
					}
					for (k = numTaps - 1; k > top; k--)
					{
						y += filter[n++] * reg[k];
					}
					filteredSignal[j] = y;
					top++;
					if (top >= numTaps)
					{
						top = 0;
					}
				}

				return filteredSignal;
			}
			
			/// <summary>
			/// Set the ds filter coefficients to those matching the gold codes
			/// </summary>
			public double[] GoldFIR_DS_Coefficients
			{
				get
				{
					return new double[] { -0.009221, -0.022961, -0.049499, -0.072660, -0.062240, 0.000000, 0.099730, 0.192645, 0.230519, 0.192645, 0.099730, 0.000000, -0.062240, -0.072660, -0.049499, -0.022961, -0.009221 };
				}
				set
				{
					this._firCoefficients_DS = new double[] { -0.009221, -0.022961, -0.049499, -0.072660, -0.062240, 0.000000, 0.099730, 0.192645, 0.230519, 0.192645, 0.099730, 0.000000, -0.062240, -0.072660, -0.049499, -0.022961, -0.009221 };
				}
			}

			/// <summary>
			/// Set the rf filter coefficients to those matching the gold codes
			/// </summary>
			public double[] GoldFIR_RF_Coefficients
			{
				get
				{
					return new double[] { -0.009168, -0.022866, -0.049365, -0.072549, -0.062206, 0.000000, 0.099806, 0.192854, 0.230794, 0.192854, 0.099806, 0.000000, -0.062206, -0.072549, -0.049365, -0.022866, -0.009168 };
				}
				set
				{
					this._firCoefficients_RF = new double[] { -0.009168, -0.022866, -0.049365, -0.072549, -0.062206, 0.000000, 0.099806, 0.192854, 0.230794, 0.192854, 0.099806, 0.000000, -0.062206, -0.072549, -0.049365, -0.022866, -0.009168 };
				}
				
			}

			/// <summary>
			/// Set the ds filter coefficients to those matching the hadamard codes
			/// </summary>
			public double[] HadamardFIR_DS_Coefficients
			{
				get
				{
					return new double[] { -0.009221, -0.022961, -0.049499, -0.072660, -0.062240, 0.000000, 0.099730, 0.192645, 0.230519, 0.192645, 0.099730, 0.000000, -0.062240, -0.072660, -0.049499, -0.022961, -0.009221 };
				}
				set
				{
					this._firCoefficients_DS = new double[] { -0.009221, -0.022961, -0.049499, -0.072660, -0.062240, 0.000000, 0.099730, 0.192645, 0.230519, 0.192645, 0.099730, 0.000000, -0.062240, -0.072660, -0.049499, -0.022961, -0.009221 };
				}
			}

			/// <summary>
			/// Set the rf filter coefficients to those matching the hadamard codes
			/// </summary>
			public double[] HadamardFIR_RF_Coefficients
			{
				get
				{
					return new double[] { -0.009165, -0.022861, -0.049358, -0.072544, -0.062204, 0.000000, 0.099810, 0.192864, 0.230808, 0.192864, 0.099810, 0.000000, -0.062204, -0.072544, -0.049358, -0.022861, -0.009165 };
				}
				set
				{
					this._firCoefficients_RF = new double[] { -0.009165, -0.022861, -0.049358, -0.072544, -0.062204, 0.000000, 0.099810, 0.192864, 0.230808, 0.192864, 0.099810, 0.000000, -0.062204, -0.072544, -0.049358, -0.022861, -0.009165 };
				}
			}

			#region old filter crap
			/*
			/// <summary>
			/// Gets the IIR coefficients for the DS filter stage or sets the internal coefficients.
			/// NOTE: Coefficients are set up for gold code usage
			/// NOTE: The value passed is not used
			/// </summary>
			public Filter.IIR_Filter.TIIRCoeff GoldIIR_DS_Coefficients
			{
				get
				{
					Filter.IIR_Filter.TIIRCoeff goldDS;
					Filter.getIIRCoefficients_Gold_DS(out goldDS);
					return goldDS;
				}
				set
				{
					Filter.getIIRCoefficients_Gold_DS(out _iirCoefficients_DS);
				}
			}

			/// <summary>
			/// Gets the IIR coefficients for the RF filter stage or sets the internal coefficients.
			/// NOTE: Coefficients are set up for gold code usage
			/// NOTE: The value passed is not used
			/// </summary>
			public Filter.IIR_Filter.TIIRCoeff GoldIIR_RF_Coefficients
			{
				get
				{
					Filter.IIR_Filter.TIIRCoeff goldRF;
					Filter.getIIRCoefficients_Gold_RF(out goldRF);
					return goldRF;
				}
				set
				{
					Filter.getIIRCoefficients_Gold_RF(out _iirCoefficients_RF);
				}
			}

			/// <summary>
			/// Gets the IIR coefficients for the DS filter stage or sets the internal coefficients.
			/// NOTE: Coefficients are set up for hadamard code usage
			/// NOTE: The value passed is not used
			/// </summary>
			public Filter.IIR_Filter.TIIRCoeff HadamardIIR_DS_Coefficients
			{
				get
				{
					Filter.IIR_Filter.TIIRCoeff goldDS;
					Filter.getIIRCoefficients_Hadamard_DS(out goldDS);
					return goldDS;
				}
				set
				{
					Filter.getIIRCoefficients_Hadamard_DS(out _iirCoefficients_DS);
				}
			}

			/// <summary>
			/// Gets the IIR coefficients for the RF filter stage or sets the internal coefficients.
			/// NOTE: Coefficients are set up for hadamard code usage
			/// NOTE: The value passed is not used
			/// </summary>
			public Filter.IIR_Filter.TIIRCoeff HadamardIIR_RF_Coefficients
			{
				get
				{
					Filter.IIR_Filter.TIIRCoeff goldRF;
					Filter.getIIRCoefficients_Hadamard_RF(out goldRF);
					return goldRF;
				}
				set
				{
					Filter.getIIRCoefficients_Hadamard_RF(out _iirCoefficients_RF);
				}
			}

			/// <summary>
			/// Set the IIR coefficients to those needed by the hadamard codes generated from matlab
			/// </summary>
			public void setIIRFiltersForHadamardCodes()
			{
				Filter.getIIRCoefficients_Hadamard_RF(out _iirCoefficients_RF);
				Filter.getIIRCoefficients_Hadamard_DS(out _iirCoefficients_DS);

			}

			/// <summary>
			/// Set the IIR coefficients to those needed by the gold codes generated from matlab
			/// </summary>
			public void setIIRFiltersForGoldCodes()
			{
				Filter.getIIRCoefficients_Gold_RF(out _iirCoefficients_RF);
				Filter.getIIRCoefficients_Gold_DS(out _iirCoefficients_DS);

			}

			/// <summary>
			/// Sets the numerator and denominator coefficients to those applicable to the gold codes
			/// </summary>
			public void setFilterCoefficients_Gold()
			{
				Filter.getDenominatorFilterCoefficients_Gold_DS(out _filterCoefficients_A_DS);
				Filter.getNumeratorFilterCoefficients_Gold_DS(out _filterCoefficients_B_DS);
				Filter.getDenominatorFilterCoefficients_Gold_RF(out _filterCoefficients_A_RF);
				Filter.getNumeratorFilterCoefficients_Gold_RF(out _filterCoefficients_B_RF);
			}

			/// <summary>
			/// Sets the numerator and denominator coefficients to those applicable to the gold codes
			/// </summary>
			public void setFilterCoefficients_Hadamard()
			{
				Filter.getDenominatorFilterCoefficients_Hadamard_DS(out _filterCoefficients_A_DS);
				Filter.getNumeratorFilterCoefficients_Hadamard_DS(out _filterCoefficients_B_DS);
				Filter.getDenominatorFilterCoefficients_Hadamard_RF(out _filterCoefficients_A_RF);
				Filter.getNumeratorFilterCoefficients_Hadamard_RF(out _filterCoefficients_B_RF);
			}

			public void setFIRFilterCoefficients_Gold()
			{
				Filter.getFIRFilterCoefficients_Gold_DS(out _firCoefficients_DS);
				Filter.getFIRFilterCoefficients_Gold_RF(out _firCoefficients_RF);
			}

			public void setFIRFilterCoefficients_Hadamard()
			{
				Filter.getFIRFilterCoefficients_Hadamard_DS(out _firCoefficients_DS);
				Filter.getFIRFilterCoefficients_Hadamard_RF(out _firCoefficients_RF);
			}
			*/
			#endregion

			#endregion

			#region Bit functions

			/// <summary>
			/// Generates numberBits of random bits.
			/// </summary>
			/// <param name="numberBits">The number of bits to generate</param>
			/// <returns>Byte array containg the generated bits</returns>
			public void GenerateRandomBits(ref byte[] bits)
			{
				//Fill with random bits
				_rng.GetBytes(bits);

				for (int i = 0; i < bits.Length; i++)
				{
					bits[i] = (byte)(bits[i] % M);
				}
			}

			/// <summary>
			/// Generates numberBits of random bits.
			/// </summary>
			/// <returns>Byte array containg the generated bits</returns>
			public byte[] GenerateRandomBits()
			{
				byte[] bits = new byte[_numberSymbolsPerWaveform];
				//Fill with random bits
				_rng.GetBytes(bits);

				for (int i = 0; i < bits.Length; i++)
				{
					bits[i] = (byte)(bits[i] % M);
				}
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
				for (int i = 0; i < inBits.Length; i++)
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

			#endregion

			#region Noise

			/// <summary>
			/// Creates an array of randomly generated white gaussian noise
			/// </summary>
			/// <param name="_numberSymbolsPerWaveform">Number of bits long the array should be</param>
			/// <param name="coefficient">Coefficient that determines the power of the noise</param>
			/// <returns>Noise array with specified power and length</returns>
			public double[] AdditiveWhiteGaussianNoise(double coefficient)
			{
				//Create an array to hold the noise sample
				double[] awgn = new double[this._numberSymbolsPerWaveform * this._samplesPerSymbol];
				double x, y;

				if (coefficient != 0)
				{
					//Loop through each sample
					for (int i = 0; i < awgn.Length; i++)
					{
						//Uses the Box-Muller transform to get a Gaussian distributed RV
						x = _rand.Value.NextDouble();
						y = _rand.Value.NextDouble();
						awgn[i] = Math.Sqrt(coefficient) * (Math.Sqrt(-2 * Math.Log(x)) * Math.Sin(2 * Math.PI * y));
					}
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
						x = _rand.Value.NextDouble();
						y = _rand.Value.NextDouble();
						awgn[i] += Math.Sqrt(coefficient) * (Math.Sqrt(-2 * Math.Log(x)) * Math.Sin(2 * Math.PI * y));
						waveform[i] += awgn[i];
						mean += awgn[i];
					}
				}
				return WaveformStatistics.PowerInWaveform(awgn, mean / waveform.Length);
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
						x = _rand.Value.NextDouble();
						y = _rand.Value.NextDouble();
						waveform[i] += Math.Sqrt(coefficient) * (Math.Sqrt(-2 * Math.Log(x)) * Math.Sin(2 * Math.PI * y));
					}
				}
			}

			#endregion
			
			#region Helper Functions

			/// <summary>
			/// Generate a time vector starting at 0
			/// </summary>
			/// <returns>Time vector</returns>
			protected double[] getTimeArray(int vecLength)
			{
				//Generate the (full length) time vector
				var time = new double[vecLength];

				//Calculate the times
				for (int i = 0; i < vecLength; i++)
				{
					time[i] = (double)i / this._sampleRate;
				}

				//Return the time array
				return time;
			}
			#endregion

			#region Code Division

			/// <summary>
			/// Initialize the spreading code instance
			/// </summary>
			/// <param name="numberChips">Number of chips per symbol</param>
			/// <param name="codeMatrix">Code matrix to use when spreading (If null, uses Gold codes of size 31)</param>
			public void InitializeCodeDivision(int numberChips=31, byte[,] codeMatrix=null)
			{
				this._codeDivision = new SpreadingCode(numberChips, _samplesPerSymbol, _numberSymbolsPerWaveform, codeMatrix);
			}


			///<summary>
			/// Applies a spreading code to a waveform
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <param name="_numberSymbolsPerWaveform">Number of symbols in the waveform</param>
			/// <returns>Original waveform with the spreading code applied.</returns>
			public void SpreadWaveform(ref double[] waveform, int user)
			{
				_codeDivision.SpreadWaveform(ref waveform, user);
			}

			///<summary>
			/// Applies a despreading code to a waveform
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <param name="_numberSymbolsPerWaveform">Number of symbols in the waveform</param>
			/// <returns>Original waveform with the despreading code applied.</returns>
			public void DespreadWaveform(ref double[] waveform, int user)
			{
				_codeDivision.DespreadWaveform(ref waveform, user);
			}

			/// <summary>
			/// Adds a seperately generated and spread waveform to the given waveform
			/// </summary>
			/// <param name="waveform">Waveform to add another spread and modulated waveform to</param>
			/// <param name="user">User number of the added signal</param>
			public void AddSpreadWaveform(ref double[] waveform, int user, double noisePower)
			{
				_codeDivision.AddSpreadWaveform(ref waveform, user, noisePower, this);
			}

			#endregion
		}
	}
}
