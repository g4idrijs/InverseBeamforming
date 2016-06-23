using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Reactive.Linq;
using System.IO;

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
		/// Describes the current state of the simulation
		/// </summary>
		public struct IntermediateSimResults
		{
			/// <summary>
			/// Total bit errors produced already
			/// </summary>
			public int TotalErrors;

			/// <summary>
			/// Number of bit errors during this iteration of the simulation
			/// </summary>
			public int NumberErrorsThisIteration;

			/// <summary>
			/// Total number of bits simulated so far
			/// </summary>
			public int TotalBitsSimulated;

			/// <summary>
			/// Bit error rate of the simulation so far
			/// </summary>
			public double BitErrorRate
			{ get { return (double)TotalErrors / TotalBitsSimulated; } }

			/// <summary>
			/// Total errors remaining before the simulation is complete
			/// </summary>
			public int TotalErrorsRemaining;

			/// <summary>
			/// Percentage of the total errors required to end the simulation
			/// </summary>
			public double PercentErrorHad
			{ get { return (double)TotalErrors / (TotalErrors + TotalErrorsRemaining); } }

			/// <summary>
			/// Constructs an instance of the IntermediateSimResults Structure
			/// </summary>
			/// <param name="totalErrors">Total number of error occured in the simulation</param>
			/// <param name="numberErrorsThisIteration">Total number of errors in this iteration</param>
			/// <param name="totalBitsSimulated">Total number of bits simulated so far</param>
			/// <param name="totalErrorsRemaining">Total number of errors remaining</param>
			public IntermediateSimResults(int totalErrors, int numberErrorsThisIteration, int totalBitsSimulated, int totalErrorsRemaining)
			{
				TotalErrors = totalErrors;
				NumberErrorsThisIteration = numberErrorsThisIteration;
				TotalBitsSimulated = totalBitsSimulated;
				TotalErrorsRemaining = totalErrorsRemaining;
			}
		}

		/// <summary>
		/// Describes the final state of the simulation
		/// </summary>
		public struct FinalSimResults
		{
			/// <summary>
			/// Total bit errors produced already
			/// </summary>
			public int TotalErrors;

			/// <summary>
			/// Total number of bits simulated so far
			/// </summary>
			public int TotalBitsSimulated;

			/// <summary>
			/// Bit error rate of the simulation so far
			/// </summary>
			public double BitErrorRate
			{ get { return (double)TotalErrors / TotalBitsSimulated; } }

			/// <summary>
			/// Constructs an instance of the IntermediateSimResults Structure
			/// </summary>
			/// <param name="totalErrors">Total number of error occured in the simulation</param>
			/// <param name="totalBitsSimulated">Total number of bits simulated so far</param>
			public FinalSimResults(int totalErrors, int totalBitsSimulated)
			{
				TotalErrors = totalErrors;
				TotalBitsSimulated = totalBitsSimulated;
			}
		}

		/// <summary>
		/// Defines an interface that allows multiple modulation types to be created
		/// </summary>
		public abstract class ModulationType : I_SpreadingCode, IObservable<IntermediateSimResults>
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
			protected double[,] _reference;

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
			/// Provides a list to hold the suscribers of simulation events
			/// </summary>
			List<IObserver<IntermediateSimResults>> observers;

			/// <summary>
			/// Creates a new instance of the ModulationType class with seed as the seed for the RNG
			/// </summary>
			/// <param name="seed">Seed to initialize the RNG</param>
			public ModulationType(int seed, double carrierFrequency, int samplingRate, int samplesPerSymbol, double signalPower, double[] firCoefficients, int numberSymbolsPerWaveform, int M)
			{
				observers = new List<IObserver<IntermediateSimResults>>();
				this._rng = new Random(seed);
				this.CarrierFrequency = carrierFrequency;
				this.SampleRate = samplingRate;
				this.SamplesPerSymbol = samplesPerSymbol;
				this.SignalPower = signalPower;
				this.FIR_Coefficients = firCoefficients;
				this.NumberSymbolsPerWaveform = numberSymbolsPerWaveform;
				this.M = M;

				//Initalize the reference array
				this._reference = new double[M, samplesPerSymbol];
			}

			/// <summary>
			/// (Needs to be implemented in inheriting class)
			/// Modulates numberBits and produces a waveform at the carrier frequency and with the given sampling rate
			/// </summary>
			/// <param name="numberBits">Number of bits to modulate</param>
			/// <returns>Waveform containg the modulated bits</returns>
			public double[] ModulateBits()
			{
				//Generate and modulate bits
				return this.ModulateBits(this.GenerateRandomBits());
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
							z[j] += waveform[i * _samplesPerSymbol + k] * _reference[j, k];
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

			/// <summary>
			/// Implements FIR filtering on the signal using the coefficients given when the instance was constructed
			/// </summary>
			/// <param name="signal">Signal to filter</param>
			/// <returns>Filtered signal</returns>
			public double[] FIR_Filter(double[] signal)
			{

				int numTaps = this._firCoefficients.Length;
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
						y += this._firCoefficients[n++] * reg[k];
					}
					for (k = numTaps - 1; k > top; k--)
					{
						y += this._firCoefficients[n++] * reg[k];
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

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Bit functions
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/// <summary>
			/// Generates numberBits of random bits.
			/// </summary>
			/// <param name="numberBits">The number of bits to generate</param>
			/// <returns>Byte array containg the generated bits</returns>
			public byte[] GenerateRandomBits()
			{
				//Create the array
				var bits = new byte[this._numberSymbolsPerWaveform];

				//Fill each bit with a random value
				for (int i = 0; i < this._numberSymbolsPerWaveform; i++)
				{
					bits[i] = (byte)this._rng.Next(0, M);
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

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Noise
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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

				//Loop through each sample
				for (int i = 0; i < awgn.Length; i++)
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
				return Math.Pow(waveform.Mean(), 2) + waveform.Variance();
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
			/// <param name="noisePower">Spectral power of the noise</param>
			/// <returns>Overall bit error rate of the simulation</returns>
			public double RunSimulationOneNoisePowerIdealFiltering(int numberToGetWrongEventually, double noisePower)
			{
				byte[] inbits;
				double[] waveform;
				byte[] outbits;

				int totalNumWrong = 0, totalBitsSimulated = 0;

				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
				{
					//Generate Random bits
					inbits = GenerateRandomBits();

					//Modulate the bits
					waveform = ModulateBits(inbits);

					//Add noise to the waveform
					AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

					//Demodulate the waveform+noise
					outbits = CorrelationReceiver(waveform);

					//Get the total number of bits that were estimated wrongly
					totalNumWrong += numberDifferentBits(inbits, outbits);

					//Update the total number of bits simulated
					totalBitsSimulated += this._numberSymbolsPerWaveform;
				}

				//Return the bit error rate
				return totalNumWrong / (double)totalBitsSimulated;
			}

			/// <summary>
			/// Runs a simulation for the given parameters and outputs a bit error rate for that simulation (Real Filters)
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bits to eventually mis-estimate</param>
			/// <param name="noisePower">Spectral power of the noise</param>
			/// <returns>Overall bit error rate of the simulation</returns>
			public double RunSimulationOneNoisePowerRealFiltering(int numberToGetWrongEventually, double noisePower)
			{
				byte[] inbits;
				double[] waveform;
				byte[] outbits;

				int totalNumWrong = 0, totalBitsSimulated = 0;

				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
				{
					//Generate Random bits
					inbits = GenerateRandomBits();

					//Modulate the bits
					waveform = ModulateBits(inbits);

					//Add noise to the waveform
					AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

					//Filter the waveform
					waveform = FIR_Filter(waveform);

					//Demodulate the waveform+noise
					outbits = CorrelationReceiver(waveform);

					//Get the total number of bits that were estimated wrongly
					totalNumWrong += numberDifferentBits(inbits, outbits);

					//Update the total number of bits simulated
					totalBitsSimulated += this._numberSymbolsPerWaveform;
				}

				//Return the bit error rate
				return totalNumWrong / (double)totalBitsSimulated;
			}

			/// <summary>
			/// Runs a simulation over many given noise powers, and returns the bit error rate from each simulation
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bits to eventually mis-estimate</param>
			/// <param name="noisePowers">Array of spectral powers of noise to simulate</param>
			/// <returns>Bit error rates of every simulation</returns>
			public double[] RunSimulationManyNoisePowersIdealFiltering(int numberToGetWrongEventually, double[] noisePowers)
			{
				double[] bers = new double[noisePowers.Length];
				//Loop through each of the noise powers
				Parallel.For(0, noisePowers.Length - 1, i =>
				  {
				  //yield return the BER from that simulation
				  bers[i] = RunSimulationOneNoisePowerIdealFiltering(numberToGetWrongEventually, noisePowers[i]);
				  });
				return bers;
			}

			/// <summary>
			/// Runs a simulation over many given noise powers, and returns the bit error rate from each simulation (Real Filter)
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bits to eventually mis-estimate</param>
			/// <param name="noisePowers">Array of spectral powers of noise to simulate</param>
			/// <returns>Bit error rates of every simulation</returns>
			public double[] RunSimulationManyNoisePowersRealFiltering(int numberToGetWrongEventually, double[] noisePowers)
			{
				double[] bers = new double[noisePowers.Length];
				//Loop through each of the noise powers
				Parallel.For(0, noisePowers.Length - 1, i =>
				{
				//yield return the BER from that simulation
				bers[i] = RunSimulationOneNoisePowerRealFiltering(numberToGetWrongEventually, noisePowers[i]);
				});
				return bers;
			}

			/// <summary>
			/// Provides a way to unsubscribe from the simulation updates
			/// </summary>
			private class Unsubscriber : IDisposable
			{
				/// <summary>
				/// List of observers of the updates
				/// </summary>
				private List<IObserver<IntermediateSimResults>> _observers;
				
				/// <summary>
				/// Observer of the updates
				/// </summary>
				private IObserver<IntermediateSimResults> _observer;

				/// <summary>
				/// Construct a new instance of the Unsubscriber class
				/// </summary>
				/// <param name="observers">List of observers of the simulation updates</param>
				/// <param name="observer">New observer of the simulation updates</param>
				public Unsubscriber(List<IObserver<IntermediateSimResults>> observers, IObserver<IntermediateSimResults> observer)
				{
					this._observers = observers;
					this._observer = observer;
				}

				/// <summary>
				/// Dispose of the observer from the list
				/// </summary>
				public void Dispose()
				{
					//If the list isn't null, remove the observer from the list
					if (! (_observer==null))
					{
						_observers.Remove(_observer);
					}
				}
			}

			/// <summary>
			/// Subscribe to the simulation updates
			/// </summary>
			/// <param name="observer">New oberver of the simulation updates</param>
			/// <returns>An interface detailing how to unsubscribe from the updates</returns>
			public IDisposable Subscribe(IObserver<IntermediateSimResults> observer)
			{
				//If the observer is not already in the list of observers, add it
				if (!observers.Contains(observer))
					observers.Add(observer);

				//Return an Unsubscriber
				return new Unsubscriber(observers, observer);
			}

			/// <summary>
			/// Runs a simulation that reports its progress
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bit errors to simulate to</param>
			/// <param name="noisePower">Power of the AWGN added to the signals</param>
			/// <returns>Final results of the simulation</returns>
			public FinalSimResults RunSimpleSimulationObservable(int numberToGetWrongEventually, double noisePower)
			{
				byte[] inbits;
				double[] waveform;
				byte[] outbits;

				int numberWrongThisIteration = 0;
				int totalNumWrong = 0, totalBitsSimulated = 0;

				foreach (var obs in observers)
				{

					(obs as IntermediateSimResultsReporter).OnStart();
				}

				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
					{
						//Generate Random bits
						inbits = GenerateRandomBits();

						//Modulate the bits
						waveform = ModulateBits(inbits);

						//Add noise to the waveform
						AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

						//Demodulate the waveform+noise
						outbits = CorrelationReceiver(waveform);

						//Get the total number of bits that were estimated wrongly
						numberWrongThisIteration= numberDifferentBits(inbits, outbits);
						totalNumWrong += numberWrongThisIteration;

						//Update the total number of bits simulated
						totalBitsSimulated += this._numberSymbolsPerWaveform;

						if (numberWrongThisIteration!=0)
						{
							foreach (var obs in observers)
							{
								obs.OnNext(new IntermediateSimResults(totalNumWrong, numberWrongThisIteration, totalBitsSimulated, numberToGetWrongEventually - totalNumWrong));
							}
						}
					}

				//Simulation is done, inform the observers
				foreach (var obs in observers)
				{
					obs.OnCompleted();
				}

				return new FinalSimResults( totalNumWrong, totalBitsSimulated);
			}

			/// <summary>
			/// Runs a simulation that reports its progress
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bit errors to simulate to</param>
			/// <param name="noisePower">Power of the AWGN added to the signals</param>
			/// <returns>Final results of the simulation</returns>
			public FinalSimResults RunFIRFilterSimulationObservable(int numberToGetWrongEventually, double noisePower)
			{
				byte[] inbits;
				double[] waveform;
				byte[] outbits;

				int numberWrongThisIteration = 0;
				int totalNumWrong = 0, totalBitsSimulated = 0;

				//Call the OnStart function for each observer
				foreach (var obs in observers)
				{
					(obs as IntermediateSimResultsReporter).OnStart();
				}
				
				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
				{
					//Generate Random bits
					inbits = GenerateRandomBits();

					//Modulate the bits
					waveform = ModulateBits(inbits);

					//Add noise to the waveform
					AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

					//Filter the waveform
					waveform = FIR_Filter(waveform);

					//Demodulate the waveform+noise
					outbits = CorrelationReceiver(waveform);

					//Get the total number of bits that were estimated wrongly
					numberWrongThisIteration = numberDifferentBits(inbits, outbits);
					totalNumWrong += numberWrongThisIteration;

					//Update the total number of bits simulated
					totalBitsSimulated += this._numberSymbolsPerWaveform;

					//If there were any bit errors this iteration
					if (numberWrongThisIteration != 0)
					{
						//Call the OnNext function for each observer
						foreach (var obs in observers)
						{
							obs.OnNext(new IntermediateSimResults(totalNumWrong, numberWrongThisIteration, totalBitsSimulated, numberToGetWrongEventually - totalNumWrong));
						}
					}
				}

				//Simulation is done, inform the observers
				foreach (var obs in observers)
				{
					obs.OnCompleted();
				}

				//Return the final simulation results
				return new FinalSimResults(totalNumWrong, totalBitsSimulated);
			}

			/// <summary>
			/// Defines a class that reports the results of the simulation in real time
			/// </summary>
			public class IntermediateSimResultsReporter : IObserver<IntermediateSimResults>
			{
				/// <summary>
				/// Starting time of the simulation
				/// </summary>
				private DateTime startTime;
				/// <summary>
				/// Ending time of the simulation
				/// </summary>
				private DateTime endTime;

				/// <summary>
				/// Hold the information needed to unsubscribe from the simulation
				/// </summary>
				private IDisposable unsubscriber;

				/// <summary>
				/// Hold a filename to write the results simultanesouly
				/// </summary>
				private string _filename;

				/// <summary>
				/// Construct a new instance of the class with a file to write the results to
				/// </summary>
				/// <param name="filename">Name of the file to write the results to</param>
				public IntermediateSimResultsReporter(string filename)
				{
					this._filename = filename;
				}

				/// <summary>
				/// Subscribe to receive the updates from the simulation
				/// </summary>
				/// <param name="provider">Simulation OBservable</param>
				public virtual void Subscribe(IObservable<IntermediateSimResults> provider)
				{
					//Provide information to unsubscribe from the observable
					unsubscriber = provider.Subscribe(this);
				}

				/// <summary>
				/// Unsubscribe
				/// </summary>
				public virtual void Unsubscribe()
				{
					unsubscriber.Dispose();
				}

				/// <summary>
				/// Sets up the file to log the simulation
				/// </summary>
				public virtual void OnStart()
				{
					//Get the starting time
					startTime = DateTime.Now;

					//Write the start time to the simulation
					File.WriteAllText(_filename, "Start of simulation. Time: " + startTime.ToLongDateString() + " " + startTime.ToLongTimeString() + Environment.NewLine);
				}

				/// <summary>
				/// Function to be run when the simulation is done
				/// </summary>
				public virtual void OnCompleted()
				{
					//Get the end time of the simulation
					endTime = DateTime.Now;

					//Append some timing information to the end
					File.AppendAllText(_filename, String.Format("Completed simulation at {0} {1}.\n", endTime.ToLongDateString(), endTime.ToLongTimeString()));
					File.AppendAllText(_filename, "The simulation took: " + (endTime - startTime).ToString());
				}

				/// <summary>
				/// Function gets called when there is an error in the simulation
				/// </summary>
				/// <param name="error"></param>
				public virtual void OnError(Exception error)
				{
					endTime = DateTime.Now;
					File.AppendAllText(_filename, String.Format("An error was encoutered at {0} {1}.\n\n{1}", endTime.ToLongDateString(), endTime.ToLongTimeString(),error.Message));
				}

				/// <summary>
				/// Function gets called when an update to the simulation is pushed
				/// </summary>
				/// <param name="isr">Struct holding the intermediate simulation results</param>
				public virtual void OnNext(IntermediateSimResults isr)
				{
					//Write the information to the log file
					File.AppendAllText(_filename, String.Format("{0} {1} Total Errors: {2, 5}, Errors this iteration: {3, 5}, Total Bits Simulated: {4, 9}, Bit Error Rate: {5: 0.00}, Percent errors found: {6: 0.00}\n", DateTime.Now.ToLongDateString(), DateTime.Now.ToLongTimeString(), isr.TotalErrors, isr.NumberErrorsThisIteration, isr.TotalBitsSimulated, isr.BitErrorRate, isr.PercentErrorHad*100));
				}
				
			}

			 

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Spreading Codes
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////

			/// <summary>
			/// Gets a full waveform length spreading code, ready to multiply with the signal. 
			/// </summary>
			/// <param name="user">Index of the spreading code to use (zero based)</param>
			/// <param name="_samplesPerSymbol">Number of samples per symbol duration</param>
			/// <param name="_numChips">Number of chips per symbol duration</param>
			/// <param name="numSymbols">Number of symbols</param>
			/// <returns>Array containing the spreading code ready to be mixed with the signal</returns>
			public byte[] GetSpreadingCode(int user, int numSymbols)
			{
				//Create a matrix for the spreading code
				byte[] spreadingCode = new byte[_samplesPerSymbol * numSymbols];
				int j = 0;

				//If the number of samples per chip is not an integer, throw an exception
				if (_samplesPerSymbol % _numChips != 0)
					throw new ArgumentException("The number of samples per chip is not an integer.", "numSamples, numChips");

				//Loop through each symbol
				for (int i = 0; i < numSymbols; i++)
				{
					//Loop through each chip
					for (int k = 0; k < _numChips; i++)
					{
						//Loop through the samples in each chip
						for (j = 0; j < _samplesPerSymbol / _numChips; i++)
						{//              | symbol offset |       chip offset         | left in chip                Chip number
							spreadingCode[i * _samplesPerSymbol + k * _samplesPerSymbol / _numChips + j] = this._codeMatrix[user, k];
						}
					}
				}

				//Return the spreading code
				return spreadingCode;
			}

			///<summary>
			/// Applies a spreading code to a waveform
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <param name="_samplesPerSymbol">Number of samples per symbol</param>
			/// <param name="_numChips">Number of chips per symbol</param>
			/// <param name="_numberSymbolsPerWaveform">Number of symbols in the waveform</param>
			/// <returns>Original waveform with the spreading code applied.</returns>
			public void SpreadWaveform(ref double[] waveform, int user)
			{
				//If the number of samples per chip is not an integer, throw an exception
				if (_samplesPerSymbol % _numChips != 0)
					throw new ArgumentException("The number of samples per chip is not an integer.", "numSamples, numChips");

				//Loop through each symbol
				for (int i = 0; i < _numberSymbolsPerWaveform; i++)
				{
					//Loop through each chip
					for (int k = 0; k < _numChips; k++)
					{
						//Loop through the samples in each chip
						for (int j = 0; j < _samplesPerSymbol / _numChips; j++)
						{//         | symbol offset |       chip offset         | left in chip										Chip number
							waveform[i * _samplesPerSymbol + k * _samplesPerSymbol / _numChips + j] *= ((double)this._codeMatrix[user, k] * 2 - 1);
						}
					}
				}
			}

			/// <summary>
			/// Applies a spreading code to a waveform
			/// </summary>
			/// <param name="waveform">Waveform to spread</param>
			/// <param name="user">Specific spreading code to use (row in the spreading code matrix)</param>
			/// <param name="_samplesPerSymbol">Number of samples per symbol</param>
			/// <param name="_numChips">Number of chips per symbol</param>
			/// <param name="_numberSymbolsPerWaveform">Number of symbols in the waveform</param>
			/// <returns>Original waveform with the spreading code applied.</returns>
			/// <remarks>In the implementation, it is functionally equivalent to SpreadWaveform. Applying the same spreading code twice gives the original signal back.</remarks>
			public void DespreadWaveform(ref double[] waveform, int user)
			{
				SpreadWaveform(ref waveform, user);
			}

			/// <summary>
			/// Initializes the components necesary to utilize the spreading functions
			/// </summary>
			/// <param name="codeMatrix">Code matrix to use to spread the waveforms</param>
			/// <param name="numChips">Number of chips per symbol</param>
			public void initializeSpreadingCodes(byte[,] codeMatrix, int numChips)
			{
				this.NumChips = numChips;
				this.CodeMatrix = codeMatrix;
			}

			/// <summary>
			/// Returns the gold codes of length 31
			/// </summary>
			/// <returns>Byte array containing the gold codes of length 31</returns>
			private byte[,] getGoldCodes()
			{
				return new byte[,] {{0,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1,0,0,0,0,1},
								{1,1,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0,0,1,1,1,1,1,0,1,1,1,0,0,0,1,0},
								{1,1,0,0,1,1,0,0,0,1,1,0,0,1,0,1,0,1,0,1,0,0,0,0,0,1,0,0,0,0,1,1},
								{0,0,1,0,0,0,1,0,0,1,0,0,1,1,1,0,0,0,0,1,0,0,1,1,0,1,1,0,0,1,0,0},
								{1,1,1,1,1,1,1,0,0,0,0,1,1,0,0,0,1,0,0,1,0,1,0,1,0,0,1,0,1,0,1,1},
								{1,1,0,0,0,1,1,0,1,0,1,1,0,1,0,1,1,0,0,1,1,0,0,1,1,0,1,1,0,1,0,0},
								{0,0,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,0,0,0,0,0,0,0,1,0,0,0,1,0,1,0},
								{1,1,0,1,0,1,0,1,0,1,0,1,1,0,1,1,1,0,1,1,0,0,1,0,1,1,1,1,0,1,1,1},
								{0,0,0,1,0,0,0,0,0,0,1,1,0,0,1,1,1,1,0,1,0,1,1,0,0,0,0,0,1,1,0,0},
								{0,0,0,1,1,0,1,0,1,1,1,0,0,0,1,1,0,0,0,1,1,1,1,1,1,1,1,1,1,0,1,1},
								{0,0,0,0,1,1,1,1,0,1,0,0,0,0,1,0,1,0,0,0,1,1,0,0,0,0,0,1,0,1,0,1},
								{0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,1,1,0,1,0,1,0,1,1,1,1,0,0,1,0,0,1},
								{1,1,1,1,0,0,1,0,1,0,0,0,0,1,1,1,1,1,1,0,0,1,0,0,0,1,1,1,0,0,0,1},
								{1,1,0,1,1,1,1,1,1,0,0,0,1,0,1,1,0,1,1,1,1,0,1,1,0,0,0,0,0,0,0,0},
								{0,0,0,0,0,1,0,1,1,0,0,1,0,0,1,0,0,1,0,0,0,1,0,1,1,1,1,0,0,0,1,0},
								{0,0,1,1,0,0,0,1,1,0,1,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,1,1},
								{1,1,0,1,1,0,0,1,1,1,0,0,0,1,0,0,1,1,0,0,0,0,1,1,1,0,1,0,1,1,0,1},
								{0,0,0,0,1,0,0,1,0,0,0,0,1,1,0,1,0,0,1,1,0,1,0,0,1,0,1,1,1,0,0,0},
								{0,0,1,0,1,0,0,0,1,0,0,1,1,1,1,0,1,1,0,1,1,0,1,0,1,0,0,1,0,0,1,1},
								{1,1,1,0,1,0,1,1,1,0,1,1,1,0,0,1,0,0,0,0,0,1,1,0,1,1,0,0,0,1,0,1},
								{1,1,1,0,1,1,0,1,1,1,1,1,0,1,1,0,1,0,1,1,1,1,1,0,0,1,1,0,1,0,0,0},
								{1,1,1,0,0,0,0,1,0,1,1,0,1,0,0,1,1,1,0,0,1,1,1,1,0,0,1,1,0,0,1,0},
								{1,1,1,1,1,0,0,0,0,1,0,1,0,1,1,1,0,0,1,0,1,1,0,1,1,0,0,0,0,1,1,0},
								{1,1,0,0,1,0,1,0,0,0,1,0,1,0,1,0,1,1,1,0,1,0,0,0,1,1,1,0,1,1,1,0},
								{0,0,1,0,1,1,1,0,1,1,0,1,0,0,0,1,0,1,1,0,0,0,1,0,0,0,1,1,1,1,1,0},
								{1,1,1,0,0,1,1,1,0,0,1,0,0,1,1,0,0,1,1,1,0,1,1,1,1,0,0,1,1,1,1,1},
								{1,1,1,1,0,1,0,0,1,1,0,0,1,0,0,0,0,1,0,1,1,1,0,0,1,1,0,1,1,1,0,0},
								{1,1,0,1,0,0,1,1,0,0,0,1,0,1,0,0,0,0,0,0,1,0,1,0,0,1,0,1,1,0,1,0},
								{0,0,0,1,1,1,0,0,1,0,1,0,1,1,0,0,1,0,1,0,0,1,1,1,0,1,0,1,0,1,1,0},
								{0,0,0,0,0,0,1,1,1,1,0,1,1,1,0,1,1,1,1,1,1,1,0,1,0,1,0,0,1,1,1,1},
								{0,0,1,1,1,1,0,1,0,0,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0,1,1,1,1,1,0,1},
								{1,1,0,0,0,0,0,0,1,1,1,1,1,0,1,0,0,0,1,0,0,0,0,1,0,0,0,1,1,0,0,1},
								{0,0,1,1,1,0,1,1,0,1,1,1,0,0,0,0,1,1,1,1,0,0,0,1,1,1,0,1,0,0,0,0}};
			}

			//////////////////////////////////////////////////////////////////////////////////////////////////////////////
			//Helper Functions
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
		}
	}
}
