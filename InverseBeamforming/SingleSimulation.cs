using System;
using System.Collections.Generic;
using System.IO;
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
	/// Holds methods that create and report on simulations of Bit error rate
	/// </summary>
	public partial class Simulations
	{
		/// <summary>
		/// Class provides methods to complete a simulation of bit error rates
		/// </summary>
		public class SingleSimulation
		{
			/// <summary>
			/// Modulation type, and most of the parameters used in the simulation having to deal with modulation
			/// </summary>
			protected Modulations.ModulationType _modOriginal;

			/// <summary>
			/// Current intermediate state of the simulation
			/// </summary>
			public IntermediateSimResults isr;

			/// <summary>
			/// Reporter used to report progress on the simulation
			/// </summary>
			protected List<SimulationReporter> _reporters;

			/// <summary>
			/// Construct a new simulation
			/// </summary>
			/// <param name="mod">Modulation class that contains most of the parameters of the simulation</param>
			/// <param name="reporter">Reporter to use to report progress on the simulation</param>
			public SingleSimulation(Modulations.ModulationType mod, List<SimulationReporter> reporters)
			{
				_modOriginal = mod;
				_reporters = reporters;
			}

			/// <summary>
			/// Runs a simulation that reports its progress
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bit errors to simulate to</param>
			/// <param name="noisePower">Power of the AWGN added to the signals</param>
			/// <returns>Final results of the simulation</returns>
			public FinalSimResults RunSimpleSimulationObservable(int numberToGetWrongEventually, double noisePower)
			{
				byte[] inbits = new byte[_modOriginal.NumberSymbolsPerWaveform];
				double[] waveform;
				byte[] outbits;

				int numberWrongThisIteration = 0;
				int totalNumWrong = 0, totalBitsSimulated = 0;

				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
				{
					//Modulate the bits
					waveform = _modOriginal.ModulateRandomBits(ref inbits);

					//Add noise to the waveform
					_modOriginal.AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

					//Demodulate the waveform+noise
					outbits = _modOriginal.CorrelationReceiver(waveform);

					//Get the total number of bits that were estimated wrongly
					numberWrongThisIteration = _modOriginal.numberDifferentBits(inbits, outbits);
					totalNumWrong += numberWrongThisIteration;

					//Update the total number of bits simulated
					totalBitsSimulated += _modOriginal.NumberSymbolsPerWaveform;

					//If there were any bit errors this iteration, report them
					if (numberWrongThisIteration != 0)
					{
						isr = new IntermediateSimResults(totalNumWrong, numberWrongThisIteration, totalBitsSimulated, numberToGetWrongEventually - totalNumWrong, noisePower.ToString());
						_reporters.OnNext_List(this);
					}
				}

				//Simulation is done, inform the observers
				_reporters.OnCompleted_List();

				//Return the final simulation results
				return new FinalSimResults(totalNumWrong, totalBitsSimulated, noisePower, _modOriginal);
			}

			/// <summary>
			/// Runs a simulation that reports its progress
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bit errors to simulate to</param>
			/// <param name="noisePower">Power of the AWGN added to the signals</param>
			/// <returns>Final results of the simulation</returns>
			public FinalSimResults RunFIRFilterSimulationObservable(int numberToGetWrongEventually, double noisePower)
			{
				byte[] inbits = new byte[_modOriginal.NumberSymbolsPerWaveform];
				double[] waveform;
				byte[] outbits;

				int numberWrongThisIteration = 0;
				int totalNumWrong = 0, totalBitsSimulated = 0;

				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
				{
					//Modulate the bits
					waveform = _modOriginal.ModulateRandomBits(ref inbits);

					//Add noise to the waveform
					_modOriginal.AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

					//Filter the waveform
					waveform = _modOriginal.FIR_Filter(waveform);

					//Demodulate the waveform+noise
					outbits = _modOriginal.CorrelationReceiver(waveform);

					//Get the total number of bits that were estimated wrongly
					numberWrongThisIteration = _modOriginal.numberDifferentBits(inbits, outbits);
					totalNumWrong += numberWrongThisIteration;

					//Update the total number of bits simulated
					totalBitsSimulated += _modOriginal.NumberSymbolsPerWaveform;

					//If there were any bit errors this iteration, report them
					if (numberWrongThisIteration != 0)
					{
						isr = new IntermediateSimResults(totalNumWrong, numberWrongThisIteration, totalBitsSimulated, numberToGetWrongEventually - totalNumWrong, noisePower.ToString());
						_reporters.OnNext_List(this);
					}
				}

				//Simulation is done, inform the observers
				_reporters.OnCompleted_List();

				//Return the final simulation results
				return new FinalSimResults(totalNumWrong, totalBitsSimulated, noisePower, _modOriginal);
			}

			/// <summary>
			/// Run a simulation using Code division that reports its progress
			/// </summary>
			/// <param name="numberToGetWrongEventually">Number of bit errors to simulate</param>
			/// <param name="noisePower">Power of the gaussian noise in the simulation</param>
			/// <param name="user">Row of the code of the user of interest</param>
			/// <returns>Simulation results</returns>
			public FinalSimResults RunCodeDivisionSimulationObservable(int numberToGetWrongEventually, double noisePower, double[] otherUserPowers, int numberChips=31, int user = 3, byte[,] codeMatrix=null)
			{
				
				//Set up code matrices and filter coefficients
				_modOriginal.GoldFIR_DS_Coefficients = null;
				_modOriginal.GoldFIR_RF_Coefficients = null;
				_modOriginal.InitializeCodeDivision(numberChips,codeMatrix);

				byte[] inbits = new byte[_modOriginal.NumberSymbolsPerWaveform];
				double[] waveform;//, other=new double[_modOriginal.NumberSymbolsPerWaveform*_modOriginal.SamplesPerSymbol];
				byte[] outbits;
				

				int numberWrongThisIteration = 0;
				int totalNumWrong = 0, totalBitsSimulated = 0;

				//Loop while there haven't been enough bit estimation errors
				while (totalNumWrong < numberToGetWrongEventually)
				{
					//Modulate the bits
					waveform = _modOriginal.ModulateRandomBits(ref inbits);

					//Spread the waveform
					_modOriginal.SpreadWaveform(ref waveform, user);

					//Add noise to the waveform
					_modOriginal.AdditiveWhiteGaussianNoiseNR(ref waveform, noisePower);

					for (int i = 0; i < otherUserPowers.Length; i++)
					{
						_modOriginal.AddSpreadWaveform(ref waveform, user + i + 1, otherUserPowers[i]);
					}

					//Waveform now holds what the receiver what actually receive
					//Everything below is in the Receiver

					//Filter the received waveform (RF Filter)
					waveform = Filter.NewIIRFilter.RF_GoldFilter(waveform);

					//Despread the filtered signal
					_modOriginal.DespreadWaveform(ref waveform, user);

					//Filter the despread signal
					waveform = Filter.NewIIRFilter.RF_GoldFilter(waveform);

					//Demodulate the waveform+noise
					outbits = _modOriginal.CorrelationReceiver(waveform);

					//Get the total number of bits that were estimated wrongly
					numberWrongThisIteration = _modOriginal.numberDifferentBits(inbits, outbits);
					totalNumWrong += numberWrongThisIteration;

					//Update the total number of bits simulated
					totalBitsSimulated += _modOriginal.NumberSymbolsPerWaveform;

					//If there were any bit errors this iteration, report them
					if (numberWrongThisIteration != 0)
					{
						isr = new IntermediateSimResults(totalNumWrong, numberWrongThisIteration, totalBitsSimulated, numberToGetWrongEventually - totalNumWrong, noisePower.ToString());
						_reporters.OnNext_List(this);
					}
				}

				//Simulation is done, inform the observers
				_reporters.OnCompleted_List();

				//Return the final simulation results
				return new FinalSimResults(totalNumWrong, totalBitsSimulated, noisePower, _modOriginal);
			}


			protected void writeToCSV<T>(T[] mat, string filename, double sampleRate = 0)
			{
				filename = filename + ".csv";
				var csv = new StringBuilder();

				if (sampleRate != 0)
					csv.AppendLine(sampleRate.ToString());

				for (int i = 0; i < mat.Length - 1; i++)
				{
					csv.AppendLine(mat[i].ToString() + ",");
				}
				csv.AppendLine(mat[mat.Length - 1].ToString());

				File.WriteAllText(filename, csv.ToString());
			}
		}
	} 
}
