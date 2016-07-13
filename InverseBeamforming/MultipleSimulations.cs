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
		/// Provides methods to run multiple simulations asynchronusly
		/// </summary>
		public class MultipleSimulations
		{
			/// <summary>
			/// Enumeration that holds the types of simulations possible
			/// </summary>
			public enum ESimulationType { Simple, FIR, CDMA };

			/// <summary>
			/// Original modulation class used in the simulations
			/// </summary>
			private Modulations.ModulationType _originalModulation;

			/// <summary>
			/// Number of bit errors to simulate before stopping
			/// </summary>
			private int _numberToGetWrongEventually;

			/// <summary>
			/// Noise powers to simulate
			/// </summary>
			private double[] _noisePowers;

			/// <summary>
			/// Base of the log file for each of the simulations. Note that the noise power will be attached to the end of each simulations actual log filename
			/// </summary>
			private string _logFilename;

			/// <summary>
			/// Type of simulation to run. If an invalid type is selected, the default is to run a simple simulation
			/// </summary>
			private ESimulationType _simulationType;

			/// <summary>
			/// Number of total users of the communications system;
			/// </summary>
			protected int _numberTotalUsers;

			/// <summary>
			/// Construct a class that can run multiple simulations simultaneously
			/// </summary>
			/// <param name="modulation">Modulation to use</param>
			/// <param name="simulationType">Simulation type to run</param>
			/// <param name="noisePowers">Powers of noise to use in each simulation</param>
			/// <param name="numberToGetWrongEventually">Number of bit errors to simulate to in every simulation</param>
			/// <param name="logFilename">Log file name base to use.
			/// Note: The noise power will be appended to the end, and each simulation will get its own log file.</param>
			public MultipleSimulations(Modulations.ModulationType modulation, ESimulationType simulationType, double[] noisePowers, int numberToGetWrongEventually, string logFilename, int numberTotalUsers = 1)
			{
				_originalModulation = modulation;
				_numberToGetWrongEventually = numberToGetWrongEventually;
				_noisePowers = noisePowers;
				_logFilename = logFilename;
				_simulationType = simulationType;
				_numberTotalUsers = numberTotalUsers;
			}

			/// <summary>
			/// Runs the simulation specified by the simulation type passed into the constructor
			/// </summary>
			/// <returns>List of Final Simulation results from each of the simulations</returns>
			public async Task<List<FinalSimResults>> RunSimulations()
			{
				if(_numberTotalUsers>1)
				{
					return await RunManyCDMASimulationsObservableAsync();
				}
				//Switch on the simulation type
				switch (_simulationType)
				{
					case ESimulationType.FIR:
						return await RunManyFIRSimulationsObservableAsync();
					case ESimulationType.Simple:
						return await RunManySimpleSimulationsObservableAsync();
					case ESimulationType.CDMA:
						return await RunManyCDMASimulationsObservableAsync();
					default:
						return await RunManySimpleSimulationsObservableAsync();
				}
			}

			/// <summary>
			/// Runs many simple simulations that each report their progress to log files
			/// </summary>
			/// <param name="_numberToGetWrongEventually">Number of bit errors to simulate before stopping</param>
			/// <param name="_noisePowers">Array of noise powers to simulate</param>
			/// <returns>List of objects that contain information about each of the simulations</returns>
			public async Task<List<FinalSimResults>> RunManySimpleSimulationsObservableAsync()
			{
				//Initialize a list to hold the results
				List<FinalSimResults> results = new List<FinalSimResults>();
				var provider = new SimulationTracker[_noisePowers.Length];
				var simulations = new SingleSimulation[_noisePowers.Length];
				Modulations.ModulationType modulation;

				List<SimulationReporter>[] reporters = new List<SimulationReporter>[_noisePowers.Length];
				for (int i = 0; i < _noisePowers.Length; i++)
				{
					provider[i] = new SimulationTracker();
					reporters[i] = new List<SimulationReporter>();
					reporters[i].Add(new SimulationLogFileReporter(Path.ChangeExtension(_logFilename, null) + "_NP_" + _noisePowers[i].ToString() + ".csv"));
					reporters[i].Add(new SimulationProgressBarReporter("Noise Power:" + _noisePowers[i].ToString()));
					reporters[i].Subscribe_List(provider[i]);
				}

				//Loop through each of the noise powers
				Parallel.For(0, _noisePowers.Length, i =>
				{
					string modType = _originalModulation.GetType().Name;
					switch (modType)
					{
						case "MPSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.MPSK_Modulation(_originalModulation as Modulations.MPSK_Modulation);
							break;
						case "BPSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.BPSK_Modulation(_originalModulation as Modulations.BPSK_Modulation);
							break;
						case "MFSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.MFSK_Modulation(_originalModulation as Modulations.MFSK_Modulation);
							break;
						default:
							modulation = (Modulations.ModulationType)new Modulations.MPSK_Modulation(_originalModulation as Modulations.MPSK_Modulation);
							break;
					}

					reporters[i].OnStart_List();

					simulations[i] = new SingleSimulation(modulation, reporters[i]);
					//Add the results of the simulation to the list
					results.Add(simulations[i].RunSimpleSimulationObservable(_numberToGetWrongEventually, _noisePowers[i]));
				});

				results.Sort();

				//Return the results
				return results;
			}

			/// <summary>
			/// Runs many simulations with FIR filtering that each report their progress to log files
			/// </summary>
			/// <param name="_numberToGetWrongEventually">Number of bit errors to simulate before stopping</param>
			/// <param name="_noisePowers">Array of noise powers to simulate</param>
			/// <returns>List of objects that contain information about each of the simulations</returns>
			public async Task<List<FinalSimResults>> RunManyFIRSimulationsObservableAsync()
			{
				//Initialize a list to hold the results
				List<FinalSimResults> results = new List<FinalSimResults>();
				var provider = new SimulationTracker[_noisePowers.Length];
				var simulations = new SingleSimulation[_noisePowers.Length];
				Modulations.ModulationType modulation;

				List<SimulationReporter>[] reporters = new List<SimulationReporter>[_noisePowers.Length];
				for (int i = 0; i < _noisePowers.Length; i++)
				{
					provider[i] = new SimulationTracker();
					reporters[i] = new List<SimulationReporter>();
					reporters[i].Add(new SimulationLogFileReporter(Path.ChangeExtension(_logFilename, null) + "_NP_" + _noisePowers[i].ToString() + ".csv"));
					reporters[i].Add(new SimulationProgressBarReporter("Noise Power:" + _noisePowers[i].ToString()));
					reporters[i].Subscribe_List(provider[i]);
				}

				//Loop through each of the noise powers
				Parallel.For(0, _noisePowers.Length, i =>
				{
					string modType = _originalModulation.GetType().Name;
					switch (modType)
					{
						case "MPSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.MPSK_Modulation(_originalModulation as Modulations.MPSK_Modulation);
							break;
						case "BPSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.BPSK_Modulation(_originalModulation as Modulations.BPSK_Modulation);
							break;
						case "MFSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.MFSK_Modulation(_originalModulation as Modulations.MFSK_Modulation);
							break;
						default:
							modulation = (Modulations.ModulationType)new Modulations.MPSK_Modulation(_originalModulation as Modulations.MPSK_Modulation);
							break;
					}

					reporters[i].OnStart_List();

					simulations[i] = new SingleSimulation(modulation, reporters[i]);
					//Add the results of the simulation to the list
					results.Add(simulations[i].RunFIRFilterSimulationObservable(_numberToGetWrongEventually, _noisePowers[i]));
				});

				results.Sort();

				//Return the results
				return results;
			}

			/// <summary>
			/// Runs many simulations with FIR filtering that each report their progress to log files
			/// </summary>
			/// <param name="_numberToGetWrongEventually">Number of bit errors to simulate before stopping</param>
			/// <param name="_noisePowers">Array of noise powers to simulate</param>
			/// <returns>List of objects that contain information about each of the simulations</returns>
			public async Task<List<FinalSimResults>> RunManyCDMASimulationsObservableAsync()
			{
				//Initialize a list to hold the results
				List<FinalSimResults> results = new List<FinalSimResults>();
				var provider = new SimulationTracker[_noisePowers.Length];
				var simulations = new SingleSimulation[_noisePowers.Length];
				var otherUsersPowers = new double[_numberTotalUsers-1];
				Modulations.ModulationType modulation;

				for(int i = 0; i<otherUsersPowers.Length; i++)
				{
					otherUsersPowers[i] = _originalModulation.SignalPower;
				}

				List<SimulationReporter>[] reporters = new List<SimulationReporter>[_noisePowers.Length];
				for (int i = 0; i < _noisePowers.Length; i++)
				{
					provider[i] = new SimulationTracker();
					reporters[i] = new List<SimulationReporter>();
					reporters[i].Add(new SimulationLogFileReporter(Path.ChangeExtension(_logFilename, null) + "_NP_" + _noisePowers[i].ToString() + ".csv"));
					reporters[i].Add(new SimulationProgressBarReporter("Noise Power:" + _noisePowers[i].ToString()));
					reporters[i].Subscribe_List(provider[i]);
				}

				//Loop through each of the noise powers
				Parallel.For(0, _noisePowers.Length, i =>
				{
					string modType = _originalModulation.GetType().Name;
					switch (modType)
					{
						case "MPSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.MPSK_Modulation(_originalModulation as Modulations.MPSK_Modulation);
							break;
						case "BPSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.BPSK_Modulation(_originalModulation as Modulations.BPSK_Modulation);
							break;
						case "MFSK_Modulation":
							modulation = (Modulations.ModulationType)new Modulations.MFSK_Modulation(_originalModulation as Modulations.MFSK_Modulation);
							break;
						default:
							modulation = (Modulations.ModulationType)new Modulations.MPSK_Modulation(_originalModulation as Modulations.MPSK_Modulation);
							break;
					}

					reporters[i].OnStart_List();

					simulations[i] = new SingleSimulation(modulation, reporters[i]);

					//Add the results of the simulation to the list
					results.Add(simulations[i].RunCodeDivisionSimulationObservable(_numberToGetWrongEventually, _noisePowers[i], otherUsersPowers));
				});

				results.Sort();

				//Return the results
				return results;
			}
		}
	}
}
