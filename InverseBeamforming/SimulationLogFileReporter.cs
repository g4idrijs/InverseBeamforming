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
		/// Class that describes the process of being subscribed to a simulation
		/// </summary>
		public class SimulationLogFileReporter : SimulationReporter
		{
			/// <summary>
			/// List of intermediate results that have not yet been written to the log file
			/// </summary>
			protected List<IntermediateSimResults> _pastISRs;

			/// <summary>
			/// Construct a new instances of the SimulationReporter class
			/// </summary>
			/// <param name="logFilename">Name of the logfile to use for this reporter</param>
			public SimulationLogFileReporter(string logFilename)
			{
				this._logFilename = logFilename;
				System.IO.Directory.CreateDirectory(Path.GetDirectoryName(logFilename));
				this._pastISRs = new List<IntermediateSimResults>();
			}

			/// <summary>
			/// Sets up the file to log the simulation
			/// </summary>
			public override void OnStart()
			{
				//Get the starting time
				startTime = DateTime.Now;
				bool notWritten = true;

				//Ensure that the starting stuff is writen to the file
				while (notWritten)
				{
					notWritten = false;
					try
					{
						using (StreamWriter file = new StreamWriter(_logFilename, false))
						{
							//Write the information to the log file
							file.WriteLine("Start of simulation. Time: " + startTime.ToLongDateString() + " " + startTime.ToLongTimeString());
						}
					}
					//An IOException means that the file is probably syncing with OneDrive or something like that
					catch (IOException)
					{
						notWritten = true;
					}
				}
			}

			/// <summary>
			/// Report the status of the simulation to the logfile
			/// </summary>
			/// <param name="sim">Simulation to report the status of</param>
			public override void OnNext(SingleSimulation sim)
			{
				//Get the intermediate results from the simulation
				IntermediateSimResults isr = sim.isr;

				//Try writing the intermediate results to the log file
				try
				{
					using (StreamWriter file = new StreamWriter(_logFilename, true))
					{
						//If there are previous ISRs that haven't been written to the file, then write them now while it is open
						foreach (var oldISR in _pastISRs.ToList())
						{
							file.WriteLine(String.Format("{1:MM/dd/yyyy HH:mm:ss:fffffff} Total Errors: {2, 5}, Errors this iteration: {3, 5}, Total Bits Simulated: {4, 9}, Bit Error Rate: {5: #.00e0}, Percent errors found: {6} %", oldISR.ToString(), DateTime.Now, oldISR.TotalErrors, oldISR.NumberErrorsThisIteration, oldISR.TotalBitsSimulated, oldISR.BitErrorRate, oldISR.PercentErrorHad.ToString("##.00").PadLeft(5)));
						}

						// Clear the list of unwritten ISRs because they just got written
						_pastISRs.Clear();

						//Write the current information to the log file
						file.WriteLine(String.Format("{1:MM/dd/yyyy HH:mm:ss:fffffff} Total Errors: {2, 5}, Errors this iteration: {3, 5}, Total Bits Simulated: {4, 9}, Bit Error Rate: {5: #.00e0}, Percent errors found: {6} %", isr.ToString(), DateTime.Now, isr.TotalErrors, isr.NumberErrorsThisIteration, isr.TotalBitsSimulated, isr.BitErrorRate, isr.PercentErrorHad.ToString("##.00").PadLeft(5)));
					}
				}
				//If there was an IOException, then the file is probably locked, so add the current ISR to the list to be added later
				catch (IOException)
				{
					_pastISRs.Add(isr);
				}
			}

			/// <summary>
			/// Report an error in the simulation to the logfile
			/// </summary>
			/// <param name="error">Exception describing the error</param>
			public override void OnError(Exception error)
			{
				bool notWritten = true;

				//Ensure that the error message is writen to the file
				while (notWritten)
				{
					try
					{
						//Set not written
						notWritten = false;

						//Write the error message and current time to the file
						using (var file = new StreamWriter(_logFilename, true))
						{
							file.WriteLineAsync(error.Message);
							file.WriteLineAsync(String.Format("An error was encoutered at {0:MM/dd/yyyy HH:mm:ss:fffff}. {1}\n", DateTime.Now, error.Message));
						}
					}
					//If there was an IOException, then the error didn't get written to the logfile
					catch (IOException)
					{
						notWritten = true;
					}
				}
			}

			/// <summary>
			/// Write completion information to the logfile and unsuscribe from the simulation
			/// </summary>
			public override void OnCompleted()
			{
				//Get the end time of the simulation
				endTime = DateTime.Now;

				bool notWritten = true;

				//Ensure that the error message is writen to the file
				while (notWritten)
				{
					try
					{
						//Set not written
						notWritten = false;

						using (StreamWriter file = new StreamWriter(_logFilename, true))
						{
							//If there are previous ISRs that haven't been written to the file, then write them now while it is open
							foreach (var oldISR in _pastISRs.ToList())
							{
								file.WriteLine(String.Format("Simulation: {0} {1:MM/dd/yyyy HH:mm:ss:fffff} Total Errors: {2, 5}, Errors this iteration: {3, 5}, Total Bits Simulated: {4, 9}, Bit Error Rate: {5: 0.00}, Percent errors found: {6: 0.00}", oldISR.ToString(), DateTime.Now, oldISR.TotalErrors, oldISR.NumberErrorsThisIteration, oldISR.TotalBitsSimulated, oldISR.BitErrorRate, oldISR.PercentErrorHad * 100));
							}

							// Clear the list of unwritten ISRs because they just got written
							_pastISRs.Clear();

							//Append some timing information to the end
							file.WriteLine(String.Format("Completed simulation at {0} {1}.\n", endTime.ToLongDateString(), endTime.ToLongTimeString()));
							file.WriteLine("The simulation took: " + (endTime - startTime).ToString());
						}
					}
					//If there was an IOException, then the error didn't get written to the logfile
					catch (IOException)
					{
						notWritten = true;
					}
				}

				//Unsubscribe from the simulation
				this.Unsubscribe();
			}
		}
	}
}
