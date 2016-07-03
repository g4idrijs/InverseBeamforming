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
	/// Holds methods that create and report on simulations of Bit error rate
	/// </summary>
	public partial class Simulations
	{
		/// <summary>
		/// Describes the current state of the simulation
		/// </summary>
		public class IntermediateSimResults
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

			public string SimulationName { get; set; }

			/// <summary>
			/// Constructs an instance of the IntermediateSimResults Structure
			/// </summary>
			/// <param name="totalErrors">Total number of error occured in the simulation</param>
			/// <param name="numberErrorsThisIteration">Total number of errors in this iteration</param>
			/// <param name="totalBitsSimulated">Total number of bits simulated so far</param>
			/// <param name="totalErrorsRemaining">Total number of errors remaining</param>
			public IntermediateSimResults(int totalErrors, int numberErrorsThisIteration, int totalBitsSimulated, int totalErrorsRemaining, string simulationName = "")
			{
				TotalErrors = totalErrors;
				NumberErrorsThisIteration = numberErrorsThisIteration;
				TotalBitsSimulated = totalBitsSimulated;
				TotalErrorsRemaining = totalErrorsRemaining;
				SimulationName = simulationName;
			}

			public override string ToString()
			{
				if (String.IsNullOrEmpty(SimulationName))
					return "";
				else
					return SimulationName;
			}
		}
	} 
}
