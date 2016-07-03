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
		/// Describes the final state of the simulation
		/// </summary>
		public struct FinalSimResults : IComparable<FinalSimResults>
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
		/// Power of the noise in this simulation
		/// </summary>
		public double NoisePower { get; }

		/// <summary>
		/// Modulation used in this simulation
		/// </summary>
		public Modulations.ModulationType Modulation;

		/// <summary>
		/// Constructs an instance of the IntermediateSimResults Structure
		/// </summary>
		/// <param name="totalErrors">Total number of error occured in the simulation</param>
		/// <param name="totalBitsSimulated">Total number of bits simulated so far</param>
		public FinalSimResults(int totalErrors, int totalBitsSimulated, double noisePower, Modulations.ModulationType modulation)
		{
			TotalErrors = totalErrors;
			TotalBitsSimulated = totalBitsSimulated;
			NoisePower = noisePower;
			Modulation = modulation;
		}

		/// <summary>
		/// Compare two instances of the FinalSimResults
		/// </summary>
		/// <param name="other">Other instance to compare this instance to</param>
		/// <returns>this.NoisePower.CompareTo(other.NoisePower)
		/// Less than 0 = This instance is less than the other value
		/// 0 = This instance is equal to value
		/// Greater than 0 = This instance is greater than the other</returns>
		public int CompareTo(FinalSimResults other)
		{
			return this.NoisePower.CompareTo(other.NoisePower);
		}
	}
	} 
}
